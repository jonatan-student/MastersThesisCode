function [a_ss, p_ss] = SteadyStateRDA(u_x, u_y, D, K, obstacle, reaction_mask, a_in, xmin, xmax)
    [ngrdy, ngrdx] = size(u_x);
    N = ngrdy * ngrdx;
    idx = @(i,j) (j-1)*ngrdy + i;

    % Allocate sparse vectors
    iA = []; jA = []; vA = []; bA = zeros(N,1);
    iB = []; jB = []; vB = []; bB = zeros(N,1);
    neigh = [1 0; -1 0; 0 1; 0 -1];

    for j = 1:ngrdx
        fprintf("\r working %d%%    ", j/ngrdx*100);
        for i = 1:ngrdy
            k = idx(i,j);

            if obstacle(i,j)
                % Obstacle: c = 0
                iA(end+1)=k; jA(end+1)=k; vA(end+1)=1; bA(k)=0;
                iB(end+1)=k; jB(end+1)=k; vB(end+1)=1; bB(k)=0;

            elseif j == 1
                % INLET Dirichlet at the *true* left boundary
                iA(end+1)=k; jA(end+1)=k;  vA(end+1)=1; bA(k)=a_in;
                iB(end+1)=k; jB(end+1)=k;  vB(end+1)=1; bB(k)=0;

            elseif j == ngrdx
                % OUTLET Neumann at the *true* right boundary
                iA(end+1)=k; jA(end+1)=idx(i,j-1); vA(end+1)=1;
                iB(end+1)=k; jB(end+1)=idx(i,j-1); vB(end+1)=1;
                % bB(k) already zero

            else
                % Interior fluid node
                ux = u_x(i,j);
                uy = u_y(i,j);

                if ux >= 0
                    % flow rightwards → take upstream value at j-1
                    % divergence term:  (u a)_i - (u a)_{i-1}  =>  -u*a(i,j) + u*a(i,j-1)
                    iA(end+1:end+2) = k;
                    jA(end+1:end+2) = [ k,   idx(i,j-1) ];
                    vA(end+1:end+2) = [ -ux,  ux       ];
                else
                    % flow leftwards → take upstream value at j+1
                    % divergence term:  (u a)_{i+1} - (u a)_i  =>   u*a(i,j+1) - u*a(i,j)
                    iA(end+1:end+2) = k;
                    jA(end+1:end+2) = [ k,   idx(i,j+1) ];
                    vA(end+1:end+2) = [  ux, -ux       ];
                end

                if uy >= 0
                    % flow upwards → take upstream at i-1
                    iA(end+1:end+2) = k;
                    jA(end+1:end+2) = [ k,   idx(max(i-1,1),j) ];
                    vA(end+1:end+2) = [ -uy,  uy            ];
                else
                    % flow downwards → take upstream at i+1
                    iA(end+1:end+2) = k;
                    jA(end+1:end+2) = [ k,   idx(min(i+1,ngrdy),j) ];
                    vA(end+1:end+2) = [  uy, -uy              ];
                end


                % --- a: diffusion + reaction ---
                Dloc = D;
                diagA = -4*Dloc - double(reaction_mask(i,j))*K;
                iA(end+1)=k; jA(end+1)=k; vA(end+1)=diagA;
                for d = 1:4
                    ii = i + neigh(d,1);
                    jj = j + neigh(d,2);
                    if ii<1 || ii>ngrdy || jj<1 || jj>ngrdx || obstacle(ii,jj)
                        iA(end+1)=k; jA(end+1)=k;    vA(end+1)=Dloc;
                    else
                        iA(end+1)=k; jA(end+1)=idx(ii,jj); vA(end+1)=Dloc;
                    end
                end

                if ux >= 0
                    iB(end+1:end+2) = k;
                    jB(end+1:end+2) = [ k, idx(i,j-1) ];
                    vB(end+1:end+2) = [ -ux,  ux ];
                else
                    iB(end+1:end+2) = k;
                    jB(end+1:end+2) = [ k, idx(i,j+1) ];
                    vB(end+1:end+2) = [  ux, -ux ];
                end

                if uy >= 0
                    iB(end+1:end+2) = k;
                    jB(end+1:end+2) = [ k, idx(max(i-1,1),j) ];
                    vB(end+1:end+2) = [ -uy,  uy ];
                else
                    iB(end+1:end+2) = k;
                    jB(end+1:end+2) = [ k, idx(min(i+1,ngrdy),j) ];
                    vB(end+1:end+2) = [  uy, -uy ];
                end

                % --- p: diffusion only ---
                iB(end+1)=k; jB(end+1)=k; vB(end+1)=-4*Dloc;
                for d = 1:4
                    ii = i + neigh(d,1);
                    jj = j + neigh(d,2);
                    if ii<1 || ii>ngrdy || jj<1 || jj>ngrdx || obstacle(ii,jj)
                        iB(end+1)=k; jB(end+1)=k;    vB(end+1)=Dloc;
                    else
                        iB(end+1)=k; jB(end+1)=idx(ii,jj); vB(end+1)=Dloc;
                    end
                end

                bB(k) = 0;
            end
        end
    end

    % Solve for a
    A = sparse(iA,jA,vA,N,N) + speye(N)*1e-12;
    a_ss = reshape(A \ bA, ngrdy, ngrdx);

    % Build p‐RHS and solve
    for j = 1:ngrdx
      for i = 1:ngrdy
        k = idx(i,j);
        if ~obstacle(i,j)
          bB(k) = double(reaction_mask(i,j)) * K * a_ss(i,j);
        end
      end
    end
    B = sparse(iB,jB,vB,N,N) + speye(N)*1e-12;
    p_ss = reshape(B \ bB, ngrdy, ngrdx);
    p_ss = p_ss .* -1; %%Because for some reason it is a negative thing
    fprintf(" bB positive entries:   min = %g,   max = %g\n", min(bB(bB~=0)), max(bB(:)));
    fprintf(" p_ss range: min = %g, max = %g\n", min(p_ss(:)), max(p_ss(:)));
end



