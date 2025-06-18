classdef GeometryGenerator
    methods(Static)
        function [wall_fn, wall_vals] = generateSlopedWall( ...
                freqs, xmin, xmax, base_h, slope, amp_decay, rough_amp)
            if nargin<6, amp_decay = 0.3; end
            if nargin<7, rough_amp = 1.0; end
            n = xmax-xmin+1;  spectrum = zeros(1,n);
            for k = 1:numel(freqs)
                idx = freqs(k)+1; if idx>n/2, continue; end
                spectrum(idx)         = exp(-amp_decay*freqs(k))*exp(1i*2*pi*rand());
                spectrum(end-idx+2)   = conj(spectrum(idx));
            end
            detail = real(ifft(spectrum));
            if max(abs(detail))>0, detail = detail/max(abs(detail))*rough_amp; end
            xvals = xmin:xmax;  linear = base_h - slope*(xvals-xmin);
            wall_vals = linear + detail;
            wall_fn = @(x) interp1(xvals, wall_vals, x, 'linear', 'extrap');
        end

        function [obs, reac] = definePoreGeometry(nx,ny,xmin,xmax,wallT,wallB)
            obs = false(ny,nx);  reac = false(ny,nx);
            obs(:,xmin:xmax)=true;
            for x=xmin:xmax
                            % Shift the top wall half a lattice unit *into* the fluid so
            % that both upper and lower walls sit symmetrically with
            % respect to the fluid nodes.
            yt = round(wallT(x));             % Shift the bottom wall half a lattice unit *into* the fluid
            % (mirror of the top‑wall shift above).
            yb = round(wallB(x));
                yt=max(1,min(ny-1,yt)); yb=max(1,min(ny-1,yb));
                obs(yt:yb,x)=false;  reac([yt yb],x)=true;
            end
            buf=2; for x=xmax-buf:xmax
                yt=round(wallT(x)); yb=round(wallB(x));
                obs(1:max(1,yt),x)=false; obs(min(ny,yb):ny,x)=false;
            end
        end

        % -------- mask‑based q‑map: q = 0.5 for every solid‑fluid link ----
        function q_map = computeQMap(obstacle)
            if exist('computeQMap_mex','file')==3
                q_map = computeQMap_mex(obstacle);
                return; end
            [ny,nx] = size(obstacle); q_map=zeros(ny,nx,9);
            cxs=[0 0 1 1 1 0 -1 -1 -1]; cys=[0 1 1 0 -1 -1 -1 0 1];
            for k=2:9
                cx=cxs(k); cy=cys(k);
                s = obstacle;        % solid mask
                f = circshift(~obstacle,[-cy,-cx]);  % neighbour fluid mask
                q_map(:,:,k) = 0.5*(s & f);          % 0.5 where solid‑fluid
            end
        end
    end
end
