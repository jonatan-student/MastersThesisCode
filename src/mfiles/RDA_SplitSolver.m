classdef RDA_SplitSolver < handle
%RDA_SPLITSOLVER  Fully‑split, unconditionally stable R–A–D integrator.
%   ½‑reaction → semi‑Lagrangian advection → PCG diffusion → ½‑reaction.
%   Tuned for Octave speed:
%     • Incomplete‑Cholesky preconditioner for PCG (builds once)
%     • Tolerance scaled with dt so you do not over‑solve when the time
%       step is large.
%     • Same PCG solve reused for A and P to trim factorisation overhead.

properties
    a, p          % fdQuant2d concentration fields
    u, v          % fdQuant2d velocity fields (lu/ts)
    D      % diffusion coef (lu²/ts)
    K      % reaction rate  (ts⁻¹)
    mask          % reaction‑active map (double or logical)
    dt     % macro time‑step (ts)

    L     % Laplacian cached once
    M1            % IC(0) preconditioner (lower‑triangular)
end

methods
    %-------------------------- constructor ---------------------------
    function obj = RDA_SplitSolver(a, p, u, v, D, K, mask, dt)
        obj.a = a; obj.p = p; obj.u = u; obj.v = v;
        obj.D = D; obj.K = K; obj.mask = mask; obj.dt = dt;

        % build sparse Laplacian once
        obj.L = a.laplaceMatrix();

        % preconditioner for PCG: IC(0) of (I - αL) with α=1 (worst‑case)
        N = numel(a.value);
        I = speye(N);
        Aic = I - obj.D*obj.dt * obj.L;   % α≈1 for safety
        obj.M1 = ichol(Aic, struct('type','ict','droptol',1e-2));
    end

    %----------------------------- step ------------------------------
    function cstep(obj)
        %% 1) first half‑reaction (exact)
        decay = exp(-0.5*obj.K*obj.dt*obj.mask);
        obj.a.value = obj.a.value .* decay;
        obj.p.value = obj.p.value + (1-decay).*obj.a.value;

        %% 2) semi‑Lagrangian advection
        obj.a.value = obj.semiLag2D(obj.a.value, obj.u.value, obj.v.value, obj.dt, obj.a.dx, obj.a.dy);
        obj.p.value = obj.semiLag2D(obj.p.value, obj.u.value, obj.v.value, obj.dt, obj.a.dx, obj.a.dy);

        %% 3) implicit diffusion via PCG  (I - α∇²) q^{n+1} = q★
        alpha = obj.D*obj.dt;
        Aop = @(x) x - alpha * (obj.L*x);
        tol = min(1e-2, 1e-4 * sqrt(obj.dt));   % always < 1  → no warning        % looser when dt is large
        maxit = 100;

        rhsA = obj.a.value(:);
        rhsP = obj.p.value(:);

        [xa,~] = pcg(Aop, rhsA, tol, maxit, obj.M1, obj.M1');
        [xp,~] = pcg(Aop, rhsP, tol, maxit, obj.M1, obj.M1');

        obj.a.value = reshape(xa, size(obj.a.value));
        obj.p.value = reshape(xp, size(obj.p.value));

        %% 4) second half‑reaction
        obj.a.value = obj.a.value .* decay;
        obj.p.value = obj.p.value + (1-decay).*obj.a.value;

        %% 5) re‑apply BCs
        obj.a.applybcs();
        obj.p.applybcs();
    end
end

methods(Static, Access = private)
    function qn = semiLag2D(q, ux, uy, dt, dx, dy)
        % First‑order semi‑Lagrangian back‑trace with bilinear interpolation.
        [ny,nx] = size(q);
        [X,Y]   = meshgrid(1:nx,1:ny);
        Xd = X - ux.*dt/dx;
        Yd = Y - uy.*dt/dy;
        Xd = min(max(Xd,1), nx);
        Yd = min(max(Yd,1), ny);
        qn = interp2(X, Y, q, Xd, Yd, 'linear', 0);
    end
end
end
