classdef fdSemiImplicit < fdIntegrator

    methods
        function this = fdSemiImplicit(dt, nquant)
            if nargin < 2, nquant = 1; end
            this.dt = dt;
            this.tnow = 0;
            this.nnow = 0;
            this.ndim = nquant;
        end

        function cphi = cstep(this, rhs, cphi, params)
            % Extract fields
            a = cphi{1}; p = cphi{2};
            ux = cphi{3}; uy = cphi{4};
            mask = cphi{5};
            D = params(1); K = params(2);

            % --- Upwind spatial derivatives for a ---
            % Compute forward and backward diffs
            dxa_f = a.calcddx('forward');
            dxa_b = a.calcddx('backward');
            dya_f = a.calcddy('forward');
            dya_b = a.calcddy('backward');
            % Choose upwind based on sign of velocity
            ddx_a = dxa_b;
            ddx_a(ux.value <= 0) = dxa_f(ux.value <= 0);
            ddy_a = dya_b;
            ddy_a(uy.value <= 0) = dya_f(uy.value <= 0);

            % --- Upwind spatial derivatives for p ---
            dxp_f = p.calcddx('forward');
            dxp_b = p.calcddx('backward');
            dyp_f = p.calcddy('forward');
            dyp_b = p.calcddy('backward');
            ddx_p = dxp_b;
            ddx_p(ux.value <= 0) = dxp_f(ux.value <= 0);
            ddy_p = dyp_b;
            ddy_p(uy.value <= 0) = dyp_f(uy.value <= 0);

            % --- Build right-hand sides with advection, reaction ---
            rhs_a = a.value + this.dt * ( ...
                -K * a.value .* mask ...        % reaction sink
                - ux.value .* ddx_a ...         % advection x
                - uy.value .* ddy_a );          % advection y

            rhs_p = p.value + this.dt * ( ...
                K * a.value .* mask ...        % reaction source
                - ux.value .* ddx_p ...         % advection x
                - uy.value .* ddy_p);          % advection y

            % --- Implicit diffusion update ---
            a_next = semi_implicit_update(a.value, rhs_a, this.dt, D, a.dx, 10);
            p_next = semi_implicit_update(p.value, rhs_p, this.dt, D, p.dx, 10);

            % Assign new values
            cphi{1}.value = a_next;
            cphi{2}.value = p_next;

            this.update();
        end
    end
end
