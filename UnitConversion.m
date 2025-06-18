classdef UnitConversion
    %UNITCONVERSION Provides mapping between LBM and physical units

    methods(Static)

        function [dx_phys, dt_phys, achieved_dp, adjusted_deltaRho, nu_lbm] = getMapping(H_phys, dp_target, nu_phys, ngrdy, tau, rho0_phys, deltaRho)
            cs2 = 1/3;
            nu_lbm = (1/3) * (tau - 0.5);
            dx_phys = H_phys / ngrdy;
            dt_phys = dx_phys^2 * nu_lbm / nu_phys;

            if isempty(deltaRho)
                deltaRho = (dp_target * dt_phys^2) / (cs2 * rho0_phys * dx_phys^2);
            end

            achieved_dp = cs2 * deltaRho * rho0_phys * dx_phys^2 / dt_phys^2;
            adjusted_deltaRho = deltaRho;

            fprintf("Mapping for H = %.1e m, dp = %.1f Pa, nu = %.2e m^2/s:\n", H_phys, dp_target, nu_phys);
            fprintf("  dx_phys          = %.3e m (%.2f nm)\n", dx_phys, dx_phys * 1e9);
            fprintf("  dt_phys          = %.3e s (%.2f ns)\n", dt_phys, dt_phys * 1e9);
            fprintf("  deltaRho used    = %.4e\n", deltaRho);
            fprintf("  achieved dp_phys = %.2f Pa\n\n", achieved_dp);

                        % 1) channel height
            assert( abs(H_phys - ngrdy*dx_phys) < 1e-12 )

            % 2) viscosity mapping
            nu_lbm   = (tau-0.5)/3;
            nu_check = nu_phys * dt_phys / dx_phys^2;
            assert( abs(nu_lbm - nu_check) < 1e-12 )

            % 3) pressure mapping
            assert( abs(dp_target - achieved_dp)/dp_target < 1e-6 )
        end

        function phys = rescaleToPhysical(lbm, dx_phys, dt_phys, rho0_phys)
            cs2 = 1/3;
            phys = struct();

            if isfield(lbm, 'ux')
                phys.ux = lbm.ux * dx_phys / dt_phys;
            end
            if isfield(lbm, 'uy')
                phys.uy = lbm.uy * dx_phys / dt_phys;
            end
            if isfield(lbm, 'nu')
                phys.nu = lbm.nu * dx_phys^2 / dt_phys;
            end
            if isfield(lbm, 'deltaRho')
                phys.pressure = cs2 * lbm.deltaRho * rho0_phys * dx_phys^2 / dt_phys^2;
            end
            if isfield(lbm, 'rho')
                phys.pressure_field = cs2 * lbm.rho * rho0_phys * dx_phys^2 / dt_phys^2;
            end
            if isfield(lbm, 'Umax') && isfield(lbm, 'L') && isfield(lbm, 'nu')
                phys.Re = lbm.Umax * lbm.L * dx_phys / (lbm.nu * dx_phys^2 / dt_phys);
            end
        end

        function [k_phys, D_phys] = convertReactionDiffusionParams(k_lbm, D_lbm, dx_phys, dt_phys)
            % Converts LBM reaction rate and diffusion coefficient to physical units
            % Inputs:
            %   k_lbm  : LBM reaction rate (1/time_step)
            %   D_lbm  : LBM diffusion coefficient (lu^2/ts)
            %   dx_phys: physical lattice spacing (m)
            %   dt_phys: physical timestep (s)
            % Outputs:
            %   k_phys : reaction rate in 1/s
            %   D_phys : diffusion coefficient in m^2/s

            k_phys = k_lbm / dt_phys;
            D_phys = D_lbm * dx_phys^2 / dt_phys;
        end

    end
end
