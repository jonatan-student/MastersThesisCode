
classdef LBMSimulator %LBMSIMULATOR Encapsulates LBM setup and simulation steps

    properties
        tau
        Ls
        obstacle
        domainSize
        f
        integrator
        weights
        slipLength
        rhoIn
        rhoOut
        lbLoops
        output_dir
    end

    methods
        function obj = LBMSimulator(tau, Ls, obstacle, domainSize, rhoIn, rhoOut, lbLoops)
            obj.tau = tau;
            obj.Ls = Ls;
            obj.obstacle = obstacle;
            obj.domainSize = domainSize;
            obj.rhoIn = rhoIn;
            obj.rhoOut = rhoOut;
            obj.lbLoops = lbLoops;

            obj.f = fdQuant2d();
            obj.integrator = fdlb(tau, Ls, obstacle);
            obj.f.lbset(domainSize);
            obj.weights = obj.integrator.weights;
            obj.slipLength = obj.integrator.slipLength;
        end

        function obj = initializeGradient(obj)
            obj.f.value = LBMSimulator.initializePressureGradient(obj.f.value, obj.rhoIn, obj.rhoOut, obj.weights);
        end

        function [obj, results] = run(obj, xmin, xmax, yup, ydown, dx_phys, dt_phys, rho0_phys, nu_lbm, showevery, visualizeBool, output_dir, flowPars)
            FinalResults = struct();
            for n = 1:obj.lbLoops
                fprintf("\r Done %.1f %%  ", n / obj.lbLoops * 100);

                obj.f.value = LBMSimulator.setEquilibriumPressure(obj.f.value, 1:2, obj.rhoIn, obj.weights);
                %obj.f.value = LBMSimulator.setEquilibriumPressure(obj.f.value, obj.domainSize(1)-1:obj.domainSize(1), obj.rhoOut, obj.weights);

                [ux, uy, rho, obj.f] = obj.integrator.step(obj.f, obj.obstacle);

%{
                numeric_umax_lbm = max( ux(:,ceil((xmin+xmax)/2)) );
                drho = mean(rho(yup(1):ydown(1), xmin)) - mean(rho(yup(end):ydown(end), xmax));
                H_lbm = ydown((xmax -xmin)/2) - yup((xmax -xmin)/2);
                u_max_ana_lbm = (((1/3) *drho) * (H_lbm^2))/(8*nu_lbm * (xmax-xmin));
                fprintf("\r Lattice‐unit u_max = %.6g  vs. analytic = %.6g  ⇒ ratio = %.3g", ...
                        numeric_umax_lbm, u_max_ana_lbm, numeric_umax_lbm / u_max_ana_lbm);
%}

                if rem(n, showevery) == 0;
                    results.ux = ux;
                    results.uy = uy;
                    results.rho = rho;
                    results.deltaRho = obj.rhoIn - obj.rhoOut;
                    results.nu = (1/3)*(obj.tau - 0.5);
                    xmid   = ceil((xmax - xmin)/2);        % absolute column index in the domain
                    results.U_center = ux(:, xmid);
                    results.Umax_center = max(results.U_center);
                    results.Umean_center = mean(results.U_center);
                    results.L = xmax - xmin;

                    phys = UnitConversion.rescaleToPhysical(results, dx_phys, dt_phys, rho0_phys);
                        %--- after you have phys = UnitConversion… ---
                    %--- compute the correct mid‐column indices ---
                    ncols  = xmax - xmin + 1;
                    midRel = ceil(ncols/2);       % 1-based index inside the pore
                    xmid   = xmin + midRel - 1;   % absolute column in the full grid

                    % top/bottom rows of the pore at that column:
                    y_top = yup(midRel);
                    y_bot = ydown(midRel);

                    % slice the physical velocity between those rows:
                    u_slice = phys.ux(y_top:y_bot, xmid);

                    % compute centre‐line Re:
                    Phys_Re = max(u_slice) * ((y_bot - y_top) * dx_phys) / phys.nu;

                    % store for output
                    results.U_center    = u_slice;
                    results.Umax_center = max(u_slice);
                    results.L           = (y_bot - y_top)*dx_phys;  % pore height

                    if visualizeBool
                        visCache = Visualizer.visualizeFlow(phys.ux, phys.uy, yup, ydown, xmin, xmax, dx_phys, phys.nu, showevery ,rho0_phys ,obj.obstacle, obj.slipLength, phys.pressure_field, Phys_Re, false, false, output_dir, flowPars);
                    end   
                    if n == obj.lbLoops
                        visCache = Visualizer.visualizeFlow(phys.ux, phys.uy, yup, ydown, xmin, xmax, dx_phys, phys.nu, showevery ,rho0_phys ,obj.obstacle, obj.slipLength, phys.pressure_field, Phys_Re, true, true, output_dir, flowPars);
                        FinalResults.Visuals = visCache;
                    end

                end
            end


            FinalResults.ux = ux;
            FinalResults.uy = uy;
            FinalResults.Phys_ux = phys.ux;
            disp(size(phys.ux));
            disp(size(FinalResults.Phys_ux));
            FinalResults.Phys_uy = phys.uy;
            FinalResults.Phys_Re_center = Phys_Re;

            dp_Pore_real = mean(phys.pressure_field(round(yup(1)):round(ydown(1)), xmin)) - mean(phys.pressure_field(round(yup(end)):round(ydown(end)), xmax));
            disp("\nReal pressure across pore");
            disp(dp_Pore_real);
            FinalResults.pressure_field = phys.pressure_field;
            FinalResults.Pore_dP = dp_Pore_real;
            flowVarstring = flowPars.flowVarstring;
            flowVar = flowPars.flowVar;
            flowVarStr = strrep(num2str(flowVar, '%.5g'), '.', '_');
            fld = sprintf('%s/%s=%s', output_dir, flowVarstring, flowVarStr)

            save(fullfile(fld, "tmpvelresults.mat"), 'FinalResults');
        end
    end

    methods(Static)
        function f = setEquilibriumPressure(f, region_indices, rho_target, weights)
            % Octave currently has a bug when you try to assign into a 3D array
            % with a *vector* slice in one dimension.  To work around it, loop:
            for k = 1:numel(weights)
                for xi = region_indices
                    f(:, xi, k) = rho_target * weights(k);
                end
            end
        end

        function f = initializePressureGradient(f, rho_left, rho_right, weights)
            [ny, nx, ~] = size(f);
            rho_profile = linspace(rho_left, rho_right, nx);

            for x = 1:nx
                rho_x = rho_profile(x);
                for k = 1:9
                    f(:, x, k) = rho_x * weights(k);
                end
            end
        end
    end
end
