classdef Visualizer
    %VISUALIZER Static visualization utilities for flow and concentration fields
    properties (Access=private, Static)
        % cache to store data from the last LBM iteration
        cache = struct();
    end

    methods(Static)

        function VisualizerCache = visualizeFlow(ux, uy, yup, ydown, xmin, xmax, dx_phys, nu_phys, showevery, rho0_phys, obstacle, slipLength, pressure_field, reynolds, lastImage, saveFig, output_dir, flowPars)
            [ngrdy, ngrdx] = size(ux);
            xmid = round((xmin + xmax) / 2);
            ymiddle = round(ngrdy/2);
            % original shallow cutoff for plotting
            plot_cutoff = round((xmax-xmin)/10)*1+1;
            % deeper cutoff for line-fit calculation
            calc_cutoff = plot_cutoff * 3;

            clf;

            % initialize global cache
            global VisualizerCache;
            if isempty(VisualizerCache)
                VisualizerCache = struct();
            end
    
            %% quiver subsampling
            ux_red = ux(1:4:end, 1:10:end);
            uy_red = uy(1:4:end, 1:10:end);
            x = linspace(1, ngrdx, size(ux_red, 2));
            y = linspace(1, ngrdy, size(ux_red, 1));
            [X, Y] = meshgrid(x, y);

            %% compute true inlet/outlet pressures
            P_in  = mean(pressure_field(round(yup(1)):round(ydown(1)), xmin+1));
            P_out = mean(pressure_field(round(yup(end)):round(ydown(end)), xmax-1));
            dp    = P_in - P_out;

            %% build a shifted pressure for plotting
            pressure_vis = pressure_field;
            pressure_vis(obstacle == 1) = P_out;
            pressure_vis = pressure_vis - P_out;

            %% first panel: shifted pressure + quiver
            subplot(2,2,1);
            imagesc(pressure_vis);
            xlim([1, ngrdx-1]); colormap('jet'); colorbar();
            title("Velocity field due to pressure difference");
            hold on; quiver(X, Y, ux_red, uy_red, 'k'); hold off;

            %% prepare arrays for velocity and dpdx
            height = max(ydown - yup);
            ncols  = xmax - xmin;
            u_anal  = zeros(height, ncols);
            ux_profile = zeros(height, ncols);
            adjacent_to_wall_collection = false(height, ncols);
            dpdx    = zeros(1, ncols);

            %% loop: compute dpdx and analytic profiles
            for x_pos = xmin:xmax-1
                col = x_pos - xmin + 1;
                fluid_mask = (obstacle(:, x_pos) == 0);
                y_vals = find(fluid_mask);
                ux_profile(1:numel(y_vals), col) = ux(y_vals, x_pos);

                % wall-adjacency detection
                adjacent_to_wall = false(size(y_vals));
                for idx = 1:length(y_vals)
                    yv = y_vals(idx);
                    if (yv>1 && obstacle(yv-1, x_pos)) || (yv<ngrdy && obstacle(yv+1, x_pos))
                        adjacent_to_wall(idx) = true;
                    end
                end
                adjacent_to_wall_collection(1:numel(y_vals), col) = adjacent_to_wall;

                % dp/dx calculation
                P_here = mean(pressure_vis(ymiddle-5:ymiddle+5, x_pos));
                P_next = mean(pressure_vis(ymiddle-5:ymiddle+5, x_pos+1));
                dpdx(col) = (P_here - P_next) / dx_phys;

                % analytic velocity profile
                G = dpdx(col)/rho0_phys;
                y_wall_top = max(min(y_vals)-1,1);
                H = ((max(y_vals)-min(y_vals))+2)*dx_phys;
                y_rel = ((y_vals) - y_wall_top)*dx_phys;
                Ls_phys_top = slipLength(y_wall_top, x_pos)*dx_phys;
                u_anal(1:numel(y_vals), col) = (G/(2*nu_phys)) * (y_rel .* (H - y_rel) + Ls_phys_top * H);
            end

            %% second panel: center-line velocity profile
            mid = round(ncols/2);
            y_mid = find(obstacle(:, xmin+mid-1)==0);
            prof  = ux_profile(1:numel(y_mid), mid);
            subplot(2,2,2); hold on;
                plot(prof, y_mid, 'b-', 'LineWidth',2);
                plot(prof(adjacent_to_wall_collection(1:numel(y_mid), mid)), y_mid(adjacent_to_wall_collection(1:numel(y_mid), mid)), 'ro', 'MarkerFaceColor','r');
                plot(u_anal(1:numel(y_mid), mid), y_mid, 'k--', 'LineWidth',1.5);
            set(gca,'YDir','reverse'); xlabel('u_x [m/s]'); ylabel('y-index');
            title('Center-line velocity profile'); grid on; axis tight; hold off;

            %% third panel: dpdx numeric vs analytic (with two cutoffs)
            n = numel(dpdx);
            dpdx_plot = dpdx(plot_cutoff : n-plot_cutoff);
            np = numel(dpdx_plot);
            x_plot = (0:np-1) * dx_phys;
            % compute analytic straight line from deep cutoff
            alpha = dpdx(calc_cutoff);
            beta  = -1* (dpdx(calc_cutoff) - dpdx(n-calc_cutoff)) / (n - 2*calc_cutoff);
            dpdx_fn = @(x) alpha + beta * x;
            % evaluate analytic at shifted indices for plotting
            shift = plot_cutoff - calc_cutoff;
            dpdx_anal_plot = dpdx_fn((0:np-1) + shift);

            subplot(2,2,3); hold on;
                plot(x_plot, dpdx_plot,      'b-', 'LineWidth',1.5);
                plot(x_plot, dpdx_anal_plot, 'k--','LineWidth',1.5);
            hold off;
            title(sprintf('Pressure drop over pore: %.3g Pa', dp));
            xlabel('x [m]'); ylabel('dp/dx [Pa/m]'); xlim([0, x_plot(end)]); grid on;
            % compute entranceHeight & slope from yup/ydown
            H_in = (ydown(calc_cutoff) - yup(calc_cutoff) + 1) * dx_phys;
            H_out = (ydown(end-calc_cutoff) - yup(end-calc_cutoff) + 1) * dx_phys;
            SLOPE = (H_out - H_in) / ((xmax - xmin) * dx_phys);
            H_x = @(x) H_in + (SLOPE*x);

            %% Panel 4: total flux comparison (using plot_cutoff)
            y_phys = (0:ngrdy-1) * dx_phys;
            ix_range = xmin+plot_cutoff : xmax-plot_cutoff;
            nsh = numel(ix_range);
            S_x_sh = nan(1, nsh);
            S_ana_sh = nan(1, nsh);
            S_ana_approx = nan(1,nsh);
            for k = 1:nsh
                ix = ix_range(k);
                fluid = find(obstacle(:, ix)==0);
                if isempty(fluid), continue; end
                h_sz = H_x(k*dx_phys);
                Ls_top = slipLength(max(min(fluid)-1,1), ix) * dx_phys;
                Ls_bot = slipLength(min(max(fluid)+1,ngrdy), ix) * dx_phys;
                Ls_avg = 0.5 * (Ls_top + Ls_bot);
                S_x_sh(k)   = 1/h_sz * trapz(y_phys(fluid), ux(fluid, ix));
                mu = nu_phys * rho0_phys;
                S_0_ana = ((alpha * H_in^2)+(alpha*6*Ls_avg*H_in));
                S_1_ana = ((beta*H_in^2) + (6*beta*H_in*Ls_avg)- (alpha *H_in*SLOPE)+(alpha*6*Ls_avg));
                S_ana_approx(k) = (1/(12*mu))*(S_0_ana + (S_1_ana*(k*dx_phys)));
                S_ana_sh(k) = (dpdx_anal_plot(k) / (12 * mu)) * (h_sz^2 + 6 * Ls_avg * h_sz);
            end

            subplot(2,2,4);
            hold on;
            plot(x_plot, S_x_sh, 'b-', 'LineWidth',1.5);
            #plot(x_plot, S_ana_sh, 'k--','LineWidth',1.5);
            plot(x_plot, S_ana_approx, 'k--', 'LineWidth', 1.5);
            hold off;
            xlabel('x [m]'); ylabel('S(x)');
            legend('numeric','analytic');
            title(sprintf('Total flux  Q_{num}=%.3g   Q_{ana}=%.3g m^3/s', trapz(x_plot,S_x_sh), trapz(x_plot,S_ana_sh)));
            grid on; axis tight;
            %ylim([1e-8, max(prof)*2])
            %pause(0.0001)

            % cache relevant data
            global VisualizerCache;
            VisualizerCache.dpdx            = dpdx;
            VisualizerCache.dpdx_anal_plot = dpdx_anal_plot;
            VisualizerCache.x_plot         = x_plot;
            VisualizerCache.S_x_sh         = S_x_sh;
            VisualizerCache.S_ana_sh       = S_ana_sh;
            VisualizerCache.y_phys         = y_phys;
            VisualizerCache.xmin           = xmin;
            VisualizerCache.xmax           = xmax;
            VisualizerCache.dx_phys        = dx_phys;
            VisualizerCache.nu_phys        = nu_phys;
            VisualizerCache.rho0_phys      = rho0_phys;
            VisualizerCache.slipLength     = slipLength;    % input scalar
            VisualizerCache.entranceHeight = H_in;
            VisualizerCache.slope = SLOPE;
            VisualizerCache.reynolds = reynolds;
            VisualizerCache.obstacle = obstacle;
            VisualizerCache.plot_cutoff = plot_cutoff;
            VisualizerCache.calc_cutoff = calc_cutoff;
            cache = VisualizerCache;

            % save to folder if requested
            if lastImage && saveFig
                flowVarstring = flowPars.flowVarstring;
                flowVar = flowPars.flowVar;
                flowVarStr = strrep(num2str(flowVar, '%.5g'), '.', '_');
                fld = sprintf('%s/%s=%s', output_dir, flowVarstring, flowVarStr);
                if ~exist(fld,'dir'), mkdir(fld); end
                save(fullfile(fld,'flowCache.mat'),'VisualizerCache');
                set(gcf, 'PaperPositionMode', 'auto');
                % export at 300 dpi:
                print(gcf, fld, '-dpng', '-r300');

                fprintf('\n[Visualizer] Saved flow data to:\n%s \nand figure to:\n%s\n', fld, fld);
            end
        end

        function visualizeConcentrations(af, pf, reaction_mask, K_phys, saveFig, output_dir ,params, flowPars, Qx, Rp)
            % VISUALIZECONCENTRATIONS Plot concentrations using cached flow data
            if nargin < 5, saveFig = false; end
            if nargin < 6, params = struct(); end

            global VisualizerCache;
            if isempty(VisualizerCache) || ~isfield(VisualizerCache,'x_plot')
                error('Missing flow cache: run visualizeFlow(...,lastImage,saveFig) first.');
            end
            c = VisualizerCache;
            c.Qx = Qx;
            c.Rp = Rp;
            Pe = params.Pe;
            Da = params.Da;
            % geometry and cutoffs
            RDAPLOTcutoff = c.plot_cutoff;
            y_phys = c.y_phys;
            x_rel  = ((c.xmin+RDAPLOTcutoff : c.xmax-RDAPLOTcutoff) - (c.xmin+RDAPLOTcutoff)) * c.dx_phys;
            ncol   = numel(x_rel);
            L_phys = x_rel(end);
            mu = c.nu_phys * c.rho0_phys;

            % flux arrays
            S_num = c.S_x_sh;
            S_lub = c.S_ana_sh;

            % numeric concentration sum
            a_num = sum(af(:, c.xmin+RDAPLOTcutoff : c.xmax-RDAPLOTcutoff), 1);
            a0    = a_num(1);

            % effective rate constant
            k_eff = nan(1, ncol);
            for j = 1:ncol
                ix = c.xmin + j - 1;
                fluid = find(c.obstacle(:,ix)==0);
                if isempty(fluid), continue; end
                areaFluid  = trapz(y_phys(fluid), double(~c.obstacle(fluid,ix)));
                perimReact = trapz(y_phys(fluid), double(reaction_mask(fluid,ix)));
                k_eff(j)   = 18 * K_phys * (perimReact / areaFluid);
            end

            % plug-flow analytic
            a_semiAna = nan(1, ncol);
            a_ana     = nan(1, ncol);
            a_semiAna(1) = a0;
            a_ana(1)     = a0;
            D_phys = c.nu_phys;
            for j = 2:ncol
                S_here = [S_num(j); S_lub(j)];
                for type = 1:2
                    c1 = S_here(type) / D_phys;
                    c2 = k_eff(j) / D_phys;
                    disc = sqrt(c1^2 + 4*c2);
                    r1 = (c1 + disc)/2;
                    r2 = (c1 - disc)/2;
                    denom = r1*exp(r1*L_phys) - r2*exp(r2*L_phys);
                    Acoef = -a0 * (r2*exp(r2*L_phys)) / denom;
                    Bcoef =  a0 * (r1*exp(r1*L_phys)) / denom;
                    if type == 1
                        a_semiAna(j) = Acoef*exp(r1*x_rel(j)) + Bcoef*exp(r2*x_rel(j));
                    else
                        a_ana(j)     = Acoef*exp(r1*x_rel(j)) + Bcoef*exp(r2*x_rel(j));
                    end
                end
            end


            % plotting
%{
            figure;
            subplot(2,2,1);
            colormap('jet'); imagesc(af); axis image; title('Steady A'); colorbar;
            subplot(2,2,2); hold on;
            plot(x_rel, a_num,     'r-',  'LineWidth',1.5);
            plot(x_rel, a_semiAna, 'g-.', 'LineWidth',1.5);
            plot(x_rel, a_ana,     'b--', 'LineWidth',1.5);
            xlabel('x [m]'); ylabel('<A>'); legend('numeric','semi','analytic'); grid on; hold off;
            subplot(2,2,3);
            colormap('jet'); imagesc(pf); axis image; title('Steady P'); colorbar;
%}
            % save if requested
            if saveFig
                flowVarstring = flowPars.flowVarstring;
                flowVar = flowPars.flowVar;
                flowVarStr = strrep(num2str(flowVar, '%.5g'), '.', '_');
                fld = sprintf('%s/%s=%s', output_dir, flowVarstring, flowVarStr);
                if ~exist(fld,'dir'), mkdir(fld); end
                Pe_string = strrep(num2str(Pe, '%.5g'), '.', '_');
                Da_string = strrep(num2str(Da, '%.5g'), '.', '_');
                dataFile = fullfile(fld, sprintf('conc_Pe%s_Da%s.mat', Pe_string, Da_string));
                %figFile  = fullfile(fld, sprintf('concPlot_Pe%s_Da%s.png',  Pe_string, Da_string));
                save(dataFile,'af','pf','c','params');
                %set(gcf, 'PaperPositionMode', 'auto');
                % export at 300 dpi:
                %print(gcf, figFile, '-dpng', '-r300');
                %fprintf('[Visualizer] Saved concentration data to %s and figure to %s\n', dataFile, figFile);
            end
        end

    end
end
