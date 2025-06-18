function [these_stats, Rp, Rescaler, xScalar] = plot_specific(stringInputs, chosen_plots, opts)
    if nargin<3
        opts = struct('fontsize', 16);
    end
    FONTSIZE = opts.fontsize;
    %%%Loading inputs
    basePath   = stringInputs{1};
    Folders = stringInputs{2};
    Target_test = stringInputs{3};

    %%%%loading specified data for plot displaying names for debug Purposes
    results = load_data(basePath, Folders, Target_test);
    keys = fieldnames(results);
    VelRes = results.(keys{1}).Velocity.tmpvelresults.FinalResults;
    concRes = results.(keys{1}).Concentration.(Target_test(1:end-4));
%{
    disp("VelRes Fieldnames")
    disp(fieldnames(VelRes))
    disp("concRes Fieldnames")
    disp(fieldnames(concRes))
    disp("concRes.c Fieldnames")
    disp(fieldnames(concRes.c))
%}

    %%Loading important values from results
    Phys_ux = VelRes.Phys_ux;        Phys_uy = VelRes.Phys_uy;
    rho0_phys = concRes.c.rho0_phys; nu_phys = concRes.c.nu_phys;

    mu_phys = rho0_phys*nu_phys;

    dpdx = concRes.c.dpdx;           dx_phys = concRes.c.dx_phys
    obstacle = concRes.c.obstacle; xmin = concRes.c.xmin;  xmax = concRes.c.xmax;
    af = concRes.af; pf = concRes.pf; slipLength = concRes.c.slipLength;
    Ls = mean(slipLength(find(slipLength ~=0)));
    if Ls ~= 0 && any(Ls)
        Ls_phys = Ls*dx_phys;
    else
        Ls_phys =0;
        Ls = 0;
    end

    plot_cutoff = concRes.c.plot_cutoff;

    Umax = max(max(Phys_ux(:, [xmin:xmax])));
    Pe = concRes.params.Pe; Da = concRes.params.Da;


    %%important stuff to define pressure used in other analytical solutions so needed to be a bit more global
    [ngrdx, ngrdy] = size(Phys_ux);
    xspace = linspace(0, ngrdx, ngrdx);
    cutoffin = ceil(ngrdx/2)-20;
    cutoffout = ceil(ngrdx/2)+20;
    xScalar = dx_phys/1e-9%%gives a scale to multiply a nondimensional xspace to nanometers
    Heights = zeros(1, ngrdx);
    for i = 1:ngrdx
        fluid_here = (obstacle(:, i) == 0);
        H = sum(fluid_here);
        H = H*dx_phys;
        Heights(i) = H;
    end

    %%%__________Pressure calculation______________
    slope = 0; (Heights(xmin)-Heights(xmax-5))/((xmax-xmin-5)*dx_phys)
    
    if slope ~= 0
        cutoffin = 1;
        cutoffout = 60;
    end

    len = (cutoffout-cutoffin);
    inlet_dpdx = dpdx(cutoffin);
    outlet_dpdx= dpdx(cutoffout);
    [dpdx_ana_fn, alpha, beta]= analytical_dpdx(inlet_dpdx, outlet_dpdx, len, slope);
    dpdx_ana = zeros(1, ngrdx);
    dpdx_num = nan(1, ngrdx);
    for i = 1:ngrdx
        xpos = i-xmin;
        dpdx_ana(i) = dpdx_ana_fn(xpos);
        if xpos > 0 && xpos <= numel(dpdx)
            dpdx_num(i) = dpdx(xpos);
        end
    end


    xpos = floor(ngrdx/2);
    test = dpdx_num(xpos)
    test= dpdx_ana(xpos)
    [Ux_Prof, Ux_profile_ana, num_fit, Ux_fit ,Ux_stats, yspace, y_doublespace, yFitspace] = specific_Ux_profile(ngrdy, xpos, Phys_ux, obstacle, dpdx_ana, Ls_phys, dx_phys, rho0_phys, nu_phys);

    %%%%____________S(X): AVERAGE VELOCITY IN X_________________
    S= zeros(1, ngrdx);
    Q_2d = zeros(1, ngrdx);
    for i = 1:ngrdx
        mask = (obstacle(:, i) == 0);
        y_phys = (find(mask))*dx_phys;
        S(i) = 1/Heights(i) * trapz(y_phys, Phys_ux(mask, i));
        Q_2d(i) = trapz(y_phys, Phys_ux(mask,i));
    end

    dw = 1; % unit length
    Q = Q_2d(xmax-5)*dw; %m^3/s

    [sx_ana_fn, S_0, S_1] = analytical_sx(slope, Heights(xmin), alpha, beta, mu_phys, Ls_phys, dx_phys);
    S_ana = zeros(1, ngrdx);
    for i = 1:ngrdx
        xpos = (i-xmin);
        S_ana(i) = sx_ana_fn(xpos);
    end
    S_0 = S_ana(xmin);

    %%%______CONCENTRATION A PROFILE__________
    Cref =1; %mol/m3
    a_prof = zeros(1,ngrdx);
    for i = 1:ngrdx
        mask = (obstacle(:, i)==0);
        y_phys = (find(mask))*dx_phys;
        a_prof(i) = trapz(y_phys, af(mask, i));
    end
    cutout = 0
    a0 = a_prof(xmin+cutout);
    aL = a_prof(xmax-5);
    deltaA = a0-aL;
    Rp = Q * (deltaA)* Cref * dw; %if Cref and dw dimensionless then in units of m^2/s multiply by a concentration Cref(mol/m^3) and dw(m) to return mol/s;

    fprintf("\n\n___________\nQ: %d, Delta a: %d, Rp: %d\n____________\n\n", Q, deltaA, Rp);

    H_phys_0 = Heights(xmin+5);
    K_phys = Da * Umax /H_phys_0;
    L = (xmax-xmin-(cutout+5));
    Rescaler = S_0 / Umax * ((L*dx_phys)/H);

    effectivity = (2*(dx_phys*L))/(H_phys_0*L)
    K_eff = K_phys * effectivity

    D_phys = Umax * H_phys_0/Pe;
    fprintf("_____________\n\n Pe: %d -> D_phys: %d\nDa: %d -> K_phys: %d\n\n______________\n", Pe, D_phys, Da, K_phys)
    plotFit = true;
    if plotFit
        fitPar = fminsearch(@(par) ConcFit(a_prof, ngrdx, xmin, xmax, cutout, dx_phys, a0, K_eff, D_phys, S_0, S_1, par), 2.0, optimset('Display','iter','TolX',1e-6));
        fprintf('\nOptimal Param = %g\n', fitPar)
    else 
        fitPar = 1;
    end
    ax_fn = analytical_Conc(ngrdx, xmin, xmax, cutout, dx_phys, a0, K_eff, D_phys, S_0, S_1, fitPar);
    ax_ana = zeros(1, ngrdx);

    for i = 1:ngrdx
        xpos = (i-xmin-cutout)/ngrdx;
        ax_ana(i) = ax_fn(xpos);
    end

    stats_A_prof = validate_profiles(a_prof([xmin:xmax]), ax_ana([xmin:xmax]));

    % ---------- gather validation numbers for this folder -------------------
    these_stats = collect_metrics( ...
        num_fit, Ux_fit, ...
        S([xmin:xmax-5]),       S_ana([xmin:xmax-5]), ...
        dpdx_num([xmin+20:xmax-20]),dpdx_ana([xmin+20:xmax-20]), ...
        a_prof([(xmin+cutout):(xmax-5)]),  ax_ana([(xmin+cutout):(xmax-5)]) ...
    );
    %disp(these_stats)

    if ismember("all", chosen_plots);
        chosen_plots = {"Velocity Field", "Velocity Profile", "Pressure in x", "S of x", "Concentrations", "Concentration a Profile"};
    end

    for i = 1:length(chosen_plots)
        plot_type = chosen_plots{i};
        switch plot_type
            case "Velocity Field"
            %%%%________________Velocity Field___________________
                %% quiver subsampling
                ux_red = Phys_ux(1:8:end, 1:8:end);
                uy_red = Phys_uy(1:8:end, 1:8:end);
                x = linspace(1, ngrdx, size(ux_red, 2));
                y = linspace(1, ngrdy, size(ux_red, 1));
                [X, Y] = meshgrid(x, y);

                mag = sqrt(Phys_ux.^2 + Phys_uy.^2);
                figure;
                hold on;
                colormap("jet")
                imagesc(mag);
                quiver(X, Y, ux_red, uy_red, 'k')
                hold off;
                drawnow("expose")
                %caxis([0, 1]);
                colorbar();
                axis image;
                waitforbuttonpress;

            case "Velocity Profile"
                clf;
                figure('position', [1 1, 1460, 840]);
                subplot(1,2,1);
                hold on;
                plot(yspace.*xScalar, Ux_Prof, 'b', 'LineWidth', 1.5);
                plot(y_doublespace.*xScalar, Ux_profile_ana,'Color', 'm' ,'LineStyle', '--', 'LineWidth', 1.5);
                line([-Ls numel(Ux_Prof)+Ls], [0 0], 'Color','k','LineStyle',':','LineWidth',1);
                line([1.5*xScalar 1.5*xScalar], [0 max(Phys_ux(:))], 'Color','k','LineStyle',':','LineWidth',1);
                line([(numel(Ux_Prof)-0.5)*xScalar numel(Ux_Prof)-0.5*xScalar], [0 max(Phys_ux(:))], 'Color','k','LineStyle',':','LineWidth',1);
                
                axis tight;
                for ax = findobj(gcf,'Type','Axes')'
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                title("Velocity Profiles", 'Fontsize', FONTSIZE, 'interpreter', 'latex');
                legend("Computational", "Analytical", "Channel Walls", 'Fontsize', FONTSIZE, 'interpreter', 'latex')
                xlabel('y (nm)', 'fontsize', FONTSIZE, 'interpreter', 'latex');
                ylabel('$u_x$ (m/s)','fontsize', FONTSIZE, 'interpreter', 'latex');
                ylim([0, max(Phys_ux(:))]);
                xlim([-Ls*xScalar, (numel(Ux_Prof)+Ls)*xScalar]);
                grid on;
                

                subplot(1,2,2);
                hold on;
                plot(yFitspace*xScalar, Ux_stats.relError,'color', 'b', 'LineWidth', 1.5, 'linestyle', ':')
                for ax = findobj(gcf,'Type','Axes')'
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                title("Relative Error Profile", 'Fontsize', FONTSIZE, 'interpreter', 'latex');
                xlabel('y (nm)', 'fontsize', FONTSIZE, 'interpreter', 'latex');
                ylabel('Relative Error','fontsize', FONTSIZE,'interpreter', 'latex');
                xlim([0 max(yFitspace)*xScalar])
                grid on;
                

                waitforbuttonpress;
            case "Velocity validation"
                MAPE_array = zeros(1, xmax-xmin-5);
                MAR_array = zeros(1, xmax-xmin-5);
                idx = 0
                for i = xmin+1:xmax-4
                    idx = idx +1;
                    [U_Prof, U_profile_ana, Numerical_fit ,U_fit, U_stats, yspace, y_doublespace, YFit] = specific_Ux_profile(ngrdy, i, Phys_ux, obstacle, dpdx_ana, Ls_phys, dx_phys, rho0_phys, nu_phys)
                    %disp(U_stats)
                    MAPE_array(idx) = U_stats.MAPE;
                    MAR_array(idx) = U_stats.MaxAbsRel;
                end
                
                disp(MAPE_array)
                %figure('Name', 'U_x profile Metrics', 'Position', [100 100 800 900]);

                figure('position', [1 1, 1460, 840]);
                subplot(1,2,1);
                hold on;
                sp = linspace(0, numel(MAPE_array), numel(MAPE_array));
                plot(sp*xScalar, MAPE_array, 'r', 'LineWidth', 2, 'Linestyle', ':');
                for ax = findobj(gcf,'Type','Axes')'
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                title("MAPE Distribution", 'Fontsize', FONTSIZE, 'interpreter', 'latex');
                xlabel('x (nm)', 'fontsize', FONTSIZE, 'interpreter', 'latex');
                ylabel('MAPE (\%)','fontsize', FONTSIZE,'interpreter', 'latex');
                grid on;
                xlim([0 max(sp)*xScalar]);
                %ylim([0 1.2])
                %xlim([ numel(R2_array)]);

                subplot(1,2,2);
                hold on;
                plot(sp*xScalar, MAR_array(:), 'r', 'LineWidth', 2, 'Linestyle', ':');
                for ax = findobj(gcf,'Type','Axes')'
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                title("MARE Distribution", 'Fontsize', FONTSIZE, 'interpreter', 'latex');
                xlabel('x (nm)', 'fontsize', FONTSIZE, 'interpreter', 'latex');
                ylabel('MARE ($Rel_{err}$)','fontsize', FONTSIZE,'interpreter', 'latex');
                grid on;
                xlim([0 max(sp)*xScalar]);
                %axis tight;
                drawnow("expose");
                waitforbuttonpress;

            case "Pressure in x"
                %%%%______________Pressure in X___________________

                clf;
                figure('position', [1 1, 1460, 840]);
                subplot(1, 2, 1);
                hold on;
                plot(xspace*xScalar, dpdx_num, "b", 'LineWidth',1.5 )
                plot(xspace*xScalar, dpdx_ana, "k--", 'LineWidth',1.5)
                line([xmin xmin], [min(dpdx_num) max(dpdx_num)], 'Color','r','LineStyle',':','LineWidth',1)
                line([xmax xmax], [min(dpdx_num) max(dpdx_num)], 'Color','r','LineStyle',':','LineWidth',1)
                for ax = findobj(gcf,'Type','Axes')'
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                legend('Computational', 'Analytical', 'Boundary', 'Fontsize', FONTSIZE, 'interpreter', 'latex');
                xlim([0, ngrdx]);
                ylabel('dpdx', "Fontsize", FONTSIZE, 'interpreter', 'latex')
                xlabel('x (nm)', "Fontsize", FONTSIZE, 'interpreter', 'latex');
                title("Pressure gradient in pore","Fontsize", FONTSIZE, 'interpreter', 'latex');
                grid on;
            
                subplot(1, 2, 2);
                n_x = numel(these_stats.relErr_dp);
                p_sp = linspace(0, n_x, n_x);
                hold on;
                plot(p_sp*xScalar, these_stats.relErr_dp, "b", 'LineStyle', ':', 'LineWidth', 1.5)
                for ax = findobj(gcf,'Type','Axes')
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                xlabel('x (nm)',"Fontsize", FONTSIZE, 'interpreter', 'latex');
                ylabel('$Rel_{err}$',"Fontsize", FONTSIZE, 'interpreter', 'latex')
                title("Relative Error Profile", "Fontsize", FONTSIZE, 'interpreter', 'latex')
                xlim([0, max(p_sp)*xScalar])
                grid on;
                waitforbuttonpress;

            case "S of x"

                clf;
                figure('position', [1 1, 1460, 840]);
                subplot(1, 2, 1)
                hold on;
                plot(xspace*xScalar, S, 'b', "LineWidth", 1.5)
                plot(xspace*xScalar, S_ana, 'm--', "LineWidth", 1.5)
                for ax = findobj(gcf,'Type','Axes')
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                xlabel('x (nm)',"Fontsize", FONTSIZE, 'interpreter', 'latex');
                ylabel('$S(x)$ (m/s)',"Fontsize", FONTSIZE, 'interpreter', 'latex')
                title("$S(x)$ Profile", "Fontsize", FONTSIZE, 'interpreter', 'latex')
                ylim([0, max(S_ana)+0.00005])
                grid on;

                subplot(1,2,2)
                n_s = numel(these_stats.relErr_S);
                s_sp = linspace(0, n_s, n_s);
                hold on;
                plot(s_sp, these_stats.relErr_S, 'b', 'LineStyle', ':', 'LineWidth', 1.5)
                for ax = findobj(gcf,'Type','Axes')
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                xlabel('x (nm)',"Fontsize", FONTSIZE, 'interpreter', 'latex');
                ylabel('$Rel_{err}$',"Fontsize", FONTSIZE, 'interpreter', 'latex')
                title("Relative Error Profile", "Fontsize", FONTSIZE, 'interpreter', 'latex')
                xlim([0, max(s_sp)*xScalar])
                grid on;

                waitforbuttonpress;

            case "Concentrations"
                c_bo = 3;
                af_bo = af(:, 1:end-c_bo);
                pf_bo = pf(:, 1:end-c_bo);
                [nx, ny] = size(af_bo)
                x_sp = linspace(0, nx, nx)*xScalar;
                y_sp = linspace(0, ny, ny)*xScalar;
                [X, Y] = meshgrid(y_sp, x_sp);
                clf;
                figure('position', [1, 1, 1460, 840]);
                colormap('jet');
                subplot(1,2,1);
                hold on;
                contourf(X, Y, af_bo ./ (max(max(af_bo))));
                hold off;
                for ax = findobj(gcf,'Type','Axes')'
                    set(ax,'fontsize', FONTSIZE-4);
                end
                xlabel('x (nm)', 'fontsize', FONTSIZE, 'interpreter', 'latex');
                ylabel('y (nm)','fontsize', FONTSIZE, 'interpreter', 'latex');
                title('Concentration of $a$', 'fontsize', FONTSIZE, 'interpreter', 'latex');
                axis tight;

                subplot(1, 2, 2);
                hold on;
                contourf(X, Y, pf_bo./ (max(max(af_bo))));
                hold off;
                for ax = findobj(gcf,'Type','Axes')
                    set(ax,'fontsize', FONTSIZE-4);
                end
                xlabel('x (nm)', 'fontsize', FONTSIZE, 'interpreter', 'latex');
                ylabel('y (nm)','fontsize', FONTSIZE, 'interpreter', 'latex');
                title('Concentration of $b$', 'fontsize', FONTSIZE, 'interpreter', 'latex')
                axis tight;
                waitforbuttonpress;
                
            case "Concentration a Profile"
                clf;
                figure('position', [1, 1, 1460, 840]);
                subplot(1,2,1)
                sp = linspace(0, numel(a_prof), numel(a_prof));
                hold on;
                plot(sp*xScalar, a_prof, 'b', 'LineWidth', 1.5);
                plot(sp*xScalar, ax_ana, 'm--', 'LineWidth' , 1.5);
                %line([0 ngrdx], [a0 a0], 'Color','k','LineStyle',':','LineWidth', 1)
                for ax = findobj(gcf,'Type','Axes')
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                xlim([(xmin-cutout)*xScalar, (xmax-4)*xScalar]);
                ylim([a_prof(xmax-5)-(2*dx_phys), a_prof(xmin+cutout)]);
                title("Concentration a(x)", 'Fontsize', FONTSIZE, "interpreter", "latex")
                xlabel("x (nm)", 'Fontsize', FONTSIZE, 'interpreter', 'latex')
                ylabel("a(x) (1/m)",'Fontsize',FONTSIZE , 'margin', 1, 'interpreter', 'latex')
                grid on;

                subplot(1,2,2)
                e_sp = linspace(0, numel(these_stats.relErr_a),  numel(these_stats.relErr_a));
                hold on;
                plot(e_sp*xScalar, these_stats.relErr_a, "b", 'LineStyle', ':', 'LineWidth', 1.5);
                for ax = findobj(gcf,'Type','Axes')
                    set(ax,'fontsize', FONTSIZE-6);
                end
                hold off;
                xlim([0,  (numel(these_stats.relErr_a)-5)*xScalar])
                title("Relative error distribution", "Fontsize", FONTSIZE, "interpreter", "latex")
                xlabel("x (nm)", 'Fontsize', FONTSIZE, 'interpreter', 'latex')
                grid on;
                waitforbuttonpress;
        end

    end
end