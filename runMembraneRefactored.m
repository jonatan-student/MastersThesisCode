function runMembraneRefactored(var, Varstring, output_dir)
    flowPars = struct();
    flowPars.flowVarstring = Varstring;
    flowPars.flowVar = var;
    fprintf("\n\n FlowvarString:%s, Current value: %d\n\n",  flowPars.flowVarstring, flowPars.flowVar)

    % --- Physical and simulation DEFAULT parameters ---
    H_phys = 100e-9;        % simulation grid total height default at 500 nm
    HMult = 1;               %multiplyier for height of pore
    nu_phys = 1e-6;         % Meters^2 / second
    dp_phys = 1000;            % Pascals
    rho0_phys = 1000;       % Kg/M^3
    tau = 0.995;
    Ls = 0; m = 0;
    Pe_range = linspace(0.0001, 1, 100); Da_range = linspace(0.0001, 1, 100);
    a_initial= 1.0;
    ngrdx = 100; ngrdy = 100;
    lbLoops = 20000;
    nloops = 10000; #this is redundant if mode = "steady"
    dx = 1; dt = 1; rho0 = 1; rdaSpeedup = 1;
    pure_advection = false; mode = 'steady';  % 'transient' or 'steady'
    showevery = 100; ShowVisualization = false;  %applys to LBM currently
    runLBM = true; runRDA = true;

    switch flowPars.flowVarstring
        case "Pressures"
            dp_phys = flowPars.flowVar;
        case "SlipLengths"
            Ls = flowPars.flowVar;
        case "Height_mult"
            HMult = flowPars.flowVar;
        case "Slopes"
            m = flowPars.flowVar;
    end
    H_phys = 100e-9 *HMult;        % simulation grid total height default at 500 nm



    % --- Geometry setup ---
    xmin = round(ngrdx/10); xmax = round(ngrdx-(ngrdx/5));
    ymid = round(ngrdy/2);
    freqs_top = [12, 34];
    freqs_bottom = [17, 43];
    H_b0 = ymid + ((ngrdy/10)); H_t0 = ymid - ((ngrdy/10));
    slope = 0.01*m; ampdecay = 0.3; ampScalar = 0.0;

    [wall_bot_fn, y_bot_vals] = GeometryGenerator.generateSlopedWall(freqs_bottom, xmin, xmax, H_b0, slope, ampdecay, ampScalar);
    [wall_top_fn, y_top_vals] = GeometryGenerator.generateSlopedWall(freqs_top, xmin, xmax, H_t0, -1*slope, ampdecay, ampScalar);
    [obstacle, reaction_mask] = GeometryGenerator.definePoreGeometry(ngrdx, ngrdy, xmin, xmax, wall_top_fn, wall_bot_fn);
    %q_map = GeometryGenerator.computeQMap(obstacle);

    fluid_mask = zeros([ngrdy, ngrdx]);
    fluid_mask(obstacle(:, :) == 0) = 1;

    if pure_advection
        K_phys = 0;
        D = 0;
        reaction_mask = reaction_mask* 0.0;
    end
    % --- Unit conversion ---
    [dx_phys, dt_phys, dp, deltaRho, nu_lbm] = UnitConversion.getMapping(H_phys, dp_phys, nu_phys, ngrdy, tau, rho0_phys, []);
    rho_in = rho0 + deltaRho / 2;
    rho_out = rho0 - deltaRho / 2;

    if runLBM
        % --- LBM Simulation ---
        lbm = LBMSimulator(tau, Ls, obstacle, [ngrdx, ngrdy], rho_in, rho_out, lbLoops);
        lbm = lbm.initializeGradient();
        [lbm, lbm_result] = lbm.run(xmin, xmax, y_top_vals, y_bot_vals, dx_phys, dt_phys, rho0_phys, nu_lbm ,showevery, ShowVisualization, output_dir, flowPars);
    end

    for Pe = Pe_range
        for Da = Da_range
            clf;
            flowVarstring = flowPars.flowVarstring;
            flowVar = flowPars.flowVar;
            flowVarStr = strrep(num2str(flowVar, '%.5g'), '.', '_');
            fld = sprintf('%s/%s=%s', output_dir, flowVarstring, flowVarStr);
            % --- Load velocity field ---
            out = load(fullfile(fld, "tmpvelresults.mat"));
            FinalResults = out.FinalResults;
            ux = FinalResults.ux;
            uy = FinalResults.uy;
            ux_phys = FinalResults.Phys_ux;
            uy_phys = FinalResults.Phys_uy;

            % --- Reaction-Diffusion Setup ---
            Height_phys =  (mean(H_b0) - mean(H_t0))*dx_phys;

            U_max = max(max(abs(ux_phys(round(mean(y_top_vals(:))):round(mean(y_bot_vals)), xmin:xmax))));     % get mean speed [m/s]
            #Pe     = U_max * Height_phys / D_phys;   % Peclet
            #Da     = K_phys * Height_phys / U_max;   % Damköhler
            K_phys = Da * U_max / Height_phys;
            D_phys = U_max * Height_phys / Pe;

            % Build dimensionless inputs for solver:
            u_x_nd = ux_phys ./ U_max;            % now O(1)
            u_y_nd = uy_phys ./ U_max;
            D_nd   = 1 / Pe;                      % dimensionless diffusion
            K_nd   = Da;                          % dimensionless reaction

            fprintf("\nPe: %d, Da: %d\nD_Phys: %d, K_Phys: %d\nU_max: %d\nChannel Height_phys: %d\n\n ", Pe, Da, D_phys, K_phys, U_max, Height_phys);

            % --- Reaction–diffusion setup (main.m) -----------------
            a = fdQuant2d([ngrdx, ngrdy], [dx, dx], "nnnn");
            p = fdQuant2d([ngrdx, ngrdy], [dx, dx], "nnnn");

            % NEW: give A and P explicit boundary conditions
            a.setbc('left' , 'dirichlet', a_initial);   % 0.8 at x = xmin-1
            a.setbc('right', 'neumann'  , 0        );   % ∂A/∂x = 0 at outlet
            p.setbc('left' , 'dirichlet', 0        );   % no P fed in
            p.setbc('right', 'neumann'  , 0        );   % ∂P/∂x = 0
            rda_dt=dt;

            u = fdQuant2d(); u.value = ux;
            v = fdQuant2d(); v.value = uy;
            
            if runRDA
                switch mode
                    case 'transient'
                    % --- Run Reaction-Diffusion Simulation ---
                        rdSim = ReactionDiffusionSimulator(a, p, K, D, u, v, reaction_mask, rda_dt);
                        rdSim = rdSim.evolve(nloops, 200, a_initial ,xmin, xmax, FinalResults.Phys_ux, dx_phys, D_phys, K_phys, y_top_vals, y_bot_vals, obstacle, output_dir);

                        % --- Make movie ---
                        system(sprintf("ffmpeg -r 30 -i %s/image-%%05d.jpg %s/test.mpg", output_dir, output_dir));
                        system(sprintf("ffmpeg -i %s/test.mpg -acodec copy -vcodec copy -f mp4 %s/test.mp4", output_dir, output_dir));
                

                        % --- return concentration fields ---
                        af = rdSim.getFinalA();
                        pf = rdSim.getFinalP();
                    case 'steady'
                        [af, pf] = SteadyStateRDA(u_x_nd, u_y_nd, D_nd, K_nd, obstacle, reaction_mask, a_initial, xmin, xmax);
                        % Qx already computed from LBM
                end

                % --- finding stuff --- redundant this is done better in Plotting_finals.m
                deltaA = sum(af(y_top_vals(1):y_bot_vals(1), xmin)) - sum(af(y_top_vals(end):y_bot_vals(end), xmax));
                deltaP = sum(pf(y_top_vals(end):y_bot_vals(end), xmax)) - sum(pf(y_top_vals(1):y_bot_vals(1), xmin));

                if deltaA + deltaP ~= 0.0
                    fprintf("\n \Delta_A: %d, \Delta_P: %d \n", deltaA, deltaP);
                end

                Qx = sum(sum(FinalResults.Phys_ux(:, xmin:xmax)));

                Rp = Qx * deltaA;

                % --- Storing FinalResults ---
                FinalResults.af = af;
                FinalResults.pf = pf;
                FinalResults.deltaA = deltaA;
                FinalResults.deltaP = deltaP;
                FinalResults.Qx = Qx;
                FinalResults.Rp = Rp;

                fprintf("\n Qx = %d, Rp = %d \n", Qx, Rp)


                #save("-mat7-binary", fullfile(output_dir, "Concfield.mat"), "af", "pf");
                flowVarstring = flowPars.flowVarstring;
                flowVar = flowPars.flowVar;
                flowVarStr = strrep(num2str(flowVar, '%.5g'), '.', '_');
                fld = sprintf('%s/%s=%s', output_dir, flowVarstring, flowVarStr);
                save(fullfile(fld, "FinalResults.mat"), "FinalResults");
            end

            load(fullfile(fld, "FinalResults.mat"));
            pars = struct();
            pars.Pe = Pe;
            pars.Da = Da;
            Visualizer.visualizeConcentrations( FinalResults.af, ...
                                            FinalResults.pf, ...
                                            reaction_mask, ...
                                            K_phys, true, output_dir ,pars, flowPars, Qx, Rp);
        end
    end
end