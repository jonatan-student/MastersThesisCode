classdef ReactionDiffusionSimulator
    %REACTIONDIFFUSIONSIMULATOR Handles reaction-advection-diffusion evolution

    properties
        a
        p
        u
        v
        K
        D
        reactionMask
        integrator
        dt
    end

    methods
        function obj = ReactionDiffusionSimulator(a_init, p_init, K, D, ux, uy, reaction_mask, dt)
            obj.a = a_init;
            obj.p = p_init;
            obj.K = K;
            obj.D = D;
            obj.u = ux;
            obj.v = uy;
            obj.reactionMask = reaction_mask;
            obj.dt = dt;
            obj.integrator = fdSemiImplicit(dt, 2)#RDA_SplitSolver(obj.a, obj.p, obj.u, obj.v, obj.D, obj.K, obj.reactionMask, obj.dt);;
        end

        function obj = evolve(obj, nsteps, visualize_every, a_initial ,xmin, xmax, ux_phys, dx_phys, D_phys, k_phys, yup, ydown, obstacle, output_dir)
            system("rm -f *jpg");
                            
            obj.a.applybcs();
            obj.p.applybcs();

            for n = 1:nsteps
                if isa(obj.integrator, 'fdSemiImplicit')
                    if n==1
                        fprintf("\r using fdSemiImplicit \n")
                    end
                    obj.integrator.cstep([], {obj.a, obj.p, obj.u, obj.v, obj.reactionMask}, [obj.D, obj.K]);
                elseif isa(obj.integrator, 'RDA_SplitSolver')
                    if n==1
                        fprintf("\r using RDA_SplitSolver \n")
                    end
                    obj.integrator.cstep();
                else
                    if n==1
                        fprintf("\r using euler integration \n")
                    end
                    obj.integrator.cstep(@ReactionDiffusionSimulator.dyneq, {obj.a, obj.p, obj.u, obj.v, obj.reactionMask}, [obj.D, obj.K]);
                end

                #obj.a.value(:, 1:2) = 0.8;
                #obj.p.value(:, 1:xmin) = 0.0;
                #obj.a.value(:, end) = obj.a.value(:, end-1);
                #obj.p.value(:, end) = obj.p.value(:, end-1);

                obstacle_r = fdObstacle2d(size(obj.a.value));
                obstacle_r.value = obstacle;
                obj.a.value = obstacle_r.correct(obj.a.value);
                obj.p.value = obstacle_r.correct(obj.p.value);

                obj.a.applybcs();
                obj.p.applybcs();
                % --- after obstacle correction -----------------
                obj.a.value(:, 1:xmin) = a_initial;   % inject A at the first fluid column
                #obj.p.value(:, 1:xmin) = 0.0;
                

                if rem(n, visualize_every) == 0
                    Visualizer.visualizeConcentrations(n, obj.a, obj.p,yup, ydown,ux_phys,dx_phys, D_phys, k_phys, xmin, xmax, visualize_every, obstacle, obj.reactionMask, output_dir);
                end

                if rem(n, 100) == 0
                    fprintf("\r Done %.1f %%  ", n/nsteps*100);
                    fflush(stdout);
                end
            end
        end

        function afinal = getFinalA(obj)
            afinal = obj.a.value;
        end

        function pfinal = getFinalP(obj)
            pfinal = obj.p.value;
        end
    end

    methods(Static)
        function r = reaction(a, k, reaction_mask)
            r = k * a.value .* reaction_mask;
        end

        function dadt = dyneq(t, quantities, params)
            a = quantities{1}; p = quantities{2};
            ux = quantities{3}; uy = quantities{4};
            reaction_mask = quantities{5};
            nu = params(1); k = params(2);

            r = ReactionDiffusionSimulator.reaction(a, k, reaction_mask);

            dadt1 = -r - ux.value .* a.calcddx('forward') - uy.value .* a.calcddy('forward') + nu * a.laplace();
            dadt2 =  r - ux.value .* p.calcddx('forward') - uy.value .* p.calcddy('forward') + nu * p.laplace();

            dadt = {dadt1, dadt2};
        end
    end
end
