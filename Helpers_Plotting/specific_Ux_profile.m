function [Ux_Prof, Ux_profile_ana, Num_fitSpace, Ux_fit, Ux_stats, yspace, y_doublespace, yFitspace] = specific_Ux_profile(ngrdy, xpos, Phys_ux, obstacle, dp ,Ls_phys, dx_phys, rho0_phys, nu_phys)
    x_pos = xpos;  
    fluid = (obstacle(:, x_pos) == 0);
    H = sum(fluid);
    Half_H = ceil(H/2);
    y_mid = ceil(ngrdy/2);
    Ls = Ls_phys / dx_phys;

    Ux_Prof = Phys_ux([y_mid-Half_H:y_mid+Half_H], x_pos);
    
    yspace = linspace(1, numel(Ux_Prof), numel(Ux_Prof));

    yFitspace = linspace(1, numel(Ux_Prof(2:end-1)), numel(Ux_Prof(2:end-1)));


    y_doublespace = linspace(-Ls, numel(Ux_Prof)+Ls+1.5, numel(Ux_Prof)*2);
    Ux_profile_ana_fn = analytical_ux(dp(x_pos), rho0_phys, nu_phys, (H)*dx_phys, Ls_phys);
    Ux_profile_ana = zeros(1, numel(Ux_Prof)*2);
    Ux_fit = zeros(1, numel(Ux_Prof(2:end-1)));

    for i = 1:numel(y_doublespace)
        Ux_profile_ana(i) = Ux_profile_ana_fn((y_doublespace(i)-1.5)*dx_phys);
    end

    f = Ux_Prof(2:end-1);
    Num_fitSpace =zeros(1, numel(f));
    for i = 1:numel(yFitspace)
        Ux_fit(i) = Ux_profile_ana_fn((yFitspace(i)-.5)*dx_phys);
        Num_fitSpace(i) = f(i);
    end

    Ux_stats = validate_profiles(Num_fitSpace, Ux_fit);
end