function Ux_profile_ana_fn = analytical_ux(dpdx_ana, rho0_phys, nu_phys, H, Ls_phys)
    mu_phys = nu_phys*rho0_phys;
    G = 1/(2*mu_phys) * dpdx_ana;
    Ux_profile_ana_fn = @(y) G * (y * (H - y) + Ls_phys * H);
end