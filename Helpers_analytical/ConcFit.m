function Rsse = ConcFit(a_prof, ngrdx, xmin, xmax, cutout, dx_phys, a0, K_eff, D_phys, S_0, S_1, par)

    ax_fn = analytical_Conc(ngrdx, xmin, xmax, cutout, dx_phys, a0, K_eff, D_phys, S_0, S_1, par);
    ax_ana = zeros(1, ngrdx);
    for i = 1:ngrdx
        xpos = (i-xmin-cutout)/ngrdx;
        ax_ana(i) = ax_fn(xpos);
    end
    ana_slice = ax_ana((xmin+cutout):(xmax-5));
    prof_slice = a_prof(((xmin+cutout):(xmax-5)));
    Rsse = sum(((prof_slice - ana_slice)./(ana_slice)).^2);
end