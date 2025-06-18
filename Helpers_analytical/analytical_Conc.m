function ax_fn = analytical_Conc(ngrdx, xmin, xmax, cutoffin, dx_phys, a0, K_eff, D_phys, S_0, S_1, fitpar)
    K_eff = K_eff*(fitpar)
    daleth_0 = S_0/D_phys;
    daleth_1 = S_1/D_phys;
    lambda = K_eff/D_phys;
    L  = ((xmax-5)-(xmin+cutoffin))*dx_phys;
    Pe0 = daleth_0*L;
    Pe1 = daleth_1 *L^2;
    Da = lambda * L^2;

    if abs (Pe1) < 1        % ---------- constant-velocity branch ----------
        gamma = sqrt (Da + Pe0^2 / 4);
        r1 = -0.5*Pe0 + gamma;
        r2 = -0.5*Pe0 - gamma;
    
        % constants from a(0)=a0,  a'(L)=0   (L = xvec(end))
        C1 =  a0 * r2 * exp(r2) / (r2*exp(r2) - r1*exp(r1));
        C2 =  -a0* r1 * exp(r1) / (r2*exp(r2) - r1*exp(r1)); 
    
        ax_fn = @(xhat) ( C1 .* exp(r1*xhat) + C2 .* exp(r2*xhat));
        return;
    endif



    gamma = -1* lambda/daleth_1;
    z_fn = @(x) abs(daleth_1)^(1/2) .* (x + (daleth_0/daleth_1));
    z = zeros(1, ngrdx);
    D_gamma_Z_array = zeros(1, ngrdx);
    D_gamma_minus_z_array = zeros(1, ngrdx);

    for i = 1:ngrdx
        xpos = i - cutoffin;
        z(i) = z_fn(xpos);
        [D_z_tmp, _] = pcd_safe(gamma, z(i));
        [D_minusz_tmp, _] = pcd_safe(gamma, -z(i));
        D_gamma_Z_array(i) = D_z_tmp;
        D_gamma_minus_z_array(i) = D_minusz_tmp;
    end

    [D_gamma_Z0 ,D_gamma_Z0Prime] = pcd_safe(gamma, z(cutoffin));
    [D_gamma_minusZ0 ,D_gamma_minusZ0Prime] = pcd_safe(gamma, -z(cutoffin));

    [D_gamma_ZL ,D_gamma_ZLPrime] = pcd_safe(gamma, z(xmax-5));
    [D_gamma_minusZL, D_gamma_minusZLPrime] = pcd_safe(gamma, -z(xmax-5));


    aleph = (daleth_0/2) + (daleth_1/2)*(xmax-cutoffin)*dx_phys;
    F1 = aleph*D_gamma_ZL + abs(daleth_1)^(0.5)* D_gamma_ZLPrime;
    F2 = aleph*D_gamma_minusZL - abs(daleth_1)^(0.5)* D_gamma_minusZLPrime;
    bet = D_gamma_Z0*F2 - D_gamma_minusZ0*F1;
    C_1 = a0 * (F2/bet);
    C_2 = -a0 * (F1/bet);

    Aeff  =  ( K_eff / ( D_phys*daleth_0*sqrt(abs(daleth_1)) ) )^(1/4);
    A = K_eff;
    

    mu_fn = @(x) exp((-A)*(1/2*daleth_0*x + 1/4*daleth_1*(x^2)));
    ax_fn = @(x) a0* mu_fn(x); %*(C_1 * D_gamma_Z_array(round(x/dx_phys)+cutoffin)+C_2*D_gamma_minus_z_array(round(x/dx_phys)+cutoffin));
end