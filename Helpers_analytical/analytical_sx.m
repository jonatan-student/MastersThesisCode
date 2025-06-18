function [sx_ana_fn, S_0, S_1] = analytical_sx(slope, H_0, alpha, beta, mu_phys, Ls_phys, dx_phys)
    H_0 = H_0;
    if slope ~= 0
        H_0 = H_0-dx_phys;
    end
    G = 1/(12*mu_phys);
    alpha = alpha;
    beta = beta;
    m = slope;
    %disp(m) here to debug the slope

    S_0 =G*alpha*(H_0^2 + 6*Ls_phys*H_0);
    S_1 =(1/2)*G*(-alpha*(H_0*m - 6*Ls_phys*m)+beta*(H_0^2 + 6*Ls_phys*H_0));
    sx_ana_fn = @(x) (S_0 + (S_1)*x);
end