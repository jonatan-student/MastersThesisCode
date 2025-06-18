function [dpdx_ana_fn, alpha, beta]= analytical_dpdx(inlet, outlet, len, slope)
    alpha = inlet;
    beta = (outlet-inlet)/(len);
    if abs(slope) <= eps
        beta = 0;
    end
    %fprintf("Alpha: %d, Beta: %d\n", alpha, beta)
    dpdx_ana_fn = @(x) alpha + beta*x;
end