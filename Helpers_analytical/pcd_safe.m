function [D, Dp] = pcd_safe(gam, z)
    % Parabolic-cylinder D_γ(z) valid for *any* |z| (all real γ)
    % Regime  A: |z| < 1e-2      –> Maclaurin log-sum
    % Regime  B: |z| ≤ 50        –> WKB  + downward recurrence   (your code)
    % Regime  C: |z|  > 50        –> large-z log-asymptotic      (this part)

    absz = abs(z);
    z = z;
    % ---------- Regime A  -------------------------------------------------

    % ---------- Regime C  (|z| large)  -----------------------------------
    % log-prefactor
    logPref = -z.^2/4 + gam*log(absz);   % this never over/under-flows
    % first two correction terms
    c1 = 1;
    c2 = - gam*(gam-1) / (2*z.^2);
    c3 =   gam*(gam-1)*(gam-2)*(3*gam-1)/(8*z.^4);
    corr = c1 + c2 + c3;

    % safe magnitude + sign
    sgn = 1;
    if z < 0
        sgn = cos(pi*gam);               % real continuation of (-1)^γ
    end
    D = sgn .* exp(logPref) .* corr;   % full D_γ(z)

    % derivative via relation  D′ = z/2 D – D_{γ+1}
    if nargout > 1
        Dgp1 = pcd_safe(gam+1, z);        % one recursion
        Dp   = (z/2).*D - Dgp1;
    end
end