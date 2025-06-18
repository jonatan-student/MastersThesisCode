function stats = validate_profiles(num, ana)
    assert(numel(num)==numel(ana), ...
            'NUM and ANA must have the same length.');

    % ----- turn into column vectors & drop NaNs/Infs ---------------------
    y  = num(:);
    yhat  = ana(:);
    good = isfinite(y) & isfinite(yhat);
    y  = y(good);
    yhat  = yhat(good);

    N = numel(y);

    % ----- error vector --------------------------------------------------
    stats.err = y - yhat;
    safe_yhat = yhat;
    safe_yhat(safe_yhat == 0) = eps;
    stats.relError = stats.err./yhat;


    % ----- metrics -------------------------------------------------------
    stats.N      = N;
    stats.MSE    = mean(stats.err.^2);
    stats.RMSE   = sqrt(stats.MSE);
    stats.MAE    = mean(abs(stats.err));
    stats.MaxAbs = max(abs(stats.err));
    stats.MaxAbsRel = max(abs(stats.relError));
    stats.MAPE = 100*mean(abs(stats.relError));

    % coefficient of determination
    y_mean = mean(y);

    ss_tot = sqrt(sum( (y - y_mean).^2 ));
    ss_res = sum( stats.err.^2 );
    stats.R2 = (1 - ss_res/ss_tot)^2;
    % disp(stats.R2)
end