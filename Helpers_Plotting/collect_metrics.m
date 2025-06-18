function M = collect_metrics(ux_mid, ux_ana, S_num,  S_ana,  dpdx_num, dpdx_ana,a_num,  a_ana)
    f = @(n,a) validate_profiles(n(:), a(:));   % anonymous shortcut
    M = struct();
    m_ux       = f(ux_mid,  ux_ana);
    m_S        = f(S_num,   S_ana);
    m_dp       = f(dpdx_num,dpdx_ana);
    m_a        = f(a_num,   a_ana);

    M.R2_ux    = m_ux.R2;     M.RMSE_ux = m_ux.RMSE;    M.MAPE_ux = m_ux.MAPE;      M.relErr_ux = m_ux.relError;
    M.R2_S     = m_S.R2;      M.RMSE_S  = m_S.RMSE;     M.MAPE_S  = m_S.MAPE;       M.relErr_S  = m_S.relError;
    M.R2_dpdx  = m_dp.R2;     M.RMSE_dp = m_dp.RMSE;    M.MAPE_dp = m_dp.MAPE;      M.relErr_dp = m_dp.relError;
    M.R2_a     = m_a.R2;      M.RMSE_a  = m_a.RMSE;     M.MAPE_a  = m_a.MAPE;       M.relErr_a  = m_a.relError;

    M.MAR_ux = m_ux.MaxAbsRel;
    M.MAR_S  = m_S.MaxAbsRel;
    M.MAR_dp = m_dp.MaxAbsRel;
    M.MAR_a  = m_a.MaxAbsRel;
end