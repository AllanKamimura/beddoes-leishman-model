function RMS_error_calculator(c_n,c_m,c_c,cn_interp,cm_interp,cc_interp,N)

    # Sum of squares of the errors
    eps_cn = sum((c_n - vec(cn_interp)).^2);
    eps_cm = sum((c_m - vec(cm_interp)).^2);
    eps_cc = sum((c_c - vec(cc_interp)).^2);

    # Root mean square error normalized by the mean
    cn_NRMSE = sqrt(eps_cn/N)/abs(maximum(cn_interp)-minimum(cn_interp));
    cm_NRMSE = sqrt(eps_cm/N)/abs(maximum(cm_interp)-minimum(cm_interp));
    cc_NRMSE = sqrt(eps_cc/N)/abs(maximum(cc_interp)-minimum(cc_interp));
    return cn_NRMSE,cm_NRMSE,cc_NRMSE
end