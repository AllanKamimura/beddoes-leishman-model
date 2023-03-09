function interp_data(c_n,c_m,c_c,t,t_cycle,n_discard,data,authors,a_0,a_1,INPUT)

    # Find time indices
    iti = findnext(t .> n_discard * t_cycle, 1);   # First index after discarded time
    itt = findnext(t .> t[end] - 5/4 * t_cycle, 1);  # Index of begin of cycle
    ittf = findnext(t .> t[end] - 1/4 * t_cycle, 1); # Index of end of cycle
    # Spline interpolation for c_n and c_c from c_l and c_d
    if length(data) > 1 && (cmp(authors[1],"McAlister et al. (1982)") == 1 || cmp(authors[1],"GU") == 1)
        if authors[1] == "McAlister et al. (1982)"
            time_interp,alpha_interp,cn_interp,cm_interp,cc_interp = interp_coefs(data,a_0,a_1,ittf-itt+1);
            # NRMSE
            cn_NRMSE,cm_NRMSE,cc_NRMSE = RMS_error_calculator(c_n[itt:ittf],c_m[itt:ittf],c_c[itt:ittf],cn_interp,cm_interp,cc_interp,ittf-itt+1);
        else
            time_interp,alpha_interp,cn_interp,cm_interp,cc_interp,itt,ittf = interp_coefs_GUD(data,a_0,a_1,t,t_cycle);
            # NRMSE
            cn_NRMSE,cm_NRMSE,cc_NRMSE = RMS_error_calculator(c_n[itt:ittf],c_m[itt:ittf],c_c[itt:ittf],cn_interp,cm_interp,cc_interp,ittf-itt+1);
        end
    else
        time_interp = NaN;
        alpha_interp = NaN;
        cn_interp = NaN;
        cm_interp = NaN;
        cc_interp = NaN;
        cn_NRMSE = NaN;
        cm_NRMSE = NaN;
        cc_NRMSE = NaN;
    end

    return iti,itt,ittf,time_interp,alpha_interp,cn_interp,cm_interp,cc_interp,cn_NRMSE,cm_NRMSE,cc_NRMSE
end