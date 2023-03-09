function read_data(input)
    ## Handle inputs
    if length(input) == 1
        if input > 1e6
            GUD = input;
        elseif input > 1000
            frame = input;
        else
            run = input;
        end
    else
        test = input;
    end

    ## Load data from frame, run or desired test condition
    if @isdefined frame
        b,a_inf,airfoil,M,alpha_0,delta_alpha,k,alpha_exp_cl,cl_exp,alpha_exp_cm,cm_exp,alpha_exp_cd,cd_exp,alpha_exp_cn,cn_exp,alpha_exp_cc,cc_exp,time_exp_cl,clt_exp,time_exp_cm,cmt_exp,time_exp_cd,cdt_exp,authors = load_frame(frame);
        time_exp_cn = NaN; cnt_exp = NaN; time_exp_cc = NaN; cct_exp = NaN;
        alpha_mod_cl = NaN; cl_mod = NaN; alpha_mod_cm = NaN; cm_mod = NaN; alpha_mod_cd = NaN; cd_mod = NaN; alpha_mod_cn = NaN; cn_mod = NaN; alpha_mod_cc = NaN; cc_mod = NaN;
    elseif @isdefined GUD
        b,a_inf,airfoil,M,alpha_0,delta_alpha,k,alpha_exp_cl,cl_exp,alpha_exp_cm,cm_exp,alpha_exp_cd,cd_exp,alpha_exp_cn,cn_exp,alpha_exp_cc,cc_exp,time_exp_cl,clt_exp,time_exp_cm,cmt_exp,time_exp_cd,cdt_exp,authors,time_exp_cn,cnt_exp,time_exp_cc,cct_exp = load_GUD(GUD);
        alpha_mod_cl = NaN; cl_mod = NaN; alpha_mod_cm = NaN; cm_mod = NaN; alpha_mod_cd = NaN; cd_mod = NaN; alpha_mod_cn = NaN; cn_mod = NaN; alpha_mod_cc = NaN; cc_mod = NaN;
    elseif @isdefined run
        b,a_inf,airfoil,M,alpha_0,delta_alpha,k,alpha_exp_cl,cl_exp,alpha_exp_cn,cn_exp,alpha_exp_cm,cm_exp,alpha_exp_cd,cd_exp,alpha_exp_cc,cc_exp,alpha_mod_cl,cl_mod,alpha_mod_cn,cn_mod,alpha_mod_cm,cm_mod,alpha_mod_cd,cd_mod,alpha_mod_cc,cc_mod,authors = load_run(run);
        time_exp_cl = NaN; clt_exp = NaN; time_exp_cm = NaN; cmt_exp = NaN; time_exp_cd = NaN; cdt_exp = NaN; time_exp_cn = NaN; cnt_exp = NaN; time_exp_cc = NaN; cct_exp = NaN;
    end

    ## Gather data
    data = [alpha_exp_cl, cl_exp, alpha_exp_cm, cm_exp, alpha_exp_cd, cd_exp, time_exp_cl, clt_exp, time_exp_cm, cmt_exp, time_exp_cd, cdt_exp, alpha_exp_cn, cn_exp, alpha_exp_cc, cc_exp, alpha_mod_cl, cl_mod, alpha_mod_cm, cm_mod, alpha_mod_cd, cd_mod, alpha_mod_cn, cn_mod, alpha_mod_cc, cc_mod, time_exp_cn, cnt_exp, time_exp_cc, cct_exp];

    ## Flow properties
    U = M*a_inf;
    beta = sqrt(1-M^2);

    ## Airfoil parameters 

    return authors,data,U,beta,M,k,alpha_0,delta_alpha,b,a_inf,airfoil
end