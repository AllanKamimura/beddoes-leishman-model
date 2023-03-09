function load_GUD(GUD)

    authors = ["GU"];
    b = 0.55/2;     # Semi-chord
    a_inf = 340;    # Assumed MSL conditions
    if GUD >= 11E6 && GUD <=12E6
    airfoil = "NACA0012";
    end

    name = Printf.@sprintf("GUD_%d.mat", 11012812)
    vars = MAT.matread(name)

    M = vars["M"]
    time = vars["time"]
    alpha = vars["alpha"]
    k = vars["k"]
    delta_alpha = vars["delta_alpha"]
    cn = vars["cn"]
    alpha_0 = vars["alpha_0"]
    cc = vars["cc"]
    cm = vars["cm"]

    alpha_exp_cl = alpha;
    alpha_exp_cn = alpha;
    alpha_exp_cm = alpha;
    alpha_exp_cc = alpha;
    alpha_exp_cd = alpha;
    cl_exp = cn.*cosd.(alpha)+cc.*sind.(alpha);
    cn_exp = cn;
    cm_exp = cm;
    cc_exp = cc;
    cd_exp = cn.*sind.(alpha)-cc.*cosd.(alpha);
    time_exp_cl = time;
    time_exp_cn = time;
    time_exp_cm = time;
    time_exp_cc = time;
    time_exp_cd = time;
    clt_exp = cl_exp;
    cnt_exp = cn_exp;
    cmt_exp = cm_exp;
    cct_exp = cc_exp;
    cdt_exp = cd_exp;

    return b,a_inf,airfoil,M,alpha_0,delta_alpha,k,alpha_exp_cl,cl_exp,alpha_exp_cm,cm_exp,alpha_exp_cd,cd_exp,alpha_exp_cn,cn_exp,alpha_exp_cc,cc_exp,time_exp_cl,clt_exp,time_exp_cm,cmt_exp,time_exp_cd,cdt_exp,authors,time_exp_cn,cnt_exp,time_exp_cc,cct_exp
end