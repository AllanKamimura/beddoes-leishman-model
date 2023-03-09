function load_frame(frame::Int64)
    b = 0.61/2;     # McAlister (1982) - page 14
    a_inf = 340;    # Assumed MSL conditions
    if frame >= 7019 && frame <= 14220
        airfoil = "NACA0012";
    elseif frame >=24022 && frame <= 31310
        airfoil = "AMES-01";
    elseif frame >=67000
        airfoil = "NLR-7301";
    end
    name = Printf.@sprintf "./example/frame_%d.mat" frame    
    vars = MAT.matread(name)

    M = vars["M"]
    alpha_exp_cl = vars["alpha_exp_cl"]
    alpha_exp_cd = vars["alpha_exp_cd"]
    time_exp_cd = vars["time_exp_cd"]
    cmt_exp = vars["cmt_exp"]
    k = vars["k"]
    time_exp_cl = vars["time_exp_cl"]
    clt_exp = vars["clt_exp"]
    delta_alpha = vars["delta_alpha"]
    alpha_exp_cm = vars["alpha_exp_cm"]
    cd_exp = vars["cd_exp"]
    cl_exp = vars["cl_exp"]
    cdt_exp = vars["cdt_exp"]
    time_exp_cm = vars["time_exp_cm"]
    cm_exp = vars["cm_exp"]
    alpha_0 = vars["alpha_0"]
    authors = vars["authors"]

    alpha_exp_cn = NaN;
    cn_exp = NaN;
    alpha_exp_cc = NaN;
    cc_exp = NaN;

    if ~@isdefined time_exp_cl
        time_exp_cl = NaN;
        clt_exp = NaN;
        time_exp_cm = NaN;
        cmt_exp = NaN;
        time_exp_cd = NaN;
        cdt_exp = NaN;
    end

    return b,a_inf,airfoil,M,alpha_0,delta_alpha,k,alpha_exp_cl,cl_exp,alpha_exp_cm,cm_exp,alpha_exp_cd,cd_exp,alpha_exp_cn,cn_exp,alpha_exp_cc,cc_exp,time_exp_cl,clt_exp,time_exp_cm,cmt_exp,time_exp_cd,cdt_exp,authors
end