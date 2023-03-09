function BL_model(INPUT,new_df,disp_info)
    # Modified BL model

    ## Choose experimental data to load from McAlister's frame, GU's data, some other experimental run or a desired test condition
    test = [0.12; 0.01; 10*pi/180; 20*pi/180; "S809"]; # M; k; a_0; a_1
    authors,data,U,beta,M,k,a_0,a_1,b,a_inf,airfoil = read_data(INPUT);
    params = airfoil_parameters(airfoil,M,U,b,new_df);
    
    ## Initialization 
    # Time span
    n_discard = 2;              # Number of cycles to discard for Plots.plots
    n_cycles = n_discard+2;     # Number of cycles to run
    t_cycle = 2*pi*b/(U*k);     # Time of one cycle
    tf = t_cycle*n_cycles;      # Test time
    tspan = [0 tf];
    # Initial conditions: {x0} = {0}
    N = 14;                     # Number of states
    x0 = zeros(N,1);
    y0 = zeros(35,1)
    y0[1:4] = [a_0; 0; 0; 0];
    y0[5:35] .= NaN;

    ## Indicial parameters
    params = read_indicial_params(M,beta,b,a_inf,airfoil,params);

    ## State space matrices
    A,B = get_SS_matrices(U,b,beta,M,a_inf,airfoil,params);

    ## ODE solver
    dt = 1e-4; options = ["hlim",dt]; # Maximum time step
    t,x,y,xdot = BL_RKF45(tspan,x0,y0,A,B,U,M,b,a_0,a_1,k,beta,params,options);
    
    ## Output variables
    alpha,c_n,c_m,c_c,c_l,c_d,c_nf,c_nI,c_mf,c_mI,c_nv,
    c_mv,c_cv,c_lv,f,f2prime,tau_v,dCP,so_lim,qR,alpha1,
    K_f,c_nC,Tf_n,dalpha1,fprime,fprime_cm,fprime_cc,q,dalpha1_cm,
    dalpha1_cc,R_dot,Tf_m,Tf_c,theta,theta_min,theta_max,P,S,R,
    RD_theta,alpha_E = BL_output_vars(x,y);

    ## Interpolate data and find normalized error coefficients
    iti,itt,ittf,time_interp,alpha_interp,cn_interp,cm_interp,cc_interp,cn_NRMSE,cm_NRMSE,cc_NRMSE = interp_data(c_n,c_m,c_c,t,t_cycle,n_discard,data,authors,a_0,a_1,INPUT);
    if disp_info
        println(string(airfoil, ": ", INPUT," (", authors[1], ")", ", NMRS errors: cn = ", cn_NRMSE, ", cm = ", cm_NRMSE, ", cc = ", cc_NRMSE));
        cl = Plots.plot(alpha[iti:end].*180/pi, c_l[iti:end], label = "modelo", ylabel = "\n\nCL", xlabel = "α [deg]")
        cl = Plots.plot!(data[1], data[2], label = "experimento", xlabel = "α [deg]")
        cl = Plots.plot!(legend = placelegend())
        cm = Plots.plot(alpha[iti:end].*180/pi, c_m[iti:end], label = "modelo", ylabel = "\n\n\nCM", xlabel = "α [deg]")
        cm = Plots.plot!(data[3], data[4], label = "experimento", xlabel = "α [deg]")
        cm = Plots.plot!(legend = placelegend()) 
        cd = Plots.plot(alpha[iti:end].*180/pi, c_d[iti:end], label = "modelo", ylabel = "\n\nCD", xlabel = "α [deg]")
        cd = Plots.plot!(data[5], data[6], label = "experimento", xlabel = "α [deg]")
        cd = Plots.plot!(legend = placelegend()) 
        cn = Plots.plot(alpha[iti:end].*180/pi, c_n[iti:end], label = "modelo", ylabel = "\n\n\nCN", xlabel = "α [deg]")
        cn = Plots.plot!(alpha_interp.*180/pi, cn_interp', label = "experimento", xlabel = "α [deg]")
        cn = Plots.plot!(legend = placelegend()) 
        cc = Plots.plot(alpha[iti:end].*180/pi, c_c[iti:end], label = "modelo", ylabel = "\n\nCC", xlabel = "α [deg]")
        cc = Plots.plot!(alpha_interp.*180/pi, cc_interp', label = "experimento", xlabel = "α [deg]")
        cc = Plots.plot!(legend = placelegend()) 
        alpha_plot = Plots.plot(cl, cm, cd, cn, cc, layout = (5, 1))

        cl = Plots.plot((t[itt:ittf] .-t[itt])./ t_cycle*360 .-90, c_l[itt:ittf], label = "modelo", xlabel = "ωt [deg]", ylabel = " ")
        cl = Plots.plot!(data[7], data[8], label = "experimento", xlabel = "ωt [deg]", ylabel = " ")
        cl = Plots.plot!(legend = placelegend())   
        cm = Plots.plot((t[itt:ittf] .-t[itt])./ t_cycle*360 .-90, c_m[itt:ittf], label = "modelo", xlabel = "ωt [deg]", ylabel = " ")
        cm = Plots.plot!(data[9], data[10], label = "experimento", xlabel = "ωt [deg]", ylabel = " ")
        cm = Plots.plot!(legend = placelegend()) 
        cd = Plots.plot((t[itt:ittf] .-t[itt])./ t_cycle*360 .-90, c_d[itt:ittf], label = "modelo", xlabel = "ωt [deg]", ylabel = " ")
        cd = Plots.plot!(data[11], data[12], label = "experimento", xlabel = "ωt [deg]", ylabel = " ")
        cd = Plots.plot!(legend = placelegend()) 
        cn = Plots.plot((t[itt:ittf] .-t[itt])./ t_cycle*360 .-90, c_n[itt:ittf], label = "modelo", xlabel = "ωt [deg]", ylabel = " ")
        cn = Plots.plot!(time_interp, cn_interp', label = "experimento", xlabel = "ωt [deg]", ylabel = " ")
        cn = Plots.plot!(legend = placelegend()) 
        cc = Plots.plot((t[itt:ittf] .-t[itt])./ t_cycle*360 .-90, c_c[itt:ittf], label = "modelo", xlabel = "ωt [deg]", ylabel = " ")
        cc = Plots.plot!(time_interp, cc_interp', label = "experimento", xlabel = "ωt [deg]", ylabel = " ")
        cc = Plots.plot!(legend = placelegend()) 
        
        time_plot = Plots.plot(cl, cm, cd, cn, cc, layout = (5, 1))

        this_plot = Plots.plot(alpha_plot, time_plot, size = (1000, 2000))

        Plots.savefig(this_plot, string(INPUT) * ".png")
        this_plot
    else
        return [cn_NRMSE, cm_NRMSE, cc_NRMSE]
    end
    # ## Plots
    # axes_size,lw,ms,xlim_vec,figure1,tabgp = init_Plots.plotter(authors,INPUT,a_0,a_1);
    # 
    # # Coefficients x alpha
    # alpha_plotter(alpha,c_l,alpha_interp,NaN,data{1,1},data{2,1},data{17,1},data{18,1},'Lift coefficient','c_l',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,xlim_vec,iti)
    # alpha_plotter(alpha,c_m,alpha_interp,NaN,data{3,1},data{4,1},data{19,1},data{20,1},'Moment coefficient','c_m',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,xlim_vec,iti)
    # alpha_plotter(alpha,c_d,alpha_interp,NaN,data{5,1},data{6,1},data{21,1},data{22,1},'Drag coefficient','c_d',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,xlim_vec,iti)
    # alpha_plotter(alpha,c_n,alpha_interp,cn_interp,data{13,1},data{14,1},data{23,1},data{24,1},'Normal coefficient','c_n',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,xlim_vec,iti)
    # alpha_plotter(alpha,c_c,alpha_interp,cc_interp,data{15,1},data{16,1},data{25,1},data{26,1},'Chordwise coefficient','c_c',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,xlim_vec,iti)
    # 
    # # Coefficients x time
    # time_plotter(t,c_l,time_interp,NaN,data{7,1},data{8,1},'Lift coefficient','c_l x t',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,t_cycle,itt,ittf)
    # time_plotter(t,c_m,time_interp,NaN,data{9,1},data{10,1},'Moment coefficient','c_m x t',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,t_cycle,itt,ittf)
    # time_plotter(t,c_d,time_interp,NaN,data{11,1},data{12,1},'Drag coefficient','c_d x t',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,t_cycle,itt,ittf)
    # time_plotter(t,c_n,time_interp,cn_interp,data{27,1},data{28,1},'Normal coefficient','c_n x t',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,t_cycle,itt,ittf)
    # time_plotter(t,c_c,time_interp,cc_interp,data{29,1},data{30,1},'Chordwise coefficient','c_c x t',tabgp,authors,frame,GUD,M,k,a_0,a_1,axes_size,lw,ms,t_cycle,itt,ittf)
    # 
    # # Separation points, angle offsets and time delay constants
    # if ~isnan(alpha_interp)
    #     f_da_Tf_Plots.plotter(alpha_interp,cn_interp,cc_interp,cm_interp,x,y,tabgp,INPUT,frame,GUD,axes_size,lw,xlim_vec,iti,itt,ittf)
    # end

    # Other variables
    # figure;plot(alpha(iti:end)*180/pi,theta(iti:end));grid  

end