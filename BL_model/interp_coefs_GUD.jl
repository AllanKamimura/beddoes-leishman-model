function interp_coefs_GUD(data,alpha_0,delta_alpha,t,t_cycle)

    alpha_exp = vec(data[1]);
    time_exp = vec(data[9]);
    cm_exp = vec(data[10]);
    cn_exp = vec(data[28]);
    cc_exp = vec(data[30]);

    alpha_init = alpha_exp[1]*pi/180;
    t_init = asin((alpha_init-alpha_0)/delta_alpha)*180/pi; # cycle angle [deg] at which the data starts
    delta_t = -t_init/360;     

    # according time shift (fraction of t_cycle)
    itt2 = findnext(t .> t[end] - (1+delta_t)*t_cycle, 1);
    ittf2 = findnext(t .> t[end] - (delta_t)*t_cycle,1);
    N = ittf2-itt2+1;

    time = range(-90,stop = 270, length = N);
    alpha = alpha_0 .+ delta_alpha*sind.(time);
    cn_interp = Dierckx.Spline1D(time_exp,cn_exp)(time);
    cm_interp = Dierckx.Spline1D(time_exp,cm_exp)(time);
    cc_interp = Dierckx.Spline1D(time_exp,cc_exp)(time);

    return time,alpha,cn_interp,cm_interp,cc_interp,itt2,ittf2
end