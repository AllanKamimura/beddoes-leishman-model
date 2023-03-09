function interp_coefs(data,alpha_0,delta_alpha,N)

    time_cl_exp = vec(data[7]);
    cl_exp = vec(data[8]);
    time_cm_exp = vec(data[9]);
    cm_exp = vec(data[10]);
    time_cd_exp = vec(data[11]);
    cd_exp = vec(data[12]);

    time = range(-90, stop = 270, length = N);
    alpha = alpha_0 .+ delta_alpha*sind.(time);

    cl_interp = Dierckx.Spline1D(time_cl_exp,cl_exp; k = 3)(time);
    cm_interp = Dierckx.Spline1D(time_cm_exp,cm_exp; k = 3)(time);
    cd_interp = Dierckx.Spline1D(time_cd_exp,cd_exp; k = 3)(time);
    cn_interp = zeros(1,N);
    cc_interp = zeros(1,N);
    for i=1:N
        # c_n and c_c
        A = [cos(alpha[i]) sin(alpha[i]);
            sin(alpha[i]) -cos(alpha[i])];
        b = [cl_interp[i]; cd_interp[i]];
        x = A*b;
        cn_interp[i] = x[1];
        cc_interp[i] = x[2];
    end

    return time, alpha, cn_interp, cm_interp, cc_interp
end