function BL_RKF45(tspan,x0,y0,A,B,U,M,b,a_0,a_1,k,beta,params,options)

    ## Handle the inputs
    if ~isempty(options)
        if options[1] == "hlim"; ind=1; hlim = options[ind+1]; end
        if options[1] == "RKFtol"; ind=1; RKFtol = options[ind+1]; end
        if options[1] == "b_max_it"; ind=1; b_max_it = options[ind+1]; end
        if options[1] == "RKF_it_max"; ind=1; RKF_it_max = options[ind+1]; end
    end        

    if ~@isdefined hlim; hlim = 5e-4; end                               # Default maximum timestep
    if ~@isdefined RKFtol; RKFtol = 1e-8; end                           # Default RFK tolerance
    if ~@isdefined b_max_it; b_max_it = 10; end                         # Default maximum number of iterations on boundary
    if ~@isdefined RKF_it_max; RKF_it_max = 10; end                     # Default maximum number of iterations for RKF approximations

    ## Setup the algorithm
    # Boundary
    delta_b = -1e-16;       # Boundary tolerance
    boundary = zeros(9,1);  # Initialize boundary 
    # Time variables
    ti = tspan[1];          # Initial time
    tf = tspan[end];        # Final time
    tc = ti;                # Current time
    # Pre-allocate
    N = length(x0);
    p_ratio = 1+sqrt(tf-ti);           # Pre-allocation ratio relative to constant time step - my simple formula
    s0 = convert(Int64, round(p_ratio*(tf-ti)/hlim));  # Initial size of the vectors
    tp = zeros(s0); xp = zeros(N,s0); yp = zeros(length(y0),s0); xdotp = zeros(N,s0); 
    # Set initial conditions
    x_i = copy(x0);
    # Initialize outputs
    tp[1] = copy(ti);  
    xp[:,1] = copy(x0); 
    yp[:,1] = copy(y0);
    xdot = copy(xdotp[:,1]);
    # Unpack parameters
    alpha_0L,alpha_ds0,alpha_ss,gamma_LS,gamma_TvL,delta_alpha_0,delta_alpha_1,kappa,nu_1,nu_2,c_d0,c_m0,c_n_alpha,d_cc,d_cm,df0_c,E0,E1,f0,fb,fS2_c_down,fS2_c_up,fS2_n,fSig_1n,fSig2_n,fSS,g_v,K0,K1,K2,K3,r0,S1,S2,Ta,Tf0,TvL,Vm,Vn1,Vn2,x_ac,z_cc,z_cm,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I = upck_params(params);
    # Initialize other variables
    tv0 = 0; so_im1 = 0; so_lim_im1 = 1; so_i = 0; so_lim_i = 1; RD_tv0 = 0; f_diff_tv0 = 0; TvL_tv0 = -1e4; theta_tv0 = 0; RD_tv0_2 = 1; f_diff_tv0_2 = 0; TvL_tv0_2 = -1e4; theta_min = 0; theta_max = 1; RD_m = 1; V2F = 0;

    ## Solve the ODEs
    i = 1;
    while (tp[i] < tf)
        hc = hlim;                # Reset current timestep to maximum allowable 
        eps = 10*RKFtol;          # Reset RKF epsilon
        b_it = 0;                 # Reset boundary iteration
        RKF_it = 0;               # Reset RKF approximations iteration
        x_ip1 = zeros(size(x_i))
        while (eps > RKFtol || any(boundary .< delta_b))
            if tc+hc > tf; hc = tf-tc; end
            k1,tv0,so_i,so_lim_i,alpha1n_i,alpha1m_i,alpha1c_i,alpha_i,q_i = BL_dxdt(tc,x_i,xdot,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,B,U,b,a_0,a_1,k,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,f0,fb,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS,fSig_1n,fSig2_n,fS2_n,fS2_c_up,fS2_c_down,gamma_TvL,theta_max,theta_min,RD_m); 
            k2,tv0 = BL_dxdt(tc+hc/4,x_i + k1*hc/4,xdot,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,B,U,b,a_0,a_1,k,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,f0,fb,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS,fSig_1n,fSig2_n,fS2_n,fS2_c_up,fS2_c_down,gamma_TvL,theta_max,theta_min,RD_m);
            k3,tv0 = BL_dxdt(tc+3/8*hc,x_i + (3/32*k1+9/32*k2)*hc,xdot,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,B,U,b,a_0,a_1,k,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,f0,fb,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS,fSig_1n,fSig2_n,fS2_n,fS2_c_up,fS2_c_down,gamma_TvL,theta_max,theta_min,RD_m);
            k4,tv0 = BL_dxdt(tc+12/13*hc,x_i + (1932*k1-7200*k2+7296*k3)/2197*hc,xdot,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,B,U,b,a_0,a_1,k,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,f0,fb,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS,fSig_1n,fSig2_n,fS2_n,fS2_c_up,fS2_c_down,gamma_TvL,theta_max,theta_min,RD_m);
            k5,tv0,so_ip1,so_lim_ip1,alpha1n_ip1,alpha1m_ip1,alpha1c_ip1,alpha_ip1,q_ip1,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F = BL_dxdt(tc+hc,x_i + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*hc,xdot,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,B,U,b,a_0,a_1,k,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,f0,fb,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS,fSig_1n,fSig2_n,fS2_n,fS2_c_up,fS2_c_down,gamma_TvL,theta_max,theta_min,RD_m);
            k6,tv0 = BL_dxdt(tc+hc/2,x_i + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*hc,xdot,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,B,U,b,a_0,a_1,k,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,f0,fb,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS,fSig_1n,fSig2_n,fS2_n,fS2_c_up,fS2_c_down,gamma_TvL,theta_max,theta_min,RD_m);
            xdot = (16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6);
            x_ip1 = x_i + xdot*hc
            eps = LinearAlgebra.norm(1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)*hc; # Error between 5th and 4th order approximations        
            boundary = BL_boundaries(boundary,tc+hc,tc,tv0,so_i,so_lim_i,alpha_i,q_i,alpha1n_i,alpha1m_i,alpha1c_i,so_ip1,so_lim_ip1,alpha_ip1,q_ip1,alpha1n_ip1,alpha1m_ip1,alpha1c_ip1,alpha_ss,TvL,r0);
            # Check minima and maxima of theta
            theta_min,theta_max,RD_m = theta_extremes(q_i,so_im1,so_i,so_ip1,so_lim_im1,so_lim_i,so_lim_ip1,x_i[13],theta_min,theta_max,RD_m);
            # Check iterations on boundary
            if any(boundary .< delta_b) 
                if boundary[9] < delta_b && abs(so_ip1/so_lim_ip1) > 1 # Check for maximum theta at begin of downstroke as well - Stalled conditions are garanteed even though the actual maximum theta has not been reached yet
                    theta_max = so_ip1/so_lim_ip1;
                end
                if boundary[2] < delta_b && abs(so_ip1/so_lim_ip1) > 0 # Zero theta_max when theta = 0 and theta is growing
                    theta_max = 0;
                end
                b_it = b_it+1;
                if b_it == b_max_it      
                    break
                end
                hc = hc/2;
            # Check RKF iterations    
            elseif eps > RKFtol
                RKF_it = RKF_it+1;
                if RKF_it == RKF_it_max
                    break
                end
                hc = hc/2;
            end
        end
        # Update outputs
        tp[i+1] = tc+hc;      
        xp[:,i+1] = copy(x_ip1)
        xdotp[:,i] = copy(xdot)
        yp[:,i] = BL_outputs(tc,x_i,tv0,U,b,beta,k,a_0,a_1,M,x_ac,K0,K1,K2,K3,kappa,K_a,K_aM,K_q,K_qM,T_I,b1,b2,b3,b4,b5,A1,A2,A3,A4,alpha_0L,E0,E1,c_m0,c_n_alpha,r0,alpha_ds0,alpha_ss,S1,S2,Tf0,f0,fb,Vn1,Vn2,Vm,g_v,c_d0,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,gamma_LS,df0_c,fSig_1n,fSig2_n,fS2_n,fS2_c_up,fS2_c_down,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,theta_max,theta_min,RD_m,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL);
        # Setup next time step
        tc = tc+hc
        x_i = copy(x_ip1)
        so_im1 = so_i
        so_lim_im1 = so_lim_i
        i = i+1
        if rem(i,1e4) == 0
    #         disp(["RKF45 progress: ",num2str((tc-ti)/(tf-ti)*100,"#10.2f") "#"])
        end
    end

    # Assume last time step derivatives and outputs equal to previous
    xdotp[:,i] = xdotp[:,i-1];
    yp[:,i] = yp[:,i-1];

    # Truncate pre-allocated vectors
    tp = tp[1:i]; xp = xp[:,1:i]; yp = yp[:,1:i]; xdotp = xdotp[:,1:i];
    return tp,xp,yp,xdotp
end