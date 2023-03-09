function read_indicial_params(M,beta,b,a_inf,airfoil,params)

    params.A1 = 0.3;
    params.A2 = 0.7;
    params.A3 = 1.5;
    params.A4 = -0.5;
    b1_0 = 0.14;
    b2_0 = 0.53;
    params.b3 = 0.25;
    params.b4 = 0.1;
    params.b5 = 0.5;

    # Adjustments for NACA0012
    if airfoil == "NACA0012"
        params.b1 = b1_0*2.5;
        params.b2 = b2_0*0.8;
    elseif airfoil == "AMES-01"
        params.b1 = b1_0*3.0;
        params.b2 = b2_0*0.6;
    elseif airfoil == "NLR-7301"
        params.b1 = b1_0*2.5;
        params.b2 = b2_0*0.8;
    elseif airfoil == "S809"
        params.b1 = b1_0*2.5;
        params.b2 = b2_0*0.8;    
    end

    # Leishman's data constants - Principles of Helicopter
    # Aerodynamics, table 8.12, pag. 473
    # A1 = 0.482;
    # A2 = 0.518;
    # b1 = 0.684;
    # b2 = 0.235;

    params.K_a = 1/(1-M+pi*beta*M^2*(params.A1*b1_0+params.A2*b2_0)); 
    params.K_q = 1/(1-M+2*pi*beta*M^2*(params.A1*b1_0+params.A2*b2_0)); 
    params.K_aM = (params.A3*params.b4+params.A4*params.b3)/(params.b3*params.b4*(1-M)); 
    params.K_qM = 7/(15*(1-M)+3*pi*beta*M^2*params.b5); 
    params.T_I = 2*b/a_inf;
    # Leishman and Nguyen 1989 say they have reduced the K constants by 25# to
    # match experimental data (see Conclusions of their paper).
    if M > 0.07
        fac = 0.75;
    else
        fac = 1-0.25*(M/0.07)^2;
    end
    params.K_a = params.K_a*fac; 
    params.K_q = params.K_q*fac;
    if airfoil == "NACA0012"
        params.K_aM = params.K_aM*1.0; 
        params.K_qM = params.K_qM*1.0;
    elseif airfoil == "AMES-01"
        params.K_aM = params.K_aM*1.25; 
        params.K_qM = params.K_qM*1.25;
    elseif airfoil == "NLR-7301"
        params.K_aM = params.K_aM*1.0; 
        params.K_qM = params.K_qM*1.0;
    elseif airfoil == "S809"
        params.K_aM = params.K_aM*1.0; 
        params.K_qM = params.K_qM*1.0;    
    end

    return params
end