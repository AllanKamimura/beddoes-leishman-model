function get_SS_matrices(U,b,beta,M,a_inf,airfoil,params)

    ## Indicial parameters
    params = read_indicial_params(M,beta,b,a_inf,airfoil,params);

    ## SS matrices
    # BL potential flow states
    a = [-U/b*beta^2*params.b1; -U/b*beta^2*params.b2; -1/(params.K_a*params.T_I); -1/(params.K_q*params.T_I); -1/(params.b3*params.K_aM*params.T_I); -1/(params.b4*params.K_aM*params.T_I); -params.b5*U/b*beta^2; -1/(params.K_qM*params.T_I)];
    A = LinearAlgebra.diagm(a);
    B = [1 1/2; 1 1/2; 1 0; 0 1; 1 0; 1 0; 0 1; 0 1];

    return A,B
end