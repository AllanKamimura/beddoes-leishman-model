function cc_coeff(alpha,f2prime_c,theta,c_n_alpha,alpha_E,E0,E1,c_d0,alpha_0L,RD)

    if abs(theta) < 1
        c_c = -c_d0*cos(alpha)+c_n_alpha*(alpha_E-alpha_0L)*sin((alpha_E-0))^1*(f2prime_c^(abs(theta))-E0*RD^2); # RD^2 is to ignore the effect of -E0 at low pitch rates 
    else
        zeta = minimum([1/RD, 20]); # Factor defining how sharp is the drop in c_c, limited to 20
        c_c = -c_d0*cos(alpha)+c_n_alpha*(alpha_E-alpha_0L)*sin((alpha_E-0))^1*(f2prime_c^(abs(theta))-E0*(RD^2+E1*(1-exp(-zeta*(abs(theta)-1))))); # E1 is a constant defining the limit rate of c_c drop
    end

    return c_c
end