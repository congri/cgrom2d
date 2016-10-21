function [log_q, d_log_q, Tc] = log_q_i(Xi, Tf_i_minus_mu, theta_cf, theta_c, Phi,  domainc)

conductivity = exp(Xi);
if any(conductivity < 1e-6)
    %lower bound on conductivity for stability
    conductivity(conductivity < 1e-6) = 1e-6;
end
if any(conductivity > 1e6)
    %upper bound on conductivity for stability
    conductivity(conductivity > 1e6) = 1e6;
end

[lg_p_c, d_lg_p_c] = log_p_c(Xi, Phi, theta_c.theta, theta_c.sigma);
[lg_p_cf, d_lg_p_cf, Tc] = log_p_cf(Tf_i_minus_mu, domainc, conductivity, theta_cf);

if(~isscalar(lg_p_c))
    lg_p_c
end
assert(isscalar(lg_p_c), 'error: lg_p_c not scalar')
assert(isscalar(lg_p_cf), 'error: lg_p_cf not scalar')

log_q = lg_p_cf + lg_p_c;

d_log_q = d_lg_p_c + d_lg_p_cf;
if any(imag(d_log_q))
    Xi
    theta_c.theta
    theta_c.sigma
    Phi
    d_log_q
    d_lg_p_c
    d_lg_p_cf
    error('Imaginary gradient')
end


%Finite difference gradient check
FDcheck = false;
if FDcheck
    disp('Gradient check log q_i')
    d = 1e-5;
    gradFD = zeros(domainc.nEl, 1);
    for i = 1:domainc.nEl
        dXi = zeros(domainc.nEl, 1);
        dXi(i) = d;
        conductivityFD = conductivity + conductivity.*dXi;
        
        [lg_p_c, ~] = log_p_c(Xi + dXi, Phi, theta_c.theta, theta_c.sigma);
        [lg_p_cf, ~] = log_p_cf(Tf_i_minus_mu, domainc, conductivityFD, theta_cf);
        
        log_qFD = lg_p_cf + lg_p_c;
        gradFD(i) = (log_qFD - log_q)/d;
    end
    
    relgrad = gradFD./d_log_q
    if(any(abs(relgrad - 1) > .1))
        %Note: if d_log_q << d_log_p_c, d_log_p_cf, then this might be due to numerical issues, i.e.
        %FD gradient is unprecise
        %for small log q, it is possible that the FD gradient is unprecise
        conductivity
        conductivityFD
        Xi
        XiFD = Xi + dXi
        d_log_q
        d_lg_p_c
        d_lg_p_cf
        pause 
    end
    
end

end

