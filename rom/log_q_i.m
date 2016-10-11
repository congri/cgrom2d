function [log_q, Tc] = log_q_i(Xi, Tf_i, theta_cf, theta_c, Phi, dCoarse, physC, W)


[lg_p_c, d_lg_p_c] = log_p_c(Xi, Phi, theta_c.theta, theta_c.sigma);
[lg_p_cf, Tc] = log_p_cf(Tf_i, dCoarse, physC, W, theta_cf.S, Xi);

log_q = lg_p_cf + lg_p_c;

% d_log_q = d_lg_p_c + d_lg_p_cf;


%Finite difference gradient check
% FDcheck = false;
% if FDcheck
%     d = 1e-4;
%     gradFD = zeros(Cmesh.N_el, 1);
%     for i = 1:Cmesh.N_el
%         dXi = zeros(Cmesh.N_el, 1);
%         dXi(i) = d;
%         Cmesh.conductivity = Xi + dXi;
%         
%         [lg_p_c, d_lg_p_c] = log_p_c(Xi + dXi, Phi, theta_c.theta, theta_c.sigma);
%         [lg_p_cf, d_lg_p_cf] = log_p_cf(Tf_i, Cmesh, heatSource, boundary, W, theta_cf.S);
%         
%         log_qFD = lg_p_cf + lg_p_c;
%         gradFD(i) = (log_qFD - log_q)/d;
%     end
%     
%     relgrad = gradFD./d_log_q
%     
% end

end

