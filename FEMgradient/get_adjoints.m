function [adjoints] = get_adjoints(K, theta_cf, Tc, domain, Tf_i)
%Compute adjoints for gradient computation

%Tc is the full temperature vector here, including essential nodes
d_log_p_cf = theta_cf.W'*(theta_cf.Sinv*(Tf_i - theta_cf.mu - theta_cf.W*Tc));
%Gradient with respect to natural nodes; take out derivatives w.r.t. essential nodes
d_log_p_cf(~isnan(domain.essentialTemperatures)) = [];

adjoints = K\d_log_p_cf;

end

