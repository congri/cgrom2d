function [adjoints] = get_adjoints(K, theta_cf, Y, domain, Tf_i)
%Compute adjoints for gradient computation

%Y is the full temperature vector here, including essential nodes
d_log_p_cf = theta_cf.W'*(theta_cf.S\(Tf_i - theta_cf.mu - theta_cf.W*Y));
%Gradient with respect to natural nodes; take out derivatives w.r.t. essential nodes
d_log_p_cf(~isnan(domain.nodalCoordinates(4, :))) = [];

adjoints = K\d_log_p_cf;
end

