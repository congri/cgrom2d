function [log_p, d_log_p, Tc] = log_p_cf(Tf_i, domainc, physicalc, conductivity, theta_cf)
%Coarse-to-fine map
%ignore constant prefactor
%log_p = -.5*logdet(S, 'chol') - .5*(Tf - mu)'*(S\(Tf - mu));
%diagonal S

W = theta_cf.W;
S = theta_cf.S;
D = zeros(2, 2, domainc.nEl);
control.plt = false;
%Conductivity matrix D, only consider isotropic materials here
for j = 1:domainc.nEl
    D(:, :, j) =  conductivity(j)*eye(2);
end
FEMout = heat2d(domainc, physicalc, control, D);

WTc = W*FEMout.Tff(:);
%only for diagonal S!
assert(isdiag(S), 'Error: matrix S not diagonal');
Sinv = diag(1./diag(S));
log_p = -.5*sum(log(diag(S))) - .5*(Tf_i - WTc)'*(Sinv*(Tf_i - WTc));

if nargout > 1
    %Gradient of FEM equation system w.r.t. conductivities
    d_r = FEMgrad(FEMout, domainc, physicalc, conductivity);
    %We need gradient of r w.r.t. log conductivities X, multiply each row with resp. conductivity
    d_rx = diag(conductivity)*d_r;
    Tc = FEMout.Tff';
    Tc = Tc(:);
    adjoints = get_adjoints(FEMout.globalStiffness, theta_cf, Tc, domainc);
    d_log_p = - d_rx*adjoints;
end






% [Tc, d_r, K] = FEMmain(domainc, heatSource, boundary);
% 
% WTc = W*Tc;
% 
% log_p = -.5*sum(log(diag(S))) - .5*(Tf_i - WTc)'*(Sinv*(Tf_i - WTc));
% 
% if nargout > 1  %gradient computation
%     %d_r is derivative w.r.t. conductivity lambda. we want derivative w.r.t. x = log(lambda).
%     d_rx = d_r*diag(domainc.conductivity);
%     lambda = get_adjoints(K, W, Tc, Tf_i, 0*Tf_i, S, domainc);
%     d_log_p = - d_rx'*lambda;
%     
%     %Finite difference gradient check
%     FDcheck = false;
%     d = 1e-3;
%     if FDcheck
%         d_log_pFD = zeros(domainc.N_el, 1);
%         CmeshFD = domainc;
%         for i = 1:domainc.N_el
%             dXq = zeros(domainc.N_el, 1);
%             dXq(i) = d;
%             CmeshFD.conductivity = domainc.conductivity + domainc.conductivity.*dXq;
%             TcFD = FEMmain(CmeshFD, heatSource, boundary);
%             muFD = W*TcFD;
%             log_pFD = -.5*sum(log(diag(S))) - .5*(Tf_i - muFD)'*(Sinv*(Tf_i - muFD));
%             d_log_pFD(i) = (log_pFD - log_p)/d;
%         end 
%         
%         relGrad = d_log_pFD./d_log_p
% %         d_log_p
% %         d_log_pFD
% 
%         if any(relGrad > 1.2) || any(relGrad < .8)
%             relGrad
%             pause
%         end
%         
%     end
% end

end

