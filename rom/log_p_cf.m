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

Tc = FEMout.Tff';
Tc = Tc(:);
WTc = W*Tc;
%only for diagonal S!
assert(isdiag(S), 'Error: matrix S not diagonal');
Sinv = diag(1./diag(S));
log_p = -.5*sum(log(diag(S))) - .5*(Tf_i - theta_cf.mu - WTc)'*(Sinv*(Tf_i - theta_cf.mu - WTc));

if nargout > 1
    %Gradient of FEM equation system w.r.t. conductivities
    d_r = FEMgrad(FEMout, domainc, conductivity);
    %We need gradient of r w.r.t. log conductivities X, multiply each row with resp. conductivity
    d_rx = diag(conductivity)*d_r;
    adjoints = get_adjoints(FEMout.globalStiffness, theta_cf, Tc, domainc, Tf_i);
    d_log_p = - d_rx*adjoints;

    
    %Finite difference gradient check
    FDcheck = false;
    if FDcheck
        disp('Gradient check log p_cf')
        d = 1e-4;
        FDgrad = zeros(domainc.nEl, 1);
        for e = 1:domainc.nEl
            conductivityFD = conductivity;
            conductivityFD(e) = conductivityFD(e) + d;
            
            DFD = zeros(2, 2, domainc.nEl);
            for j = 1:domainc.nEl
                DFD(:, :, j) =  conductivityFD(j)*eye(2);
            end
            FEMoutFD = heat2d(domainc, physicalc, control, DFD);
            TcFD = FEMoutFD.Tff';
            TcFD = TcFD(:);
            
            WTcFD = W*TcFD;
            log_pFD = -.5*sum(log(diag(S))) - .5*(Tf_i - theta_cf.mu - WTcFD)'*(Sinv*(Tf_i - theta_cf.mu - WTcFD));
            FDgrad(e) = conductivity(e)*(log_pFD - log_p)/d;
        end
        relgrad = FDgrad./d_log_p
%         d_r
%         d_rx
%         adjoints
%         d_log_p
%         FDgrad
%         conductivity
% 
%         conductivityFDcheck = conductivity + .001*(FDgrad./conductivity);
%         DFDcheck = zeros(2, 2, domainc.nEl);
%         for j = 1:domainc.nEl
%             DFDcheck(:, :, j) =  conductivityFDcheck(j)*eye(2);
%         end
%         FEMoutFDcheck = heat2d(domainc, physicalc, control, DFDcheck);
%         TcFDcheck = FEMoutFDcheck.Tff';
%         TcFDcheck = TcFDcheck(:);
%         WTcFDcheck = W*TcFDcheck;
%         log_pFDcheck = -.5*sum(log(diag(S))) - .5*(Tf_i - WTcFDcheck)'*(Sinv*(Tf_i - WTcFDcheck));
%         checkFD = log_pFDcheck - log_p
%         
%         conductivitycheck = conductivity + .001*(d_log_p./conductivity);
%         Dcheck = zeros(2, 2, domainc.nEl);
%         for j = 1:domainc.nEl
%             Dcheck(:, :, j) =  conductivitycheck(j)*eye(2);
%         end
%         FEMoutcheck = heat2d(domainc, physicalc, control, Dcheck);
%         Tccheck = FEMoutcheck.Tff';
%         Tccheck = Tccheck(:);
%         WTccheck = W*Tccheck;
%         log_pcheck = -.5*sum(log(diag(S))) - .5*(Tf_i - WTccheck)'*(Sinv*(Tf_i - WTccheck));
%         check = log_pcheck - log_p
    end %FD gradient check
    
    
end

end

