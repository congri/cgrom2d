function [d_r] = FEMgrad(FEMout, domain, physical, conductivity)
%Compute derivatives of FEM equation system r = K*Y - F w.r.t. Lambda_e
%ONLY VALID FOR ISOTROPIC HEAT CONDUCTIVITY MATRIX D!!!


% (d/d Lambda_e) k^(e) = (1/Lambda_e) k^(e)     as k^(e) linear in Lambda_e
d_r = zeros(domain.nEl, domain.nEq);
for e = 1:domain.nEl
    gradLocStiff = zeros(4, 4, domain.nEl);
    gradLocStiff(:, :, e) = FEMout.localStiffness(:, :, e)/conductivity(e);     %gradient of local stiffnesses
    gradK = get_glob_stiff2(domain, gradLocStiff);
    gradF = get_glob_force(domain, physical, gradLocStiff);
    
    d_r(e, :) = (gradK*FEMout.naturalTemperatures - gradF)';
end


end

