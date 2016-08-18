function [x, PhiArray, Tf] = genFineData(dFine, dCoarse, fineCond, phi, physF)
%Generating full order data

%generating conductivity from lognormal distribution
x = genConductivity(fineCond, dFine.nEl);

%Construct design matrices for each data point
PhiArray = zeros(dCoarse.nEl, size(phi, 1), size(x, 2));
for i = 1:size(x, 2)
    PhiArray(:,:,i) = designMatrix(phi, x(:, i), dFine, dCoarse);
end

Tf = zeros(dFine.nElX + 1, dFine.nElY + 1, fineCond.nSamples);
control.plt = 0;
for i = 1:fineCond.nSamples
    Dfine = zeros(2, 2, dFine.nEl);
    for j = 1:dFine.nEl
        Dfine(:,:,j) = x(j, i)*eye(2); %only isotropic material
    end
    [fineOut] = heat2d(dFine, physF, control, Dfine);
    Tf(:, :, i) = fineOut.Tff;
end

end

