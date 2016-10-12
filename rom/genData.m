function [x, Tf] = genData(domain, physical, fineData)
%Generating full order data

%% Draw conductivity/ log conductivity
if strcmp(fineData.dist, 'uniform')
    %conductivity uniformly distributed between lo and up
    cond = (fineData.up - fineData.lo)*rand(domain.nEl, fineData.nSamples) + fineData.lo;
    x = log(cond);
elseif strcmp(fineData.dist, 'gaussian')
    %log conductivity gaussian distributed
    x = normrnd(fineData.mu, fineData.sigma, domain.nEl, fineData.nSamples);
    cond = exp(x);
elseif strcmp(fineData.dist, 'binary')
    %binary distribution of conductivity (Bernoulli)
    r = rand(domain.nEl, fineData.nSamples);
    cond = fineData.lo*ones(domain.nEl, fineData.nSamples);
    cond(r > fineData.p_lo) = fineData.up;
    x = log(cond);
else
    error('unknown FOM conductivity distribution');
end

%% Compute output data (finescale nodal temperatures)
Tf = zeros(domain.nNodes, fineData.nSamples);
D = zeros(2, 2, domain.nEl);
for i = 1:fineData.nSamples
    control.plt = false;
    %Conductivity matrix D, only consider isotropic materials here
    for j = 1:domain.nEl
        D(:, :, j) =  cond(j, i)*eye(2);
    end
    FEMout = heat2d(domain, physical, control, D);
    %Store fine temperatures as a vector Tf. Use reshape(Tf(:, i), domain.nElX + 1, domain.nElY + 1)
    %and then transpose result to reconvert it to original temperature field
    Ttemp = FEMout.Tff';
    Tf(:, i) = Ttemp(:);
end

end

