function [cond, Tf] = genData(domain, fineData, FD)
%Generating full order data

fineData
FD
pause

%% Draw conductivity/ log conductivity
if strcmp(FD.distributionType, 'uniform')
    %conductivity uniformly distributed between lo and up
    cond = (fineData.up - fineData.lo)*rand(domain.nEl, fineData.nSamples) + fineData.lo;
elseif strcmp(FD.distributionType, 'gaussian')
    %log conductivity gaussian distributed
    x = normrnd(fineData.mu, fineData.sigma, domain.nEl, fineData.nSamples);
    cond = exp(x);
elseif strcmp(FD.distributionType, 'binary')
    %binary distribution of conductivity (Bernoulli)
    r = rand(domain.nEl, fineData.nSamples);
    cond = fineData.lo*ones(domain.nEl, fineData.nSamples);
    cond(r > fineData.p_lo) = fineData.up;
elseif strcmp(FD.distributionType, 'correlated_binary')
    %Compute coordinates of element centers
    x = (domain.lElX/2):domain.lElX:(1 - (domain.lElX/2));
    y = (domain.lElY/2):domain.lElY:(1 - (domain.lElY/2));
    [X, Y] = meshgrid(x, y);
    %directly clear potentially large arrays
    clear y;
    x = [X(:) Y(:)]';
    clear X Y;
    nBochnerBasis = 1e3;    %Number of cosine basis functions
    parPoolInit();
    cond{1} = zeros(domain.nEl, 1);
    cond = repmat(cond, 1, (fineData.nSamples + fineData.nTest));
    disp('Generating conductivity samples...')
    for i = 1:(fineData.nSamples + fineData.nTest)
        p{i} = genBochnerSamples(fineData.lx, fineData.sigma_f2, nBochnerBasis);
    end
    parfor i = 1:(fineData.nSamples + fineData.nTest)
        %use for-loop instead of vectorization to save memory
        p_lo = norminv(1 - fineData.theoreticalVolumeFraction, 0, fineData.sigma_f2);
        for j = 1:domain.nEl
            ps = p{i}(x(:, j));
            cond{i}(j) = fineData.up*(ps > fineData.p_lo) +...
                fineData.lo*(ps <= fineData.p_lo);
        end
    end

else
    error('unknown FOM conductivity distribution');
end
t = toc

%% Compute output data (finescale nodal temperatures)
Tf = zeros(domain.nNodes, (fineData.nSamples + fineData.nTest));
D{1} = zeros(2, 2, domain.nEl);
D = repmat(D, (fineData.nSamples + fineData.nTest), 1);
disp('Solving finite element system...')
parfor i = 1:(fineData.nSamples + fineData.nTest)
    %Conductivity matrix D, only consider isotropic materials here
    for j = 1:domain.nEl
        D{i}(:, :, j) =  cond{i}(j)*eye(2);
    end
    FEMout = heat2d(domain, D{i});
    %Store fine temperatures as a vector Tf. Use reshape(Tf(:, i), domain.nElX + 1, domain.nElY + 1)
    %and then transpose result to reconvert it to original temperature field
    Ttemp = FEMout.Tff';
    Tf(:, i) = Ttemp(:);
end
disp('FEM systems solved...')
t = toc



end

