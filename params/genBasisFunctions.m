%% Some predefined basis functions for linear model p_c
phi_1 = @(lambda) size(lambda, 1)/sum(1./lambda);
phi_2 = @(lambda) mean(lambda);
phi_3 = @(lambda) exp(mean(log(lambda)));
phi_4 = @(lambda) 1;

dLinPathMax = 0;
dLinPathMin = 0;
dLinPathIncr = 0;
i = 1;
nElc = [domainc.nElX domainc.nElY];
nElf = [domainf.nElX domainf.nElY];
for d = dLinPathMin:dLinPathIncr:dLinPathMax
    linPath{i} = @(lambda) .5*linealPath(lambda, d, 'x', 2, fineData, nElc, nElf) +...
    .5*linealPath(lambda, d, 'y', 2, fineData, nElc, nElf);
    i = i + 1;
end


d2pointCorrMax = 0;
d2pointCorrMin = 0;
d2pointCorrIncr = 0;
i = 1;
for d = d2pointCorrMin:d2pointCorrIncr:d2pointCorrMax
    twoPointCorr{i} = @(lambda) log(.5*twoPointCorrelation(lambda, d, 'x', 2, fineData, nElc, nElf) +...
        .5*twoPointCorrelation(lambda, d, 'y', 2, fineData, nElc, nElf));
    i = i + 1;
end

pathLengths = (0:4:60)';
% lpa = @(lambda) linPathParams(lambda, pathLengths, fineData, domainc, domainf, 'a');
lpb = @(lambda) abs(linPathParams(lambda, pathLengths, fineData, domainc, domainf, 'b'));

mps = @(lambda) meanPoreSize(lambda, 2, fineData, nElc, nElf, 'mean')/(nElf(1)/2);
vps = @(lambda) sqrt(meanPoreSize(lambda, 2, fineData, nElc, nElf, 'var'))/(nElf(1)/2);

% phi = linPath;
% phi{1} = lpa;
% phi{end + 1} = lpb;
% phi{end + 1} = mps;
% phi{end + 1} = vps;
% phi{end + 1} = phi_1;
% phi{end + 1} = phi_2;
% phi{end + 1} = phi_3;

%harmonic, arhithmetic and geometric exponentiated
% minexp = 0;
% maxexp = 5;
% expincr = .5;
% exponents = minexp:expincr:maxexp;
% nExponents = length(exponents);
% i = 1;
% for ex = exponents
%     phi{i} = @(lambda) phi_2(lambda)^ex;
%     phi{i + nExponents} = @(lambda) phi_1(lambda)^ex;
%     phi{i + 2*nExponents} = @(lambda) phi_3(lambda)^ex;
%     i = i + 1;
% end

%generalized means
min_z = -2;
max_z = 2;
z_incr = .1;
z = min_z:z_incr:max_z;
i = 1;
for zz = z
    phi{i} = @(lambda) generalizedMean(lambda, zz)/max(lambda);
    i = i + 1;
end
phi{end + 1} = mps;
phi{end + 1} = vps;
phi{end + 1} = @(lambda) specificSurface(lambda, 2, fineData, nElc, nElf);
phi{end + 1} = lpb;



%Pixel grid
%Fine elements per coarse element in x and y directions
% xc = nElf(1)/nElc(1);
% yc = nElf(2)/nElc(2);
% firstRow = 33;
% lastRow = yc - 32;
% firstCol = 33;
% lastCol = xc - 32;
% pixelIncr = 64;
% for row = firstRow:pixelIncr:lastRow
%     for col = firstCol:pixelIncr:lastCol
%         phi{end + 1} = @(lambda) lambda(col + xc*(row - 1));
%     end
% end
% 
nBasis = numel(phi);