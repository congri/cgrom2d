%% Some predefined basis functions for linear model p_c
phi_1 = @(lambda) log(size(lambda, 1)/sum(1./lambda));
phi_2 = @(lambda) log(mean(lambda));
phi_3 = @(lambda) mean(log(lambda));
phi_4 = @(lambda) 1;

dLinPathMax = 21;
dLinPathMin = 0;
dLinPathIncr = 1;
i = 1;
nElc = [domainc.nElX domainc.nElY];
nElf = [domainf.nElX domainf.nElY];
for d = dLinPathMin:dLinPathIncr:dLinPathMax
    linPath{i} = @(lambda) log(.5*linealPath(lambda, d, 'x', 2, fineData, nElc, nElf) +...
    .5*linealPath(lambda, d, 'y', 2, fineData, nElc, nElf));
    i = i + 1;
end


d2pointCorrMax = 21;
d2pointCorrMin = 2;
d2pointCorrIncr = 1;
i = 1;
for d = d2pointCorrMin:d2pointCorrIncr:d2pointCorrMax
    twoPointCorr{i} = @(lambda) .5*twoPointCorrelation(lambda, d, 'x', 2, fineData, nElc, nElf) +...
        .5*twoPointCorrelation(lambda, d, 'y', 2, fineData, nElc, nElf);
    i = i + 1;
end

maxPathLength = 70;
lpa = @(lambda) log(linPathParams(lambda, maxPathLength, fineData, domainc, domainf, 'a'));
lpb = @(lambda) log(abs(linPathParams(lambda, maxPathLength, fineData, domainc, domainf, 'b')));

phi = {phi_2};
nBasis = numel(phi);