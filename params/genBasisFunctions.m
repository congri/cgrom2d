%% Some predefined basis functions for linear model p_c
phi_1 = @(lambda) log(size(lambda, 1)/sum(1./lambda));
phi_2 = @(lambda) log(mean(lambda));
phi_3 = @(lambda) mean(log(lambda));
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
    twoPointCorr{i} = @(lambda) .5*twoPointCorrelation(lambda, d, 'x', 2, fineData, nElc, nElf) +...
        .5*twoPointCorrelation(lambda, d, 'y', 2, fineData, nElc, nElf);
    i = i + 1;
end

pathLengths = (0:4:60)';
lpa = @(lambda) log(linPathParams(lambda, pathLengths, fineData, domainc, domainf, 'a'));
lpb = @(lambda) abs(linPathParams(lambda, pathLengths, fineData, domainc, domainf, 'b'));

phi = {phi_1 phi_2 phi_3 lpb};
nBasis = numel(phi);