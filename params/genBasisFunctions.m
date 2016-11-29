%% Some predefined basis functions for linear model p_c
phi_1 = @(lambda) log(size(lambda, 1)/sum(1./lambda));
phi_2 = @(lambda) log(mean(lambda));
phi_3 = @(lambda) mean(log(lambda));
phi_4 = @(lambda) 1;

dLinPathMax = 100;
dLinPathMin = 0;
dLinPathIncr = 5;
i = 1;
nElc = [domainc.nElX domainc.nElY];
nElf = [domainf.nElX domainf.nElY];
for d = dLinPathMin:dLinPathIncr:dLinPathMax
    linPath{i} = @(lambda) log(.5*linealPath(lambda, d, 'x', 2, fineData, nElc, nElf) +...
    .5*linealPath(lambda, d, 'y', 2, fineData, nElc, nElf));
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
lpa = @(lambda) log(linPathParams(lambda, pathLengths, fineData, domainc, domainf, 'a'));
lpb = @(lambda) log(abs(linPathParams(lambda, pathLengths, fineData, domainc, domainf, 'b')));

mps = @(lambda) meanPoreSize(lambda, 2, fineData, nElc, nElf);

phi = linPath;
phi{end + 1} = lpa;
phi{end + 1} = lpb;
phi{end + 1} = @(lambda) specificSurface(lambda, 2, fineData, nElc, nElf);
phi{end + 1} = mps;
nBasis = numel(phi);