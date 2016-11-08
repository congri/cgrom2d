%% Some predefined basis functions for linear model p_c
phi_1 = @(lambda) log(size(lambda, 1)/sum(1./lambda));
phi_2 = @(lambda) log(mean(lambda));
phi_3 = @(lambda) mean(log(lambda));
dLinPathMax = 21;
dLinPathMin = 0;
dLinPathIncr = 1;
i = 1;
for d = dLinPathMin:dLinPathIncr:dLinPathMax
    linPath{i} = @(lambda) .5*linealPath(lambda, d, 'x', 2, fineData, domainc, domainf) +...
    .5*linealPath(lambda, d, 'y', 2, fineData, domainc, domainf);
    i = i + 1;
end

d2pointCorrMax = 21;
d2pointCorrMin = 2;
d2pointCorrIncr = 1;
i = 1;
for d = d2pointCorrMin:d2pointCorrIncr:d2pointCorrMax
    twoPointCorr{i} = @(lambda) .5*twoPointCorrelation(lambda, d, 'x', 2, fineData, domainc, domainf) +...
    .5*twoPointCorrelation(lambda, d, 'y', 2, fineData, domainc, domainf);
    i = i + 1;
end

phi = {linPath{:}, twoPointCorr{:}};
nBasis = numel(phi);