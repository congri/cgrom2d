function [variationalParams] = varInf2(trueLogDist, initialVariationalParams, params)
%Variational inference


meanGradFun = @(samples, variationalParams) ELBOgrad(samples, variationalParams, trueLogDist, params);
samplesFun = @(variationalParams, nSamples) normrnd(0, 1, nSamples, size(variationalParams, 1)/2);

[variationalParams] = oBFGS(meanGradFun, samplesFun, .3, 100, .1, 1, 1e-10, initialVariationalParams, params.nSamples);


end

