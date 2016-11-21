%Test script to test variational inference

clear all
%Set VI params
dim = 1;        %ATTENTION: set params accordingly
params.family = 'diagonalGaussian';
params.initialParams{1} = zeros(1, dim);
params.nSamples = 10;
params.RMstepWidth = .1;
params.RMoff = 100;
params.robbinsMonro.relXtol = 1e-4;
params.BFGS.c = .1;
params.BFGS.lambda = 1;
params.BFGS.epsilon = 1e-10;


if strcmp(params.family, 'diagonalGaussian')
    %Take a diagonal Gaussian as the true distribution
    params.initialParams{2} = zeros(1, dim);
    trueMu = [3.5];
    trueCovarDiag = [4];
    dim = length(trueMu);
    trueCovarInvDiagMat = sparse(1:dim, 1:dim, 1./trueCovarDiag);
    logTrueDist = @(x) -.5*dim*log(2*pi) - .5*sum(log(trueCovarDiag)) - .5*(x - trueMu)*trueCovarInvDiagMat*(x - trueMu)';
end



[optVarDist] = variationalInference(logTrueDist, params)