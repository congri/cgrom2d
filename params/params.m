%main parameter file for 2d coarse-graining
%CHANGE JOBFILE IF YOU CHANGE LINE NUMBERS!
%Number of training data samples
nTrain = 16;
%jobname for saving data
dt = datestr(now, 'mmmddHHMMSS');
jobname = 'trainModel*nTrain=16*contrast=10';
jobname = strcat(dt, jobname)

%% Initialize coarse domain
nc = 2;
domainc = Domain(nc, nc, 1, 1);
domainc = setBoundaries(domainc, [2:(4*nc)], Tb, qb);           %ATTENTION: natural nodes have to be set manually
                                                                %and consistently in domainc and domainf
% domainc = setHeatSource(domainc, zeros(domainc.nEl, 1));
%% Some predefined basis functions for linear model p_c
phi_1 = @(lambda) log(size(lambda, 1)/sum(1./lambda));
phi_2 = @(lambda) log(mean(lambda));
phi_3 = @(lambda) mean(log(lambda));
dLinPathMax = 160;
dLinPathMin = 0;
dLinPathIncr = 8;
i = 1;
for d = dLinPathMin:dLinPathIncr:dLinPathMax
    linPath{i} = @(lambda) .5*linealPath(lambda, d, 'x', 2, fineData, domainc, domainf) +...
    .5*linealPath(lambda, d, 'y', 2, fineData, domainc, domainf);
    i = i + 1;
end

d2pointCorrMax = 160;
d2pointCorrMin = 8;
d2pointCorrIncr = 8;
i = 1;
for d = d2pointCorrMin:d2pointCorrIncr:d2pointCorrMax
    twoPointCorr{i} = @(lambda) .5*twoPointCorrelation(lambda, d, 'x', 2, fineData, domainc, domainf) +...
    .5*twoPointCorrelation(lambda, d, 'y', 2, fineData, domainc, domainf);
    i = i + 1;
end

phi = {linPath{:}, twoPointCorr{:}};
nBasis = numel(phi);

%% EM params
basisFunctionUpdates = 0;
basisUpdateGap = 80;
maxIterations = (basisFunctionUpdates + 1)*basisUpdateGap - 1;

%% Start value of model parameters
%Shape function interpolate in W
theta_cf.W = shapeInterp(domainc, domainf);
theta_cf.S = 1*ones(domainf.nNodes, 1);
theta_cf.Sinv = sparse(1:domainf.nNodes, 1:domainf.nNodes, 1./theta_cf.S);
%precomputation to save resources
theta_cf.WTSinv = theta_cf.W'*theta_cf.Sinv;
theta_cf.mu = zeros(domainf.nNodes, 1);
% theta_c.theta = (1/size(phi, 1))*ones(size(phi, 1), 1);
theta_c.theta = .1*ones(nBasis, 1);
theta_c.sigma = 1;


%what kind of prior for theta_c
prior_type = 'hierarchical_gamma';                  %hierarchical_gamma, hierarchical_laplace, laplace, gaussian or none
%prior hyperparams; obsolete for no prior
% prior_hyperparam = 100*eye(size(phi, 1));         %variance of prior gaussian
% prior_hyperparam = 1;                             %Exponential decay parameter for laplace
prior_hyperparam = [0; 1];                        %parameters a, b of gamma hyperprior, a, b > 0, but not too small.
                                                    %The smaller b, the more aggressive sparsity

%% MCMC options
MCMC.method = 'MALA';                                %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 10;
MCMC.nThermalization = 0;                            %thermalization steps
nSamplesBeginning = [50];
MCMC.nSamples = 100;                                 %number of samples
MCMC.nGap = 200;                                     %decorrelation gap
MCMC.Xi_start = 20*ones(domainc.nEl, 1);
%only for random walk
MCMC.MALA.stepWidth = .005;
stepWidth = 2e-0;
MCMC.randomWalk.proposalCov = stepWidth*eye(domainc.nEl);   %random walk proposal covariance
MCMC = repmat(MCMC, nTrain, 1);

%% MCMC options for test chain to find step width
MCMCstepWidth = MCMC;
for i = 1:fineData.nSamples
    MCMCstepWidth(i).nSamples = 2;
    MCMCstepWidth(i).nGap = 100;
end

%% Control convergence velocity - take weighted mean of adjacent parameter estimates
mix_sigma = 0;
mix_S = 0.5;
mix_W = 0;
mix_theta = 0;



