%main parameter file for 2d coarse-graining
%CHANGE JOBFILE IF YOU CHANGE LINE NUMBERS!
%Number of training data samples
nTrain = 16;
dt = datestr(now, 'mmddHHMMSS')
jobname = 'trainModel_nTrain16contrast100'

%% Initialize coarse domain
genCoarseDomain;
                                                                
%% Generate basis function for p_c
genBasisFunctions;

%% EM params
basisFunctionUpdates = 0;
basisUpdateGap = 80;
maxIterations = (basisFunctionUpdates + 1)*basisUpdateGap - 1;

%% Start value of model parameters
%Shape function interpolate in W
theta_cf.W = shapeInterp(domainc, domainf);
theta_cf.S = (1e-5)*ones(domainf.nNodes, 1);
theta_cf.Sinv = sparse(1:domainf.nNodes, 1:domainf.nNodes, 1./theta_cf.S);
%precomputation to save resources
theta_cf.WTSinv = theta_cf.W'*theta_cf.Sinv;
theta_cf.mu = zeros(domainf.nNodes, 1);
% theta_c.theta = (1/size(phi, 1))*ones(size(phi, 1), 1);
% theta_c.theta = (1/nBasis)*ones(nBasis, 1);
theta_c.theta = zeros(nBasis, 1);
% theta_c.theta(end) = 1;
% theta_c.theta = 0;
theta_c.sigma = .1;


%what kind of prior for theta_c
theta_prior_type = 'hierarchical_gamma';                  %hierarchical_gamma, hierarchical_laplace, laplace, gaussian or none
sigma_prior_type = 'exponential_sigma';
%prior hyperparams; obsolete for no prior
theta_prior_hyperparam = [0, 1e-8];                   %a and b params for Gamma hyperprior
sigma_prior_hyperparam = 1e7;

%% MCMC options
MCMC.method = 'MALA';                                %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 10;
MCMC.nThermalization = 0;                            %thermalization steps
nSamplesBeginning = [40];
MCMC.nSamples = 40;                                 %number of samples
MCMC.nGap = 40;                                     %decorrelation gap
MCMC.Xi_start = 20*ones(domainc.nEl, 1);
%only for random walk
MCMC.MALA.stepWidth = .01;
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
mix_S = 0;
mix_W = 0;
mix_theta = 0;

%% Variational inference params
dim = domainc.nEl;
VIparams.family = 'diagonalGaussian';
initialParamsArray{1} = [ones(1, domainc.nEl) 10*ones(1, domainc.nEl)];
initialParamsArray = repmat(initialParamsArray, nTrain, 1);
VIparams.nSamples = 50;    %Gradient samples per iteration
VIparams.inferenceSamples = 100;
VIparams.optParams.optType = 'adam';
VIparams.optParams.dim = domainc.nEl;
VIparams.optParams.stepWidth = .01;
VIparams.optParams.XWindow = 20;    %Averages dependent variable over last iterations
VIparams.optParams.offset = 10000;  %Robbins-Monro offset
VIparams.optParams.relXtol = 1e-12;
VIparams.optParams.maxIterations = 30;
VIparams.optParams.meanGradNormTol = 10;    %Converged if norm of mean of grad over last k iterations is smaller
VIparams.optParams.gradNormTol = 10;    %Converged if average norm of gradient in last gradNormWindow iterations is below
VIparams.optParams.gradNormWindow = 10;  %gradNormTol
VIparams.optParams.decayParam = .7;   %only works for diagonal Gaussian
VIparams.optParams.adam.beta1 = .9;     %The higher this parameter, the more gradient information from previous steps is retained
VIparams.optParams.adam.beta2 = .999;



