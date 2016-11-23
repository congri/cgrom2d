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
basisUpdateGap = 40;
maxIterations = (basisFunctionUpdates + 1)*basisUpdateGap - 1;

%% Start value of model parameters
%Shape function interpolate in W
theta_cf.W = shapeInterp(domainc, domainf);
theta_cf.S = 2*ones(domainf.nNodes, 1);
theta_cf.Sinv = sparse(1:domainf.nNodes, 1:domainf.nNodes, 1./theta_cf.S);
%precomputation to save resources
theta_cf.WTSinv = theta_cf.W'*theta_cf.Sinv;
theta_cf.mu = zeros(domainf.nNodes, 1);
% theta_c.theta = (1/size(phi, 1))*ones(size(phi, 1), 1);
% theta_c.theta = (1/nBasis)*ones(nBasis, 1);
theta_c.theta = .5*ones(nBasis, 1);
% theta_c.theta = 0;
theta_c.sigma = 1;


%what kind of prior for theta_c
prior_type = 'spikeAndSlab';                  %hierarchical_gamma, hierarchical_laplace, laplace, gaussian or none
%prior hyperparams; obsolete for no prior
% prior_hyperparam = 100*eye(size(phi, 1));         %variance of prior gaussian
% prior_hyperparamArray = .1;                          %For hierarchical Laplace: the bigger, the sharper
% prior_hyperparam = [0; 1];                        %parameters a, b of gamma hyperprior, a, b > 0, but not too small.
                                                    %The smaller b, the more aggressive sparsity
% prior_hyperparamArray = [1e-10, 1e-10];
prior_hyperparamArray = [.9 1e-4 1e1];              %Spike and slab: spike weight, spike var, slab var

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
mix_S = 0.5;
mix_W = 0;
mix_theta = 0;

%% Variational inference params
dim = domainc.nEl;
VIparams.family = 'diagonalGaussian';
initialParamsArray{1} = [2*ones(1, domainc.nEl) 10*ones(1, domainc.nEl)];
initialParamsArray = repmat(initialParamsArray, nTrain, 1);
VIparams.nSamples = 10;    %Gradient samples per iteration
VIparams.inferenceSamples = 10;
VIparams.optParams.optType = 'adam';
VIparams.optParams.dim = domainc.nEl;
VIparams.optParams.stepWidth = .03;
VIparams.optParams.XWindow = 20;    %Averages dependent variable over last iterations
VIparams.optParams.offset = 10000;
VIparams.optParams.relXtol = 1e-12;
VIparams.optParams.maxIterations = 200;
VIparams.optParams.meanGradNormTol = 10;    %Converged if norm of mean of grad over last k iterations is smaller
VIparams.optParams.gradNormTol = 10;    %Converged if average norm of gradient in last gradNormWindow iterations is below
VIparams.optParams.gradNormWindow = 10;  %gradNormTol
VIparams.optParams.decayParam = .7;   %only works for diagonal Gaussian
VIparams.optParams.adam.beta1 = .5;     %The higher this parameter, the more gradient information from previous steps is retained
VIparams.optParams.adam.beta2 = .999;



