%main parameter file for 2d coarse-graining
%CHANGE JOBFILE IF YOU CHANGE LINE NUMBERS!
%Number of training data samples
nTrain = 16;
dt = datestr(now, 'mmddHHMMSS')
jobname = 'trainModel_nTrain16contrast100'
% jobname = strcat(dt, jobname)     %Better to give datestring to jobfolder name, see jobfile

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
theta_c.theta = zeros(nBasis, 1);
% theta_c.theta = 0;
theta_c.sigma = 1;


%what kind of prior for theta_c
prior_type = 'hierarchical_gamma';                  %hierarchical_gamma, hierarchical_laplace, laplace, gaussian or none
%prior hyperparams; obsolete for no prior
% prior_hyperparam = 100*eye(size(phi, 1));         %variance of prior gaussian
% prior_hyperparam = 1;                             %Exponential decay parameter for laplace
% prior_hyperparam = [0; 1e-10];                        %parameters a, b of gamma hyperprior, a, b > 0, but not too small.
                                                    %The smaller b, the more aggressive sparsity
prior_hyperparamArray = [zeros(16, 1), (logspace(5, -10, 16))'];

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
VIparams.initialParams{1} = 2*ones(1, dim);
VIparams.RMtype = 'adam';
VIparams.initialParams{2} = .3*ones(1, dim);
VIparams.nSamples = 100;
VIparams.robbinsMonro.stepWidth = .1;
VIparams.robbinsMonro.offset = 100;
VIparams.robbinsMonro.relXtol = 1e-3;
VIparams.decayParam = .9;   %only works for diagonal Gaussian
VIparams.adam.beta1 = .9;
VIparams.adam.beta2 = .999;



