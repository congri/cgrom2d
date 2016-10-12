%main parameter file for 2d coarse-graining

%% Initialize domain
nf = 8;
nc = 2;
assert(mod(nf, nc) == 0, 'error: nF not divisible by nC')
domainf = Domain(nf, nf, 1, 1);
domainc = Domain(nc, nc, 1, 1);

%% Boundary conditions
boundaryConditions;

%% Finescale data params
fineData.genData = true;    %generate new dataset?
fineData.dist = 'binary';   %uniform, gaussian or binary (dist of log conductivity)
fineData.nSamples = 16;
if strcmp(fineData.dist, 'gaussian')
    fineData.mu = 1.2;      %mean of log of lambda
    fineData.sigma = .3;    %sigma of log of lambda
elseif (strcmp(fineData.dist, 'uniform') || strcmp(fineData.dist, 'binary'))
    %for uniform & binary
    fineData.lo = 2;
    fineData.up = 100;
    contrast = fineData.up/fineData.lo;
    %for binary
    if strcmp(fineData.dist, 'binary')
        fineData.p_lo = .4;
    end
else
    error('unknown fineCond distribution');
end


%% Some predefined basis functions for linear model p_c
phi_1 = @(lambda) log(size(lambda, 1)/sum(1./lambda));
phi_2 = @(lambda) log(mean(lambda));
phi = {phi_1; phi_2};
nBasis = numel(phi);

%% Object containing EM optimization params and stats
EM = EMstats;
basisUpdateGap = 30;        %After this number of iterations, include new basis function in p_c
EM = EM.setMaxIterations(1*basisUpdateGap - 1);
EM = EM.prealloc(fineData, domainf, domainc, nBasis);           %preallocation of data arrays

%% Start value of model parameters
Winterp = true;
if Winterp
    %Shape function interpolate in W
    theta_cf.W = shapeInterp(domainc, domainf);
else
    %random initialization of W
    theta_cf.W = (2/domainc.nEq)*rand(domainf.nNodes, domainc.nNodes);
end
theta_cf.S = 1*eye(domainf.nNodes);
theta_cf.mu = zeros(domainf.nNodes, 1);
theta_c.theta = (1/size(phi, 1))*ones(size(phi, 1), 1);
theta_c.sigma = 2;

%what kind of prior for theta_c
prior_type = 'none';                  %hierarchical_gamma, hierarchical_laplace, laplace, gaussian or none
%prior hyperparams; obsolete for no prior
% prior_hyperparam = 100*eye(size(phi, 1));         %variance of prior gaussian
% prior_hyperparam = 1;                             %Exponential decay parameter for laplace
prior_hyperparam = [0; .01];                      %parameters a, b of gamma hyperprior, a, b > 0, but not too small.
                                                    %The smaller b, the more aggressive sparsity

%% MCMC options
MCMC.method = 'MALA';                               %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 6;
MCMC.nThermalization = 0;                           %thermalization steps
MCMC.nSamples = 50;                                 %number of samples
MCMC.nGap = 100;                                    %decorrelation gap
MCMC.Xi_start = 20*ones(domainc.nEl, 1);
%only for random walk
MCMC.MALA.stepWidth = .1;
stepWidth = 2e-0;
MCMC.randomWalk.proposalCov = stepWidth*eye(domainc.nEl);   %random walk proposal covariance
MCMC = repmat(MCMC, fineData.nSamples, 1);

%% MCMC options for test chain to find step width
MCMCstepWidth = MCMC;
for i = 1:fineData.nSamples
    MCMCstepWidth(i).nSamples = 2;
    MCMCstepWidth(i).nGap = 500;
end

%% Control convergence velocity - take weighted mean of adjacent parameter estimates
mix_sigma = 0;
mix_S = 0.2;
mix_W = 0;
mix_theta = 0;



