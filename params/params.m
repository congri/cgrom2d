%main parameter file for 2d coarse-graining

%initialize domain
nF = 6;
nC = 3;
%conductivity in unit cell
lambdaUnit = [1 3; 2 4];
assert(mod(nF, 2) == 0, 'error: nF not divisible by 2')
dFine = Domain(nF, nF, 1, 1);
dCoarse = Domain(nC, nC, 1, 1);
%boundary conditions
boundaryConditions;
dFine = setNodalCoordinates(dFine, physical);
dFine = setBvec(dFine);

dCoarse = setNodalCoordinates(dCoarse, physical);
dCoarse = setBvec(dCoarse);

%physical params, boundary conditions etc.
physC = physicalParams(dCoarse, physical);
physF = physicalParams(dFine, physical);

%Finescale conductivity params
fineCond.mu = 3;    %mean of log of lambda
fineCond.sigma = .5; %sigma of log of lambda
fineCond.nSamples = 14;

%Define basis functions for p_c here
phi_1 = @(x) size(x, 1)/sum(1./x);
phi_2 = @(x) mean(x);
phi = {phi_1, phi_2};

%maximum iterations in EM
maxIterations = 4000;

%start value of model parameters
theta_cf.W = (2/dCoarse.nEq)*rand(dFine.nEq, dCoarse.nEq);
theta_cf.S = 30*eye(dFine.nEq);
theta_cf.mu = zeros(dFine.nEq, 1);
theta_c.theta = [.5, .5];
theta_c.sigma = 10;

%MCMC options
MCMC.method = 'randomWalk';                             %proposal type: randomWalk, nonlocal or MALA
MCMC.seed = 6;
MCMC.nThermalization = 50;                              %thermalization steps
MCMC.nSamples = 40;                                    %number of samples
MCMC.nGap = 100;
MCMC.Xi_start = 20*ones(dCoarse.nEl, 1);
%only for random walk
MCMC.MALA.stepWidth = .01;
stepWidth = 2e-0;
MCMC.randomWalk.proposalCov = stepWidth*eye(dCoarse.nEl);   %random walk proposal covariance
MCMC = repmat(MCMC, fineCond.nSamples, 1);

%Control convergence velocity - take weighted mean of adjacent parameter estimates
mix_sigma = 0;
mix_S = 0;
mix_W = 0;
mix_theta = 0;



