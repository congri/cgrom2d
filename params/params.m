%main parameter file for 2d coarse-graining

%% Initialize coarse domain
nc = 2;
domainc = Domain(nc, nc, 1, 1);
domainc = setBoundaries(domainc, [2:(4*nc)], Tb, qb);           %ATTENTION: natural nodes have to be set manually
                                                                %and consistently in domainc and domainf
domainc = setNodalCoordinates(domainc);
domainc = setBvec(domainc);
domainc = setHeatSource(domainc, zeros(domainc.nEl, 1));

%% Some predefined basis functions for linear model p_c
phi_1 = @(lambda) log(size(lambda, 1)/sum(1./lambda));
phi_2 = @(lambda) log(mean(lambda));
phi_3 = @(lambda) mean(log(lambda));
linpathphase1l2 = @(lambda) .5*linealPath(lambda, 2, 'x', 1, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 2, 'y', 1, fineData, domainc, domainf);
linpathphase1l3 = @(lambda) .5*linealPath(lambda, 3, 'x', 1, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 3, 'y', 1, fineData, domainc, domainf);
linpathphase2l2 = @(lambda) .5*linealPath(lambda, 2, 'x', 2, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 2, 'y', 2, fineData, domainc, domainf);
linpathphase2l3 = @(lambda) .5*linealPath(lambda, 3, 'x', 2, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 3, 'y', 2, fineData, domainc, domainf);
linpathphase2l6 = @(lambda) .5*linealPath(lambda, 6, 'x', 2, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 6, 'y', 2, fineData, domainc, domainf);
linpathphase2l9 = @(lambda) .5*linealPath(lambda, 9, 'x', 2, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 9, 'y', 2, fineData, domainc, domainf);
linpathphase1l1 = @(lambda) .5*linealPath(lambda, 1, 'x', 1, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 1, 'y', 1, fineData, domainc, domainf);
linpathphase2l1 = @(lambda) .5*linealPath(lambda, 1, 'x', 2, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 1, 'y', 2, fineData, domainc, domainf);
volFrac1 = @(lambda) .5*linealPath(lambda, 0, 'x', 1, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 0, 'y', 1, fineData, domainc, domainf);
volFrac2 = @(lambda) .5*linealPath(lambda, 0, 'x', 2, fineData, domainc, domainf) +...
    .5*linealPath(lambda, 0, 'y', 2, fineData, domainc, domainf);

phi = {phi_3; volFrac2; linpathphase2l3; linpathphase2l6; linpathphase2l9};
nBasis = numel(phi);

%% EM params
basisFunctionUpdates = 0;
basisUpdateGap = 50;
maxIterations = (basisFunctionUpdates + 1)*basisUpdateGap;
%jobname for saving data
jobname = datestr(now, 'mmmddHHMMSS');
mkdir(strcat('./data/', jobname));

%% Start value of model parameters
Winterp = true;
if Winterp
    %Shape function interpolate in W
    theta_cf.W = shapeInterp(domainc, domainf);
else
    %random initialization of W
    theta_cf.W = (2/domainc.nEq)*rand(domainf.nNodes, domainc.nNodes);
end

theta_cf.S = 1*ones(domainf.nNodes, 1);
theta_cf.Sinv = sparse(1:domainf.nNodes, 1:domainf.nNodes, 1./theta_cf.S);
%precomputation to save resources
theta_cf.WTSinv = theta_cf.W'*theta_cf.Sinv;
theta_cf.mu = zeros(domainf.nNodes, 1);
% theta_c.theta = (1/size(phi, 1))*ones(size(phi, 1), 1);
theta_c.theta = 3*ones(nBasis, 1);
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
MCMC.nSamples = 50;                                 %number of samples
MCMC.nGap = 100;                                     %decorrelation gap
MCMC.Xi_start = 20*ones(domainc.nEl, 1);
%only for random walk
MCMC.MALA.stepWidth = .005;
stepWidth = 2e-0;
MCMC.randomWalk.proposalCov = stepWidth*eye(domainc.nEl);   %random walk proposal covariance
MCMC = repmat(MCMC, fineData.nSamples, 1);

%% MCMC options for test chain to find step width
MCMCstepWidth = MCMC;
for i = 1:fineData.nSamples
    MCMCstepWidth(i).nSamples = 2;
    MCMCstepWidth(i).nGap = 200;
end

%% Control convergence velocity - take weighted mean of adjacent parameter estimates
mix_sigma = 0;
mix_S = 0.6;
mix_W = 0;
mix_theta = 0;



