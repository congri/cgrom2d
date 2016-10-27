%% Main script for 2d coarse-graining
%% Preamble
tic;    %measure runtime
clear all;

addpath('./params')
addpath('./heatFEM')
addpath('./rom')
addpath('./computation')
addpath('./FEMgradient')
addpath('./MCMCsampler')
addpath('./optimization')
addpath('./genConductivity')

rng('shuffle')  %system time seed

%Open parallel pool
parPoolInit();

%Load training data
loadTrainingData;
%Get model and training parameters
params;
Tf = Tffile.Tf(:, 1:nTrain);        %Finescale temperatures - load partially to save memory

% Compute design matrix for each data point
PhiArray = designMatrix(phi, domainf, domainc, Tffile, nTrain);
% Compute sum_i Phi^T(x_i)^Phi(x_i)
sumPhiSq = zeros(numel(phi), numel(phi));
for i = 1:nTrain
    sumPhiSq = sumPhiSq + PhiArray(:,:,i)'*PhiArray(:,:,i);
end
%shrink finescale domain object to save memory
domainf = domainf.shrink();

%% EM optimization - main body
k = 1;          %EM iteration index
collectData;    %Write initial parametrizations to disk

for k = 2:(maxIterations + 1)
    for i = 1:nTrain
        %take MCMC initializations at mode of p_c
        MCMC(i).Xi_start = PhiArray(:, :, i)*theta_c.theta;
    end
    %% Test run for step sizes
    disp('test sampling...')
    parfor i = 1:nTrain
        Tf_i_minus_mu = Tf(:, i) - theta_cf.mu;
        log_qi{i} = @(Xi) log_q_i(Xi, Tf_i_minus_mu, theta_cf, theta_c,...
            PhiArray(:, :, i), domainc);
        %find maximum of qi for thermalization
        %start value has some randomness to drive transitions between local optima
        X_start{i} = normrnd(MCMC(i).Xi_start, .01);
        Xmax{i} = max_qi(log_qi{i}, X_start{i});
        
        %sample from every q_i
        outStepWidth(i) = MCMCsampler(log_qi{i}, Xmax{i}, MCMCstepWidth(i));
        while (outStepWidth(i).acceptance < .5 || outStepWidth(i).acceptance > .9)
            outStepWidth(i) = MCMCsampler(log_qi{i}, Xmax{i}, MCMCstepWidth(i));
            MCMCstepWidth(i).Xi_start = outStepWidth(i).samples(:, end);
            if strcmp(MCMCstepWidth(i).method, 'MALA')
                MCMCstepWidth(i).MALA.stepWidth = (1/.7)*(outStepWidth(i).acceptance + (1 - outStepWidth(i).acceptance)*.1)*...
                    MCMCstepWidth(i).MALA.stepWidth;
            elseif strcmp(MCMCstepWidth(i).method, 'randomWalk')
                MCMCstepWidth(i).randomWalk.proposalCov = (1/.7)*(outStepWidth(i).acceptance + (1 - outStepWidth(i).acceptance)*.1)*MCMCstepWidth(i).randomWalk.proposalCov;
            else
                error('Unknown MCMC method')
            end
        end
        %Set step widths and start values
        if strcmp(MCMCstepWidth(i).method, 'MALA')
            MCMC(i).MALA.stepWidth = MCMCstepWidth(i).MALA.stepWidth;
        elseif strcmp(MCMCstepWidth(i).method, 'randomWalk')
            MCMC(i).randomWalk.proposalCov = MCMCstepWidth(i).randomWalk.proposalCov;
        else
            error('Unknown MCMC method')
        end
        MCMC(i).Xi_start = MCMCstepWidth(i).Xi_start;
    end
    
    for i = 1:nTrain
        if(k - 1 <= length(nSamplesBeginning))
            %less samples at the beginning
            MCMC(i).nSamples = nSamplesBeginning(k - 1);
        end
    end
    
    disp('actual sampling...')
    %% Generate samples from every q_i
    parfor i = 1:nTrain
        Tf_i_minus_mu = Tf(:, i) - theta_cf.mu;
        log_qi{i} = @(Xi) log_q_i(Xi, Tf_i_minus_mu, theta_cf, theta_c,...
            PhiArray(:, :, i), domainc);
        %sample from every q_i
        out(i) = MCMCsampler(log_qi{i}, Xmax{i}, MCMC(i));
        %avoid very low acceptances
        while out(i).acceptance < .1
            out(i) = MCMCsampler(log_qi{i}, Xmax{i}, MCMC(i));
            %if there is a second loop iteration, take last sample as initial position
            MCMC(i).Xi_start = out(i).samples(:,end);
            if strcmp(MCMC(i).method, 'MALA')
                MCMC(i).MALA.stepWidth = (1/.9)*(out(i).acceptance + (1 - out(i).acceptance)*.1)*MCMC(i).MALA.stepWidth;
            elseif strcmp(MCMC(i).method, 'randomWalk')
                MCMC(i).randomWalk.proposalCov = .2*MCMC(i).randomWalk.proposalCov;
                                MCMC(i).randomWalk.proposalCov = (1/.7)*(out(i).acceptance + (1 - out(i).acceptance)*.1)*MCMC(i).randomWalk.proposalCov;
            else
                error('Unknown MCMC method')
            end
            warning('Acceptance ratio below .1')
        end
               
        %Refine step width
        if strcmp(MCMC(i).method, 'MALA')
            MCMC(i).MALA.stepWidth = (1/.7)*out(i).acceptance*MCMC(i).MALA.stepWidth;
        elseif strcmp(MCMC(i).method, 'randomWalk')
            MCMC(i).randomWalk.proposalCov = (1/.7)*out(i).acceptance*MCMC(i).randomWalk.proposalCov;
        else
        end
        
        XMean(:, i) = mean(out(i).samples, 2);
        XNormSqMean(i) = mean(sum(out(i).samples.^2));
        
        %for S
        %Tc_samples(:,:,i) contains coarse nodal temperature samples (1 sample == 1 column) for full order data
        %sample i
        Tc_samples(:, :, i) = reshape(cell2mat(out(i).data), domainc.nNodes, MCMC(i).nSamples);
        %only valid for diagonal S here!
        p_cf_exponent(:, i) = mean((repmat(Tf_i_minus_mu, 1, MCMC(i).nSamples) - theta_cf.W*Tc_samples(:, :, i)).^2, 2);
        
    end
    clear Tc_samples;
    
    %% M-step: determine optimal parameters given the sample set
    disp('M-step: find optimal params')
    %Optimal S (decelerated convergence)
    lowerBoundS = 1e-6;
    theta_cf.S = (1 - mix_S)*mean(p_cf_exponent, 2)...
        + mix_S*theta_cf.S + lowerBoundS*ones(domainf.nNodes, 1);
    clear p_cf_exponent;
    theta_cf.Sinv = sparse(1:domainf.nNodes, 1:domainf.nNodes, 1./theta_cf.S);
    theta_cf.WTSinv = theta_cf.W'*theta_cf.Sinv;        %Precomputation for efficiency

    %optimal theta_c and sigma
    %sum_i Phi_i^T <X^i>_qi
    sumPhiTXmean = zeros(numel(phi), 1);
    for i = 1:nTrain
        sumPhiTXmean = sumPhiTXmean + PhiArray(:,:,i)'*XMean(:,i);
    end

    [theta_c] = optTheta_c(theta_c, nTrain, domainc.nEl, XNormSqMean,...
        sumPhiTXmean, sumPhiSq, prior_type, prior_hyperparam);
    disp('M-step done, current params:')
    curr_theta = theta_c.theta
    curr_sigma = theta_c.sigma
    mean_S = mean(theta_cf.S)
    
    if(mod(k - 1, basisUpdateGap) == 0)
        %this still needs to be generalized for 2d!
        [phi, PhiArray, sumPhiSq] = updateBasisFunction(XMean, x, theta, PhiArray, nFine, nCoarse, phi);
    end
    %collect data and write it to disk periodically to save memory
    collectData;
end
%tidy up
clear i j k m Wa Wa_mean Tc_dyadic_mean log_qi p_cf_exponent curr_theta XMean XNormSqMean;
runtime = toc

































