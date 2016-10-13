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

rng('shuffle')  %system time seed

%% Load params
params;

%% Generate finescale dataset
if fineData.genData
    [cond, Tf] = genData(domainf, physicalf, fineData);
    % Compute and store design matrix for each data point
    PhiArray = zeros(domainc.nEl, numel(phi), fineData.nSamples);
    for i = 1:fineData.nSamples
        PhiArray(:, :, i) = designMatrix(phi, cond(:, i), domainf, domainc);
    end
    % Compute inverse of sum_i Phi^T(x_i)^Phi(x_i)
    sumPhiSq = zeros(size(phi, 1), size(phi, 1));
    for i = 1:fineData.nSamples
        sumPhiSq = sumPhiSq + PhiArray(:,:,i)'*PhiArray(:,:,i);
    end
    sumPhiSqInv = inv(sumPhiSq);
    % save data
    save('./data/fineData/fineData', 'cond', 'Tf', 'PhiArray', 'sumPhiSqInv');
else
    load('./data/fineData/fineData')
    sumPhiSq = inv(sumPhiSqInv);
end
for i = 1:fineData.nSamples
    %take MCMC initializations at mode of p_c
    MCMC(i).Xi_start = PhiArray(:,:,i)*theta_c.theta;
end

%% Open parallel pool
parPoolInit(fineData.nSamples);

%% EM optimization - main body
%store handle to every q_i in a cell array lq
log_qi = cell(fineData.nSamples, 1);

%collect data in data arrays
k = 1;  %EM iteration index
zOptVec = [];
collectData;

for k = 2:(EM.maxIterations + 1)
    %% Test run for step sizes
    disp('test sampling...')
    parfor i = 1:fineData.nSamples
        log_qi{i} = @(Xi) log_q_i(Xi, Tf(:, i), theta_cf, theta_c,...
            PhiArray(:, :, i), domainf, domainc, physicalc);
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
    
    disp('actual sampling...')
    %% Generate samples from every q_i
    parfor i = 1:fineData.nSamples
        log_qi{i} = @(Xi) log_q_i(Xi, Tf(:, i), theta_cf, theta_c,...
            PhiArray(:, :, i), domainf, domainc, physicalc);
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
        p_cf_exponent(:, i) = mean((repmat(Tf(:, i) - theta_cf.mu, 1, MCMC(i).nSamples) - theta_cf.W*Tc_samples(:, :, i)).^2, 2);
        
        if(~Winterp)
            %still has to be generalized for 2d
            %First factor for matrix W
            Wa(:, :, i) = (Tf(:, i) - theta_cf.mu)*mean(Tc_samples(:, :, i), 2)';
        end
    end
    
    %% Compute params of p_cf
    if(~Winterp)
        %this still has to be generalized to 2d
        Tc_dyadic_mean = TcDyadicMean(Tc_samples, fineData.nSamples, MCMC);
        Wa_mean = mean(Wa, 3);
        theta_cf.W = compW(Tc_dyadic_mean,  Wa_mean,...
            inv(theta_cf.S), theta_cf.W, paramIndices, constIndices);
    end
    
    %decelerate convergence of S
    lowerBoundS = 1e-10;
    theta_cf.S = (1 - mix_S)*diag(mean(p_cf_exponent, 2)) + mix_S*theta_cf.S + lowerBoundS*eye(size(theta_cf.S, 1));
    
    %% Compute theta_c and sigma if there is a prior on theta_c, sigma
    
    %sum_i Phi_i^T <X^i>_qi
    sumPhiTXmean = zeros(size(phi, 1), 1);
    for i = 1:fineData.nSamples
        sumPhiTXmean = sumPhiTXmean + PhiArray(:,:,i)'*XMean(:,i);
    end

    disp('M-step: find optimal params')
    [theta_c] = optTheta_c(theta_c, fineData, domainc.nEl, XNormSqMean,...
        sumPhiTXmean, sumPhiSq, sumPhiSqInv, mix_sigma, prior_type, prior_hyperparam);
    disp('M-step done, current theta:')
        
    curr_theta = theta_c.theta
    curr_sigma = theta_c.sigma
    
    %Start next chain at mean of p_c
    for i = 1:fineData.nSamples
        MCMC(i).Xi_start = PhiArray(:,:,i)*theta_c.theta;
        % MCMC(i).Xi_start = out(i).samples(:, end);
    end
    
    
    S = diag(theta_cf.S)'
    
    if(mod(k - 1, basisUpdateGap) == 0)
        %this still needs to be generalized for 2d
        disp('Updating basis functions phi in p_c...')

%         [thetaTildeOpt, zOpt] = optNewPhi(XMean, x, theta_c.theta, PhiArray, nFine, nCoarse);
%         zOptVec = [zOptVec zOpt]
%         if zOpt >= 30
%             phi{end + 1, 1} = @(x) log(max(x));
%         elseif zOpt <= -30
%             phi{end + 1, 1} = @(x) log(min(x));
%         else
%             phi{end + 1, 1} = @(x) (1/zOpt)*log((1/FperC)*sum(x.^zOpt));
%         end
%         theta_c.theta(end + 1, 1) = thetaTildeOpt;
%         EM.theta = [EM.theta; zeros(1, size(EM.theta, 2))];
        
        if size(phi, 1) == 1
            phi{2, 1} = phi_2;
            theta_c.theta = [theta_c.theta; 0];
            EM.theta = [EM.theta; zeros(1, size(EM.theta, 2))];
        elseif size(phi, 1) == 2
            phi{3, 1} = phi_1;
            theta_c.theta = [theta_c.theta; 0];
            EM.theta = [EM.theta; zeros(1, size(EM.theta, 2))];
        elseif size(phi, 1) == 3
            phi{4, 1} = phi_4;
            theta_c.theta = [theta_c.theta; 0];
            EM.theta = [EM.theta; zeros(1, size(EM.theta, 2))];
        else
            error('Which basis function to add?')
        end
        
        % Compute and store design matrix for each data point
        PhiArray = zeros(domainc.nEl, size(phi, 1), fineData.nSamples);
        for i = 1:size(cond, 2)
            PhiArray(:,:,i) = designMatrix(phi, cond, domainf, domainc);
        end
        % Compute inverse of sum_i Phi^T(x_i)^Phi(x_i)
        sumPhiSq = zeros(size(phi, 1), size(phi, 1));
        for i = 1:fineData.nSamples
            sumPhiSq = sumPhiSq + PhiArray(:,:,i)'*PhiArray(:,:,i);
        end
        sumPhiSqInv = inv(sumPhiSq);
    end
    % collect data in data arrays
    collectData;
end
clear i j k m Wa Wa_mean Tc_dyadic_mean log_qi p_cf_exponent curr_theta XMean XNormSqMean;























runtime = toc












