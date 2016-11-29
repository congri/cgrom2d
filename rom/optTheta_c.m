function [theta_c] = optTheta_c(theta_c, nTrain, nCoarse, XNormSqMean,...
    sumPhiTXmean, sumPhiSq, theta_prior_type, theta_prior_hyperparam,...
    sigma_prior_type, sigma_prior_hyperparam)
%% Find optimal theta_c and sigma

% %compute gradients of posterior lower bound F
% %sigmaTheta = [sigma^-2; theta_c]
% gradHess = @(theta) dF(theta, theta_c, prior_type, prior_hyperparam, nTrain,...
%     XNormSqMean, sumPhiTXmean, sumPhiSq, nCoarse);
% 
% %Newton-Raphson maximization
% startValue = [theta_c.theta];
% Xtol = 1e-8;
% provide_objective = false;
% debugNRmax = false;
% theta_opt = newtonRaphsonMaximization(gradHess, startValue,...
%     Xtol, provide_objective, debugNRmax);
% sigma2 = (1/(nCoarse*nTrain))*(sum(XNormSqMean) - 2*theta_opt'*sumPhiTXmean...
%     + theta_opt'*sumPhiSq*theta_opt);
% theta_c.sigma = sqrt(sigma2);
% theta_c.theta = theta_opt;


%Solve self-consistently: compute optimal sigma2, then theta, then sigma2 again and so on
theta = theta_c.theta;
sigma2 = theta_c.sigma^2;
iter = 0;
converged = false;
while(~converged)

%     sigma2 = (1/(nCoarse*nTrain))*(sum(XNormSqMean) - 2*theta'*sumPhiTXmean...
%         + theta'*sumPhiSq*theta)
    
    
    %Newton-Raphson maximization
    startValueTheta = theta;
    Xtol = 1e-8;
    provide_objective = false;
    debugNRmax = false;
    gradHessTheta = @(theta) dF2(theta, sigma2, theta_c, theta_prior_type, theta_prior_hyperparam, nTrain,...
    sumPhiTXmean, sumPhiSq);
    theta_old = theta;
    stepSizeTheta = .8;
    theta = newtonRaphsonMaximization(gradHessTheta, startValueTheta,...
        Xtol, provide_objective, stepSizeTheta, debugNRmax);
    
    gradHessLogSigmaMinus2 = @(logSigmaMinus2) dFlogSigmaMinus2(logSigmaMinus2, theta, nCoarse, nTrain, XNormSqMean,...
    sumPhiTXmean, sumPhiSq, sigma_prior_type, sigma_prior_hyperparam);
    startValueLogSigmaMinus2 = -log(sigma2);
    stepSizeSigma = .1;
    logSigmaMinus2 = newtonRaphsonMaximization(gradHessLogSigmaMinus2, startValueLogSigmaMinus2,...
        Xtol, provide_objective, stepSizeSigma, debugNRmax);
    
    sigmaMinus2 = exp(logSigmaMinus2);
    sigma2 = 1/sigmaMinus2;
    
%     sigma2 = (1/(nCoarse*nTrain))*(sum(XNormSqMean) - 2*theta'*sumPhiTXmean...
%         + theta'*sumPhiSq*theta) - (2/(nCoarse*nTrain))*...
%         d_log_prior_sigma(sqrt(sigma2), sigma_prior_type, sigma_prior_hyperparam);
    
    iter = iter + 1;
    if(iter > 30 || norm(theta_old - theta)/norm(theta) < 1e-4)
        converged = true;
    end
    
end
theta_c.theta = theta;
theta_c.sigma = sqrt(sigma2);

end

