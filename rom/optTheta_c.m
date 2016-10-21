function [theta_c] = optTheta_c(theta_c, fineData, nCoarse, XNormSqMean,...
    sumPhiTXmean, sumPhiSq, prior_type, prior_hyperparam)
%% Find optimal theta_c and sigma

%compute gradients of posterior lower bound F
%sigmaTheta = [sigma^-2; theta_c]
gradHess = @(theta) dF(theta, theta_c, prior_type, prior_hyperparam, fineData,...
    XNormSqMean, sumPhiTXmean, sumPhiSq, nCoarse);

%Newton-Raphson maximization
startValue = [theta_c.theta];
Xtol = 1e-8;
provide_objective = false;
debugNRmax = false;
theta_opt = newtonRaphsonMaximization(gradHess, startValue,...
    Xtol, provide_objective, debugNRmax);
sigma2 = (1/(nCoarse*fineData.nSamples))*(sum(XNormSqMean) - 2*theta_opt'*sumPhiTXmean...
    + theta_opt'*sumPhiSq*theta_opt);
theta_c.sigma = sqrt(sigma2);
theta_c.theta = theta_opt;

end

