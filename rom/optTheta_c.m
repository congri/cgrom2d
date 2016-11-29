function [theta_c] = optTheta_c(theta_c, nTrain, nCoarse, XNormSqMean,...
    sumPhiTXmean, sumPhiSq, prior_type, prior_hyperparam)
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
    startValue = theta;
    Xtol = 1e-8;
    provide_objective = false;
    debugNRmax = false;
    gradHess = @(theta) dF2(theta, sigma2, theta_c, prior_type, prior_hyperparam, nTrain,...
    sumPhiTXmean, sumPhiSq);
    theta_old = theta;
    theta = newtonRaphsonMaximization(gradHess, startValue,...
        Xtol, provide_objective, debugNRmax);
    
    sigma2 = (1/(nCoarse*nTrain))*(sum(XNormSqMean) - 2*theta'*sumPhiTXmean...
        + theta'*sumPhiSq*theta);
    
    iter = iter + 1;
    if(iter > 30 || norm(theta_old -theta) < 1e-4)
        converged = true;
    end
    
end
theta_c.theta = theta;
theta_c.sigma = sqrt(sigma2);

end

