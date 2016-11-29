function [dF_dlogSigmaMinus2, d2F_dlogSigmaMinus22] = dFlogSigmaMinus2(logSigmaMinus2, theta, nCoarse, nTrain, XNormSqMean,...
    sumPhiTXmean, sumPhiSq, prior_type, prior_param)
%Gradient and Hessian of EM lower bound w.r.t. log sigma^-2

sigmaMinus2 = exp(logSigmaMinus2);
if strcmp(prior_type, 'none')
    dlogprior_dsigmaMinus2 = 0;
    d2logprior_dsigmaMinus22 = 0;
else
    [~, dlogprior_dsigmaMinus2, d2logprior_dsigmaMinus22] = log_prior_sigma((1/sigmaMinus2), prior_type, prior_param);
end

dF_dsigmaMinus2 = .5*nCoarse*nTrain/sigmaMinus2 - .5*(sum(XNormSqMean) - 2*theta'*sumPhiTXmean...
        + theta'*sumPhiSq*theta) + dlogprior_dsigmaMinus2;
d2F_dsigmaMinus22 = -.5*nCoarse*nTrain/(sigmaMinus2^2) + d2logprior_dsigmaMinus22;

dF_dlogSigmaMinus2 = sigmaMinus2*dF_dsigmaMinus2;
d2F_dlogSigmaMinus22 = dF_dlogSigmaMinus2 + sigmaMinus2*d2F_dsigmaMinus22;

end

