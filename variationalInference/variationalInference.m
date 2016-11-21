function [optVarDist, RMsteps] = variationalInference(logTrueCondDist, params)
%Function for variational inference
%Input:
%   trueDist:           handle to true conditional distribution p(hidden| observable)
%   family:             family of variational distributions
%   initialParams:      parameter initialization for variational family





if strcmp(params.family, 'isotropicGaussian')
    variationalGaussMean = params.initialParams{1};
    dim = length(params.initialParams{1});
    varVariance = params.initialParams{2};
    mlog_varVariance = -log(varVariance);
    GMean = zeros(1, dim);
    GVar = 0;
    
    converged = false;
    RMsteps = 0;
    RMsamples_mu = zeros(params.nSamples, dim);
    RMsamples_var = zeros(params.nSamples, 1);
    while ~converged
        varMean = repmat(variationalGaussMean, params.nSamples, 1);
        %Draw samples from variational distribution
        samples = normrnd(varMean, varVariance); %1 sample = 1 row
        
        %Gradient of variational distribution w.r.t. params
        samplesMinusMu = (samples - varMean);
        samplesMinusMuNormSq = sum(samplesMinusMu.^2, 2);
        d_muSamples = (1/varVariance)*samplesMinusMu;
        d_sigmaMinus2Samples = .5*dim*varVariance - .5*samplesMinusMuNormSq;
        
        for i = 1:params.nSamples
            rightFactor = logTrueCondDist(samples(i, :)) + .5*dim*log(2*pi) +...
               .5*dim*log(varVariance) + .5*(1/varVariance)*samplesMinusMuNormSq(i, :);
            RMsamples_mu(i, :) = d_muSamples(i, :)*rightFactor;
            %factor 1/varVariance due to log transformation
            RMsamples_var(i) = (1/varVariance)*d_sigmaMinus2Samples(i)*rightFactor;
        end
        
        %Robbins-Monro
        variationalGaussMean = mean(varMean);
        oldParams = [variationalGaussMean varVariance];
        meanGrad = mean(RMsamples_mu);
        GMean = GMean + meanGrad.^2;    %AdaGrad
        rhoMean = params.robbinsMonro.stepWidth*GMean.^(-.5);
%         varMeanMean = varMeanMean + (params.robbinsMonro.stepWidth/(params.robbinsMonro.offset + RMsteps))*...
%             mean(RMsamples_mu);
        variationalGaussMean = variationalGaussMean + rhoMean.*meanGrad;
        varGrad = mean(RMsamples_var);
        GVar = GVar + varGrad^2;    %AdaGrad
        rhoVar = params.robbinsMonro.stepWidth*GVar^(-.5);
        mlog_varVariance = mlog_varVariance + rhoVar*varGrad;
        varVariance = 1/(exp(mlog_varVariance));
        %converged?
        if (norm(oldParams - [variationalGaussMean varVariance])/norm([variationalGaussMean varVariance]) < params.robbinsMonro.relXtol)
            converged = true;
        end
        RMsteps = RMsteps + 1;
    end
    
    %Construct function handle to optimal variational distribution
    optVarDist.dist = @(x) normpdf(x, variationalGaussMean, sqrt(varVariance));
    optVarDist.params{1} = variationalGaussMean;
    optVarDist.params{2} = varVariance;
    
elseif strcmp(params.family, 'diagonalGaussian')

    variationalGaussMean = params.initialParams{1};
    dim = length(variationalGaussMean);
    variationalGaussVariance = params.initialParams{2}; %Vector of variances in each dimension
    variationalNegLogVariance = - log(variationalGaussVariance);
    variationalParameters = [variationalGaussMean variationalNegLogVariance];
    GMean = zeros(1, dim);
    deltaMeanSq = ones(1, dim);
    deltaVarSq = ones(1, dim);
    GVar = zeros(1, dim);
    offset = 1e-3;
    if strcmp(params.RMtype, 'adam')
       momentum = zeros(1, 2*dim);
       uncenteredParameterVariance = zeros(1, 2*dim);
    end

    
    converged = false;
    RMsteps = 0;
    RMsamples_mu = zeros(params.nSamples, dim);
    RMsamples_var = zeros(params.nSamples, dim);
    d_logSigmaMinus2 = RMsamples_var;
    while ~converged

        samples = normrnd(0, 1, params.nSamples, dim); %1 sample = 1 row
        grad = ELBOgrad(samples, variationalParameters, logTrueCondDist, params)
        meanGrad = grad(1:dim);
        varGrad = grad((dim + 1):end);
        
        %Robbins-Monro
        oldParams = variationalParameters;
           
        if strcmp(params.RMtype, 'adaGrad')
            GMean = GMean + meanGrad.^2;    %AdaGrad
            rhoMean = params.robbinsMonro.stepWidth*GMean.^(-.5);
            variationalGaussMean = variationalGaussMean + rhoMean.*meanGrad
            GVar = GVar + varGrad.^2;
            rhoVar = params.robbinsMonro.stepWidth*GVar.^(-.5);
            variationalNegLogVariance = variationalNegLogVariance + rhoVar.*varGrad;
            variationalGaussVariance = 1./(exp(variationalNegLogVariance))
        elseif strcmp(params.RMtype, 'adaDelta')
            GMean = params.decayParam*GMean + (1 - params.decayParam)*meanGrad.^2;
            GVar = params.decayParam*GVar +(1 - params.decayParam)*varGrad.^2;
            currDeltaMean = sqrt(deltaMeanSq./(GMean + offset)).*meanGrad;
            currDeltaVar = sqrt(deltaVarSq./(GVar + offset)).*varGrad;
            variationalGaussMean = variationalGaussMean + currDeltaMean;
            variationalNegLogVariance = variationalNegLogVariance + currDeltaVar;
            deltaMeanSq = params.decayParam*deltaMeanSq + (1 - params.decayParam)*currDeltaMean.^2;
            deltaVarSq = params.decayParam*deltaVarSq + (1 - params.decayParam)*currDeltaVar.^2;
            variationalGaussVariance = 1./(exp(variationalNegLogVariance));
        elseif strcmp(params.RMtype, 'rmsProp')
            GMean = params.decayParam*GMean + (1 - params.decayParam)*meanGrad.^2;
            GVar = params.decayParam*GVar +(1 - params.decayParam)*varGrad.^2;
            variationalGaussMean = variationalGaussMean + (params.robbinsMonro.stepWidth/(params.robbinsMonro.offset + RMsteps))*(1./sqrt(GMean)).*meanGrad;
            variationalNegLogVariance = variationalNegLogVariance + (params.robbinsMonro.stepWidth/(params.robbinsMonro.offset + RMsteps))*(1./sqrt(GVar)).*varGrad;
            variationalGaussVariance = 1./(exp(variationalNegLogVariance));
        elseif strcmp(params.RMtype, 'adam')
            momentum = params.adam.beta1*momentum + (1 - params.adam.beta1)*grad;
            uncenteredParameterVariance = params.adam.beta2*uncenteredParameterVariance...
                + (1 - params.adam.beta2)*grad.^2;
            
            %parameter updates
            %First half is means, second half is log sigma^-2
            variationalParameters = variationalParameters +...
                (params.robbinsMonro.stepWidth*params.robbinsMonro.offset/(params.robbinsMonro.offset + RMsteps))*...
                (1./(sqrt(uncenteredParameterVariance) + 1e-8)).*momentum
        end
        
        %converged?
        if (norm(oldParams - variationalParameters)/norm(variationalParameters) < params.robbinsMonro.relXtol)
            converged = true;
        end
        RMsteps = RMsteps + 1;
    end
    
    %Construct function handle to optimal variational distribution
    variationalGaussMean = variationalParameters(1:dim);
    variationalGaussVariance = 1./(exp(variationalNegLogVariance));    
    optVarDist.dist = @(x) normpdf(x, variationalGaussMean, sqrt(variationalGaussVariance));
    %variationalParameters = [mean, log(sigma^-2)]
    optVarDist.params = variationalParameters;
    
    elseif strcmp(params.family, 'fullRankGaussian')
    varMean = params.initialParams{1};
    dim = length(params.initialParams{1});
    varCovar = params.initialParams{2};
    varPrecision = inv(varCovar);
    L = chol(varPrecision);
    GMean = zeros(1, dim);
    GL = zeros(dim);
    logpi = .5*dim*log(2*pi);
    
    converged = false;
    RMsteps = 0;
    RMsamples_mu = zeros(params.nSamples, dim);
    RMsamples_L = zeros(dim, dim, params.nSamples);
    samplesMinusMu = zeros(params.nSamples, dim);
    d_LSamples = zeros(dim, dim, params.nSamples);
    while ~converged
        %Draw samples from variational distribution
        samples = mvnrnd(varMean, varCovar, params.nSamples); %1 sample = 1 row
        
        %Gradient of variational distribution w.r.t. params
        for i = 1:params.nSamples
            samplesMinusMu(i, :) = (samples(i, :) - varMean);
        end
        d_muSamples = samplesMinusMu*varPrecision;
        %L is an upper triangle Cholesky factorization matrix of varPrecision
        invLT = inv(L)';
        for i = 1:params.nSamples
            d_LSamples(:, :, i) = invLT - (samplesMinusMu(i, :)'*samplesMinusMu(i, :))*L;
        end
        
        logdetL = sum(log(diag(abs(L))));
        for i = 1:params.nSamples
            %The right factor in eq. 3 in Ranganath, Gerrish, Blei
            rightFactor = logTrueCondDist(samples(i, :)) + logpi - logdetL...
                + .5*samplesMinusMu(i, :)*varPrecision*samplesMinusMu(i, :)';
            RMsamples_mu(i, :) = d_muSamples(i, :)*rightFactor;
            RMsamples_L(:, :, i) = d_LSamples(:, :, i)*rightFactor;
        end
        
        %Robbins-Monro
        oldParams = [varMean(:); varCovar(:)];
        meanGrad = mean(RMsamples_mu);
        GMean = GMean + meanGrad.^2;    %AdaGrad
        rhoMean = params.robbinsMonro.stepWidth*GMean.^(-.5);
        varMean = varMean + rhoMean.*meanGrad;
        LGrad = mean(RMsamples_L, 3);
        GL = GL + LGrad.^2;
        rhoL = params.robbinsMonro.stepWidth*GL.^(-.5);
        %Add stability term to avoid singular matrix L
        L = L + rhoL.*LGrad + (1e-10)*eye(dim);
        varPrecision = L'*L;
        varCovar = inv(varPrecision);
        
        %converged?
        if (norm(oldParams - [varMean(:); varCovar(:)])/norm([varMean(:); varCovar(:)]) < params.robbinsMonro.relXtol)
            converged = true;
        end
        RMsteps = RMsteps + 1;
    end
    
    %Construct function handle to optimal variational distribution
    optVarDist.dist = @(x) mvnpdf(x, varMean, varCovar);
    optVarDist.params{1} = varMean;
    optVarDist.params{2} = varCovar;
    
else
    error('Unknown parametric family for variational inference')
end

end

