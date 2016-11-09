function [optVarDist] = varInf(logTrueCondDist, params)
%Function for variational inference
%Input:
%   trueDist:           handle to true conditional distribution p(hidden| observable)
%   family:             family of variational distributions
%   initialParams:      parameter initialization for variational family
%   dim:                Distribution input dimension

if strcmp(params.family, 'isotropicGaussian')
    varMeanMean = params.initialParams{1};
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
        varMean = repmat(varMeanMean, params.nSamples, 1);
        %Draw samples from variational distribution
        varSamples = normrnd(varMean, varVariance); %1 sample = 1 row
        
        %Gradient of variational distribution w.r.t. params
        samplesMinusMu = (varSamples - varMean);
        samplesMinusMuNormSq = sum(samplesMinusMu.^2, 2);
        d_muSamples = (1/varVariance)*samplesMinusMu;
        d_sigmaMinus2Samples = .5*dim*varVariance - .5*samplesMinusMuNormSq;
        
        for i = 1:params.nSamples
            rightFactor = logTrueCondDist(varSamples(i, :)) + .5*dim*log(2*pi) +...
               .5*dim*log(varVariance) + .5*(1/varVariance)*samplesMinusMuNormSq(i, :);
            RMsamples_mu(i, :) = d_muSamples(i, :)*rightFactor;
            %factor 1/varVariance due to log transformation
            RMsamples_var(i) = (1/varVariance)*d_sigmaMinus2Samples(i)*rightFactor;
        end
        
        %Robbins-Monro
        varMeanMean = mean(varMean);
        oldParams = [varMeanMean varVariance];
        meanGrad = mean(RMsamples_mu);
        GMean = GMean + meanGrad.^2;    %AdaGrad
        rhoMean = params.robbinsMonro.stepWidth*GMean.^(-.5);
%         varMeanMean = varMeanMean + (params.robbinsMonro.stepWidth/(params.robbinsMonro.offset + RMsteps))*...
%             mean(RMsamples_mu);
        varMeanMean = varMeanMean + rhoMean.*meanGrad;
        varGrad = mean(RMsamples_var);
        GVar = GVar + varGrad^2;    %AdaGrad
        rhoVar = params.robbinsMonro.stepWidth*GVar^(-.5);
        mlog_varVariance = mlog_varVariance + rhoVar*varGrad;
        varVariance = 1/(exp(mlog_varVariance));
        %converged?
        if (norm(oldParams - [varMeanMean varVariance])/norm([varMeanMean varVariance]) < params.robbinsMonro.relXtol)
            converged = true;
        end
        RMsteps = RMsteps + 1;
    end
    
    %Construct function handle to optimal variational distribution
    optVarDist.dist = @(x) normpdf(x, varMeanMean, sqrt(varVariance));
    optVarDist.params{1} = varMeanMean;
    optVarDist.params{2} = varVariance;
    
elseif strcmp(params.family, 'diagonalGaussian')
    varMeanMean = params.initialParams{1};
    dim = length(params.initialParams{1});
    varCovar_diagMean = params.initialParams{2};
    mlog_varCovarDiag = - log(varCovar_diagMean);
    GMean = zeros(1, dim);
    GVar = zeros(1, dim);
    
    converged = false;
    RMsteps = 0;
    RMsamples_mu = zeros(params.nSamples, dim);
    RMsamples_var = zeros(params.nSamples, dim);
    while ~converged
        varMean = repmat(varMeanMean, params.nSamples, 1);
        varCovar_diag = repmat(varCovar_diagMean, params.nSamples, 1);
        %Draw samples from variational distribution
        varSamples = normrnd(varMean, varCovar_diag); %1 sample = 1 row
        
        %Gradient of variational distribution w.r.t. params
        samplesMinusMu = (varSamples - varMean);
        samplesMinusMuSq = samplesMinusMu.^2;
        d_muSamples = (1./varCovar_diag).*samplesMinusMu;
        d_sigmaMinus2Samples = .5*varCovar_diag - .5*samplesMinusMuSq;
        
        for i = 1:params.nSamples
            rightFactor = logTrueCondDist(varSamples(i, :)) + .5*dim*log(2*pi) +...
               .5*sum(log(varCovar_diagMean)) + .5*sum((1./varCovar_diagMean).*samplesMinusMuSq(i, :));
            RMsamples_mu(i, :) = d_muSamples(i, :)*rightFactor;
            %factor 1/varVariance due to log transformation
            RMsamples_var(i, :) = (1./varCovar_diagMean).*d_sigmaMinus2Samples(i, :)*rightFactor;
        end
        
        %Robbins-Monro
        varMeanMean = mean(varMean)
        oldParams = [varMeanMean varCovar_diagMean];
        meanGrad = mean(RMsamples_mu);
        GMean = GMean + meanGrad.^2;    %AdaGrad
        rhoMean = params.robbinsMonro.stepWidth*GMean.^(-.5);
        varMeanMean = varMeanMean + rhoMean.*meanGrad
        varGrad = mean(RMsamples_var);
        GVar = GVar + varGrad.^2;
        rhoVar = params.robbinsMonro.stepWidth*GVar.^(-.5);
        mlog_varCovarDiag = mlog_varCovarDiag + rhoVar.*varGrad;
        varCovar_diagMean = 1./(exp(mlog_varCovarDiag))
        pause
        
        %converged?
        if (norm(oldParams - [varMeanMean varCovar_diagMean])/norm([varMeanMean varCovar_diagMean]) < params.robbinsMonro.relXtol)
            converged = true;
        end
        RMsteps = RMsteps + 1;
    end
    
    %Construct function handle to optimal variational distribution
    optVarDist.dist = @(x) normpdf(x, varMeanMean, sqrt(varCovar_diagMean));
    optVarDist.params{1} = varMeanMean;
    optVarDist.params{2} = varCovar_diagMean;
    
else
    error('Unknown parametric family for variational inference')
end

end

