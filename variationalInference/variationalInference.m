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
    negLogVariance = - log(variationalGaussVariance);
    GMean = zeros(1, dim);
    deltaMeanSq = ones(1, dim);
    deltaVarSq = ones(1, dim);
    GVar = zeros(1, dim);
    offset = 1e-3;
    if strcmp(params.RMtype, 'adam')
       momentumMean = zeros(1, dim);
       unCenVarMean = momentumMean;
       momentumVar = zeros(1, dim);
       unCenVarVar = momentumVar;
    end

    
    converged = false;
    RMsteps = 0;
    RMsamples_mu = zeros(params.nSamples, dim);
    RMsamples_var = zeros(params.nSamples, dim);
    d_logSigmaMinus2 = RMsamples_var;
    reparamTrick = true;
    while ~converged

        varMean = repmat(variationalGaussMean, params.nSamples, 1);
        varCovar_diag = repmat(variationalGaussVariance, params.nSamples, 1);
        if reparamTrick
            varSamples = normrnd(zeros(params.nSamples, dim), ones(params.nSamples, dim)); %1 sample = 1 row
            zSamples = varMean + sqrt(varCovar_diag).*varSamples;
            meanGrad = zeros(1, dim);
            varGrad = zeros(1, dim);
            meanGradSq = zeros(1, dim);
            varGradSq = zeros(1, dim);
            for i = 1:params.nSamples
                [~, trueGrad] = logTrueCondDist(zSamples(i, :));
                trueGrad = trueGrad';
                variationGrad = variationalGaussMean - zSamples(i, :);
                meanGrad = meanGrad + (trueGrad - variationGrad).*ones(1, dim);
                meanGradSq = meanGradSq + ((trueGrad - variationGrad).*ones(1, dim)).^2;
                varGrad = varGrad + (trueGrad - variationGrad).*varSamples(i, :);
                varGradSq = varGradSq + ((trueGrad - variationGrad).*varSamples(i, :)).^2;
            end
            meanGrad = meanGrad/params.nSamples;
            MeanGradErr = sqrt(meanGradSq/params.nSamples - meanGrad.^2);
            varGrad = -.5*(varGrad.*sqrt(variationalGaussVariance))/params.nSamples;
            VarGradErr = sqrt(varGradSq/params.nSamples - varGrad.^2);

            
        else
            %Draw samples from variational distribution
            varSamples = normrnd(varMean, sqrt(varCovar_diag)); %1 sample = 1 row
            %Gradient of variational distribution w.r.t. params
            samplesMinusMu = (varSamples - varMean);
            samplesMinusMuSq = samplesMinusMu.^2;
            d_muSamples = (1./varCovar_diag).*samplesMinusMu;
            d_sigmaMinus2Samples = .5*varCovar_diag - .5*samplesMinusMuSq;
            
            for i = 1:params.nSamples
                rightFactor = logTrueCondDist(varSamples(i, :)) + .5*dim*log(2*pi) +...
                    .5*sum(log(variationalGaussVariance)) + .5*sum((1./variationalGaussVariance).*samplesMinusMuSq(i, :));
                RMsamples_mu(i, :) = d_muSamples(i, :)*rightFactor;
                d_logSigmaMinus2(i, :) = (1./variationalGaussVariance).*d_sigmaMinus2Samples(i, :);
                %factor 1/varVariance due to log transformation
                RMsamples_var(i, :) = d_logSigmaMinus2(i, :)*rightFactor;
            end
            
            meanGrad = mean(RMsamples_mu);
            varGrad = mean(RMsamples_var);
            
            %control variates
            a_mu = zeros(1, dim);
            a_var = zeros(1, dim);
            for i = 1:dim
                a_mu(i) = (mean(RMsamples_mu(:, i).*d_muSamples(:, i)) - meanGrad(i)*mean(d_muSamples(:, i)))/var(d_muSamples(:, i));
                a_var(i) = (mean(RMsamples_var(:, i).*d_logSigmaMinus2(:, i)) - varGrad(i)*mean(d_logSigmaMinus2(:, i)))/var(d_logSigmaMinus2(:, i));
            end
            for i = 1:dim
                meanGrad(i) = meanGrad(i) - a_mu(i)*mean(d_muSamples(:, i));
                varGrad(i) = varGrad(i) - a_var(i)*mean(d_logSigmaMinus2(:, i));
            end
        end
        
        %Robbins-Monro
        variationalGaussMean = mean(varMean);
        oldParams = [variationalGaussMean variationalGaussVariance];
        
        
        
        if strcmp(params.RMtype, 'adaGrad')
            GMean = GMean + meanGrad.^2;    %AdaGrad
            rhoMean = params.robbinsMonro.stepWidth*GMean.^(-.5);
            variationalGaussMean = variationalGaussMean + rhoMean.*meanGrad
            GVar = GVar + varGrad.^2;
            rhoVar = params.robbinsMonro.stepWidth*GVar.^(-.5);
            negLogVariance = negLogVariance + rhoVar.*varGrad;
            variationalGaussVariance = 1./(exp(negLogVariance))
        elseif strcmp(params.RMtype, 'adaDelta')
            GMean = params.decayParam*GMean + (1 - params.decayParam)*meanGrad.^2;
            GVar = params.decayParam*GVar +(1 - params.decayParam)*varGrad.^2;
            currDeltaMean = sqrt(deltaMeanSq./(GMean + offset)).*meanGrad;
            currDeltaVar = sqrt(deltaVarSq./(GVar + offset)).*varGrad;
            variationalGaussMean = variationalGaussMean + currDeltaMean;
            negLogVariance = negLogVariance + currDeltaVar;
            deltaMeanSq = params.decayParam*deltaMeanSq + (1 - params.decayParam)*currDeltaMean.^2;
            deltaVarSq = params.decayParam*deltaVarSq + (1 - params.decayParam)*currDeltaVar.^2;
            variationalGaussVariance = 1./(exp(negLogVariance));
        elseif strcmp(params.RMtype, 'rmsProp')
            GMean = params.decayParam*GMean + (1 - params.decayParam)*meanGrad.^2;
            GVar = params.decayParam*GVar +(1 - params.decayParam)*varGrad.^2;
            variationalGaussMean = variationalGaussMean + (params.robbinsMonro.stepWidth/(params.robbinsMonro.offset + RMsteps))*(1./sqrt(GMean)).*meanGrad;
            negLogVariance = negLogVariance + (params.robbinsMonro.stepWidth/(params.robbinsMonro.offset + RMsteps))*(1./sqrt(GVar)).*varGrad;
            variationalGaussVariance = 1./(exp(negLogVariance));
        elseif strcmp(params.RMtype, 'adam')
            momentumMean = params.adam.beta1*momentumMean + (1 - params.adam.beta1)*meanGrad;
%             momentumMean = momentumMean/(1 - params.adam.beta1);
            unCenVarMean = params.adam.beta2*unCenVarMean + (1 - params.adam.beta2)*meanGrad.^2;
%             unCenVarMean = unCenVarMean/(1 - params.adam.beta2);
            momentumVar = params.adam.beta1*momentumVar + (1 - params.adam.beta1)*varGrad;
%             momentumVar = momentumVar/(1 - params.adam.beta1);
            unCenVarVar = params.adam.beta2*unCenVarVar + (1 - params.adam.beta2)*varGrad.^2;
%             unCenVarVar = unCenVarVar/(1 - params.adam.beta2);
            
            %parameter updates
            variationalGaussMean = variationalGaussMean +...
                (params.robbinsMonro.stepWidth*params.robbinsMonro.offset/(params.robbinsMonro.offset + RMsteps))*...
                (1./(sqrt(unCenVarMean) + 1e-8)).*momentumMean
            negLogVariance = negLogVariance +...
                (params.robbinsMonro.stepWidth*params.robbinsMonro.offset/(params.robbinsMonro.offset + RMsteps))*...
                (1./(sqrt(unCenVarVar) + 1e-8)).*momentumVar;
            variationalGaussVariance = 1./(exp(negLogVariance))
        end
        
        %converged?
        if (norm(oldParams - [variationalGaussMean variationalGaussVariance])/norm([variationalGaussMean variationalGaussVariance]) < params.robbinsMonro.relXtol)
            converged = true;
        end
        RMsteps = RMsteps + 1;
    end
    
    %Construct function handle to optimal variational distribution
    optVarDist.dist = @(x) normpdf(x, variationalGaussMean, sqrt(variationalGaussVariance));
    optVarDist.params{1} = variationalGaussMean;
    optVarDist.params{2} = variationalGaussVariance;
    
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
        varSamples = mvnrnd(varMean, varCovar, params.nSamples); %1 sample = 1 row
        
        %Gradient of variational distribution w.r.t. params
        for i = 1:params.nSamples
            samplesMinusMu(i, :) = (varSamples(i, :) - varMean);
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
            rightFactor = logTrueCondDist(varSamples(i, :)) + logpi - logdetL...
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

