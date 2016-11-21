function [optVarDist] = variationalInference(log_trueDist, params)
%Variational inference using oBFGS

    function [g] = meanGrad(samples, theta)
        %Gives back mean gradient of Blei eq. 3 w.r.t. mu,  log sigma^-2
        if strcmp(params.family, 'diagonalGaussian')
            m = theta(1:(length(theta)/2));
            v = exp(-theta((length(theta)/2 + 1):end));
            vinv = 1./v;
            samplesMinusMu = samples - repmat(m, 1, size(samples, 2));
            samplesMinusMuSq = samplesMinusMu.^2;
            
            %Gradient samples
            dMu_samples = repmat(vinv, 1, size(samples, 2)).*samplesMinusMu;
            dLogSigmaMinus2_samples = .5 - .5*samplesMinusMuSq.*repmat(vinv, 1, size(samples, 2));
            
            %Multiplication with right factor in Blei eq. 3
            RMsamples_mu = zeros(dim, size(samples, 2));
            RMsamples_logSigmaMinus2 = RMsamples_mu;
            for i = 1:size(samples, 2)
                rightFactor = log_trueDist(samples(:, i)') + .5*dim*log(2*pi) +...
                    .5*sum(log(v)) + .5*sum((vinv).*samplesMinusMuSq(:, i))
                RMsamples_mu(:, i) = dMu_samples(:, i)*rightFactor;
                RMsamples_logSigmaMinus2(:, i) = dLogSigmaMinus2_samples(:, i)*rightFactor;
            end
            g = [mean(RMsamples_mu, 2); mean(RMsamples_logSigmaMinus2, 2)];
        end
    end
meanGradFun = @(samples, theta) meanGrad(samples, theta);

    function [samples] = samplesFun(theta, nSamples)
        if strcmp(params.family, 'diagonalGaussian')
            m = theta(1:(length(theta)/2));
            v = exp(-theta((length(theta)/2 + 1):end));
            samples = normrnd(repmat(m, 1, nSamples), repmat(sqrt(v), 1, nSamples));
        end
    end
sampFun = @(theta, nSamp) samplesFun(theta, nSamp);
dim = length(params.initialParams{1});


[theta] = oBFGS(meanGradFun, sampFun, params.RMstepWidth, params.RMoff, params.BFGS.c,...
    params.BFGS.lambda, params.BFGS.epsilon, [params.initialParams{1}'; params.initialParams{2}'], params.nSamples)
optVarDist.theta = theta;
end

