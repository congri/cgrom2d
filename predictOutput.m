function [Tf_meanArray, Tf_varArray, meanMahaErr, meanSqDist] = predictOutput(nSamples_p_c, testSample_lo, testSample_up, testFileName, modelParamsFolder, compMahalanobis)
%Function to predict finescale output from generative model

%Load test file
Tffile = matfile(strcat('~/matlab/data/fineData/', testFileName));
[theta_c, theta_cf, domainc, domainf, phi] = loadTrainedParams(modelParamsFolder);

addpath('./rom')
%design Matrix for p_c
PhiArray = designMatrix(phi, domainf, domainc, Tffile, testSample_lo, testSample_up);

tic;
%% Sample from p_c
disp('Sampling from p_c...')
Xsamples = zeros(domainc.nEl, nSamples_p_c, size(PhiArray, 3));
LambdaSamples = Xsamples;
for i = 1:size(PhiArray, 3)
    Xsamples(:, :, i) = mvnrnd(PhiArray(:, :, i)*theta_c.theta, (theta_c.sigma^2)*eye(domainc.nEl), nSamples_p_c)';
    LambdaSamples(:, :, i) = exp(Xsamples(:, :, i));
end
disp('done')

%% Run coarse model and sample from p_cf
disp('Solving coarse model and sample from p_cf...')
addpath('./heatFEM')
meanMahaErr = 0;
meanSqDist = 0;
Tf_meanArray = zeros(domainf.nNodes, size(PhiArray, 3));
Tf_varArray = Tf_meanArray;
for j = 1:size(PhiArray, 3)
    Tc = zeros(domainc.nNodes, nSamples_p_c);
    Tf_mean = zeros(domainf.nNodes, 1);
    Tf_sq_mean = zeros(domainf.nNodes, 1);
    for i = 1:nSamples_p_c
        D = zeros(2, 2, domainc.nEl);
        for e = 1:domainc.nEl
            D(:, :, e) = LambdaSamples(e, i)*eye(2);
        end
        FEMout = heat2d(domainc, D);
        Tctemp = FEMout.Tff';
        Tc(:, i) = Tctemp(:);
        
        %sample from p_cf
        mu_cf = theta_cf.mu + theta_cf.W*Tc(:, i);
        %only for diagonal S!!
        %Sequentially compute mean and <Tf^2> to save memory
        Tf_temp = normrnd(mu_cf, theta_cf.S);
        Tf_mean = ((i - 1)/i)*Tf_mean + (1/i)*Tf_temp;
        Tf_sq_mean = ((i - 1)/i)*Tf_sq_mean + (1/i)*(Tf_temp.^2);
    end
    disp('done')
    Tf_var = Tf_sq_mean - Tf_mean.^2;
    Tf_meanArray(:, j) = Tf_mean;
    Tf_varArray(:, j) = Tf_var;
    pure_prediction_time = toc
    
    if compMahalanobis
        Tf = Tffile.Tf(:, testSample_lo + j - 1);
        meanMahaErrTemp = mean(sqrt((.5./(Tf_var)).*(Tf - Tf_mean).^2));
        meanSqDistTemp = mean(((Tf - Tf_mean)./Tf).^2);
    end
    meanMahaErr = ((j- 1)/j)*meanMahaErr + (1/j)*meanMahaErrTemp;
    meanSqDist = ((j - 1)/j)*meanSqDist + (1/j)*meanSqDistTemp;
end
rmpath('./rom')
rmpath('./heatFEM')

end