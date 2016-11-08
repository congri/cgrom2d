%Mahalanobis and relative error computation script

%Get directory name of optimal parameters
temp = dir('./data');
fdr = temp(3).name;
[Tf_meanArray, Tf_varArray, meanMahaErr, meanSqDist] = predictOutput(1000, 1, 56, 'test_nf=512_contrast=100_samples=56.mat',...
    fdr, true);
save('errComputation.mat', 'meanMahaErr', 'meanSqDist');