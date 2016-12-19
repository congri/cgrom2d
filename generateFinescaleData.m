function [] = generateFinescaleData()   %Void function to keep workspace clean
%Function to generate and save finescale data
%CHANGE JOBFILE IF YOU CHANGE LINE NUMBERS!
tic;
addpath('./params')
addpath('./heatFEM')
addpath('./rom')
addpath('./genConductivity')
addpath('./computation')

%% Temperature field and gradient generating the boundary conditions
boundaryConditions;

%% Generate finescale domain
nf = 8;       %Finescale mesh size, should be 2^n
disp('Generate finescale domain...')
domainf = Domain(nf, nf, 1, 1);
domainf = setBoundaries(domainf, 2:(4*nf), Tb, qb);       %Only fix lower left corner as essential node
disp('done')
toc

%% Finescale data params
disp('Setting up finescale data parameters...')
FD = FinescaleData(1, 100);
FD.nSets = 2;
FD.nSamples = [1 2 3];
FD.distributionType = 'correlated_binary';
FD.distributionParams = {.3 [5 5] 1};
FD.nSets = 3;

if strcmp(FD.distributionType, 'correlated_binary')
    %Folder where finescale data is saved
    fineDataPath = '~/matlab/data/fineData/';
    %System size
    fineDataPath = strcat(fineDataPath, 'systemSize=', num2str(domainf.nElX), 'x', num2str(domainf.nElY), '/');
    %Type of conductivity distribution
    fineDataPath = strcat(fineDataPath, FD.distributionType, '/',...
        'IsoSEcov/', 'l=', num2str(FD.distributionParams{2}(1)/domainf.lElX),...
        '_sigmafSq=', num2str(FD.distributionParams{3}), '/volumeFraction=',...
        num2str(FD.distributionParams{1}), '/', 'locond=', num2str(FD.loCond),...
        '_upcond=', num2str(FD.upCond), '/', 'BCcoeffs=', mat2str(boundaryCoeffs), '/');
else
    error('No savepath set for other distribution than correlated_binary')
end

if ~exist(fineDataPath, 'dir')
    mkdir(fineDataPath);
end

%Generate finescale conductivity samples and solve FEM
for i = 1:FD.nSets
    filename = strcat(fineDataPath, 'set', num2str(i), '-samples=', num2str(FD.nSamples(i)));
    FD.solveFEM(domainf, i, filename);
end

%save params
save(strcat(fineDataPath, 'params.mat'), 'domainf', 'FD', 'boundaryCoeffs', 'Tb', 'qb', 'nf');

% %% Generate finescale data
% [condAll, TfAll] = genData(domainf, fineData, FD);
% 
% error
% 
% 
% 
% 
% 
% Tf = TfAll(:, 1:fineData.nSamples);
% TfTest = TfAll(:, (fineData.nSamples + 1):end);
% clear TfAll;
% condAll = cell2mat(condAll)';
% cond = condAll(1:fineData.nSamples, :);
% condTest = condAll((fineData.nSamples + 1):end, :);
% clear condAll;
% 
% %% Save data
% %Name of training data file
% trainFileName = strcat('set1-samples=', num2str(fineData.nSamples), '.mat');
% %Name of test data file
% testFileName = strcat('set2-samples=', num2str(fineData.nTest), '.mat');
% %Name of parameter file
% paramFileName = strcat('params','.mat');
% %Folder where finescale data is saved
% fineDataPath = '/home/constantin/matlab/data/fineData/';
% %System size
% fineDataPath = strcat(fineDataPath, 'systemSize=', num2str(domainf.nElX), 'x', num2str(domainf.nElY), '/');
% %Type of conductivity distribution
% fineDataPath = strcat(fineDataPath, fineData.dist, '/', 'IsoSEcov/', 'l=',...
%     num2str(fineData.lx/domainf.lElX), '_sigmafSq=', num2str(fineData.sigma_f2), '/volumeFraction=',...
%     num2str(fineData.theoreticalVolumeFraction), '/', 'locond=', num2str(fineData.lo),...
%     '_hicond=', num2str(fineData.up), '/', 'BCcoeffs=', mat2str(boundaryCoeffs), '/');
% 
% if ~exist(fineDataPath, 'dir')
%     mkdir(fineDataPath);
% end
% 
% disp('Finescale data generated. Saving and quitting...')
% %Directly save to disc and load where needed. This saves memory.
% save(strcat(fineDataPath, trainFileName), 'cond', 'Tf', '-v7.3')    %partial loading only for -v7.3
% %Save test data variables under the same name
% cond = condTest;
% Tf = TfTest;
% save(strcat(fineDataPath, testFileName), 'cond', 'Tf', '-v7.3')     %partial loading only for -v7.3
% 
% %save params
% save(strcat(fineDataPath, paramFileName), 'domainf', 'fineData', 'boundaryCoeffs', 'Tb', 'qb', 'nf');
% time_for_data_generation = toc

end