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
nf = 128;       %Should be 2^n
disp('Generate domain object...')
domainf = Domain(nf, nf, 1, 1);
domainf = setBoundaries(domainf, [2:(4*nf)], Tb, qb);
disp('done')
toc

%% Finescale data params
disp('Setting up finescale data parameters...')
fineData.genData = true;    %generate new dataset?
fineData.dist = 'correlated_binary';   %uniform, gaussian or binary (dist of log conductivity)
fineData.nSamples = 15;
fineData.nTest = 1;
if strcmp(fineData.dist, 'gaussian')
    fineData.mu = 1.2;      %mean of log of lambda
    fineData.sigma = .3;    %sigma of log of lambda
elseif (strcmp(fineData.dist, 'uniform') || strcmp(fineData.dist, 'binary')...
        || strcmp(fineData.dist, 'predefined_binary') || strcmp(fineData.dist, 'correlated_binary'))
    %for uniform & binary
    fineData.lo = 1;
    fineData.up = 100;   %Change jobfile if you change this line number!
    contrast = fineData.up/fineData.lo;
    %for binary
    if (strcmp(fineData.dist, 'binary') || strcmp(fineData.dist, 'correlated_binary'))
        fineData.p_lo = .5;  %this controls volume fraction
        if strcmp(fineData.dist, 'correlated_binary')
            fineData.lx = 10*domainf.lElX;
            fineData.ly = 10*domainf.lElY;
            fineData.sigma_f2 = 1; %GP variance parameter. Has this parameter any impact?
        end
    end
else
    error('unknown fineCond distribution');
end
disp('done')
toc



%% Generate finescale data
[condAll, TfAll] = genData(domainf, fineData);
Tf = TfAll(:, 1:fineData.nSamples);
TfTest = TfAll(:, (fineData.nSamples + 1):end);
clear TfAll;
condAll = cell2mat(condAll)';
cond = condAll(1:fineData.nSamples, :);
condTest = condAll((fineData.nSamples + 1):end, :);
clear condAll;

%% Save data
%Name of training data file
trainFileName = strcat('train_', 'nf=', num2str(nf), '_locond=', num2str(fineData.lo),...
    '_hicond=', num2str(fineData.up), '_samples=', num2str(fineData.nSamples),...
    '_corrlength=', num2str(fineData.lx/domainf.lElX), '_volfrac=', num2str(fineData.p_lo));
%Name of test data file
testFileName = strcat('test_', 'nf=', num2str(nf), '_locond=', num2str(fineData.lo),...
    '_hicond=', num2str(fineData.up), '_samples=', num2str(fineData.nSamples),...
    '_corrlength=', num2str(fineData.lx/domainf.lElX), '_volfrac=', num2str(fineData.p_lo));
%Name of parameter file
paramFileName = strcat('param_', 'nf=', num2str(nf), '_locond=', num2str(fineData.lo),...
    '_hicond=', num2str(fineData.up), '_samples=', num2str(fineData.nSamples),...
    '_corrlength=', num2str(fineData.lx/domainf.lElX), '_volfrac=', num2str(fineData.p_lo));
%Folder where finescale data is saved
fineDataPath = '/home/constantin/matlab/data/fineData/';

disp('Finescale data generated. Saving and quitting...')
%Directly save to disc and load where needed. This saves memory.
save(strcat(fineDataPath, trainFileName), 'cond', 'Tf', '-v7.3')    %partial loading only for -v7.3
%Save test data variables under the same name
cond = condTest;
Tf = TfTest;
save(strcat(fineDataPath, testFileName), 'cond', 'Tf', '-v7.3')     %partial loading only for -v7.3

%save params
save(strcat(fineDataPath, paramFileName), 'domainf', 'fineData', 'Tb', 'qb', 'nf');
time_for_data_generation = toc

end