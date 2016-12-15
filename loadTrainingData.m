%Script that loads finescale data

%Which data to load? We will get error if data doesn't exist
nf = 512;
fineData.lo = 1;
fineData.up = 100;
nSamples = 1024;
corrlength = '10';
volfrac = '0.3';  %high conducting phase volume fraction
sigma_f2 = '1';
cond_distribution = 'correlated_binary';
bc = '[-5 3 0 0]';


%Folder where finescale data is saved
fineDataPath = '/home/constantin/matlab/data/fineData/';
%System size
fineDataPath = strcat(fineDataPath, 'systemSize=', num2str(nf), 'x', num2str(nf), '/');
%Type of conductivity distribution
fineDataPath = strcat(fineDataPath, cond_distribution, '/', 'IsoSEcov/', 'l=',...
    corrlength, '_sigmafSq=', sigma_f2, '/volumeFraction=',...
    volfrac, '/', 'locond=', num2str(fineData.lo),...
    '_hicond=', num2str(fineData.up), '/', 'BCcoeffs=', bc, '/');


%Name of training data file
trainFileName = strcat('set1-samples=', num2str(nSamples), '.mat');
%Name of parameter file
paramFileName = strcat('params','.mat');

%load data params
load(strcat(fineDataPath, paramFileName));

%load finescale temperatures partially
Tffile = matfile(strcat(fineDataPath, trainFileName));