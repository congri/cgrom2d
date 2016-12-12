%Script that loads finescale data

%Which data to load? We will get error if data doesn't exist
nf = 512;
fineData.lo = 1;
fineData.up = 10000;
nSamples = 1024;
corrlength = 20;
volfrac = 0.8;
%Folder where finescale data is saved
fineDataPath = strcat('/home/constantin/matlab/data/fineData/');
%Name of training data file
trainFileName = strcat('train_', 'nf=', num2str(nf), '_locond=', num2str(fineData.lo),...
    '_hicond=', num2str(fineData.up), '_samples=', num2str(nSamples),...
    '_corrlength=', num2str(20), '_volfrac=', num2str(volfrac), '.mat');
%Name of parameter file
paramFileName = strcat('param_', 'nf=', num2str(nf), '_locond=', num2str(fineData.lo),...
    '_hicond=', num2str(fineData.up), '_samples=', num2str(nSamples),...
    '_corrlength=', num2str(20), '_volfrac=', num2str(volfrac), '.mat');

%load data params
load(strcat(fineDataPath, paramFileName));

%load finescale temperatures partially
Tffile = matfile(strcat(fineDataPath, trainFileName));