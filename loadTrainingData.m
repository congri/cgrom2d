%Script that loads finescale data

%Which data to load? We will get error if data doesn't exist
nf = 512;
fineData.lo = 10;
fineData.up = 100;
nSamples = 96;
corrlength = 20;
%Folder where finescale data is saved
fineDataPath = strcat('/home/constantin/matlab/data/fineData/');
%Name of training data file
trainFileName = strcat('train_', 'nf=', num2str(nf), '_locond=', num2str(fineData.lo),...
    '_hicond=', num2str(fineData.up), '_samples=', num2str(nSamples),...
    '_corrlength=', num2str(20), '_volfrac=', num2str(0));
%Name of parameter file
paramFileName = strcat('param_', 'nf=', num2str(nf), '_locond=', num2str(fineData.lo),...
    '_hicond=', num2str(fineData.up), '_samples=', num2str(nSamples),...
    '_corrlength=', num2str(20), '_volfrac=', num2str(0));

%load data params
load(strcat(fineDataPath, paramFileName));

%load finescale temperatures partially
Tffile = matfile(strcat(fineDataPath, trainFileName));