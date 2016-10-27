%Script that loads finescale data

%Which data to load? We will get error if data doesn't exist
nf = 512;
contrast = 2;
nSamples = 255;
%Folder where finescale data is saved
fineDataPath = './data/fineData/';
%Name of training data file
trainFileName = strcat('train_', 'nf=', num2str(nf), '_contrast=', num2str(contrast), '_samples=',...
    num2str(nSamples));
%Name of parameter file
paramFileName = strcat('param_', 'nf=', num2str(nf), '_contrast=', num2str(contrast));

%load data params
load(strcat(fineDataPath, paramFileName));

%load finescale temperatures partially
Tffile = matfile(strcat(fineDataPath, trainFileName));