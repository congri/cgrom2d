%Script that generates a specific job folder
%This file is typically changed from the jobfile. Thus modifications have to be translated
%jobname for saving data
dt = datestr(now, 'mmmddHHMMSS');
jobname = 'trainModel*nTrain=16*contrast=2';
jobname = strcat(dt, jobname)
%create job directory
sv = true;
if (~exist(strcat('./data/', jobname), 'dir') && sv)
    mkdir(strcat('./data/', jobname));
end

%Copy parameter files to folder
copyfile('./params/params.m', strcat('./data/', jobname, '/params', jobname, '.m'));
copyfile('./params/genBasisFunctions.m', strcat('./data/', jobname, '/genBasisFunctions', jobname, '.m'));