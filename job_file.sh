#PBS -N trainModel*nTrain=16*contrast=2
#PBS -l nodes=1:ppn=16,walltime=120:00:00
#PBS -o /home/constantin/OEfiles
#PBS -e /home/constantin/OEfiles
#PBS -m abe
#PBS -M mailscluster@gmail.com

#Switch to project directory
cd ~/matlab/projects/cgrom2d
#Set parameters
sed -i "4s/.*/nTrain = 16;/" ./params/params.m
sed -i "7s/.*/jobname = 'trainModel*nTrain=16*contrast=2';/" ./params/params.m
sed -i "5s/contrast = 2;/" ./loadTrainingData.m


#Run Matlab
/home/constantin/Software/matlab2016b/bin/matlab -nodesktop -nodisplay -nosplash -r "trainModel ; quit;"