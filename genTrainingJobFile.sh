PROJECTDIR="~/matlab/projects/cgrom2d"
NTRAIN=16
CONTRAST=2
JOBNAME="trainModel*nTrain=$NTRAIN*contrast=$CONTRAST"

rm job_file.sh

#write job file
printf "#PBS -N $JOBNAME
#PBS -l nodes=1:ppn=16,walltime=120:00:00
#PBS -o /home/constantin/OEfiles
#PBS -e /home/constantin/OEfiles
#PBS -m abe
#PBS -M mailscluster@gmail.com

#Switch to project directory
cd $PROJECTDIR
#Set parameters
sed -i \"4s/.*/nTrain = $NTRAIN;/\" ./params/params.m
sed -i \"7s/.*/jobname = '$JOBNAME';/\" ./params/params.m
sed -i \"5s/contrast = $CONTRAST;/\" ./loadTrainingData.m


#Run Matlab
/home/constantin/Software/matlab2016b/bin/matlab -nodesktop -nodisplay -nosplash -r \"trainModel ; quit;\"" >> job_file.sh

#directly submit job file
qsub job_file.sh

