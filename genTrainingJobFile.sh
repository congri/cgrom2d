NTRAIN=96
CONTRAST=100

DATESTR=`date +%m-%d-%H:%M:%S`	#datestring for jobfolder name
PROJECTDIR="/home/constantin/matlab/projects/cgrom2d"
JOBNAME="trainModel_nTrain${NTRAIN}contrast$CONTRAST"
JOBDIR="/home/constantin/matlab/data/$DATESTR$JOBNAME"

#Create job directory and copy source code
mkdir $JOBDIR
cp -r $PROJECTDIR/* $JOBDIR
#Change directory to job directory; completely independent from project directory
cd $JOBDIR
rm job_file.sh

#write job file
printf "#PBS -N $JOBNAME
#PBS -l nodes=1:ppn=16,walltime=120:00:00
#PBS -o $JOBDIR
#PBS -e $JOBDIR
#PBS -m abe
#PBS -M mailscluster@gmail.com

#Switch to job directory
cd $JOBDIR
#Set parameters
sed -i \"4s/.*/nTrain = $NTRAIN;/\" ./params/params.m
sed -i \"6s/.*/jobname = '$JOBNAME';/\" ./params/params.m
sed -i \"5s/.*/contrast = $CONTRAST;/\" ./loadTrainingData.m


#Run Matlab
/home/constantin/Software/matlab2016b/bin/matlab -nodesktop -nodisplay -nosplash -r \"trainModel ; quit;\"" >> job_file.sh

chmod +x job_file.sh
#directly submit job file
qsub job_file.sh
#./job_file.sh

