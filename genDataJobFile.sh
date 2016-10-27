NF=512
CONTRAST=100
PROJECTDIR="~/matlab/projects/cgrom2d"
NTRAIN=200
NTEST=56
JOBNAME="genData*nf=$NF*nTrain=$NTRAIN*contrast=$CONTRAST"

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
sed -i \"35s/.*/    fineData.up = $CONTRAST;   %%Change jobfile if you change this line number!/\" ./generateFinescaleData.m
sed -i \"15s/.*/nf = $NF;       %%Should be 2^n/\" ./generateFinescaleData.m
sed -i \"26s/.*/fineData.nSamples = $NTRAIN;/\" ./generateFinescaleData.m
sed -i \"27s/.*/fineData.nTest = $NTEST;/\" ./generateFinescaleData.m

#Run Matlab
/home/constantin/Software/matlab2016b/bin/matlab -nodesktop -nodisplay -nosplash -r \"generateFinescaleData ; quit;\"" >> job_file.sh

#directly submit job file
qsub job_file.sh

