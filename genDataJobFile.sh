NF=512
CONTRAST=2
NTRAIN=200
NTEST=56

#Set up file paths
PROJECTDIR="/home/constantin/matlab/projects/cgrom2d"
JOBNAME="genDataNf${NF}contrast${CONTRAST}"
JOBDIR="/home/constantin/matlab/data/$JOBNAME"

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
sed -i \"35s/.*/    fineData.up = $CONTRAST;   %%Change jobfile if you change this line number!/\" ./generateFinescaleData.m
sed -i \"15s/.*/nf = $NF;       %%Should be 2^n/\" ./generateFinescaleData.m
sed -i \"26s/.*/fineData.nSamples = $NTRAIN;/\" ./generateFinescaleData.m
sed -i \"27s/.*/fineData.nTest = $NTEST;/\" ./generateFinescaleData.m

#Run Matlab
/home/constantin/Software/matlab2016b/bin/matlab -nodesktop -nodisplay -nosplash -r \"generateFinescaleData ; quit;\"" >> job_file.sh

chmod +x job_file.sh
#directly submit job file
qsub job_file.sh


