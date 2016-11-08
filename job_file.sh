#PBS -N genDataNf512nTrain200contrast2
#PBS -l nodes=1:ppn=16,walltime=120:00:00
#PBS -o /home/constantin/OEfiles
#PBS -e /home/constantin/OEfiles
#PBS -m abe
#PBS -M mailscluster@gmail.com

#Switch to project directory
cd /home/constantin/matlab/data/genDataNf512nTrain200contrast2
#Set parameters
sed -i "35s/.*/    fineData.up = 2;   %Change jobfile if you change this line number!/" ./generateFinescaleData.m
sed -i "15s/.*/nf = 512;       %Should be 2^n/" ./generateFinescaleData.m
sed -i "26s/.*/fineData.nSamples = 200;/" ./generateFinescaleData.m
sed -i "27s/.*/fineData.nTest = 56;/" ./generateFinescaleData.m

#Run Matlab
/home/constantin/Software/matlab2016b/bin/matlab -nodesktop -nodisplay -nosplash -r "generateFinescaleData ; quit;"