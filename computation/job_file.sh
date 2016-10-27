#PBS -N train5
#PBS -l nodes=1:ppn=16,walltime=120:00:00
#PBS -o /home/constantin/OEfiles
#PBS -e /home/constantin/OEfiles
#PBS -m abe
#PBS -M mailscluster@gmail.com
cd ~/matlab/projects/cgrom2d
/home/constantin/Software/matlab2016b/bin/matlab -nodesktop -nodisplay -nosplash -r "main ; quit;"
