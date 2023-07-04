#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=256M
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="setting6"

for i in {1..200}

do
     ID=${i}
     RUNSCRIPTS=~/projects/def-gcohenfr/ntamvada/july3_sims_beluga/runscripts/setting6
     LOGS=~/projects/def-gcohenfr/ntamvada/july3_sims_beluga/logs/setting6
echo "#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=200GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=${LOGS}/${ID}.std_out.2.txt
#SBATCH --job-name="setting6"

module load StdEnv/2020 r/4.2.2

Rscript setting6.R

">$RUNSCRIPTS/${ID}.run.job
echo "Submitting $RUNSCRIPTS/${ID}.run.job to the cluster"
sbatch $RUNSCRIPTS/${ID}.run.job
done
