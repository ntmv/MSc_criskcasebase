#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=256M
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="setting1"

for i in {1..100}

do
     ID=${i}
     RUNSCRIPTS=~/projects/def-gcohenfr/
     LOGS=~/projects/def-gcohenfr/
echo "#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --output=${LOGS}/${ID}.std_out.2.txt
#SBATCH --job-name="setting1"

module load StdEnv/2020 r/4.2.2

Rscript example_script_simulation_AR.R

">$RUNSCRIPTS/${ID}.run.job
echo "Submitting $RUNSCRIPTS/${ID}.run.job to the cluster"
sbatch $RUNSCRIPTS/${ID}.run.job
done
