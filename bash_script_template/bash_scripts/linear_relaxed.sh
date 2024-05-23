#!/bin/bash

RUNSCRIPTS=~/projects/def-gcohenfr/ntamvada/simulations/runscripts/linear_relaxed
LOGS=~/projects/def-gcohenfr/ntamvada/simulations/logs/linear_relaxed

#SBATCH --time=2:00:00
#SBATCH --mem=128GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=<aromanus@gmail.com>
#SBATCH --output=${LOGS}/linear_relaxed.txt
#SBATCH --job-name="linear_relaxed"

module load StdEnv/2020 r/4.2.2

Rscript relaxed_simulation_script.R

">$RUNSCRIPTS/run.job
echo "Submitting $RUNSCRIPTS/${ID}.run.job to the cluster"
sbatch $RUNSCRIPTS/${ID}.run.job
done