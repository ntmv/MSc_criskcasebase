#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=256M
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="batch_n400iid_midsparse"

for i in {1..1000}

do
     ID=${i}
     RUNSCRIPTS=~/projects/def-gcohenfr/ntamvada/simulations/runscripts/n400_iid_midsparse
     LOGS=~/projects/def-gcohenfr/ntamvada/simulations/logs/n400_iid_midsparse
echo "#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=128GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=${LOGS}/${ID}.std_out.2.txt
#SBATCH --job-name="n400_iid_midsparse"

module load StdEnv/2020 r/4.2.2

Rscript var_sel_iidmidsparse.R

">$RUNSCRIPTS/${ID}.run.job
echo "Submitting $RUNSCRIPTS/${ID}.run.job to the cluster"
sbatch $RUNSCRIPTS/${ID}.run.job
done