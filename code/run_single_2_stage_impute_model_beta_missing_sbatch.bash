#!/bin/bash

N_SAMPLE=$1
FRAC=$2
SIM_NUMBER=$3
MISSING_RATE=$4
#MEM=${4:-4g}

#SBATCH --account=cgillies1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=4g
#SBATCH --mail-user=cgillies@med.umich.edu
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j

echo "Sample Size: $N_SAMPLE"
echo "Fraction of sig metabolites: $FRAC"
echo "Simulation number: $SIM_NUMBER"
echo "Average missing rate: $MISSING_RATE"
#echo "Memory: $MEM"



# I recommend using the following lines to write output to indicate your script is working
if [[ $SLURM_JOB_NODELIST ]] ; then
   echo "Running on"
   scontrol show hostnames $SLURM_JOB_NODELIST
fi

module load gcc/8.2.0 R/3.5.1
# install.packages("rstan", repos="https://repo.miserver.it.umich.edu", type = "source")
# install.packages("mvtnorm", repos="https://repo.miserver.it.umich.edu", type = "source")

#PROG=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/run_sparse_combo.R
#PROG=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/run_sparse_combo_vector.R
#PROG=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/run_sparse_combo_vector_impute.R
#PROG=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/run_sparse_combo_vector_impute_beta_missing.R
PROG=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/run_single_2_stage_impute_model_beta_missing.R

CMD="Rscript --vanilla $PROG $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE"

echo $CMD

$CMD

