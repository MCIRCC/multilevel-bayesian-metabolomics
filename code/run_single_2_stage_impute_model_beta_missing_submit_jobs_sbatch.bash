#!/bin/bash

OUT_DIR=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200305_2_stage_impute_beta_missing/
PROG_DIR=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/


cd $OUT_DIR

N_SIM=200
FRAC=0.4
N_SAMPLE=100
for SIM_NUMBER in $(eval echo {1..$N_SIM})
do
echo "SIM_NUMBER: $SIM_NUMBER"
MISSING_RATE=0.01
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_beta_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.05
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_beta_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.1
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_beta_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.15
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_beta_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.2
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_beta_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_beta_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.3
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_beta_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE


done