#!/bin/bash

OUT_DIR=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200319_2_stage_impute_beta_no_missing/
PROG_DIR=/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/


cd $OUT_DIR

N_SIM=100
MISSING_RATE=0
for SIM_NUMBER in $(eval echo {1..$N_SIM})
do
echo "SIM_NUMBER: $SIM_NUMBER"

FRAC=0.1
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

FRAC=0.2
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

FRAC=0.25
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

FRAC=0.3
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

FRAC=0.4
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

FRAC=0.5
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

FRAC=0.6
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

FRAC=0.7
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

FRAC=0.8
N_SAMPLE=25
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=50
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=75
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=100
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE
N_SAMPLE=150
sbatch --account=cgillies1 --time=23:00:00 --mem-per-cpu=4g --cpus-per-task=4  $PROG_DIR/run_single_2_stage_impute_model_nu_1_beta_no_missing_sbatch.bash $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE


done
