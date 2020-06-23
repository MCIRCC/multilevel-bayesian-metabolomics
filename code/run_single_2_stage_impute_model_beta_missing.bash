#conda activate R
cd /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/
# bash run_single_2_stage_impute_model_beta_missing.bash

N_SIM=2
FRAC=0.4
N_SAMPLE=100
for SIM_NUMBER in $(eval echo {1..$N_SIM})
do
echo "SIM_NUMBER: $SIM_NUMBER"
MISSING_RATE=0.01
#Rscript --vanilla run_single_2_stage_impute_model_beta_missing.R $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.05
#Rscript --vanilla run_single_2_stage_impute_model_beta_missing.R $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.1
#Rscript --vanilla run_single_2_stage_impute_model_beta_missing.R $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.15
#Rscript --vanilla run_single_2_stage_impute_model_beta_missing.R $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.2
#Rscript --vanilla run_single_2_stage_impute_model_beta_missing.R $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.25
#Rscript --vanilla run_single_2_stage_impute_model_beta_missing.R $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

MISSING_RATE=0.3
Rscript --vanilla run_single_2_stage_impute_model_beta_missing.R $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

done

