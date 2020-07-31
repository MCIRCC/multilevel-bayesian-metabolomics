#conda activate R
cd /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/
# bash run_single_2_stage_impute_model_beta_missing.bash

#cp /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200305_2_stage_impute_beta_missing/*.stan /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200315_2_stage_impute_beta_no_missing/

cp /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200319_2_stage_impute_beta_missing/*.stan /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200319_2_stage_impute_beta_no_missing/

N_SIM=2
FRAC=0.4
N_SAMPLE=100
MISSING_RATE=0
for SIM_NUMBER in $(eval echo {1..$N_SIM})
do
echo "FRAC: $FRAC"
echo "N_SAMPLE: $N_SAMPLE"
echo "SIM_NUMBER: $SIM_NUMBER"

Rscript --vanilla run_single_2_stage_impute_model_nu_1_beta_no_missing.R $N_SAMPLE $FRAC $SIM_NUMBER $MISSING_RATE

done

