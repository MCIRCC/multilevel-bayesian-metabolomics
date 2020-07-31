# File Notes

### Diagrams

* https://www.lucidchart.com/documents/edit/4fc037b3-ee2f-44f3-ac08-3acc1528d2d8/Pn7HC6StM9yT

## Final models

* The final models were all done using the joint 2-stage model

```
2_stage_model_simulation_test.ipynb
run_single_2_stage_impute_model_beta_missing.R
run_2_stage_util.R
```

## Result Data For Teddy

#### Missing + Imputation

```
analyze_run_single_results_jointly.ipynb
```

```
/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis//20200305_2_stage_impute_beta_missing
```

#### No missing + No Imputation
```
analyze_run_single_results_jointly_no_missing.ipynb
```

```
/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis//20200315_2_stage_impute_beta_no_missing
```



#### copy data to box

```

lftp
open ftps://cgillies@umich.edu@ftp.box.com

mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200315_2_stage_impute_beta_no_missing/*.csv -O ./
mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200315_2_stage_impute_beta_no_missing/*.pdf -O ./

mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200305_2_stage_impute_beta_missing/*.csv -O ./
mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200305_2_stage_impute_beta_missing/*.pdf -O /MIDAS/results_data/20200305_2_stage_impute_beta_missing/



mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200305_2_stage_impute_beta_missing/RAW_SIM_NUMBER_1_N_SAMPLE_100_FRAC_SIG_0.4* -O /MIDAS/results_data/20200305_2_stage_impute_beta_missing/

```

#### copy data to box review response

```
lftp
open ftps://cgillies@umich.edu@ftp.box.com
mkdir /MIDAS/results_data/20200319_2_stage_impute_beta_no_missing/
cd /MIDAS/results_data/20200319_2_stage_impute_beta_no_missing/

mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200319_2_stage_impute_beta_no_missing/*.csv -O ./
mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200319_2_stage_impute_beta_no_missing/*.pdf -O ./
```


```
lftp
open ftps://cgillies@umich.edu@ftp.box.com
mkdir /MIDAS/results_data/20200728_fixed_prior_no_missing/
cd /MIDAS/results_data/20200728_fixed_prior_no_missing/

mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200728_fixed_prior_no_missing/*.csv -O ./
mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200728_fixed_prior_no_missing/*.pdf -O ./
```


#### Imputation quality

* 2_stage_model_analyze_imputation_quality.ipynb

#### Real Data Analysis

* run_2_stage_real_data.ipynb

```
lftp
open ftps://cgillies@umich.edu@ftp.box.com

mput /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200317_real_data_analysis/*.csv -O /MIDAS/results_data/20200317_real_data_analysis/


```

## OLD

## Result folders

* Normal distribution on betas
```
/home/cgillies/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200129
```

* Student t-distribution on betas
```
/home/cgillies/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200220
```

* Student t-distribution on betas and imputation results
```
/home/cgillies/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200224
```

* Imputation analysis at different missing rates
```
/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200303_beta_missing/
```

* Imputation analysis at different missing rates WITH HIGHER ADAPT DELTA AND TREE-DEPTH
    * This should improve model fitting at the cost of time
    * a more robust model is needed
```
/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200304_beta_missing/
```


* GLM only models to investigate missing rates
* Setting the seed makes the data sets consistent across runs as we increase the missing rates
```
/home/cgillies/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200228
```

* Run GLM and Imputtation model for new missing rate scenario
```
/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200302/
```

* Copy model
```
cp /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200224/imputation_model_code.stan /nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis/20200302/
```

* Two stage model results are in the following folder
  * I noticed the true effect sizes were calculated incorrectly. Sampled 10000 class 0 and 100 class 1 examples
  
```
/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis//20200305_2_stage_impute_beta_missing
```

## Ideas that seem to work
* Two stage model seems to fit better and run faster than the 'BIG' model


## Failed Ideas

* Plot the shrinkage amount against the fraction of significant metabolites
* Would laplace prior be better?
    * the answer is no
       * sparse_model_vector_impl_simulation_testing_laplace.ipynb
       
#### Ideas to test the imputation model with more missingness

* We will try to make the plot for up to 0.2 missing rate
* It looks like initialization helped with the accuracy of the fit

#### 3-04-2020 Imputation results do not look great

* The model takes a long time to fit
* And the false positive rate goes up as the missing rate goes up
* I was hoping this problem would be fixed, but it was not
* What is the benefit of imputing if it cannot fix the issue?

* Let us look at it in terms of how well the model does in terms of imputation accuracy.
    * If the imputation error was x, what would the performance be?
    * To do this, we need a model that does not impute, but opporates off of means and variances.

### 2-24-2020 Imputation model result analysis

* It looks like the imputation performance was better than the imputation model does not actually hurt performance, but this is mostly due to noise.
* We need to figure out when the power reduces for the the standard logistic regression mode
* What happens when the missingness is higher?
* What is the current missingness rate?
* Is it too low for a noticeable performance difference

## Exploration of improved performance for sparse models

* sparse_model_vector_impl.ipynb
    * This file defines an improved logistic regression model that uses a a student t-distribution for betas to allow for larger effect sizes
    * This model uses a vector implementation and is more efficient

* sparse_model_simulation_testing.ipynb
    * base simulation code for this vector implementation

* run_sparse_combo_vector.R
    * Run the vector model in Rscript

* 2_stage_model_simulation_test.ipynb
    * test code for the 2-stage imputation model

## Initial model development

* sparse_model_analysis.ipynb
    * This file defines the initial logistic regression model and fits the model in a long format

* basic_simulation.ipynb
    * This model puts together the initial code for the simulation

* run_sparse_combo.R
    * This is the code to run the simulation mode

* run_sparse_combo.bash
    * This is the driver scripe for run_sparse_combo.R

* imputation_model.ipynb
    * This file defines the