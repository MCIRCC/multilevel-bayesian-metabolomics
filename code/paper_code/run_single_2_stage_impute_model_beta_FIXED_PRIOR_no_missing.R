source('/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/stringer-race-caps/sparsity_analysis/run_2_stage_util.R')
require(tidyverse)
require(rstan)
require(ggplot2)
require(mvtnorm)
args = commandArgs(trailingOnly=TRUE)
# ADDED THIS TO SEE IF THIS WILL MAKE THE SIMULATIONS CONSISTENT ACROSS RUNS

if(length(args) != 4) {
    stop("Please supply 4 args")
}

#n_sample = 50
#frac_sig = 0.25

n_sim_per_chunk = 1

missing_val = -100
adapt_delta = 0.8
max_treedepth = 10
n_sample = as.numeric(args[1])
frac_sig = as.numeric(args[2])
sim_number = as.numeric(args[3])
MISSING_RATE = as.numeric(args[4])

max_missing = 0.4
cores=chains=4
n_iter=2500

CENSOR = ifelse(MISSING_RATE == 0, FALSE, TRUE)
print(paste0('CENSOR: ', CENSOR, ' TYPE ', class(CENSOR)))
CENSOR=as.logical(CENSOR)[1]

#OUT_FOLDER = '/20200305_2_stage_impute_beta_missing'
#OUT_FOLDER = '/20200315_2_stage_impute_beta_no_missing'
OUT_FOLDER = '/20200728_fixed_prior_no_missing'

######
# SET THE SEED SO THAT THE SIM NUMBER WILL BE THE SAME ACROSS MISSING_RATES
######
set.seed(sim_number)



#out_file = as.numeric(args[4])
#out_file = paste0( "SIM",n_sim,"_N_SAMPLE_",n_sample,"_FRAC_SIG_",frac_sig,'_MISSING_RATE_',MISSING_RATE,'.csv' )

#print(out_file)

out_file_raw = paste0("RAW_", "SIM_NUMBER_",sim_number,"_N_SAMPLE_",n_sample,"_FRAC_SIG_",frac_sig,'_MISSING_RATE_',MISSING_RATE,'.csv' )

print(out_file_raw)

path = file.path('/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/analysis', OUT_FOLDER)
dir.create(path)
setwd(path)

rstan_options(auto_write = TRUE)
stage_1_impute_code_file = 'stage_1_impute.stan'
stage_1_impute_model = stan_model(file = stage_1_impute_code_file, verbose = FALSE)

stage_2_regress_code_file = 'stage_2_regress.stan'
stage_2_regress_model = stan_model(file = stage_2_regress_code_file, verbose = FALSE)

df_nmr = read_csv('/nfs/turbo/umms-cgillies/cgillies/RACE_CAPS/11_20_2019_even_odds_interaction_prior_1/nmr_metabolites_scaled_matrix.csv')


died = df_nmr$died_90_day
sigma_0 = cov(df_nmr[died == 0,c(1:27)])
mu_0 = colMeans(df_nmr[died == 0,c(1:27)]) 
sigma_1 = cov(df_nmr[died == 1,c(1:27)])
mu_1 = colMeans(df_nmr[died == 1,c(1:27)]) 

metabolites = colnames(df_nmr)[1:27]
metabolites

# COMPUTE MISSING RATE
missing_rate = df_nmr[,metabolites] %>% summarise_all(function(x) { mean(x == min(x))  })
missing_rate_numeric = missing_rate %>% t() %>% as.numeric
mu = mean( missing_rate_numeric )
sigma_2 = var(missing_rate_numeric)





######
# Compute the true effect sizes
######
sim_data = make_df(10000,mu_0,sigma_0,10000,mu_1,sigma_1,frac_sig=1,censor=FALSE)
df_sample = sim_data$df_sample
truth = sim_data$truth
df_eff_sim = sapply(colnames(df_sample)[1:27],FUN=function( x ) { apply_logistic_sim(x,df_sample) } )
true_effs = as.numeric(df_eff_sim['Estimate',])



beta_params = get_alpha_beta_for_target_rate(df_nmr,metabolites,TARGET_MISSING_RATE = MISSING_RATE)
alpha = beta_params$alpha
beta = beta_params$beta

res = run_sim(censor=CENSOR, n_sim=n_sim_per_chunk,n_0=n_sample,n_1=n_sample, frac_sig=frac_sig,metabolites=metabolites,max_missing=max_missing,iter=n_iter ,sim_number=sim_number, MISSING_RATE=MISSING_RATE, alpha=alpha, beta=beta)


res$N_SIM = sim_number
res$N_SAMPLE = n_sample
res$FRAC_SIG = frac_sig
res$MISSING_RATE = MISSING_RATE

write_csv(res, out_file_raw)


