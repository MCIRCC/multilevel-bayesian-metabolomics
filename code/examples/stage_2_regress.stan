data {
    int<lower=0> N;
    int<lower=0> N_var;
    int<lower=0> N_miss;

    matrix[ N, N_var ] x_raw;

    int<lower=0,upper=1> y[N];

    //Threshold for imputation
    //Separate threshold per metabolite
    vector[N_var] threshold;
    
    // This vector contains the missing row index, column index
    int x_missing[ N_miss, 2 ];
    real x_missing_mean[ N_miss ];
    real x_missing_sd[ N_miss ];
    real<lower=0> sigma_beta_x_lower_bound;
}
parameters {
    // Needed for imputation
    vector<upper=0>[ N_miss ] x_impute_raw;
    
    //regression on y per metabolite in x_raw columns
    vector[N_var] beta_x;
    vector[N_var] alpha_x;

    real<lower=sigma_beta_x_lower_bound> sigma_beta_x;
    real<lower=1> nu_x;

}
transformed parameters {
    // https://mc-stan.org/docs/2_18/stan-users-guide/vectors-with-varying-bounds.html
    vector<upper=0>[ N_miss ] x_impute;
    x_impute = x_impute_raw;
    for(i in 1:N_miss) {
        int metab_row = x_missing[i,1];
        int metab_col = x_missing[i,2];
        x_impute[i] = x_impute[i] + threshold[metab_col];
    }

}
model {
    matrix[N, N_var] x_merge;
    //copy matrix
    x_merge = x_raw;
    
    //fill in the values we are going to impute
    for(i in 1:N_miss) {
        int metab_row = x_missing[i,1];
        int metab_col = x_missing[i,2];
        x_merge[metab_row,metab_col] = x_impute[ i ];
    }
    
    x_impute ~ normal(x_missing_mean, x_missing_sd);
    
    sigma_beta_x ~ cauchy(0,1);
    nu_x ~ gamma(2,0.1);
    beta_x ~  student_t(nu_x,0,sigma_beta_x);
    alpha_x ~ normal(0,5);
    
    for(j in 1:N_var) {
        vector[N] x_col = x_merge[ , j];
        real beta = beta_x[j];
        real alpha = alpha_x[j];
        y ~ bernoulli_logit(alpha + beta * x_col);
    }
}
