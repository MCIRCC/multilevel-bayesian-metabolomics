functions {
    real normal_ub_rng(real mu, real sigma, real ub) {
      real p_ub = normal_cdf(ub, mu, sigma);
      real u = uniform_rng(0, p_ub);
      real y = mu + sigma * inv_Phi(u);
      return y;
    }
}
data {
    int<lower=0> N;
    int<lower=0> N_var;
    int<lower=0> N_miss;
    int<lower=0> N_use_impute;

    matrix[ N, N_var ] x_naive_impute;
    matrix[ N, N_var ] x_is_missing;

    //contains the indicies of x_naive_imp to use for imputation
    int x_use_impute[N_use_impute, N_var ];

    //Threshold for imputation
    //Separate threshold per metabolite
    vector[N_var] threshold;
    

    // This vector contains the missing row index, column index
    int x_missing[ N_miss, 2 ];
    
}
parameters {

    //betas for each imputed value
    matrix[ N_use_impute, N_var ] beta_impute_raw;
    vector[ N_var ] alpha_impute;
    vector<lower=0>[ N_var ] sigma_impute;
    real<lower=0> sigma_beta_impute_x;

}
transformed parameters {
    matrix[ N_use_impute, N_var ] beta_impute;
    for(j in 1:N_var) {
        beta_impute[,j] = beta_impute_raw[,j] * sigma_beta_impute_x;
    }
}
model {
    sigma_beta_impute_x ~ std_normal();
    // FORGOT TO SET PRIOR ON alpha_x_pred and sigma_impute

    for(j in 1:N_var) {
        real y[N];
        vector[N_use_impute] beta_x_pred;
        int x_use_ind[N_use_impute];
        matrix[N, N_use_impute] x_pred_use;
        real alpha_x_pred;
        real sigma_x_pred;
        //predicted target
        vector[N] xj_pred;
        //target
        vector[N] xj;
        int observed[N];
        
        //ADD TO LIKELIHOOD
        beta_impute_raw[, j] ~ std_normal();
    
        //ASSIGN VALS
        xj = x_naive_impute[, j];
        //dont make a prediction using itself
        //this should be copying the values according to ref
        beta_x_pred = beta_impute[, j];
        alpha_x_pred = alpha_impute[j];
        x_use_ind = x_use_impute[, j];
        sigma_x_pred = sigma_impute[j];
        x_pred_use = x_naive_impute[, x_use_ind];
        
        //Make pred
        xj_pred = x_pred_use * beta_x_pred + alpha_x_pred;
        
        
        //observed = (x_is_missing[, j] == 0);
        
        // target +=  normal_lpdf(xj | xj_pred, sigma_x_pred);

        for(i in 1:N) {
            real mu_x_pred;
            mu_x_pred = xj_pred[i];
            if( x_is_missing[i,j] == 0) {
                // if observed use pdf for observed value
                target +=  normal_lpdf(xj[i] | mu_x_pred, sigma_x_pred);
            } else {
                // if censored, use normal CDF
                target +=  normal_lcdf(threshold[j] | mu_x_pred, sigma_x_pred);
            }
        }
    }
}
generated quantities {
    real x_impute_mean[ N_miss ];
    real x_impute[ N_miss ];

    for(i in 1:N_miss) {
        //These correspond to the rows and columns in the original data frame
        int metab_row;
        int metab_col;
        row_vector[N_use_impute] x_pred_use;
        vector[N_use_impute] beta_x_pred;
        real alpha_x_pred;
        int x_use_ind[N_use_impute];
        real mu_x;
        real sd_x;
        real ub_x;

        metab_row = x_missing[i,1];
        metab_col = x_missing[i,2];
        x_use_ind = x_use_impute[, metab_col ];
        beta_x_pred = beta_impute[, metab_col ];
        alpha_x_pred = alpha_impute[ metab_col ];
        x_pred_use = x_naive_impute[metab_row, x_use_ind];

        mu_x = x_pred_use * beta_x_pred + alpha_x_pred;
        sd_x =  sigma_impute[ metab_col ];
        ub_x = threshold[ metab_col ];
        x_impute_mean[ i ] = mu_x;
        x_impute[i] = normal_ub_rng(mu_x, sd_x, ub_x);
    }
}
