missing_val = -9999.9999

get_alpha_beta_for_target_rate = function(df,metabolites,TARGET_MISSING_RATE) {
    missing_rate = df_nmr[,metabolites] %>% summarise_all(function(x) { mean(x == min(x))  })
    missing_rate_numeric = missing_rate %>% t() %>% as.numeric
    mu = mean( missing_rate_numeric )
    sigma_2 = var(missing_rate_numeric)
    RATE = TARGET_MISSING_RATE/mu
    mu = RATE * mu
    sigma_2 = RATE * sigma_2
    mu
    sigma_2
    nu = mu * (1 - mu) / sigma_2 - 1
    alpha = mu * nu
    beta = (1 - mu) * nu
    list(alpha=alpha,beta=beta,mu=mu,sigma_2=sigma_2)
}

get_alpha_beta_for_target_rate_2 = function(missing_rate_numeric,TARGET_MISSING_RATE) {
    #missing_rate = df_nmr[,metabolites] %>% summarise_all(function(x) { mean(x == min(x))  })
    #missing_rate_numeric = missing_rate %>% t() %>% as.numeric
    mu = mean( missing_rate_numeric )
    sigma_2 = var(missing_rate_numeric)
    RATE = TARGET_MISSING_RATE/mu
    mu = RATE * mu
    sigma_2 = RATE * sigma_2
    mu
    sigma_2
    nu = mu * (1 - mu) / sigma_2 - 1
    alpha = mu * nu
    beta = (1 - mu) * nu
    list(alpha=alpha,beta=beta,mu=mu,sigma_2=sigma_2)
}

make_summary_table_beta = function(fit,metabolites) {
    df_summary = data.frame(summary(fit)$summary)
    post = rstan::extract(fit)
    #compute prob > 0 or prob < 0
    n_sample = dim(post$beta_x)[1]
    g0 = colSums( (post$beta_x > 0) ) / n_sample 
    l0 = colSums( (post$beta_x < 0) ) / n_sample
    df_prob = data.frame(G0=g0,L0=l0) %>% rowwise() %>% mutate(Max=max(G0,L0))
    beta_names = c()
    for(i in seq(0,dim(post$beta_x)[2])) {
        beta_names=c(beta_names,paste0('beta_x[',i,']'))
    }
    df_summary_beta = df_summary %>% mutate(Param = rownames(.)) %>% 
    filter(Param %in% beta_names) %>% cbind(Metabolite = metabolites, P_GT_LT_0=df_prob$Max) %>%
    mutate(Z=mean/sd) %>% arrange(-P_GT_LT_0, -abs(Z)) %>% select(Metabolite,Param,mean,sd,X2.5.,X50.,X97.5.,Rhat,Z,P_GT_LT_0)
    return(df_summary_beta)   
}


stage_1_impute_get_init_vals = function(data) {
    n_var = data$N_var
    N_use_impute = data$N_use_impute
    # matrix[ N_use_impute, N_var ] beta_impute_raw;
    beta_impute_raw = matrix( rep(0, data$N_use_impute * n_var), ncol=n_var)
    alpha_impute = rep(0,n_var)
    x_use_impute = data$x_use_impute
    sigma_impute = rep(1, n_var)
    
    for(i in 1:n_var) {
        x_imp_vals_use = data$x_naive_impute[, x_use_impute[,i]   ]
        y_out_imp = data$x_naive_impute[,i]
        m_out = glm(y_out_imp ~ x_imp_vals_use,family="gaussian")
        coefs_imp = summary(m_out)$coef
        alpha_impute[i] = coefs_imp[1,1]
        for(k in 1:N_use_impute) {
            beta_impute_raw[k,i] = coefs_imp[1 + k,1]
        }
        sigma_impute[i] = sum(m_out$residuals^2) / m_out$df.residual
    }
    
    sigma_beta_impute_x = sqrt( sum( beta_impute_raw^2 ) / (prod(dim(beta_impute_raw)) - 1) ) 
    
    inits = list(
        sigma_impute = sigma_impute,
        sigma_beta_impute_x = sigma_beta_impute_x,
        alpha_impute = alpha_impute,
        beta_impute_raw = beta_impute_raw / sigma_beta_impute_x
    )
    
    inits
    
}


stage_1_impute_make_data = function(x_censored, y, threshold=NULL,
                                         N_use_impute=8) {
    
    if(is.null(threshold)) {
        threshold = apply(x_censored,MARGIN=2,FUN=function(x) { min(x,na.rm = T)})
    }
    
    naive_impute_val = log( exp(threshold) / 2)
    
    n = dim(x_censored)[1]
    N_var = dim(x_censored)[2]
    N_miss = sum(is.na(x_censored))
    
    x_naive_impute = x_censored %>% select_all()
    x_is_missing = apply( is.na(x_censored), FUN=as.numeric, MARGIN=c(1,2))             
    for(j in 1:N_var) {
        x_naive_impute[is.na(x_naive_impute[,j]),j] = naive_impute_val[j]
    }            
    # create a matrix that contains the N_use_impute indicies of the most correlated
    # variables with the metabolite of interest
    # each column contains a metabolite corresponding with x_censored
    # each row corresponds to the index of the metabolite most correlated with it up to N_use_impute
    # rows
    var_cors = cor(x_naive_impute,use = "complete")
    diag(var_cors) = 0
    get_inds_to_use_per_metabolite = function(x) { which( rank(-abs(x)) <= N_use_impute)   }
    x_use_impute = apply(var_cors,MARGIN=1,FUN=get_inds_to_use_per_metabolite )

    x_censored[is.na(x_censored)] = missing_val

    x_censored = as.matrix(x_censored)
    x_missing = matrix(nrow = N_miss, ncol = 2)
    x_missing_init = rep(0,N_miss)
    x_missing_init_raw = rep(0,N_miss)
    
    # make a table were the first column is the row number in x_censored
    # second column contains the metabolite column number in x_censored
    col1 = c()
    col2 = c()
    for(j in 1:dim(x_censored)[2]) {
        missing_in_col = which(x_censored[,j] == missing_val)
        col1 = c(col1, missing_in_col)
        col2 = c(col2, rep(j,length(missing_in_col)))
    }
    x_missing[,1] = col1
    x_missing[,2] = col2
    
    #print(x_missing)
    #print(threshold)
    # make a vector with naive imputation values
    if( N_miss > 0 ) {
        for(i in 1:dim(x_missing)[1]) {
            row = x_missing[i,1]
            col = x_missing[i,2]
            x_missing_init[i] = x_naive_impute[row,col]
            #print(row)
            #print(col)
            #print(x_missing_init[i])
            #print(threshold[col])
            x_missing_init_raw[i] = x_missing_init[i] - threshold[col]
        }
    }

    data = list(
            N = n,
            N_miss = N_miss,
            N_var = N_var,
            x_raw = x_censored,
            x_naive_impute = as.matrix(x_naive_impute),
            N_use_impute = N_use_impute,
            x_use_impute = x_use_impute,
            x_missing = x_missing,
            x_missing_init = x_missing_init,
            x_missing_init_raw = x_missing_init_raw,
            x_is_missing = x_is_missing,
            threshold = threshold,
            y=y
        )

    #str(data)
    #mean(x_censored == missing_val)
    #data$N_miss / (data$N * data$N_var)

    data   
}

replicate_init = function(init_data,chains=4) {
    init = list()
    for(i in 1:chains){
        init[paste0(i)] = list(init_data)
    }
    init
} 

default_impute = function(x_censored, x_missing) {
    N_miss = dim(x_missing)[1]
    #print(N_miss)
    x_missing_mean = rep(0,N_miss)
    x_missing_sd = rep(1,N_miss)
    
    list(
         x_missing_mean=x_missing_mean,
         x_missing_sd=x_missing_sd
        )
}

stage_2_regress_make_data = function(x_censored, y, IMPUTE_FUNC=default_impute, threshold=NULL, sigma_beta_x_lower_bound=0,...) {
    
    if(is.null(threshold)) {
        threshold = apply(x_censored,MARGIN=2,FUN=function(x) { min(x,na.rm = T)})
    }
    
    naive_impute_val = log( exp(threshold) / 2)
    
    n = dim(x_censored)[1]
    N_var = dim(x_censored)[2]
    N_miss = sum(is.na(x_censored))
    
    x_naive_impute = x_censored %>% select_all()
                 
    for(j in 1:N_var) {
        x_naive_impute[is.na(x_naive_impute[,j]),j] = naive_impute_val[j]
    }            


    x_censored[is.na(x_censored)] = missing_val

    x_censored = as.matrix(x_censored)
    x_missing = matrix(nrow = N_miss, ncol = 2)
    x_missing_init = rep(0,N_miss)
    x_missing_init_raw = rep(0,N_miss)
    x_missing_mean = rep(0,N_miss)
    x_missing_sd = rep(1,N_miss)
    
    # make a table were the first column is the row number in x_censored
    # second column contains the metabolite column number in x_censored
    col1 = c()
    col2 = c()
    for(j in 1:dim(x_censored)[2]) {
        missing_in_col = which(x_censored[,j] == missing_val)
        col1 = c(col1, missing_in_col)
        col2 = c(col2, rep(j,length(missing_in_col)))
    }
    x_missing[,1] = col1
    x_missing[,2] = col2
    
    # make a vector with naive imputation values
    if( N_miss > 0 ) {
        for(i in 1:dim(x_missing)[1]) {
            row = x_missing[i,1]
            col = x_missing[i,2]
            x_missing_init[i] = x_naive_impute[row,col]
            x_missing_init_raw[i] = x_missing_init[i] - threshold[col]
        }

        mean_sds = IMPUTE_FUNC(x_censored,x_missing,...)
        x_missing_mean = mean_sds$x_missing_mean
        x_missing_sd = mean_sds$x_missing_sd
    }
    
    data = list(
            N = n,
            N_miss = N_miss,
            N_var = N_var,
            x_raw = x_censored,
            x_naive_impute = as.matrix(x_naive_impute),
            x_missing = x_missing,
            x_missing_init = x_missing_init,
            x_missing_init_raw = x_missing_init_raw,
            x_missing_mean=x_missing_mean,
            x_missing_sd=x_missing_sd,
            threshold = threshold,
            y=y,
            sigma_beta_x_lower_bound=sigma_beta_x_lower_bound
        )

    #str(data)
    #mean(x_censored == missing_val)
    #data$N_miss / (data$N * data$N_var)

    data   
}

stage_1_impute = function(x_censored, x_missing, stage_1_post=NULL) {
    N_miss = dim(x_missing)[1]
    print(N_miss)
    x_missing_mean = apply( stage_1_post$x_impute_mean, FUN=mean, MARGIN=2 )
    x_missing_sd = apply( stage_1_post$x_impute, FUN=sd, MARGIN=2 )
            
    list(
         x_missing_mean=x_missing_mean,
         x_missing_sd=x_missing_sd
        )
}


run_bayes_model = function(df_censored,metabolites,cores=4,chains=4,iter=2000, 
                           adapt_delta=0.8,
                           max_treedepth=10, y_column='died_90_day', sigma_beta_x_lower_bound=0) {
    x_censored = df_censored[, metabolites]
    y = as.numeric( df_censored[,y_column] )
    
    stage_1_data = stage_1_impute_make_data(x_censored, y, N_use_impute=8)
    stage_1_init_data = replicate_init(
        stage_1_impute_get_init_vals(stage_1_data), chains=4
    )
    
    control = list(adapt_delta=adapt_delta,max_treedepth=max_treedepth)
    print("CONTROL")
    print(control)
    
    N_miss = stage_1_data$N_miss
    stage_1_fit = NULL
    stage_1_post = NULL
    if(N_miss > 0) {
        print("IMPUTING")
        stage_1_fit = sampling(stage_1_impute_model,  data = stage_1_data, init=stage_1_init_data, cores=cores, chains=chains, iter = iter,
                             control = control )

        stage_1_post = extract(stage_1_fit)
    } else {
        print("NO MISING DATA SKIPPING IMPUTATION")
    }
    
    stage_2_data = stage_2_regress_make_data(x_censored, y, IMPUTE_FUNC = stage_1_impute, stage_1_post=stage_1_post, sigma_beta_x_lower_bound=sigma_beta_x_lower_bound)
    stage_2_init = list(
        x_impute_raw = stage_2_data$x_missing_init_raw
    )
    
    stage_2_init = replicate_init(
        list(
            x_impute_raw = stage_2_data$x_missing_init_raw
        ), chains=4
    )
    
    
    stage_2_fit = sampling(stage_2_regress_model,  data = stage_2_data, init=stage_2_init, cores=cores, chains=chains, iter = iter, control = control )

    
    sum_table = make_summary_table_beta(stage_2_fit,metabolites)
    rownames(sum_table) = sum_table$Metabolite
    list( sum_table=sum_table, stage_1_fit=stage_1_fit, stage_2_fit=stage_2_fit)
}





make_df = function(n_0,mu_0,sigma_0,n_1,mu_1,sigma_1,frac_sig=0.5,censor=TRUE,max_missing=0.4,alpha=-1,beta=-1) {
    stopifnot(max_missing <= 0.4)
    
    diffs = abs( mu_1 - mu_0) / sqrt( (diag(sigma_0)^2)/n_0 + (diag(sigma_1)^2)/n_1  )
    #print(diffs)
    m = length(diffs)
    diffs = diffs/sum(diffs)
    #print(diffs)
    ind = seq(1,m)
    set_non_zero = sample(ind,size=ceiling(frac_sig * m),prob=diffs)
    set_zero = setdiff(ind,set_non_zero)
    truth = rep(FALSE,m)
    truth[set_non_zero] = TRUE
    mu_0[set_zero] = 0
    mu_1[set_zero] = 0
    X_0 = rmvnorm(n = n_0, mean = mu_0, sigma = sigma_0) 
    X_1 = rmvnorm(n = n_1, mean = mu_1, sigma = sigma_1)
    df_sample = data.frame(rbind(X_0,X_1))
    df_sample$died_90_day = c( rep(0,n_0),rep(1,n_1)  )
    
    df_censored = df_sample
    df_naive_impute = df_sample
    n_var = length(mu_0)
    thresholds = rep(0,n_var)
    est.thresholds = rep(0,n_var)
    est.naive_impute = rep(0,n_var)
    missing_rates=rep(0,n_var)
    if(censor == TRUE) {
        
        missing_rates = rbeta(n_var,shape1 = alpha, shape2 = beta)
        # only allow a missing rate up to a certain level
        missing_rates = sapply(missing_rates, function(x) { min(max_missing,x) })
        
        #missing_rates = rep(MISSING_RATE, n_var)
        
        for(i in 1:n_var) {
            thresholds[i] = quantile( df_censored[, i], missing_rates[i] )
            df_censored[ (df_censored[, i] < thresholds[i]), i ] = NA
            est.thresholds[i] = min( df_censored[ , i ], na.rm = T)
            est.naive_impute[i] = log( min( exp(df_censored[ , i ]), na.rm = T) / 2 )
            df_naive_impute[ (df_naive_impute[, i] < thresholds[i]), i ] = est.naive_impute[i]
        }
    }
    
    
    list(df_sample=df_sample,truth=truth, df_censored=df_censored, df_naive_impute=df_naive_impute, missing_rates=missing_rates, thresholds=thresholds, est.thresholds=est.thresholds, est.naive_impute=est.naive_impute)
}

apply_logistic_sim = function(met,df,y_column='died_90_day') {
    m = glm(I(df[,y_column]) ~ I(df[,met]) ,data=df,family="binomial")
    summary(m)$coef[2,]
}

compute_stats = function(sig, truth, label, est_effs, true_effs ) {
    tp = sum(   (sig == 1) & (truth == 1)  )
    fn = sum(  (sig == 0) & (truth == 1) )
    fp = sum(  (sig == 1) & (truth == 0)   )
    tn =  sum(  (sig == 0) & (truth == 0)  )
    avg_mag_error = mean( abs(est_effs[sig & truth]/true_effs[sig & truth]) )
    num_sign_error = sum( (sign(est_effs) != sign(true_effs))[sig & truth] )
    n_sig = sum(sig)
    n_true = sum(truth)
    
    if(n_sig == 0) {
        avg_mag_error = 0
        num_sign_error = 0
    }
    return(list(label=label,tp=tp,fn=fn,fp=fp,tn=tn,avg_mag_error=avg_mag_error, num_sign_error=num_sign_error, n_sig=n_sig, n_true=n_true))
}


run_sim = function(n_sim=100,n_0=100,n_1=100,frac_sig=0.4,
                   metabolites=NULL,censor=TRUE,max_missing=0.4,iter=2000,sim_number=1,MISSING_RATE=NULL,alpha=NULL,beta=NULL) {
    res = data.frame(stringsAsFactors = F)
    for(i in 1:n_sim) {
        print('ENTER_LOOP')
        sim_data = make_df(n_0,mu_0,sigma_0,n_1,mu_1,sigma_1,frac_sig=frac_sig,censor=censor,max_missing=max_missing,alpha = alpha, beta=beta)
        df_censored = sim_data$df_censored
        df_naive_impute = sim_data$df_naive_impute
        truth = sim_data$truth
        
        out_sim = paste0("RAW_", "SIM_NUMBER_",sim_number,"_N_SAMPLE_",n_sample,"_FRAC_SIG_",frac_sig,'_MISSING_RATE_',MISSING_RATE,'_SIM_DATA','.rds' )
        print("WRITNG OUT...")
        print(out_sim)
        write_rds(sim_data, out_sim, compress = 'gz')
        
        df_eff_sim = sapply(metabolites,FUN=function( x ) { apply_logistic_sim(x,df_naive_impute) } )
        sig_raw = as.vector( df_eff_sim['Pr(>|z|)',] < 0.05)
        eff_raw = as.vector( df_eff_sim['Estimate',])
        
        sig_bon = as.vector( p.adjust( df_eff_sim['Pr(>|z|)',],method='bonferroni') < 0.05)
        eff_bon = as.vector( df_eff_sim['Estimate',])
        
        sig_bh = as.vector( p.adjust( df_eff_sim['Pr(>|z|)',],method='BH') < 0.05)
        eff_bh = as.vector( df_eff_sim['Estimate',])
        
        bayes_model_res = run_bayes_model(df_censored,metabolites,cores=cores,chains=chains,iter=iter, 
                           adapt_delta=adapt_delta,
                           max_treedepth=max_treedepth)
        
        
        sum_table = bayes_model_res$sum_table
        stage_1_fit = bayes_model_res$stage_1_fit
        stage_2_fit = bayes_model_res$stage_2_fit
        
        out_fit = paste0("RAW_", "SIM_NUMBER_",sim_number,"_N_SAMPLE_",n_sample,"_FRAC_SIG_",frac_sig,'_MISSING_RATE_',MISSING_RATE,'_FIT_DATA','.rds' )
        print("WRITNG OUT...")
        print(out_fit)
        write_rds(bayes_model_res, out_fit, compress = 'gz')

        sig_bayes = as.vector( sum_table[metabolites,'P_GT_LT_0'] > 0.975  )
        eff_bayes = as.vector( sum_table[metabolites,'mean'] )
        
        res_bayes = compute_stats(sig_bayes,truth,'bayes', eff_bayes, true_effs )
        res = rbind(res,res_bayes,stringsAsFactors = F)
        
        res_raw = compute_stats(sig_raw,truth,'raw', eff_raw, true_effs )
        res = rbind(res,res_raw,stringsAsFactors = F)
        res_bon = compute_stats(sig_bon,truth,'bon', eff_bon, true_effs )
        res = rbind(res,res_bon,stringsAsFactors = F)
        res_bh = compute_stats(sig_bh,truth,'bh', eff_bh, true_effs )
        res = rbind(res,res_bh,stringsAsFactors = F)
        
    }
    res
}


