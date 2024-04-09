// This Stan program defines a hurdle model for sgRNA values based on 
//   total RNA values, including both a logistic & a linear component 
// This uses the previously selected logistic model, and it is written for
//   model selection of the full model (i.e., inluding the linear part) 

data {
  
  // Number of observations
  int<lower=0> N_train; // number of observations in training data
  int<lower=0> N_test; // number of observations in test data
  int<lower=0> N_test_lin; // number of observations for linear component
  int<lower=0> N_article; // for study-level hierarchy
  
  // Indicators for which predictors to consider
  int<lower=0, upper=1> dose_inc;
  int<lower=0, upper=1> dpi_inc;
  int<lower=0, upper=1> st_inc;
  int<lower=0, upper=1> sp_inc;
  int<lower=0, upper=1> tg_inc;
  int<lower=0, upper=1> age_inc;
  int<lower=0, upper=1> sex_inc;
  
  // Levels for categorical predictors
  int<lower=0> L_dpi; // day post infection
  int<lower=0> L_st; // sample type
  int<lower=0> L_sp; // species
  int<lower=0> L_tg; // target gene
  int<lower=0> L_age; // age class

  // Outcome variables
  array[N_train] int<lower=0, upper=1> sg_pos_train; // logistic model
  array[N_test] int<lower=0, upper=1> sg_pos_test; 
  array[N_train] real sg_train; // linear model
  array[N_test] real sg_test; 

  // All predictors
  array[N_train] real t_train; // primary predictor, total RNA value
  array[N_test] real t_test; 
  array[N_train] real dose_train; // log 10 pfu
  array[N_test] real dose_test; 
  array[N_train] int dpi_train; // day post infection
  array[N_test] int dpi_test;
  array[N_train] int st_train; // sample type 
  array[N_test] int st_test;
  array[N_train] int sp_train; // species 
  array[N_test] int sp_test;
  array[N_train] int tg_train; // target gene
  array[N_test] int tg_test;
  array[N_train] int age_train; // age class
  array[N_test] int age_test;
  array[N_train] int sex_train; // sex
  array[N_test] int sex_test;

  // Hierarchical error terms per study
  array[N_train] int article_train; 
  array[N_test] int article_test;
  
}

parameters {
  
  // Logistic regression model, chosen in prior x-validation step
  real gamma;
  real deltaT;
  real deltaDOSE;
  array[L_tg] real deltaTG;
  array[L_sp] real deltaSP;
  
  // Linear regression model, written to toggle predictors on/off
  real alpha;
  real<lower=0> betaT; 
  array[dose_inc ? 1 : 0] real betaDOSE_p;
  array[dpi_inc ? L_dpi : 0] real betaDPI_p;
  array[st_inc ? 1: 0] real betaST_p;
  array[sp_inc ? L_sp : 0] real betaSP_p;
  array[tg_inc ? L_tg : 0] real betaTG_p;
  array[age_inc ? L_age : 0] real betaAGE_p;
  array[sex_inc ? 1 : 0] real betaSEX_p;
  
  // Hierarchical errors for linear regression
  vector<lower=0>[N_article] sigma; // study-specific error 
  real<lower=0> sigma_bar; // average sigma
  real<lower=0> sigma_sd; // standard deviation of distribution of sigmas
  
  // Sex & age imputation
  real<lower=0, upper=1> p_adult;
  real<lower=0, upper=1> p_female;

}

transformed parameters {
  
  // To accomodate toggling predictors on/off
  real betaDOSE;
  real betaSEX;
  real betaST;
  array[L_sp] real betaSP;
  array[L_tg] real betaTG;
  array[L_dpi] real betaDPI;
  array[L_age] real betaAGE;
  
  // Continous variables
  if (dose_inc == 1) {
    betaDOSE = betaDOSE_p[1];
  }
  else if (dose_inc == 0) {
    betaDOSE = 0;
  }
  
  if (sex_inc == 1) {
    betaSEX = betaSEX_p[1];
  }
  else if (sex_inc == 0) {
    betaSEX = 0;
  }
  
  if (st_inc == 1) {
    betaST = betaST_p[1];
  }
  else if (st_inc == 0) {
    betaST = 0;
  }
  
  // Categorical variables
  if (sp_inc == 1) {
    betaSP = betaSP_p;
  }
  else if (sp_inc == 0) {
    betaSP = rep_array(0.0, L_sp);
  }
  
  if (dpi_inc == 1) {
    betaDPI = betaDPI_p;
  }
  else if (dpi_inc == 0) {
    betaDPI = rep_array(0.0, L_dpi);
  }
  
  if (tg_inc == 1) {
    betaTG = betaTG_p;
  }
  else if (tg_inc == 0) {
    betaTG = rep_array(0.0, L_tg);
  }
  
  if (age_inc == 1) {
    betaAGE = betaAGE_p;
  }
  else if (age_inc == 0) {
    betaAGE = rep_array(0.0, L_age);
  }
}


model {
  vector[4] lps; // for convenience during marginalization below

  // Logistic regression model: everything contributes
  for (n in 1:N_train) {
      if (age_train[n] != -9) {
        // demographics known, normal handling.
        target += bernoulli_logit_lpmf(sg_pos_train[n] | gamma + 
                                                         deltaT * t_train[n] +
                                                         deltaDOSE * dose_train[n] +
                                                         deltaTG[tg_train[n]] +
                                                         deltaSP[sp_train[n]]);
      }
      
      else if (age_train[n] == -9) {
        // age unknown, sex known. mixture required. update p_female.
        target += log_mix(p_adult, 
                    bernoulli_logit_lpmf(sg_pos_train[n] | gamma + 
                                                           deltaT * t_train[n] +
                                                           deltaDOSE * dose_train[n]+
                                                           deltaTG[tg_train[n]] +
                                                           deltaSP[sp_train[n]]),
                    bernoulli_logit_lpmf(sg_pos_train[n] | gamma + 
                                                           deltaT * t_train[n] +
                                                           deltaDOSE * dose_train[n]+
                                                           deltaTG[tg_train[n]]) +
                                                           deltaSP[sp_train[n]]);        
      }
  }
  
  // Linear regression model: only sgRNA+ samples contribute
  for (n in 1:N_train) {
    if (sg_pos_train[n] == 1 && sg_train[n] != -9) {
      if (age_train[n] != -9 && sex_train[n] != -9) {
        // demographics known, normal handling.
        
        target += normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                            betaDPI[dpi_train[n]] + 
                                            betaDOSE * dose_train[n] +
                                            betaST * st_train[n] + 
                                            betaAGE[age_train[n]] + 
                                            betaSEX * sex_train[n] + 
                                            betaSP[sp_train[n]] + 
                                            betaTG[tg_train[n]],
                                            sigma[article_train[n]]);

      }
      else if (age_train[n] != -9 && sex_train[n] == -9) {
        // sex unknown, age known. mixture required. 
        target += log_mix(p_female, 
                          normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                                    betaDPI[dpi_train[n]] + 
                                                    betaDOSE * dose_train[n] +
                                                    betaST * st_train[n] + 
                                                    betaAGE[age_train[n]] + 
                                                    betaSEX * 1 + 
                                                    betaSP[sp_train[n]] + 
                                                    betaTG[tg_train[n]], 
                                                    sigma[article_train[n]]),
                          normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                                    betaDPI[dpi_train[n]] + 
                                                    betaDOSE * dose_train[n] +
                                                    betaST * st_train[n] + 
                                                    betaAGE[age_train[n]] + 
                                                    betaSEX * 0 +
                                                    betaSP[sp_train[n]] + 
                                                    betaTG[tg_train[n]], 
                                                    sigma[article_train[n]]));
      }
      
      else if (age_train[n] == -9 && sex_train[n] != -9) {
        // age unknown, sex known. mixture required.
        target += log_mix(p_adult, 
                          normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                                    betaDPI[dpi_train[n]] + 
                                                    betaDOSE * dose_train[n] +
                                                    betaST * st_train[n] + 
                                                    betaAGE[2] + 
                                                    betaSEX * sex_train[n] + 
                                                    betaSP[sp_train[n]] + 
                                                    betaTG[tg_train[n]], 
                                                    sigma[article_train[n]]),
                          normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                                    betaDPI[dpi_train[n]] + 
                                                    betaDOSE * dose_train[n] +
                                                    betaST * st_train[n] + 
                                                    betaAGE[1] + 
                                                    betaSEX * sex_train[n] +
                                                    betaSP[sp_train[n]] + 
                                                    betaTG[tg_train[n]], 
                                                    sigma[article_train[n]]));
      }
      
      else if (age_train[n] == -9 && sex_train[n] == -9) {
        // age and sex unknown. large mixture required.
        lps[1] = log(p_adult) + 
                 log(p_female) + 
                 normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                           betaDPI[dpi_train[n]] +
                                           betaDOSE * dose_train[n] +
                                           betaST * st_train[n] + 
                                           betaAGE[2] + 
                                           betaSEX * 1 + 
                                           betaSP[sp_train[n]] + 
                                           betaTG[tg_train[n]], 
                                           sigma[article_train[n]]);
        lps[2] = log1m(p_adult) + 
                 log(p_female) + 
                 normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                           betaDPI[dpi_train[n]] + 
                                           betaDOSE * dose_train[n] +
                                           betaST * st_train[n] + 
                                           betaAGE[1] + 
                                           betaSEX * 1 + 
                                           betaSP[sp_train[n]] + 
                                           betaTG[tg_train[n]], 
                                           sigma[article_train[n]]);
        lps[3] = log(p_adult) + 
                 log1m(p_female) + 
                 normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                           betaDPI[dpi_train[n]] +
                                           betaDOSE * dose_train[n] +
                                           betaST * st_train[n] + 
                                           betaAGE[2] + 
                                           betaSEX * 0 + 
                                           betaSP[sp_train[n]] + 
                                           betaTG[tg_train[n]], 
                                           sigma[article_train[n]]);
        lps[4] = log1m(p_adult) + 
                 log1m(p_female) + 
                 normal_lpdf(sg_train[n] | alpha + betaT * t_train[n] + 
                                           betaDPI[dpi_train[n]] + 
                                           betaDOSE * dose_train[n] +
                                           betaST * st_train[n] + 
                                           betaAGE[1] + 
                                           betaSEX * 0 + 
                                           betaSP[sp_train[n]] + 
                                           betaTG[tg_train[n]], 
                                           sigma[article_train[n]]);
        target += log_sum_exp(lps);                                                      
      }
    }
  }
  
  // Priors: Non-Informative
  // Linear Component
  alpha ~ normal(-2, 1);
  betaT ~ gamma(2, 0.5);
  betaDOSE ~ normal(0, 1);
  betaDPI ~ normal(0, 1);
  betaAGE ~ normal(0, 1);
  betaSEX ~ normal(0, 1);
  betaST ~ normal(0, 1);
  betaSP ~ normal(0, 1);
  betaTG ~ normal(0, 1);
  sigma ~ normal(sigma_bar, sigma_sd);
  sigma_bar ~ normal(0, 1);
  sigma_sd ~ exponential(1);
  
  // Logistic  Component
  gamma ~ normal(0, 1);
  deltaT ~ normal(0, 1);
  deltaSP ~ normal(0, 1);
  deltaTG ~ normal(0, 1);
  deltaDOSE ~ normal(0, 1);

  // Priors: Informative
  // Linear Component
  //alpha ~ normal(-3, 1);
  //betaT ~ gamma(2, 0.5);
  //betaDOSE ~ normal(-0.25, 1);
  //betaDPI[1] ~ normal(-0.5, 1);
  //betaDPI[2] ~ normal(0.5, 1);
  //betaDPI[3] ~ normal(0, 1);
  //betaAGE ~ normal(0, 1);
  //betaSEX ~ normal(0, 1);
  //betaST ~ normal(0, 1);
  //betaSP ~ normal(0, 1);
  //betaTG[1] ~ normal(0.5, 1);
  //betaTG[2] ~ normal(0.5, 1);
  //betaTG[3] ~ normal(-0.5, 1);
  //betaTG[4] ~ normal(-0.5, 1);
  //sigma ~ normal(sigma_bar, sigma_sd);
  //sigma_bar ~ normal(0, 1);
  //sigma_sd ~ exponential(1);
  ////
  ////// Logistic Component
  //gamma ~ normal(-1, 2);
  //deltaT ~ normal(2, 1);
  //deltaTG[1] ~ normal(1, 1);
  //deltaTG[2] ~ normal(1, 1);
  //deltaTG[3] ~ normal(-1, 1);
  //deltaTG[4] ~ normal(-1, 1);
  //deltaDOSE ~ normal(-1, 0.5);
  //deltaSP ~ normal(0, 1);

}

generated quantities {
  
  // TEST DATA 
  
  // Log likelihoods for performance statistics
  array[N_test] real log_lik_log;
  array[N_test_lin] real log_lik_lin;
  int ii_lin = 1; // counter for indexing log_lik_lin
  
  // Predicted censoring & values
  array[N_test] real p_pos_test; 
  array[N_test] real pred_val_test; 
  
  for (ii in 1:N_test) {
    if (age_test[ii] != -9 && sex_test[ii] != -9) {
      
      // Posterior predictions on censoring & values
      p_pos_test[ii] = bernoulli_logit_rng(gamma + 
                                           deltaT * t_test[ii] +
                                           deltaDOSE * dose_test[ii]  +
                                           deltaTG[tg_test[ii]] +
                                           deltaSP[sp_test[ii]]);
      pred_val_test[ii] = normal_rng(alpha + betaT * t_test[ii] + 
                                     betaDPI[dpi_test[ii]] + 
                                     betaDOSE * dose_test[ii] +
                                     betaST * st_test[ii] + 
                                     betaAGE[age_test[ii]] + 
                                     betaSEX * sex_test[ii] + 
                                     betaSP[sp_test[ii]] + 
                                     betaTG[tg_test[ii]],
                                     sigma[discrete_range_rng(1, N_article)]);
      
      // Log likelihoods
      // All points contribute to likelihood on logistic component
      log_lik_log[ii] = bernoulli_logit_lpmf(sg_pos_test[ii] | gamma + 
                                                               deltaT * t_test[ii] +
                                                               deltaDOSE * dose_test[ii]  +
                                                               deltaTG[tg_test[ii]] +
                                                               deltaSP[sp_test[ii]]);
      
      // Only positive values contribute to likelihood on linear component
      if (sg_pos_test[ii] == 1 && sg_test[ii] != -9) {
        log_lik_lin[ii_lin] = normal_lpdf(sg_test[ii] | alpha + betaT * t_test[ii] + 
                                                        betaDPI[dpi_test[ii]] + 
                                                        betaDOSE * dose_test[ii] +
                                                        betaST * st_test[ii] + 
                                                        betaAGE[age_test[ii]] + 
                                                        betaSEX * sex_test[ii] + 
                                                        betaSP[sp_test[ii]] + 
                                                        betaTG[tg_test[ii]], 
                                                        sigma[article_test[ii]]);
       ii_lin = ii_lin + 1;                                                 
      }
    }
    
    else if (age_test[ii] == -9 && sex_test[ii] != -9) {
      
      // Posterior predictions on censoring & values     
      p_pos_test[ii] = bernoulli_logit_rng(gamma + 
                                           deltaT * t_test[ii] +
                                           deltaDOSE * dose_test[ii]  +
                                           deltaTG[tg_test[ii]] +
                                           deltaSP[sp_test[ii]]);
                                            
      pred_val_test[ii] = normal_rng(alpha + betaT * t_test[ii] + 
                                     betaDPI[dpi_test[ii]] + 
                                     betaDOSE * dose_test[ii] +
                                     betaST * st_test[ii] + 
                                     betaAGE[bernoulli_rng(0.5) + 1] + 
                                     betaSEX * sex_test[ii] + 
                                     betaSP[sp_test[ii]] + 
                                     betaTG[tg_test[ii]], 
                                     sigma[discrete_range_rng(1, 13)]);
                                     
      // Log likelihoods
      // All points contribute to likelihood on logistic component
      log_lik_log[ii] = bernoulli_logit_lpmf(sg_pos_test[ii] | gamma + 
                                                               deltaT * t_test[ii] +
                                                               deltaDOSE * dose_test[ii]  +
                                                               deltaTG[tg_test[ii]] +
                                                               deltaSP[sp_test[ii]]);        
      
      // Only positive values contribute to likelihood on linear component      
      if (sg_pos_test[ii] == 1 && sg_test[ii] != -9) {

        log_lik_lin[ii_lin] = normal_lpdf(sg_test[ii] | alpha + betaT * t_test[ii] + 
                                                        betaDPI[dpi_test[ii]] + 
                                                        betaDOSE * dose_test[ii] +
                                                        betaST * st_test[ii] + 
                                                        betaAGE[bernoulli_rng(0.5) + 1] + 
                                                        betaSEX * sex_test[ii] + 
                                                        betaSP[sp_test[ii]] + 
                                                        betaTG[tg_test[ii]], 
                                                        sigma[article_test[ii]]);
        ii_lin = ii_lin + 1;                                                           
      }
    }
    
    else if (age_test[ii] != -9 && sex_test[ii] == -9) {
      
      // Posterior predictions on censoring & values       
      p_pos_test[ii] = bernoulli_logit_rng(gamma + 
                                           deltaT * t_test[ii] +
                                           deltaDOSE * dose_test[ii]  +
                                           deltaTG[tg_test[ii]] +
                                           deltaSP[sp_test[ii]]);
                                            
      pred_val_test[ii] = normal_rng(alpha + betaT * t_test[ii] + 
                                     betaDPI[dpi_test[ii]] + 
                                     betaDOSE * dose_test[ii] +
                                     betaST * st_test[ii] + 
                                     betaAGE[age_test[ii]] + 
                                     betaSEX * bernoulli_rng(0.5) + 
                                     betaSP[sp_test[ii]] + 
                                     betaTG[tg_test[ii]], 
                                     sigma[discrete_range_rng(1, 13)]);                              
  
      // Log likelihoods
      // All points contribute to likelihood on logistic component
      log_lik_log[ii] = bernoulli_logit_lpmf(sg_pos_test[ii] | gamma + 
                                                               deltaT * t_test[ii] +
                                                               deltaDOSE * dose_test[ii]  +
                                                               deltaTG[tg_test[ii]] +
                                                               deltaSP[sp_test[ii]]);        
      
      // Only positive values contribute to likelihood on linear component         
      if (sg_pos_test[ii] == 1 && sg_test[ii] != -9) {
        log_lik_lin[ii_lin] = normal_lpdf(sg_test[ii] | alpha + betaT * t_test[ii] + 
                                                        betaDPI[dpi_test[ii]] + 
                                                        betaDOSE * dose_test[ii] +
                                                        betaST * st_test[ii] + 
                                                        betaAGE[age_test[ii]]+ 
                                                        betaSEX * bernoulli_rng(0.5) +  
                                                        betaSP[sp_test[ii]] + 
                                                        betaTG[tg_test[ii]],
                                                        sigma[article_test[ii]]);
        ii_lin = ii_lin + 1;                                                  
      }                                     
    }
    
    else if (age_test[ii] == -9 && sex_test[ii] == -9) {
      
      // Posterior predictions on censoring & values    
      p_pos_test[ii] = bernoulli_logit_rng(gamma + 
                                           deltaT * t_test[ii] +
                                           deltaDOSE * dose_test[ii] +
                                           deltaTG[tg_test[ii]] +
                                           deltaSP[sp_test[ii]]);
                                            
      pred_val_test[ii] = normal_rng(alpha + betaT * t_test[ii] + 
                                     betaDPI[dpi_test[ii]] + 
                                     betaDOSE * dose_test[ii] +
                                     betaST * st_test[ii] + 
                                     betaAGE[bernoulli_rng(0.5) + 1] + 
                                     betaSEX * bernoulli_rng(0.5) + 
                                     betaSP[sp_test[ii]] + 
                                     betaTG[tg_test[ii]], 
                                     sigma[discrete_range_rng(1, 13)]);
                                     
      // Log likelihoods
      // All points contribute to likelihood on logistic component
      log_lik_log[ii] = bernoulli_logit_lpmf(sg_pos_test[ii] | gamma + 
                                                               deltaT * t_test[ii] +
                                                               deltaDOSE * dose_test[ii]  +
                                                               deltaTG[tg_test[ii]] +
                                                               deltaSP[sp_test[ii]]);        
      // Only positive values contribute to likelihood on linear component          
      if (sg_pos_test[ii] == 1 && sg_test[ii] != -9) {
        log_lik_lin[ii_lin] = normal_lpdf(sg_test[ii] | alpha + betaT * t_test[ii] + 
                                                        betaDPI[dpi_test[ii]] + 
                                                        betaDOSE * dose_test[ii] +
                                                        betaST * st_test[ii] + 
                                                        betaAGE[bernoulli_rng(0.5) + 1] + 
                                                        betaSEX * bernoulli_rng(0.5) + 
                                                        betaSP[sp_test[ii]] + 
                                                        betaTG[tg_test[ii]], 
                                                        sigma[article_test[ii]]); 
        ii_lin = ii_lin + 1;                                                
      }        
    }
  }


  // TRAINING DATA
  
  // Posterior predictions of training data for performance statistics
  array[N_train] real p_pos_train;
  array[N_train] real pred_val_train;
  
  for (ii in 1:N_train) {
      if (age_train[ii] != -9 && sex_train[ii] != -9) {
        p_pos_train[ii] = bernoulli_logit_rng(gamma + 
                                              deltaT * t_train[ii] +
                                              deltaDOSE * dose_train[ii] +
                                              deltaTG[tg_train[ii]] +
                                              deltaSP[sp_train[ii]]);
        pred_val_train[ii] = normal_rng(alpha + betaT * t_train[ii] + 
                                        betaDPI[dpi_train[ii]] + 
                                        betaDOSE * dose_train[ii] +
                                        betaST * st_train[ii] + 
                                        betaAGE[age_train[ii]] + 
                                        betaSEX * sex_train[ii] + 
                                        betaSP[sp_train[ii]] + 
                                        betaTG[tg_train[ii]],
                                        sigma[discrete_range_rng(1, 13)]);
      }
     
     else if (age_train[ii] == -9 && sex_train[ii] != -9) {
        p_pos_train[ii] = bernoulli_logit_rng(gamma + 
                                              deltaT * t_train[ii] +
                                              deltaDOSE * dose_train[ii] +
                                              deltaTG[tg_train[ii]] + 
                                              deltaSP[sp_train[ii]]);
        pred_val_train[ii] = normal_rng(alpha + betaT * t_train[ii] + 
                                        betaDPI[dpi_train[ii]] +
                                        betaDOSE * dose_train[ii] +
                                        betaST * st_train[ii] + 
                                        betaAGE[bernoulli_rng(0.5) + 1] +
                                        betaSEX * sex_train[ii] + 
                                        betaSP[sp_train[ii]] + 
                                        betaTG[tg_train[ii]],
                                        sigma[discrete_range_rng(1, 13)]);                               
                                        
     }
      
     else if (age_train[ii] != -9 && sex_train[ii] == -9) {
        p_pos_train[ii] = bernoulli_logit_rng(gamma + 
                                              deltaT * t_train[ii] +
                                              deltaDOSE * dose_train[ii] + 
                                              deltaTG[tg_train[ii]] +
                                              deltaSP[sp_train[ii]]);
        pred_val_train[ii] = normal_rng(alpha + betaT * t_train[ii] + 
                                        betaDPI[dpi_train[ii]] + 
                                        betaDOSE * dose_train[ii] +
                                        betaST * st_train[ii] + 
                                        betaAGE[age_train[ii]] + 
                                        betaSEX * bernoulli_rng(p_female) + 
                                        betaSP[sp_train[ii]] + 
                                        betaTG[tg_train[ii]],
                                        sigma[discrete_range_rng(1, 13)]);                             
      }
      
     else if (age_train[ii] == -9 && sex_train[ii] == -9) {
        p_pos_train[ii] = bernoulli_logit_rng(gamma + 
                                              deltaT * t_train[ii] +
                                              deltaDOSE * dose_train[ii] + 
                                              deltaTG[tg_train[ii]] + 
                                              deltaSP[sp_train[ii]]);
        pred_val_train[ii] = normal_rng(alpha + betaT * t_train[ii] + 
                                        betaDPI[dpi_train[ii]] + 
                                        betaDOSE * dose_train[ii] +
                                        betaST * st_train[ii] + 
                                        betaAGE[bernoulli_rng(0.5) + 1] + 
                                        betaSEX * bernoulli_rng(p_female) + 
                                        betaSP[sp_train[ii]] + 
                                        betaTG[tg_train[ii]], 
                                        sigma[discrete_range_rng(1, 13)]);                             
        }                                                
    }
}
