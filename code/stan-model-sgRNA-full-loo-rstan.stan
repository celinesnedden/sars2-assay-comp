//
// This Stan program defines a hurdle model for sgRNA values based on 
//   total RNA values, including both a logistic & a linear component 
// This uses the previously selected logistic model, and it is written for
//   model selection of the full model (i.e., inluding the linear part) 

data {
  
  // Number of observations
  int<lower=0> N; // number of observations in training data
  int<lower=0> N_lin; // number of observations for linear component
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
  int<lower=0, upper=1> sg_pos[N]; // logistic model
  real sg[N]; // linear model

  // All predictors
  real t[N]; // primary predictor, total RNA value
  real dose[N]; // log 10 pfu
  int dpi[N]; // day post infection
  int st[N]; // sample type 
  int sp[N]; // species 
  int tg[N]; // target gene
  int age[N]; // age class
  int sex[N]; // sex

  // Hierarchical error terms per study
  int article[N]; 
  
}

parameters {
  
  // Logistic regression model, chosen in prior x-validation step
  real gamma;
  real deltaT;
  real deltaDOSE;
  real deltaTG[L_tg];
  real deltaSP[L_sp];
  
  // Linear regression model
  real alpha;
  real<lower=0> betaT; 
  real betaDOSE_p[dose_inc ? 1 : 0];
  real betaDPI_p[dpi_inc ? L_dpi : 0];
  real betaST_p[st_inc ? 1: 0];
  real betaSP_p[sp_inc ? L_sp : 0];
  real betaTG_p[tg_inc ? L_tg : 0];
  real betaAGE_p[age_inc ? L_age : 0];
  real betaSEX_p[sex_inc ? 1 : 0];
  
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
  real betaSP[L_sp];
  real betaTG[L_tg];
  real betaDPI[L_dpi];
  real betaAGE[L_age];
  
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
  for (n in 1:N) {
      if (age[n] != -9) {
        // demographics known, normal handling. 
        target += bernoulli_logit_lpmf(sg_pos[n] | gamma + 
                                                   deltaT * t[n] +
                                                   deltaDOSE * dose[n] +
                                                   deltaTG[tg[n]] +
                                                   deltaSP[sp[n]] );
      }
      
      else if (age[n] == -9) {
        // age unknown, sex known. mixture required. 
        target += log_mix(p_adult, 
                    bernoulli_logit_lpmf(sg_pos[n] | gamma + 
                                                     deltaT * t[n] +
                                                     deltaDOSE * dose[n] +
                                                     deltaTG[tg[n]] +
                                                     deltaSP[sp[n]]),
                    bernoulli_logit_lpmf(sg_pos[n] | gamma + 
                                                     deltaT * t[n] +
                                                     deltaDOSE * dose[n] +
                                                     deltaTG[tg[n]] +
                                                     deltaSP[sp[n]] ));        
      }
  }
  
  // Linear regression model: only sgRNA+ samples contribute
  for (n in 1:N) {
    if (sg_pos[n] == 1 && sg[n] != -9) {
      if (age[n] != -9 && sex[n] != -9) {
        // demographics known, normal handling. 
        target += normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                      betaDPI[dpi[n]] + 
                                      betaDOSE * dose[n] +
                                      betaST * st[n] + 
                                      betaAGE[age[n]] + 
                                      betaSEX * sex[n] + 
                                      betaSP[sp[n]] + 
                                      betaTG[tg[n]],
                                      sigma[article[n]]);

      }
      else if (age[n] != -9 && sex[n] == -9) {
        // sex unknown, age known. mixture required. 
        target += log_mix(p_female, 
                          normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                              betaDPI[dpi[n]] + 
                                              betaDOSE * dose[n] +
                                              betaST * st[n] + 
                                              betaAGE[age[n]] +  
                                              betaSEX * 1 + 
                                              betaSP[sp[n]] + 
                                              betaTG[tg[n]],
                                              sigma[article[n]]),
                          normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                              betaDPI[dpi[n]] + 
                                              betaDOSE * dose[n] +
                                              betaST * st[n] + 
                                              betaAGE[age[n]] + 
                                              betaSEX * 0 +
                                              betaSP[sp[n]] + 
                                              betaTG[tg[n]], 
                                              sigma[article[n]]));
      }
      
      else if (age[n] == -9 && sex[n] != -9) {
        // age unknown, sex known. mixture required. 
        target += log_mix(p_adult, 
                          normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                              betaDPI[dpi[n]] + 
                                              betaDOSE * dose[n] +
                                              betaST * st[n] + 
                                              betaAGE[2] + 
                                              betaSEX * sex[n] + 
                                              betaSP[sp[n]] + 
                                              betaTG[tg[n]], 
                                              sigma[article[n]]),
                          normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                              betaDPI[dpi[n]] + 
                                              betaDOSE * dose[n] +
                                              betaST * st[n] + 
                                              betaAGE[1] +
                                              betaSEX * sex[n] +
                                              betaSP[sp[n]] + 
                                              betaTG[tg[n]], 
                                              sigma[article[n]]));
      }
      
      else if (age[n] == -9 && sex[n] == -9) {
        // age and sex unknown. large mixture required.
        lps[1] = log(p_adult) + 
                 log(p_female) + 
                 normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                     betaDPI[dpi[n]] +
                                     betaDOSE * dose[n] +
                                     betaST * st[n] + 
                                     betaAGE[2] +
                                     betaSEX * 1 + 
                                     betaSP[sp[n]] + 
                                     betaTG[tg[n]], 
                                     sigma[article[n]]);
        lps[2] = log1m(p_adult) + 
                 log(p_female) + 
                 normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                     betaDPI[dpi[n]] + 
                                     betaDOSE * dose[n] +
                                     betaST * st[n] + 
                                     betaAGE[1] +
                                     betaSEX * 1 + 
                                     betaSP[sp[n]] + 
                                     betaTG[tg[n]], 
                                     sigma[article[n]]);
        lps[3] = log(p_adult) + 
                 log1m(p_female) + 
                 normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                     betaDPI[dpi[n]] +
                                     betaDOSE * dose[n] +
                                     betaST * st[n] + 
                                     betaAGE[2] + 
                                     betaSEX * 0 + 
                                     betaSP[sp[n]] + 
                                     betaTG[tg[n]], 
                                     sigma[article[n]]);
        lps[4] = log1m(p_adult) + 
                 log1m(p_female) + 
                 normal_lpdf(sg[n] | alpha + betaT * t[n] + 
                                     betaDPI[dpi[n]] + 
                                     betaDOSE * dose[n] +
                                     betaST * st[n] + 
                                     betaAGE[1] +
                                     betaSEX * 0 + 
                                     betaSP[sp[n]] + 
                                     betaTG[tg[n]], 
                                     sigma[article[n]]);
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
  //// Logistic Component
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
  
  real log_lik[N + N_lin];
  int ii_lin = 1; // 
  
  for (ii in 1:N) {
    if (age[ii] != -9 && sex[ii] != -9) {
  
      // All points contribute to likelihood on logistic component
      log_lik[ii] = bernoulli_logit_lpmf(sg_pos[ii] | gamma + 
                                                      deltaT * t[ii] +
                                                      deltaTG[tg[ii]] +
                                                      deltaDOSE * dose[ii]  +
                                                      deltaSP[sp[ii]]);
      
      // Only positive values contribute to likelihood on linear component
      if (sg_pos[ii] == 1 && sg[ii] != -9) {
        log_lik[N + ii_lin] = normal_lpdf(sg[ii] | alpha + betaT * t[ii] + 
                                                   betaDPI[dpi[ii]] + 
                                                   betaDOSE * dose[ii] +
                                                   betaST * st[ii] + 
                                                   betaAGE[age[ii]] + 
                                                   betaSEX * sex[ii] + 
                                                   betaSP[sp[ii]] + 
                                                   betaTG[tg[ii]], 
                                                   sigma[article[ii]]);
       ii_lin = ii_lin + 1;                                                 
      }
    }
    
    else if (age[ii] == -9 && sex[ii] != -9) {

      // All points contribute to likelihood on logistic component
      log_lik[ii] = bernoulli_logit_lpmf(sg_pos[ii] | gamma + 
                                                      deltaT * t[ii] +
                                                      deltaDOSE * dose[ii] +
                                                      deltaTG[tg[ii]] +
                                                      deltaSP[sp[ii]]);        
      
      // Only positive values contribute to likelihood on linear component      
      if (sg_pos[ii] == 1 && sg[ii] != -9) {
        log_lik[N + ii_lin] = normal_lpdf(sg[ii] | alpha + betaT * t[ii] + 
                                                   betaDPI[dpi[ii]] + 
                                                   betaDOSE * dose[ii] +
                                                   betaST * st[ii] + 
                                                   betaAGE[bernoulli_rng(0.5) + 1] + 
                                                   betaSEX * sex[ii] + 
                                                   betaSP[sp[ii]] + 
                                                   betaTG[tg[ii]], 
                                                   sigma[article[ii]]); 
        ii_lin = ii_lin + 1;                                                           
      }
    }
    
    else if (age[ii] != -9 && sex[ii] == -9) {

      // All points contribute to likelihood on logistic component
      log_lik[ii] = bernoulli_logit_lpmf(sg_pos[ii] | gamma + 
                                                      deltaT * t[ii] +
                                                      deltaTG[tg[ii]] +
                                                      deltaDOSE * dose[ii] +
                                                      deltaSP[sp[ii]]);        
      
      // Only positive values contribute to likelihood on linear component         
      if (sg_pos[ii] == 1 && sg[ii] != -9) {
        log_lik[N + ii_lin] = normal_lpdf(sg[ii] | alpha + betaT * t[ii] + 
                                                   betaDPI[dpi[ii]] + 
                                                   betaDOSE * dose[ii] +
                                                   betaST * st[ii] + 
                                                   betaAGE[age[ii]] + 
                                                   betaSEX * bernoulli_rng(0.5) + 
                                                   betaSP[sp[ii]] + 
                                                   betaTG[tg[ii]],
                                                   sigma[article[ii]]);
        ii_lin = ii_lin + 1;                                                  
      }                                     
    }
    
    else if (age[ii] == -9 && sex[ii] == -9) {

      // All points contribute to likelihood on logistic component
      log_lik[ii] = bernoulli_logit_lpmf(sg_pos[ii] | gamma + 
                                                      deltaT * t[ii] +
                                                      deltaDOSE * dose[ii] +
                                                      deltaTG[tg[ii]] +
                                                      deltaSP[sp[ii]]);        
      // Only positive values contribute to likelihood on linear component          
      if (sg_pos[ii] == 1 && sg[ii] != -9) {
        log_lik[N + ii_lin] = normal_lpdf(sg[ii] | alpha + betaT * t[ii] + 
                                                   betaDPI[dpi[ii]] + 
                                                   betaDOSE * dose[ii] +
                                                   betaST * st[ii] + 
                                                   betaAGE[bernoulli_rng(0.5) + 1] + 
                                                   betaSEX * bernoulli_rng(0.5) + 
                                                   betaSP[sp[ii]] + 
                                                   betaTG[tg[ii]], 
                                                   sigma[article[ii]]);  
        ii_lin = ii_lin + 1;                                                
      }        
    }
  }
}
