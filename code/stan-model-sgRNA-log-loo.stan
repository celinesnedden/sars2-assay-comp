//
// This Stan program defines a logistic regression for whether sgRNA is
//   predicted to fall above or below the LOD
// This file is written to use LOO (so does not stratify by test/train data),
//   but it is otherwise the same file as stan-model-sgRNA-log-xval.stan
// This is embedded inside the full model, but is written separately for the
//   purpose of model selection on the logistic component of the hurdle model

data {

  // Number of observations
  int<lower=0> N; // number of observations in data
  int<lower=0> N_article; // for study-level hierarchy
  
  // Indicators for which predictors to consider
  int<lower=0, upper=1> dose_inc;
  int<lower=0, upper=1> dpi_inc;
  int<lower=0, upper=1> st_inc;
  int<lower=0, upper=1> sp_inc;
  int<lower=0, upper=1> tg_inc;
  int<lower=0, upper=1> age_inc;
  int<lower=0, upper=1> sex_inc;
  int<lower=0, upper=1> article_inc;
  
  // Levels for categorical predictors
  int<lower=0> L_dpi; // day post infection
  int<lower=0> L_st; // sample type
  int<lower=0> L_sp; // species
  int<lower=0> L_tg; // target gene
  int<lower=0> L_age; // age class
  
  // Outcome variables
  array[N] int<lower=0, upper=1> sg_pos;

  // All predictors
  array[N] real t; // primary predictor, total RNA value
  array[N] real dose; // log 10 
  array[N] int dpi; // day post infection
  array[N] int st; // sample type 
  array[N] int sp; // species 
  array[N] int tg; // target gene
  array[N] int age; // age class
  array[N] int sex; // sex

  // Hierarchy
  array[N] int article; // article level hierarchy for error terms

}

parameters {
  
  real gamma; 
  real deltaT;
  
  // written to enable toggling predictors on and off
  array[dose_inc ? 1 : 0] real deltaDOSE_p;
  array[dpi_inc ? L_dpi : 0] real deltaDPI_p;
  array[st_inc ? 1: 0] real deltaST_p;
  array[sp_inc ? L_sp : 0] real deltaSP_p;
  array[tg_inc ? L_tg : 0] real deltaTG_p;
  array[age_inc ? L_age : 0] real deltaAGE_p;
  array[sex_inc ? 1 : 0] real deltaSEX_p;
  array[article_inc ? N_article : 0] real deltaARTICLE_p;
  
  // Sex & age imputation
  real<lower=0, upper=1> p_adult;
  real<lower=0, upper=1> p_female;

}

transformed parameters {
  real deltaDOSE;
  real deltaSEX;
  real deltaST;
  array[L_sp] real deltaSP;
  array[L_tg] real deltaTG;
  array[L_dpi] real deltaDPI;
  array[L_age] real deltaAGE;
  array[N_article] real deltaARTICLE;
  
  // Continous variables
  if (dose_inc == 1) {
    deltaDOSE = deltaDOSE_p[1];
  }
  else if (dose_inc == 0) {
    deltaDOSE = 0;
  }
  
  if (sex_inc == 1) {
    deltaSEX = deltaSEX_p[1];
  }
  else if (sex_inc == 0) {
    deltaSEX = 0;
  }
  
  if (st_inc == 1) {
    deltaST = deltaST_p[1];
  }
  else if (st_inc == 0) {
    deltaST = 0;
  }
  
  // Categorical variables
  if (sp_inc == 1) {
    deltaSP = deltaSP_p;
  }
  else if (sp_inc == 0) {
    deltaSP = rep_array(0.0, L_sp);
  }
  
  if (dpi_inc == 1) {
    deltaDPI = deltaDPI_p;
  }
  else if (dpi_inc == 0) {
    deltaDPI = rep_array(0.0, L_dpi);
  }
  
  if (tg_inc == 1) {
    deltaTG = deltaTG_p;
  }
  else if (tg_inc == 0) {
    deltaTG = rep_array(0.0, L_tg);
  }
  
  if (age_inc == 1) {
    deltaAGE = deltaAGE_p;
  }
  else if (age_inc == 0) {
    deltaAGE = rep_array(0.0, L_age);
  }
  
  if (article_inc == 1) {
    deltaARTICLE = deltaARTICLE_p;
  }
  else if (article_inc == 0) {
    deltaARTICLE = rep_array(0.0, N_article);
  }
}

model {
  vector[4] lps; // for convenience during marginalization below
  
  for (n in 1:N) {
      if (age[n] != -9 && sex[n] != -9) {
        // demographics known, normal handling. 
        
        target += bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                    deltaDOSE * dose[n]  +
                                                    deltaDPI[dpi[n]] + 
                                                    deltaST * st[n] + 
                                                    deltaAGE[age[n]] +
                                                    deltaSEX * sex[n] + 
                                                    deltaTG[tg[n]] +
                                                    deltaSP[sp[n]] +
                                                    deltaARTICLE[article[n]]);
      }
      
      else if (age[n] != -9 && sex[n] == -9) {
        // sex unknown, age known. mixture required. 
        
        target += log_mix(p_female, 
                    bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                      deltaDOSE * dose[n] +
                                                      deltaDPI[dpi[n]] + 
                                                      deltaST * st[n] + 
                                                      deltaAGE[age[n]] +
                                                      deltaSEX * 1 + 
                                                      deltaTG[tg[n]] +
                                                      deltaSP[sp[n]] +
                                                      deltaARTICLE[article[n]]),
                    bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                      deltaDOSE * dose[n] +
                                                      deltaDPI[dpi[n]] + 
                                                      deltaST * st[n] + 
                                                      deltaAGE[age[n]] +
                                                      deltaSEX * 0 + 
                                                      deltaTG[tg[n]] +
                                                      deltaSP[sp[n]] +
                                                      deltaARTICLE[article[n]]));
      }
      
      else if (age[n] == -9 && sex[n] != -9) {
        // age unknown, sex known. mixture required. 
        target += log_mix(p_adult, 
                    bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                      deltaDOSE * dose[n] +
                                                      deltaDPI[dpi[n]] + 
                                                      deltaST * st[n] + 
                                                      deltaAGE[2] +
                                                      deltaSEX * sex[n] + 
                                                      deltaTG[tg[n]] +
                                                      deltaSP[sp[n]] +
                                                      deltaARTICLE[article[n]]),
                    bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                      deltaDOSE * dose[n] +
                                                      deltaDPI[dpi[n]] + 
                                                      deltaST * st[n] + 
                                                      deltaAGE[1] +
                                                      deltaSEX * sex[n] + 
                                                      deltaTG[tg[n]] +
                                                      deltaSP[sp[n]] +
                                                      deltaARTICLE[article[n]]));        
      }
      
      else if (age[n] == -9 && sex[n] == -9) {
        // age and sex unknown. large mixture required.
        lps[1] = log(p_adult) + 
                 log(p_female) + 
                 bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                   deltaDOSE * dose[n] +
                                                   deltaDPI[dpi[n]] + 
                                                   deltaST * st[n] + 
                                                   deltaAGE[2] +
                                                   deltaSEX * 1 + 
                                                   deltaTG[tg[n]] +
                                                   deltaSP[sp[n]] +
                                                   deltaARTICLE[article[n]]);
        lps[2] = log1m(p_adult) + 
                 log(p_female) + 
                 bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                   deltaDOSE * dose[n] +
                                                   deltaDPI[dpi[n]] + 
                                                   deltaST * st[n] + 
                                                   deltaAGE[1] +
                                                   deltaSEX * 1 + 
                                                   deltaTG[tg[n]] +
                                                   deltaSP[sp[n]] +
                                                   deltaARTICLE[article[n]]);
        lps[3] = log(p_adult) + 
                 log1m(p_female) + 
                 bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                   deltaDOSE * dose[n] +
                                                   deltaDPI[dpi[n]] + 
                                                   deltaST * st[n] + 
                                                   deltaAGE[2] +
                                                   deltaSEX * 0 + 
                                                   deltaTG[tg[n]] +
                                                   deltaSP[sp[n]] +
                                                   deltaARTICLE[article[n]]);
        lps[4] = log1m(p_adult) + 
                 log1m(p_female) + 
                 bernoulli_logit_lpmf(sg_pos[n] | gamma + deltaT * t[n] +
                                                   deltaDOSE * dose[n] +
                                                   deltaDPI[dpi[n]] + 
                                                   deltaST * st[n] + 
                                                   deltaAGE[1] +
                                                   deltaSEX * 0 + 
                                                   deltaTG[tg[n]] +
                                                   deltaSP[sp[n]] +
                                                   deltaARTICLE[article[n]]);
        target += log_sum_exp(lps);
      }
  }

  
  // Priors: Informative
  gamma ~ normal(-1, 2);
  deltaT ~ normal(2, 1);
  deltaDPI[1] ~ normal(-1, 1);
  deltaDPI[2] ~ normal(1, 1);
  deltaDPI[3] ~ normal(0, 1);
  deltaAGE ~ normal(0, 1);
  deltaSEX ~ normal(0, 1);
  deltaST ~ normal(0, 1);
  deltaSP ~ normal(0, 1);
  deltaTG[1] ~ normal(1, 1);
  deltaTG[2] ~ normal(1, 1);
  deltaTG[3] ~ normal(-1, 1); 
  deltaTG[4] ~ normal(-1, 1);
  deltaDOSE ~ normal(-1, 0.5);
  deltaARTICLE ~ normal(0, 0.5);
 
  
  // Priors: Non-Informative
  //gamma ~ normal(0, 1);
  //deltaT ~ normal(0, 1);
  //deltaDPI ~ normal(0, 1);
  //deltaAGE ~ normal(0, 1);
  //deltaSEX ~ normal(0, 1);
  //deltaST ~ normal(0, 1);
  //deltaSP ~ normal(0, 1);
  //deltaTG ~ normal(0, 1);
  //deltaDOSE ~ normal(0, 1);
}

generated quantities {
  
  // Log likelihoods & posterior predictions
  array[N] real log_lik;
  array[N] real p_pos;
  for (ii in 1:N) {
    if (age[ii] != -9 && sex[ii] != -9) {
      log_lik[ii] = bernoulli_logit_lpmf(sg_pos[ii] | gamma + deltaT * t[ii] +
                                                       deltaDOSE * dose[ii]  +
                                                       deltaDPI[dpi[ii]] + 
                                                       deltaST * st[ii] + 
                                                       deltaAGE[age[ii]] +
                                                       deltaSEX * sex[ii] + 
                                                       deltaTG[tg[ii]] +
                                                       deltaSP[sp[ii]] +
                                                       deltaARTICLE[article[ii]]);
      p_pos[ii] = bernoulli_logit_rng(gamma + deltaT * t[ii] +
                                       deltaDOSE * dose[ii]  +
                                       deltaDPI[dpi[ii]] + 
                                       deltaST * st[ii] + 
                                       deltaAGE[age[ii]] +
                                       deltaSEX * sex[ii] + 
                                       deltaTG[tg[ii]] +
                                       deltaSP[sp[ii]] +
                                       deltaARTICLE[article[ii]]);                                                     
  
    }
    
    else if (age[ii] == -9 && sex[ii] != -9) {
      log_lik[ii] = bernoulli_logit_lpmf(sg_pos[ii] | gamma + deltaT * t[ii] +
                                                       deltaDOSE * dose[ii]  +
                                                       deltaDPI[dpi[ii]] + 
                                                       deltaST * st[ii] + 
                                                       deltaAGE[bernoulli_rng(0.5) + 1] +
                                                       deltaSEX * sex[ii] + 
                                                       deltaTG[tg[ii]] +
                                                       deltaSP[sp[ii]] +
                                                       deltaARTICLE[article[ii]]);
      p_pos[ii] = bernoulli_logit_rng(gamma + deltaT * t[ii] +
                                       deltaDOSE * dose[ii]  +
                                       deltaDPI[dpi[ii]] + 
                                       deltaST * st[ii] + 
                                       deltaAGE[bernoulli_rng(0.5) + 1] +
                                       deltaSEX * sex[ii] + 
                                       deltaTG[tg[ii]] +
                                       deltaSP[sp[ii]] +
                                       deltaARTICLE[article[ii]]);                                                       
    }
    
    else if (age[ii] != -9 && sex[ii] == -9) {
      log_lik[ii] = bernoulli_logit_lpmf(sg_pos[ii] | gamma + deltaT * t[ii] +
                                                       deltaDOSE * dose[ii]  +
                                                       deltaDPI[dpi[ii]] + 
                                                       deltaST * st[ii] + 
                                                       deltaAGE[age[ii]] +
                                                       deltaSEX * bernoulli_rng(0.5) + 
                                                       deltaTG[tg[ii]] +
                                                       deltaSP[sp[ii]] +
                                                       deltaARTICLE[article[ii]]);
      
      p_pos[ii] = bernoulli_logit_rng(gamma + deltaT * t[ii] +
                                       deltaDOSE * dose[ii]  +
                                       deltaDPI[dpi[ii]] + 
                                       deltaST * st[ii] + 
                                       deltaAGE[age[ii]] +
                                       deltaSEX * bernoulli_rng(0.5) + 
                                       deltaTG[tg[ii]] +
                                       deltaSP[sp[ii]] +
                                       deltaARTICLE[article[ii]]); 
                                                            
    }
    
    else if (age[ii] == -9 && sex[ii] == -9) {
      log_lik[ii] = bernoulli_logit_lpmf(sg_pos[ii] | gamma + deltaT * t[ii] +
                                                       deltaDOSE * dose[ii]  +
                                                       deltaDPI[dpi[ii]] + 
                                                       deltaST * st[ii] + 
                                                       deltaAGE[bernoulli_rng(0.5) + 1] +
                                                       deltaSEX * bernoulli_rng(0.5) + 
                                                       deltaTG[tg[ii]] +
                                                       deltaSP[sp[ii]] +
                                                       deltaARTICLE[article[ii]]);     
      p_pos[ii] = bernoulli_logit_rng(gamma + deltaT * t[ii] +
                                       deltaDOSE * dose[ii]  +
                                       deltaDPI[dpi[ii]] + 
                                       deltaST * st[ii] + 
                                       deltaAGE[bernoulli_rng(0.5) + 1] +
                                       deltaSEX * bernoulli_rng(0.5) + 
                                       deltaTG[tg[ii]] +
                                       deltaSP[sp[ii]] +
                                       deltaARTICLE[article[ii]]);                                                              
    }
  }
}
