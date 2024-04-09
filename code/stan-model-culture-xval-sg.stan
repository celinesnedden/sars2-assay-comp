//
// This Stan program defines a logistic regression for whether culture is 
//   predicted to be censored or not

data {

  // Number of observations
  int<lower=0> N_train; // number of observations in training data
  int<lower=0> N_test; // number of observations in test data
  int<lower=0> N_article; // for study-level hierarchy
  
  // Indicators for which predictors to consider
  int<lower=0, upper=1> t_inc;
  int<lower=0, upper=1> sg_inc;
  int<lower=0, upper=1> dose_inc;
  int<lower=0, upper=1> dpi_inc;
  int<lower=0, upper=1> st_inc;
  int<lower=0, upper=1> sp_inc;
  int<lower=0, upper=1> tg_inc;
  int<lower=0, upper=1> age_inc;
  int<lower=0, upper=1> sex_inc;
  int<lower=0, upper=1> cell_inc;
  int<lower=0, upper=1> assay_inc;
  
  // Levels for categorical predictors
  int<lower=0> L_dpi; // day post infection
  int<lower=0> L_sp; // species
  int<lower=0> L_tg; // target gene
  int<lower=0> L_age; // age class
  int<lower=0> L_cell; // cell line
  
  // Outcome variables
  array[N_train] int<lower=0, upper=1> pos_inf_train;
  array[N_test] int<lower=0, upper=1> pos_inf_test;

  // All predictors
  array[N_train] real t_train; // primary predictor, total RNA value
  array[N_test] real t_test; 
  array[N_train] real sg_train; // primary predictor, sgRNA value
  array[N_test] real sg_test; 
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
  array[N_train] int cell_train; // cell line
  array[N_test] int cell_test;
  array[N_train] int assay_train; // cell line
  array[N_test] int assay_test;

}

parameters {
  
  real gamma;
  array[t_inc ? 1 : 0] real psiT_p;
  array[sg_inc ? 1 : 0] real psiSG_p;
  array[dose_inc ? 1 : 0] real psiDOSE_p;
  array[dpi_inc ? L_dpi : 0] real psiDPI_p;
  array[st_inc ? 1: 0] real psiST_p;
  array[sp_inc ? L_sp : 0] real psiSP_p;
  array[tg_inc ? L_tg : 0] real psiTG_p;
  array[age_inc ? L_age : 0] real psiAGE_p;
  array[sex_inc ? 1 : 0] real psiSEX_p;
  array[cell_inc ? L_cell : 0] real psiCELL_p;
  array[assay_inc ? 1 : 0] real psiASSAY_p;
  
  // Sex & age imputation
  real<lower=0, upper=1> p_female; // easy two options for sex
  simplex[L_age] age_probs;
  //real<lower=0, upper=1> p_adult; // need simplex for age with 3 options

}

transformed parameters {
  real psiT;
  real psiSG;
  real psiDOSE;
  real psiSEX;
  real psiST;
  real psiASSAY;
  array[L_sp] real psiSP;
  array[L_tg] real psiTG;
  array[L_dpi] real psiDPI;
  array[L_age] real psiAGE;
  array[L_cell] real psiCELL;
  
  // Continous variables
  if (t_inc == 1) {
    psiT = psiT_p[1];
  }
  else if (t_inc == 0) {
    psiT = 0;
  }
  
  if (sg_inc == 1) {
    psiSG = psiSG_p[1];
  }
  else if (sg_inc == 0) {
    psiSG = 0;
  }
  
  if (dose_inc == 1) {
    psiDOSE = psiDOSE_p[1];
  }
  else if (dose_inc == 0) {
    psiDOSE = 0;
  }
  
  if (sex_inc == 1) {
    psiSEX = psiSEX_p[1];
  }
  else if (sex_inc == 0) {
    psiSEX = 0;
  }
  
  if (st_inc == 1) {
    psiST = psiST_p[1];
  }
  else if (st_inc == 0) {
    psiST = 0;
  }
  
  if (assay_inc == 1) {
    psiASSAY = psiASSAY_p[1];
  }
  else if (assay_inc == 0) {
    psiASSAY = 0;
  }
  
  // Categorical variables
  if (sp_inc == 1) {
    psiSP = psiSP_p;
  }
  else if (sp_inc == 0) {
    psiSP = rep_array(0.0, L_sp);
  }
  
  if (dpi_inc == 1) {
    psiDPI = psiDPI_p;
  }
  else if (dpi_inc == 0) {
    psiDPI = rep_array(0.0, L_dpi);
  }
  
  if (tg_inc == 1) {
    psiTG = psiTG_p;
  }
  else if (tg_inc == 0) {
    psiTG = rep_array(0.0, L_tg);
  }
  
  if (age_inc == 1) {
    psiAGE = psiAGE_p;
  }
  else if (age_inc == 0) {
    psiAGE = rep_array(0.0, L_age);
  }
  
  if (cell_inc == 1) {
    psiCELL = psiCELL_p;
  }
  else if (cell_inc == 0) {
    psiCELL = rep_array(0.0, L_cell);
  }
}

model {
  vector[4] lps; // for convenience during marginalization below
  vector[2] lps_small; // for convenience during marginalization below
  vector[3] theta;
  
  for (n in 1:N_train) {
      if (age_train[n] != -9 && sex_train[n] != -9) {
        // demographics known, normal handling.
        
        target += bernoulli_logit_lpmf(pos_inf_train[n] | gamma + psiT * t_train[n] +
                                                          psiSG * sg_train[n] +
                                                          psiDOSE * dose_train[n]  +
                                                          psiDPI[dpi_train[n]] + 
                                                          psiST * st_train[n] + 
                                                          psiAGE[age_train[n]] +
                                                          psiSEX * sex_train[n] + 
                                                          psiTG[tg_train[n]] +
                                                          psiSP[sp_train[n]] +
                                                          psiASSAY * assay_train[n] +
                                                          psiCELL[cell_train[n]]);
      }
      
      else if (age_train[n] != -9 && sex_train[n] == -9) {
        // sex_train unknown, age_train known. mixture required.
        
        target += log_mix(p_female, 
                    bernoulli_logit_lpmf(pos_inf_train[n] | gamma + psiT * t_train[n] + 
                                                            psiSG * sg_train[n] +
                                                            psiDOSE * dose_train[n] +
                                                            psiDPI[dpi_train[n]] + 
                                                            psiST * st_train[n] + 
                                                            psiAGE[age_train[n]] +
                                                            psiSEX * 1 + 
                                                            psiTG[tg_train[n]] +
                                                            psiSP[sp_train[n]] +
                                                            psiASSAY * assay_train[n] +
                                                            psiCELL[cell_train[n]]),
                    bernoulli_logit_lpmf(pos_inf_train[n] | gamma + psiT * t_train[n] +
                                                            psiSG * sg_train[n] +
                                                            psiDOSE * dose_train[n] +
                                                            psiDPI[dpi_train[n]] + 
                                                            psiST * st_train[n] + 
                                                            psiAGE[age_train[n]] +
                                                            psiSEX * 0 + 
                                                            psiTG[tg_train[n]] +
                                                            psiSP[sp_train[n]] +
                                                            psiASSAY * assay_train[n] +
                                                            psiCELL[cell_train[n]]));
      }
      
      else if (age_train[n] == -9 && sex_train[n] != -9) {
        // age_train unknown, sex_train known. mixture required.
        
        lps_small[1] = log(age_probs[1]) + 
                                bernoulli_logit_lpmf(pos_inf_train[n] | 
                                                     gamma + psiT * t_train[n] +
                                                     psiSG * sg_train[n] +
                                                     psiDOSE * dose_train[n] +
                                                     psiDPI[dpi_train[n]] + 
                                                     psiST * st_train[n] + 
                                                     psiAGE[1] +
                                                     psiSEX * sex_train[n] + 
                                                     psiTG[tg_train[n]] +
                                                     psiSP[sp_train[n]] +
                                                     psiASSAY * assay_train[n] +
                                                     psiCELL[cell_train[n]]);
                                                     
        lps_small[2] = log(age_probs[2]) + 
                                bernoulli_logit_lpmf(pos_inf_train[n] | 
                                                     gamma + psiT * t_train[n] +
                                                     psiSG * sg_train[n] +
                                                     psiDOSE * dose_train[n] +
                                                     psiDPI[dpi_train[n]] + 
                                                     psiST * st_train[n] + 
                                                     psiAGE[2] +
                                                     psiSEX * sex_train[n] + 
                                                     psiTG[tg_train[n]] +
                                                     psiSP[sp_train[n]] +
                                                     psiASSAY * assay_train[n] +
                                                     psiCELL[cell_train[n]]); 

                                                     
        target += log_sum_exp(lps_small);                                             

      }
      
      else if (age_train[n] == -9 && sex_train[n] == -9) {
        // age_train and sex_train unknown. large mixture required.
        lps[1] = log(age_probs[1]) + 
                 log(p_female) + 
                 bernoulli_logit_lpmf(pos_inf_train[n] | gamma + psiT * t_train[n] +
                                                         psiSG * sg_train[n] +
                                                         psiDOSE * dose_train[n] +
                                                         psiDPI[dpi_train[n]] + 
                                                         psiST * st_train[n] + 
                                                         psiAGE[1] +
                                                         psiSEX * 1 + 
                                                         psiTG[tg_train[n]] +
                                                         psiSP[sp_train[n]] +
                                                         psiASSAY * assay_train[n] +
                                                         psiCELL[cell_train[n]]);
        lps[2] = log(age_probs[2]) + 
                 log(p_female) + 
                 bernoulli_logit_lpmf(pos_inf_train[n] | gamma + psiT * t_train[n] +
                                                         psiSG * sg_train[n] +
                                                         psiDOSE * dose_train[n] +
                                                         psiDPI[dpi_train[n]] + 
                                                         psiST * st_train[n] + 
                                                         psiAGE[2] +
                                                         psiSEX * 1 + 
                                                         psiTG[tg_train[n]] +
                                                         psiSP[sp_train[n]] +
                                                         psiASSAY * assay_train[n] +
                                                         psiCELL[cell_train[n]]);
                                                      
        lps[3] = log(age_probs[1]) + 
                 log1m(p_female) + 
                 bernoulli_logit_lpmf(pos_inf_train[n] | gamma + psiT * t_train[n] +
                                                         psiSG * sg_train[n] +
                                                         psiDOSE * dose_train[n] +
                                                         psiDPI[dpi_train[n]] + 
                                                         psiST * st_train[n] + 
                                                         psiAGE[1] +
                                                         psiSEX * 0 + 
                                                         psiTG[tg_train[n]] +
                                                         psiSP[sp_train[n]] +
                                                         psiASSAY * assay_train[n] +
                                                         psiCELL[cell_train[n]]);
        lps[4] = log(age_probs[2]) + 
                 log1m(p_female) + 
                 bernoulli_logit_lpmf(pos_inf_train[n] | gamma + psiT * t_train[n] +
                                                         psiSG * sg_train[n] +
                                                         psiDOSE * dose_train[n] +
                                                         psiDPI[dpi_train[n]] + 
                                                         psiST * st_train[n] + 
                                                         psiAGE[2] +
                                                         psiSEX * 0 + 
                                                         psiTG[tg_train[n]] +
                                                         psiSP[sp_train[n]] +
                                                         psiASSAY * assay_train[n] +
                                                         psiCELL[cell_train[n]]);

        target += log_sum_exp(lps);
      }
  }

  
  // Priors: Non-Informative
  //gamma ~ normal(0, 1);
  //psiT ~ normal(0, 1);
  //psiSG ~ normal(0, 1);
  //psiDPI ~ normal(0, 1);
  //psiAGE ~ normal(0, 1);
  //psiSEX ~ normal(0, 1);
  //psiST ~ normal(0, 1);
  //psiSP ~ normal(0, 1);
  //psiTG ~ normal(0, 1);
  //psiDOSE ~ normal(0, 1);
  //psiCELL ~ normal(0, 1);
  //psiASSAY ~ normal(0, 1);
  
  // Priors: Informative
  gamma ~ normal(-1, 1); 
  psiT ~ normal(1, 0.5);
  psiSG ~ normal(1, 0.5);
  psiDPI[1] ~ normal(0, 1);
  psiDPI[2] ~ normal(0, 1); 
  psiDPI[3] ~ normal(0, 1);
  psiAGE ~ normal(0, 1);
  psiSEX ~ normal(0, 1);
  psiST ~ normal(0, 1);
  psiSP ~ normal(0, 1);
  psiTG[1] ~ normal(-0.5, 1);
  psiTG[2] ~ normal(0, 1);
  psiTG[3] ~ normal(0.5, 1);
  psiDOSE ~ normal(0, 1); 
  psiCELL[1] ~ normal(0, 1);
  psiCELL[2] ~ normal(0, 1);
  psiCELL[3] ~ normal(0.5, 1);
  psiASSAY ~ normal(-0.5, 1);

}

generated quantities {
  
  // TRAINING DATA
  
  // Log likelihoods & posterior predictions
  // array[N_train] real log_lik;
  array[N_train] real p_pos_train;
  for (ii in 1:N_train) {
    if (age_train[ii] != -9 && sex_train[ii] != -9) {

      p_pos_train[ii] = bernoulli_logit_rng(gamma + psiT * t_train[ii] +
                                            psiSG * sg_train[ii] +
                                            psiDOSE * dose_train[ii]  +
                                            psiDPI[dpi_train[ii]] + 
                                            psiST * st_train[ii] + 
                                            psiAGE[age_train[ii]] +
                                            psiSEX * sex_train[ii] + 
                                            psiTG[tg_train[ii]] +
                                            psiSP[sp_train[ii]] +
                                            psiASSAY * assay_train[ii] +
                                            psiCELL[cell_train[ii]]);                                                     
  
    }
    
    else if (age_train[ii] == -9 && sex_train[ii] != -9) {

      p_pos_train[ii] = bernoulli_logit_rng(gamma + psiT * t_train[ii] +
                                            psiSG * sg_train[ii] +
                                            psiDOSE * dose_train[ii]  +
                                            psiDPI[dpi_train[ii]] + 
                                            psiST * st_train[ii] + 
                                            psiAGE[discrete_range_rng(1, L_age)] +
                                            psiSEX * sex_train[ii] + 
                                            psiTG[tg_train[ii]] +
                                            psiSP[sp_train[ii]] +
                                            psiASSAY * assay_train[ii] +
                                            psiCELL[cell_train[ii]]);                                                       
    }
    
    else if (age_train[ii] != -9 && sex_train[ii] == -9) {

      p_pos_train[ii] = bernoulli_logit_rng(gamma + psiT * t_train[ii] +
                                            psiSG * sg_train[ii] +
                                            psiDOSE * dose_train[ii]  +
                                            psiDPI[dpi_train[ii]] + 
                                            psiST * st_train[ii] + 
                                            psiAGE[age_train[ii]] +
                                            psiSEX * bernoulli_rng(0.5) + 
                                            psiTG[tg_train[ii]] +
                                            psiSP[sp_train[ii]] +
                                            psiASSAY * assay_train[ii] +
                                            psiCELL[cell_train[ii]]); 
                                                            
    }
    
    else if (age_train[ii] == -9 && sex_train[ii] == -9) {
      
      p_pos_train[ii] = bernoulli_logit_rng(gamma + psiT * t_train[ii] +
                                            psiSG * sg_train[ii] +
                                            psiDOSE * dose_train[ii]  +
                                            psiDPI[dpi_train[ii]] + 
                                            psiST * st_train[ii] + 
                                            psiAGE[discrete_range_rng(1, L_age)] +
                                            psiSEX * bernoulli_rng(0.5) + 
                                            psiTG[tg_train[ii]] +
                                            psiSP[sp_train[ii]] +
                                            psiASSAY * assay_train[ii] +
                                            psiCELL[cell_train[ii]]);                                                              
    }
  }
  
  // TEST DATA
  array[N_test] real log_lik;
  array[N_test] real p_pos_test;
  
  for (ii in 1:N_test) {
    if (age_test[ii] != -9 && sex_test[ii] != -9) {
      log_lik[ii] = bernoulli_logit_lpmf(pos_inf_test[ii] | gamma + psiT * t_test[ii] +
                                                            psiSG * sg_test[ii] +
                                                            psiDOSE * dose_test[ii]  +
                                                            psiDPI[dpi_test[ii]] + 
                                                            psiST * st_test[ii] + 
                                                            psiAGE[age_test[ii]] +
                                                            psiSEX * sex_test[ii] + 
                                                            psiTG[tg_test[ii]] +
                                                            psiSP[sp_test[ii]] +
                                                            psiASSAY * assay_test[ii] +
                                                            psiCELL[cell_test[ii]]);
      p_pos_test[ii] = bernoulli_logit_rng(gamma + psiT * t_test[ii] +
                                           psiSG * sg_test[ii] +
                                           psiDOSE * dose_test[ii]  +
                                           psiDPI[dpi_test[ii]] + 
                                           psiST * st_test[ii] + 
                                           psiAGE[age_test[ii]] +
                                           psiSEX * sex_test[ii] + 
                                           psiTG[tg_test[ii]] +
                                           psiSP[sp_test[ii]] +
                                           psiASSAY * assay_test[ii] +
                                           psiCELL[cell_test[ii]]);                                                     
  
    }
    
    else if (age_test[ii] == -9 && sex_test[ii] != -9) {
      log_lik[ii] = bernoulli_logit_lpmf(pos_inf_test[ii] | gamma + psiT * t_test[ii] +
                                                            psiSG * sg_test[ii] +
                                                            psiDOSE * dose_test[ii]  +
                                                            psiDPI[dpi_test[ii]] + 
                                                            psiST * st_test[ii] + 
                                                            psiAGE[discrete_range_rng(1, L_age)] +
                                                            psiSEX * sex_test[ii] + 
                                                            psiTG[tg_test[ii]] +
                                                            psiSP[sp_test[ii]] +
                                                            psiASSAY * assay_test[ii] +
                                                            psiCELL[cell_test[ii]]);
      p_pos_test[ii] = bernoulli_logit_rng(gamma + psiT * t_test[ii] +
                                           psiSG * sg_test[ii] +
                                           psiDOSE * dose_test[ii]  +
                                           psiDPI[dpi_test[ii]] + 
                                           psiST * st_test[ii] + 
                                           psiAGE[discrete_range_rng(1, L_age)] +
                                           psiSEX * sex_test[ii] + 
                                           psiTG[tg_test[ii]] +
                                           psiSP[sp_test[ii]] +
                                           psiASSAY * assay_test[ii] +
                                           psiCELL[cell_test[ii]]);                                                       
    }
    
    else if (age_test[ii] != -9 && sex_test[ii] == -9) {
      log_lik[ii] = bernoulli_logit_lpmf(pos_inf_test[ii] | gamma + psiT * t_test[ii] +
                                                            psiSG * sg_test[ii] +
                                                            psiDOSE * dose_test[ii]  +
                                                            psiDPI[dpi_test[ii]] + 
                                                            psiST * st_test[ii] + 
                                                            psiAGE[age_test[ii]] +
                                                            psiSEX * bernoulli_rng(0.5) + 
                                                            psiTG[tg_test[ii]] +
                                                            psiSP[sp_test[ii]] +
                                                            psiASSAY * assay_test[ii] +
                                                            psiCELL[cell_test[ii]]);
      
      p_pos_test[ii] = bernoulli_logit_rng(gamma + psiT * t_test[ii] +
                                           psiSG * sg_test[ii] +
                                           psiDOSE * dose_test[ii]  +
                                           psiDPI[dpi_test[ii]] + 
                                           psiST * st_test[ii] + 
                                           psiAGE[age_test[ii]] +
                                           psiSEX * bernoulli_rng(0.5) + 
                                           psiTG[tg_test[ii]] +
                                           psiSP[sp_test[ii]] +
                                           psiASSAY * assay_test[ii] +
                                           psiCELL[cell_test[ii]]); 
                                                            
    }
    
    else if (age_test[ii] == -9 && sex_test[ii] == -9) {
      log_lik[ii] = bernoulli_logit_lpmf(pos_inf_test[ii] | gamma + psiT * t_test[ii] +
                                                            psiSG * sg_test[ii] +
                                                            psiDOSE * dose_test[ii]  +
                                                            psiDPI[dpi_test[ii]] + 
                                                            psiST * st_test[ii] + 
                                                            psiAGE[discrete_range_rng(1, L_age)] +
                                                            psiSEX * bernoulli_rng(0.5) + 
                                                            psiTG[tg_test[ii]] +
                                                            psiSP[sp_test[ii]] +
                                                            psiASSAY * assay_test[ii] +
                                                            psiCELL[cell_test[ii]]);     
      p_pos_test[ii] = bernoulli_logit_rng(gamma + psiT * t_test[ii] +
                                           psiSG * sg_test[ii] +
                                           psiDOSE * dose_test[ii]  +
                                           psiDPI[dpi_test[ii]] + 
                                           psiST * st_test[ii] + 
                                           psiAGE[discrete_range_rng(1, L_age)] +
                                           psiSEX * bernoulli_rng(0.5) + 
                                           psiTG[tg_test[ii]] +
                                           psiSP[sp_test[ii]] +
                                           psiASSAY * assay_test[ii] +
                                           psiCELL[cell_test[ii]]);                                                              
    }
  }
}
