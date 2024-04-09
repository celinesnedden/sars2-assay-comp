//
// This Stan program defines a logistic regression for whether culture is 
//   predicted to be censored or not

data {

  // Number of observations
  int<lower=0> N; // number of observations in training data
  
  // Levels for categorical predictors
  int<lower=0> L_dpi; // day post infection
  int<lower=0> L_sp; // species
  int<lower=0> L_tg; // target gene
  int<lower=0> L_age; // age class
  int<lower=0> L_cell; // cell line
  
  // Outcome variables
  array[N] int<lower=0, upper=1> pos_inf;

  // All predictors
  array[N] real t; // primary predictor, total RNA value
  array[N] real dose; // log 10 pfu
  array[N] int dpi; // day post infection
  array[N] int sp; // species 
  array[N] int tg; // target gene
  array[N] int age; // age class
  array[N] int cell; // cell line
  array[N] int assay; // assay type

}

parameters {
  
  real gamma;
  real psiT;
  real psiDOSE;
  array[L_dpi] real psiDPI;
  array[L_sp] real psiSP;
  array[L_tg] real psiTG;
  array[L_age] real psiAGE;
  array[L_cell] real psiCELL;
  real psiASSAY;
  
  // Sex & age imputation
  real<lower=0, upper=1> p_female; // easy two options for sex
  simplex[L_age] age_probs;

}


model {
  vector[3] lps_small; // for convenience during marginalization below
  vector[3] theta;
  
  for (n in 1:N) {
      if (age[n] != -9) {
        
        target += bernoulli_logit_lpmf(pos_inf[n] | gamma + psiT * t[n] +
                                                    psiDOSE * dose[n]  +
                                                    psiDPI[dpi[n]] + 
                                                    psiAGE[age[n]] +
                                                    psiTG[tg[n]] +
                                                    psiSP[sp[n]] +
                                                    psiASSAY * assay[n] +
                                                    psiCELL[cell[n]]);
      }

      
      
      else if (age[n] == -9) {
        
        lps_small[1] = log(age_probs[1]) + 
                                bernoulli_logit_lpmf(pos_inf[n] | 
                                                     gamma + psiT * t[n] +
                                                     psiDOSE * dose[n] +
                                                     psiDPI[dpi[n]] + 
                                                     psiAGE[1] + 
                                                     psiTG[tg[n]] +
                                                     psiSP[sp[n]] +
                                                     psiASSAY * assay[n] +
                                                     psiCELL[cell[n]]);
                                                     
        lps_small[2] = log(age_probs[2]) + 
                                bernoulli_logit_lpmf(pos_inf[n] | 
                                                     gamma + psiT * t[n] +
                                                     psiDOSE * dose[n] +
                                                     psiDPI[dpi[n]] + 
                                                     psiAGE[2] +
                                                     psiTG[tg[n]] +
                                                     psiSP[sp[n]] +
                                                     psiASSAY * assay[n] +
                                                     psiCELL[cell[n]]); 
                                                     
        lps_small[3] = log(age_probs[3]) + 
                                bernoulli_logit_lpmf(pos_inf[n] | 
                                                     gamma + psiT * t[n] +
                                                     psiDOSE * dose[n] +
                                                     psiDPI[dpi[n]] + 
                                                     psiAGE[3] + 
                                                     psiTG[tg[n]] +
                                                     psiSP[sp[n]] +
                                                     psiASSAY * assay[n] +
                                                     psiCELL[cell[n]]);   
                                                     
        target += log_sum_exp(lps_small);                                             

      }
      
  }

  
  // Priors: Non-Informative
 //gamma ~ normal(0, 1);
 //psiT ~ normal(0, 1);
 //psiDPI ~ normal(0, 1);
 //psiAGE ~ normal(0, 1);
 //psiSP ~ normal(0, 1);
 //psiTG ~ normal(0, 1);
 //psiDOSE ~ normal(0, 1);
 //psiCELL ~ normal(0, 1);
 //psiASSAY ~ normal(0, 1);
  
  // Priors: Informative
  gamma ~ normal(-1, 1); 
  psiT ~ normal(1, 0.5);
  psiDPI[1] ~ normal(0, 1);
  psiDPI[2] ~ normal(0, 1); 
  psiDPI[3] ~ normal(0, 1);
  psiAGE ~ normal(0, 1);
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
  
  // Log likelihoods & posterior predictions
  array[N] real p_pos;
  for (ii in 1:N) {
    if (age[ii] != -9) {

      p_pos[ii] = bernoulli_logit_rng(gamma + psiT * t[ii] +
                                      psiDOSE * dose[ii]  +
                                      psiDPI[dpi[ii]] + 
                                      psiAGE[age[ii]] +
                                      psiTG[tg[ii]] +
                                      psiSP[sp[ii]] +
                                      psiASSAY * assay[ii] +
                                      psiCELL[cell[ii]]);                                                     
  
    }
    
    else if (age[ii] == -9) {

      p_pos[ii] = bernoulli_logit_rng(gamma + psiT * t[ii] +
                                      psiDOSE * dose[ii]  +
                                      psiDPI[dpi[ii]] + 
                                      psiAGE[discrete_range_rng(1, 3)] +
                                      psiTG[tg[ii]] +
                                      psiSP[sp[ii]] +
                                      psiASSAY * assay[ii] +
                                      psiCELL[cell[ii]]);                                                       
    }
  }
}
