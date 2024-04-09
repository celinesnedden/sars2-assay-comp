//
// This Stan program defines the final hurdle model for sgRNA and total RNA, 
//   including both a logistic & a linear component 
// It includes only the predictors selected by the model selection procedure 


data {
  
  // MODEL FITTING
  
  // Number of observations
  int<lower=0> N; // total number of observations (when tRNA > LOD)
  int<lower=0> N_article; // for study-level hierarchy
  int<lower=0> N_lab;
  
  // Levels for categorical predictors
  int<lower=0> L_dpi; // day post infection
  int<lower=0> L_sp; // species
  int<lower=0> L_tg; // target gene

  // Outcome variables
  array[N] int<lower=0, upper=1> pos_sg; // 0: sg < LOD; 1: sg > LOD
  array[N] real val_sg; // log10 value of sg, if > LOD; -9 otherwise

  // All predictors
  array[N] real t_val; // primary predictor, total RNA value
  array[N] real dose; // log 10 pfu
  array[N] int dpi; // day post infection
  array[N] int sp; // species 
  array[N] int tg; // target gene

  // Hierarchical error terms per study
  array[N] int article; 
  array[N] int lab;

}

parameters {
  
  // Logistic regression component
  real gamma;
  real deltaT;
  real deltaDOSE;
  array[L_tg] real deltaTG;
  array[L_sp] real deltaSP;
  array[N_lab] real deltaLAB;

  // Linear regression component
  real alpha;
  real<lower=0> betaT;
  real betaDOSE; 
  array[L_tg] real betaTG; 
  array[L_sp] real betaSP;
  array[L_dpi] real betaDPI;
  array[N_lab] real betaLAB;

  // Hierarchical errors for linear regression
  vector<lower=0>[N_article] sigma; // study-specific error 
  real<lower=0> sigma_bar; // average sigma
  real<lower=0> sigma_sd; // standard deviation of distribution of sigmas
  
  // Age imputation
  real<lower=0, upper=1> p_adult;

}

model {

  // Logistic component: everything contributes
  for (n in 1:N) {
      target += bernoulli_logit_lpmf(pos_sg[n] | gamma + 
                                                 deltaT * t_val[n] +
                                                 deltaDOSE * dose[n]  +
                                                 deltaTG[tg[n]] +
                                                 deltaSP[sp[n]] +
                                                 deltaLAB[lab[n]]);
  }
  
  // Linear regression model: only sgRNA+ samples with quant data contribute
  for (n in 1:N) {
    if (pos_sg[n] == 1 && val_sg[n] != -9) {
      target += normal_lpdf(val_sg[n] | alpha + 
                                        betaT * t_val[n] + 
                                        betaTG[tg[n]] +
                                        betaSP[sp[n]] + 
                                        betaDPI[dpi[n]] + 
                                        betaDOSE * dose[n] +
                                        betaLAB[lab[n]],
                                        sigma[article[n]]);
    }
  }
  
  // Priors: Informative, same as from model w/o article effect
  // Linear Component
  alpha ~ normal(-3, 1);
  betaT ~ gamma(2, 0.5);
  //betaDOSE ~ normal(-0.5, 0.25);
  betaDOSE ~ normal(-0.25, 1);
  betaDPI[1] ~ normal(-0.5, 1);
  betaDPI[2] ~ normal(0.5, 1);
  betaDPI[3] ~ normal(0, 1);
  betaSP ~ normal(0, 1);
  betaTG[1] ~ normal(0.5, 1);
  betaTG[2] ~ normal(0.5, 1);
  betaTG[3] ~ normal(-0.5, 1);
  betaTG[4] ~ normal(-0.5, 1);
  betaLAB ~ normal(0, 0.5);
  sigma ~ normal(sigma_bar, sigma_sd);
  sigma_bar ~ normal(0, 1);
  sigma_sd ~ exponential(1);
  
  // Logistic Component
  gamma ~ normal(-1, 2);
  deltaT ~ normal(2, 1);
  deltaTG[1] ~ normal(1, 1);
  deltaTG[2] ~ normal(1, 1);
  deltaTG[3] ~ normal(-1, 1);
  deltaTG[4] ~ normal(-1, 1);
  deltaDOSE ~ normal(-1, 0.5);
  deltaSP ~ normal(0, 1);
  deltaLAB ~ normal(0, 0.5);
  
  
  
}

generated quantities {
  
  // Predicting detectability and values of known sgRNA values
  array[N] real prob_pos; 
  array[N] real pred_val;
  
  for (ii in 1:N) {
    prob_pos[ii] = bernoulli_logit_rng(gamma + 
                                       deltaT * t_val[ii] +
                                       deltaDOSE * dose[ii]  +
                                       deltaTG[tg[ii]] +
                                       deltaSP[sp[ii]] + 
                                       deltaLAB[lab[ii]]);
    
    
    pred_val[ii] = normal_rng(alpha + 
                              betaT * t_val[ii] + 
                              betaTG[tg[ii]] +
                              betaSP[sp[ii]] + 
                              betaDPI[dpi[ii]] + 
                              betaDOSE * dose[ii] + 
                              betaLAB[lab[ii]],
                              sigma[article[ii]]);
    
    }
}
