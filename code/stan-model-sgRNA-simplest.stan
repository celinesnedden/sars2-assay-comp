//
// This Stan program defines the simplest hurdle model for sgRNA and total RNA, 
//   where total RNA is the only predictor


data {
  
  // Number of observations
  int<lower=0> N; // total number of observations (when tRNA > LOD)
  int<lower=0> N_lin;

  // Outcome variables
  array[N] int<lower=0, upper=1> pos_sg; // 0: sg < LOD; 1: sg > LOD
  array[N] real val_sg; // log10 value of sg, if > LOD; -9 otherwise

  // Single predictor
  array[N] real t_val; // primary predictor, total RNA value
  
}

parameters {
  
  // Logistic regression component
  real gamma;
  real deltaT;

  // Linear regression component
  real alpha;
  real<lower=0> betaT;

  // Hierarchical errors for linear regression
  real<lower=0> sigma; 

}

model {

  // Logistic component: everything contributes
  for (n in 1:N) {
      target += bernoulli_logit_lpmf(pos_sg[n] | gamma + deltaT * t_val[n]);
      // Linear component: only positive values contribute
      if (pos_sg[n] == 1 && val_sg[n] != -9) {
          target += normal_lpdf(val_sg[n] | alpha + betaT * t_val[n], sigma);
      }
  }
  
  // Priors: Non-informative
  //gamma ~ normal(0, 1);
  //deltaT ~ normal(0, 1); 
  //alpha ~ normal(-2, 1);
  //betaT ~ gamma(2, 0.5);
  //sigma ~ normal(0, 1);
  
  // Priors: Informative
  gamma ~ normal(-1, 2);
  deltaT ~ normal(2, 1); 
  alpha ~ normal(-3, 1);
  betaT ~ gamma(2, 0.5);
  sigma ~ normal(0, 1);

}

generated quantities {
  
  // Predicting detectability and values of sgRNA 
  array[N] real prob_pos; 
  array[N] real pred_val;
  
  for (ii in 1:N) {
    
    prob_pos[ii] = bernoulli_logit_rng(gamma + deltaT * t_val[ii]);
    pred_val[ii] = normal_rng(alpha + betaT * t_val[ii], sigma); 
    
  }
    
  // Log likelihoods for performance statistics
  array[N] real log_lik_log;
  array[N_lin] real log_lik_lin;
  int ii_lin = 1;
  
  for (ii in 1:N) {
    
    log_lik_log[ii] = bernoulli_logit_lpmf(pos_sg[ii] | gamma + deltaT * t_val[ii]);
      
    // Only positive values contribute to likelihood on linear component
    if (pos_sg[ii] == 1 && val_sg[ii] != -9) {
      log_lik_lin[ii_lin] = normal_lpdf(val_sg[ii] | alpha + betaT * t_val[ii], 
                                                     sigma);
     ii_lin = ii_lin + 1;                                                 
  }
}
}
