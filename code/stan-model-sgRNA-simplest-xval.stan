// This Stan program defines the simplest hurdle model for sgRNA and total RNA, 
//   where total RNA is the only predictor, for comparison with other tested models


data {
  
  // Number of observations
  int<lower=0> N_train; // total number of observations (when tRNA > LOD)
  int<lower=0> N_test; 
  int<lower=0> N_test_lin;

  // Outcome variables
  array[N_train] int<lower=0, upper=1> pos_sg_train; // 0: sg < LOD; 1: sg > LOD
  array[N_test] int<lower=0, upper=1> pos_sg_test; // 0: sg < LOD; 1: sg > LOD

  array[N_train] real val_sg_train; // log10 value of sg, if > LOD; -9 otherwise
  array[N_test] real val_sg_test; // log10 value of sg, if > LOD; -9 otherwise

  // Single predictor
  array[N_train] real t_train; // primary predictor, total RNA value
  array[N_test] real t_test; // primary predictor, total RNA value

  
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
  for (n in 1:N_train) {
      target += bernoulli_logit_lpmf(pos_sg_train[n] | gamma + deltaT * t_train[n]);
      // Linear component: only positive values contribute
      if (pos_sg_train[n] == 1 && val_sg_train[n] != -9) {
          target += normal_lpdf(val_sg_train[n] | alpha + betaT * t_train[n], sigma);
      }
  }
  
  // Priors: Non-informative
  gamma ~ normal(0, 1);
  deltaT ~ normal(0, 1); 
  alpha ~ normal(-2, 1);
  betaT ~ gamma(2, 0.5);
  sigma ~ normal(0, 1);
  
  // Priors: Informative
  //gamma ~ normal(-1, 2);
  //deltaT ~ normal(2, 1); 
  //alpha ~ normal(-3, 1);
  //betaT ~ gamma(2, 0.5);
  //sigma ~ normal(0, 1);


}

generated quantities {
  
  // Log likelihoods for performance statistics
  array[N_test] real log_lik_log;
  array[N_test_lin] real log_lik_lin;
  int ii_lin = 1;
  
  // Predicted detectability and values of sgRNA for test data
  array[N_test] real p_pos_test; 
  array[N_test] real pred_val_test; 
  
  for (ii in 1:N_test) {
    
    p_pos_test[ii] = bernoulli_logit_rng(gamma + deltaT * t_test[ii]);
    pred_val_test[ii] = normal_rng(alpha + betaT * t_test[ii], sigma);
    
    log_lik_log[ii] = bernoulli_logit_lpmf(pos_sg_test[ii] | gamma + deltaT * t_test[ii]);
      
    // Only positive values contribute to likelihood on linear component
    if (pos_sg_test[ii] == 1 && val_sg_test[ii] != -9) {
      log_lik_lin[ii_lin] = normal_lpdf(val_sg_test[ii] | alpha + betaT * t_test[ii], 
                                                          sigma);
     ii_lin = ii_lin + 1;                                                 
    }
  }
  
  // Predicting detectability and values of sgRNA for training data
  array[N_train] real p_pos_train; 
  array[N_train] real pred_val_train;
  
  for (ii in 1:N_train) {
    p_pos_train[ii] = bernoulli_logit_rng(gamma + deltaT * t_train[ii]);
    pred_val_train[ii] = normal_rng(alpha + betaT * t_train[ii], sigma);
  }
}

