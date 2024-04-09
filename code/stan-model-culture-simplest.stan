//
// This Stan program defines a logistic regression for whether culture is 
//   predicted to be censored or not

data {

  // Number of observations
  int<lower=0> N; // number of observations in training data
  
  // Outcome variables
  array[N] int<lower=0, upper=1> pos_inf;

  // All predictors
  array[N] real t; // primary predictor, total RNA value

}

parameters {
  
  real gamma;
  real psiT;

}


model {
  
  for (n in 1:N) {
      target += bernoulli_logit_lpmf(pos_inf[n] | gamma + psiT * t[n]);
  }
  
  // Priors: Non-Informative
  //gamma ~ normal(0, 1);
  //psiT ~ normal(0, 1);
  
  // Priors: Informative
  gamma ~ normal(-1, 1); 
  psiT ~ normal(1, 0.5);

}

generated quantities {
  
  // Posterior predictions
  array[N] real p_pos;
  for (ii in 1:N) {
    p_pos[ii] = bernoulli_logit_rng(gamma + psiT * t[ii]);                                                     
    }
}
