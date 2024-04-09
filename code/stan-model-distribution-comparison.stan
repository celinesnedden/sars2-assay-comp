//

data {
  int<lower=0> N;
  int<lower=0> N_types;
  array[N] real y;
  array[N] int y_type;
}


parameters {
  array[N_types] real mu;
  array[N_types] real<lower=0> sigma;
}

model {
  for (n in 1:N) {
      target += normal_lpdf(y[n] | mu[y_type[n]],
                                   sigma[y_type[n]]);
  }
}

