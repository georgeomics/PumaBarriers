data {
  int<lower=1> N; // number of observations/individuals
  int<lower=2> K; // number of ancestral groups
  int<lower=1> PA; // number of predictor A categories/regions
  int<lower=1> PB; // number of predictor B categories/regions)
  matrix[N, PA] rivers; // rivers category for each observation (dummy categories)
  matrix[N, PB] ecoregs; // ecoregs category for each observation (dummy categories)
  simplex[K] Q[N]; // observed probs ancestral groups
}

parameters {
  vector[K-1] beta0; // intercept for K-1 groups (1st one becomes baseline)
  matrix[K-1, PA] beta_rivers; // coefs for each predictor A and group (i.e. not base)
  matrix[K-1, PB] beta_ecoregs; // coefs for each predictor B and group (i.e. not base)
  real<lower=0> phi; // precision parameter for Dirichlet
}

model {
  // Priors
  beta0 ~ normal(0, 10);
  to_vector(beta_rivers) ~ normal(0, 10);
  to_vector(beta_ecoregs) ~ normal(0, 10);
  phi ~ gamma(1, 1);  

  // Dirichlet regression
  for (n in 1:N) {
    vector[K] eta; // linear predictor (for each ancestral group)
    vector[K] mu; // expected proportions
    eta[1] = 0; // first group as baseline
    for (k in 2:K)
      eta[k] = beta0[k-1] + rivers[n] * beta_rivers[k-1]' + ecoregs[n] * beta_ecoregs[k-1]';
    mu = softmax(eta); // convert linear predictor to proportions
    Q[n] ~ dirichlet(mu * phi); // likelihood for Dirichlet for observed proportions
  }
}

generated quantities {
  vector[K] mu_pred[N];
  simplex[K] Q_pred[N];
  vector[N] log_lik;
  for (n in 1:N) {
    vector[K] eta;
    eta[1] = 0;
    for (k in 2:K)
      eta[k] = beta0[k-1] + rivers[n] * beta_rivers[k-1]' + ecoregs[n] * beta_ecoregs[k-1]';
    mu_pred[n] = softmax(eta);
    Q_pred[n] = dirichlet_rng(mu_pred[n] * phi);
    log_lik[n] = dirichlet_lpdf(Q[n] | mu_pred[n] * phi);
  }
}
