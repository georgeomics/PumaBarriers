
data {
  int<lower=0> N; // Number of observations
  int<lower=1> K; // Number of ancestral groups
  int<lower=1> E; // Number of road categories
  int<lower=1, upper=E> ecoreg[N]; // ecoreg category for each observation
  matrix[N, N] Z; // Pre-calculated Z matrix for each individual (NxN)
  simplex[K] Q[N]; // Observed probabilities for each ancestral group
}

parameters {
  matrix[K, E] beta_ecoreg; // Coefficients for each road category and ancestral group
  real<lower=0> gamma; // Scaling parameter for the Gaussian decay
  real<lower=0> epsilon;
}

transformed parameters {
  vector[N] phi; // Precision parameter for each individual
  matrix[N, K] alpha; // Dirichlet parameters
  vector[K] mu;

  // Calculate the precision parameter phi for each individual
  for (n in 1:N) {
    phi[n] = exp(gamma * sum(Z[n]) + epsilon);
  }
  // Calculate alpha for each observation
  for (n in 1:N) {
    vector[K] logit_mu = rep_vector(0.0, K);

    for (k in 1:K) {
        logit_mu[k] = beta_ecoreg[k, 1] * (ecoreg[n] == 1) + // (ecoreg[n] == r) is a logical check that evaluates to 1 if true and 0 if false
                      beta_ecoreg[k, 2] * (ecoreg[n] == 2) + 
                      beta_ecoreg[k, 3] * (ecoreg[n] == 3) + 
                      beta_ecoreg[k, 4] * (ecoreg[n] == 4) + 
                      beta_ecoreg[k, 5] * (ecoreg[n] == 5);        
    }
      
    // use multivariate logit/softmax to compute mu
    mu = softmax(logit_mu);

    // alpha for the Dirichlet distribution
    alpha[n] = to_row_vector(mu * phi[n]); 
  }
}

model {
  // Priors
  to_vector(beta_ecoreg) ~ normal(0, 1); // Prior for beta coefficients
  gamma ~ normal(0, 1); // Prior for gamma
  epsilon ~ normal(0, 1); // Prior for epsilon
  // Likelihood
  for (n in 1:N) {
    Q[n] ~ dirichlet(alpha[n]); // Dirichlet likelihood
  }
}

generated quantities {
  vector[N] log_lik;
  simplex[K] Q_pred[N];
  for (n in 1:N) {
    log_lik[n] = dirichlet_lpdf(Q[n] | alpha[n]);
    Q_pred[n] = dirichlet_rng(to_vector(alpha[n])); // Convert row_vector to vector for Dirichlet RNG
  }
}
