data {
  int<lower=0> N; // Number of observations
  int<lower=1> K; // Number of ancestral groups, K = 4 for Q1, Q2, Q3, Q4
  int<lower=1> R; // Number of rivers
  int<lower=1> I; // Number of roads
  int<lower=1> E; // Number of ecoreg
  int<lower=1, upper=I> rivers[N];
  int<lower=1, upper=I> roads[N]; // Ecoregion category for each observation
  int<lower=1, upper=E> ecoreg[N]; // River category for each observation
  matrix[N, N] Z; // Pre-calculated Z matrix for each individual (NxN)
  simplex[K] Q[N]; // Observed probabilities for each ancestral group
}

parameters {
  matrix[K, R] beta_river;
  matrix[K, I] beta_roads; // Coefficients for each road category and ancestral group
  matrix[K, E] beta_ecoreg; // Coefficients for each road category and ancestral group
  real<lower=0> gamma; // Scaling parameter for the Gaussian decay
  real<lower=0> epsilon;
}

transformed parameters {
  vector[N] phi; // Precision parameter for each individual
  vector[K] mu;
  matrix[N, K] alpha; // Dirichlet parameters

  // Calculate the precision parameter phi for each individual
  for (n in 1:N) {
    phi[n] = exp(gamma * sum(Z[n]));
  }
  // Calculate alpha for each observation
  for (n in 1:N) {
      vector[K] logit_mu = rep_vector(0.0, K);

      for (k in 1:K) {
        // Calculate logit_mu[k] based on river categories
        logit_mu[k] = beta_river[k, 1] * (rivers[n] == 1) +
                      beta_river[k, 2] * (rivers[n] == 2) +
                      beta_river[k, 3] * (rivers[n] == 3) +
                      beta_river[k, 4] * (rivers[n] == 4) +
                      beta_river[k, 5] * (rivers[n] == 5);

        // Add effect of roads categories
        logit_mu[k] += beta_roads[k, 1] * (roads[n] == 1) +
                      beta_roads[k, 2] * (roads[n] == 2) +
                      beta_roads[k, 3] * (roads[n] == 3) +
                      beta_roads[k, 4] * (roads[n] == 4) +
                      beta_roads[k, 5] * (roads[n] == 5) +
                      beta_roads[k, 6] * (roads[n] == 6);

        // Add effect of ecoregion categories
        logit_mu[k] += beta_ecoreg[k, 1] * (ecoreg[n] == 1) +
                      beta_ecoreg[k, 2] * (ecoreg[n] == 2) +
                      beta_ecoreg[k, 3] * (ecoreg[n] == 3) +
                      beta_ecoreg[k, 4] * (ecoreg[n] == 4) +
                      beta_ecoreg[k, 5] * (ecoreg[n] == 5);


        mu[k] = inv_logit(logit_mu[k]); 
      }

    /// use multivariate logit/softmax to compute mu
    mu = softmax(logit_mu);

    // alpha for the Dirichlet distribution
    alpha[n] = to_row_vector(mu * phi[n]);
  }
}

model {
  // Priors
  to_vector(beta_river) ~ normal(0, 1);
  to_vector(beta_roads) ~ normal(0, 1); // Prior for beta coefficients
  to_vector(beta_ecoreg) ~ normal(0, 1);
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







