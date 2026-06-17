functions {
  // NNGP prior on latent w (consistent with Lu Zhang's case study)
  real nngp_w_lpdf(vector w, real sigmasq, real phi,
                   matrix NN_dist, matrix NN_distM, int[,] NN_ind,
                   int N, int M) {
    vector[N] V;
    vector[N] I_Aw = w;
    int dim;
    int h;

    for (i in 2:N) {
      matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M ] iNNdistM;
      matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M ] iNNCholL;
      vector[ i < (M + 1) ? (i - 1) : M ] iNNcorr;
      vector[ i < (M + 1) ? (i - 1) : M ] v;
      row_vector[ i < (M + 1) ? (i - 1) : M ] v2;

      dim = (i < (M + 1)) ? (i - 1) : M;

      // Build neighbor covariance matrix for node i
      if (dim == 1) {
        iNNdistM[1, 1] = 1;
      } else {
        h = 0;
        for (j in 1:(dim - 1)) {
          for (k in (j + 1):dim) {
            h = h + 1;
            iNNdistM[j, k] = exp(-phi * NN_distM[(i - 1), h]);
            iNNdistM[k, j] = iNNdistM[j, k];
          }
        }
        for (j in 1:dim) {
          iNNdistM[j, j] = 1;
        }
      }
      iNNCholL = cholesky_decompose(iNNdistM);
      // Vectorized construction of iNNcorr
      iNNcorr = to_vector(exp(-phi * NN_dist[(i - 1), 1:dim]));

      v = mdivide_left_tri_low(iNNCholL, iNNcorr);
      V[i] = 1 - dot_self(v);
      v2 = mdivide_right_tri_low(v', iNNCholL);

      // Vectorized neighbor update for I_Aw[i]
      I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];
    }
    V[1] = 1;
    return -0.5 * (1 / sigmasq * dot_product(I_Aw, (I_Aw ./ V)) +
                   sum(log(V)) + N * log(sigmasq));
  }
}

data {
  int<lower=1> N; // number of observations/individuals
  int<lower=2> K; // number of ancestral groups
  int<lower=1> PA; // number of predictor A categories/regions
  int<lower=1> PB; // number of predictor B categories/regions)
  int<lower=1> PC; // number of predictor C categories/regions)
  matrix[N, PA] roads; // roads category for each observation (dummy categories)
  matrix[N, PB] rivers; // rivers category for each observation (dummy categories)
  matrix[N, PC] ecoregs; // ecoregs category for each observation (dummy categories)
  simplex[K] Q[N]; // observed probs ancestral groups
  int<lower=1> M; // number of neighbors

  int<lower=1, upper=N> NN_ind[N - 1, M];
  matrix[N - 1, M] NN_dist;
  matrix[N - 1, (M * (M - 1)) / 2] NN_distM;
}

parameters {
  vector[K-1] beta0; // intercept for K-1 groups (1st one becomes baseline)
  matrix[K-1, PA] beta_roads; // coefs for each predictor A and group (i.e. not base)
  matrix[K-1, PB] beta_rivers; // coefs for each predictor B and group (i.e. not base)
  matrix[K-1, PC] beta_ecoregs; // coefs for each predictor C and group (i.e. not base)
  real<lower=0> phi_Dirich; // precision parameter for Dirichlet
  real<lower=0> phi; // for NNGP
  real<lower=0> sigmasq; // spatial process sd (used in NNGP)
  vector[N] w; // spatial random effect (used in NNGP)
}

model {
  // Priors
  beta0 ~ normal(0, 10);
  to_vector(beta_roads) ~ normal(0, 10);
  to_vector(beta_rivers) ~ normal(0, 10);
  to_vector(beta_ecoregs) ~ normal(0, 10);
  phi ~ gamma(1, 1);  
  phi_Dirich ~ lognormal(log(10), 1);
  sigmasq ~ lognormal(log(10), 1);

  // NNGP joint prior for w
  w ~ nngp_w(sigmasq, phi, NN_dist, NN_distM, NN_ind, N, M);

  // Dirichlet regression likelihood with spatial random effect
  for (n in 1:N) {
    vector[K] eta; // linear predictor (for each ancestral group)
    vector[K] mu; // expected proportions
    eta[1] = 0; // first group as baseline
    for (k in 2:K)
      eta[k] = beta0[k-1] + roads[n] * beta_roads[k-1]' + rivers[n] * beta_rivers[k-1]' + ecoregs[n] * beta_ecoregs[k-1]'+ w[n];
    mu = softmax(eta); // convert linear predictor to proportions
    Q[n] ~ dirichlet(mu * phi_Dirich); // likelihood for Dirichlet for observed proportions
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
      eta[k] = beta0[k-1] + roads[n] * beta_roads[k-1]' + rivers[n] * beta_rivers[k-1]' + ecoregs[n] * beta_ecoregs[k-1]'+ w[n];
    mu_pred[n] = softmax(eta);
    Q_pred[n] = dirichlet_rng(mu_pred[n] * phi_Dirich);
    log_lik[n] = dirichlet_lpdf(Q[n] | mu_pred[n] * phi_Dirich);
  }
}
