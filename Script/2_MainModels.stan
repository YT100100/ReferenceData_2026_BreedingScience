functions {
  
  vector centering(vector x) {
    return x - mean(x);
  }
  
}

data {
  
  // Basic dimensions
  int<lower=1> N; // Number of observations
  int<lower=1> D; // Number of days in the weather data
  int<lower=1> K; // Number of groups
  int<lower=1> M; // Number of reference-site combinations
  
  // Response variable
  vector[N] logitCG;
  
  // Explanatory variables
  matrix[N, D] SR;
  matrix[N, D] RH;
  matrix[N, D] AT26;
  
  // Grouping information associated with each observation
  int group[N];
  int refsite[N];
  
  // Model specification
  int<lower=0, upper=1> is_interc_ran; // Whether the intercept is modeled as a random effect
  int<lower=0, upper=1> is_MET26_ran ; // Whether the MET26 effect is modeled as a random effect
  int<lower=0, upper=1> is_SR_group  ; // Whether the SR slope differs among groups
  int<lower=0, upper=1> is_RH_group  ; // Whether the RH slope differs among groups
  int<lower=0, upper=1> is_mat_ested ; // Whether the grain-filling period is estimated
  int<lower=0, upper=1> is_mat_group ; // Whether the grain-filling period differs among groups
  
}

transformed data {
  
  // Vector of days after flowering
  vector[D] DAF = linspaced_vector(D, 1, D);

}

parameters {
  
  // Vectors of ramdom effects
  vector[K] intercept;
  vector[K] eff_MET26;
  vector[K] eff_SR_group;
  vector[K] eff_RH_group;
  vector[M] eff_refsite;
  
  // Means of the ramdom effects
  real mu_intercept;
  real mu_eff_MET26;
  real mu_eff_SR;
  real mu_eff_RH;
  
  // Standard deviations of the ramdom effects
  real<lower=0> sigma_intercept;
  real<lower=0> sigma_eff_MET26;
  real<lower=0> sigma_eff_SR;
  real<lower=0> sigma_eff_RH;
  real<lower=0> sigma_eff_refsite;
  
  // Fixed effects
  real eff_SR_common;
  real eff_RH_common;
  
  // Grain-filling period
  vector<lower=0, upper=40>[K] maturity_center_group;
  real mu_maturity_center;
  real<lower=0> sigma_maturity_center;
  real<lower=0, upper=40> maturity_center_common;

  // Steepness parameter of the weighting function for weather variables
  real<lower=0, upper=20> maturity_exp;
  
  // Residual standard deviation
  real<lower=0> sigma;
  
}

transformed parameters {
  
  // Calculate group-specific weights for the weather variables
  vector[D] weight[K];
  real sum_weight[K];
  for (k in 1:K) {
    
    real maturity_center_k;
    real maturity_exp_k;
    
    // Select parameters according to the model specification
    if (is_mat_ested) {
      
      maturity_exp_k = maturity_exp;
      if (is_mat_group) {
        maturity_center_k = maturity_center_group[k];
      } else {
        maturity_center_k = maturity_center_common;
      }
      
    } else {
      
      maturity_center_k = 20.5;
      maturity_exp_k = 0.01;
      
    }
    
    // Calculate the weights
    weight[k] = 1 / (1 + exp((DAF - maturity_center_k) / maturity_exp_k));
    sum_weight[k] = sum(weight[k]);
    
  }
  
  // Calculate weighted means of the weather variables
  vector[N] MET26_wm;
  vector[N] SR_wm;
  vector[N] RH_wm;
  for (n in 1:N) {
    int k = group[n];
    MET26_wm[n] = AT26[n, ] * weight[k] / sum_weight[k];
    SR_wm   [n] = SR  [n, ] * weight[k] / sum_weight[k];
    RH_wm   [n] = RH  [n, ] * weight[k] / sum_weight[k];
  }
  
  // Center the weighted weather variables
  vector[N] MET26_wms = centering(MET26_wm);
  vector[N] SR_wms    = centering(SR_wm);
  vector[N] RH_wms    = centering(RH_wm);
  
  // Select group-specific or common effects of SR and RH
  vector[K] eff_SR;
  vector[K] eff_RH;
  for (k in 1:K) {
    eff_SR[k] = is_SR_group ? eff_SR_group[k] : eff_SR_common;
    eff_RH[k] = is_RH_group ? eff_RH_group[k] : eff_RH_common;
  }
  
  // Expected value of logitCG
  vector[N] mu;
  for (n in 1:N) {
    
    int k = group[n];   // Group associated with observation n
    int m = refsite[n]; // Reference-site combination associated with observation n
    
    // Effects of the environmental variables on logitCG
    real met_eff_n = 0;
    met_eff_n += eff_MET26[k] * MET26_wms[n];
    met_eff_n += eff_SR   [k] * SR_wms   [n];
    met_eff_n += eff_RH   [k] * RH_wms   [n];
    
    // Expected value of logitCG
    mu[n] = intercept[k] + met_eff_n + eff_refsite[m];

  }
  
}

model {
  
  // Random effects: intercept, MET26, and reference-site
  if (is_interc_ran) {
    intercept ~ normal(mu_intercept, sigma_intercept);
  } else {
    mu_intercept ~ normal(0, 1);
    sigma_intercept ~ normal(10, 1);
  }

  if (is_MET26_ran) {
    eff_MET26 ~ normal(mu_eff_MET26, sigma_eff_MET26);
  } else {
    mu_eff_MET26 ~ normal(0, 1);
    sigma_eff_MET26 ~ normal(10, 1);
  }

  eff_refsite ~ normal(0, sigma_eff_refsite);
  
  // Random effects of SR
  eff_SR_group ~ normal(mu_eff_SR, sigma_eff_SR);
  if (is_SR_group) {
    eff_SR_common ~ normal(0, 1);
  } else {
    mu_eff_SR    ~ normal(0, 1);
    sigma_eff_SR ~ normal(10, 1);
  }
  
  // Random effects of RH
  eff_RH_group ~ normal(mu_eff_RH, sigma_eff_RH);
  if (is_RH_group) {
    eff_RH_common ~ normal(0, 1);
  } else {
    mu_eff_RH    ~ normal(0, 1);
    sigma_eff_RH ~ normal(10, 1);
  }
  
  // Grain-filling period
  maturity_center_group ~ normal(mu_maturity_center, sigma_maturity_center);
  if (is_mat_group || (!is_mat_ested)) {
    maturity_center_common ~ normal(20, 1);
  }

  if (!is_mat_group) {
    mu_maturity_center    ~ normal(20, 1);
    sigma_maturity_center ~ normal(10, 1);
  }
  
  // Likelihood
  logitCG ~ normal(mu, sigma);
  
}

generated quantities {
  
  // Means of the weighted weather variables
  real mean_MET26 = mean(MET26_wm);
  real mean_SR    = mean(SR_wm);
  real mean_RH    = mean(RH_wm);
  
  // Pointwise log-likelihood
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(logitCG[n] | mu[n], sigma);
  }
  
}
