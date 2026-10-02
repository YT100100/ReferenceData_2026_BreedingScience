functions {
  
  vector centering(vector x) {
    return x - mean(x);
  }
  
}

data {
  
  // 基礎情報
  int<lower=1> N; // サンプル数
  int<lower=1> D; // 気象データの日数
  int<lower=1> K; // グループの数
  int<lower=1> M; // Reference-siteの数
  
  // 目的変数
  vector[N] logitCG;
  
  // 説明変数
  matrix[N, D] SR;
  matrix[N, D] RH;
  matrix[N, D] AT26;
  
  // 各データに付随する情報
  int group[N];
  int refsite[N];
  
  // モデル選択
  int<lower=0, upper=1> is_interc_ran; // 切片は変量効果か
  int<lower=0, upper=1> is_MET26_ran ; // MET26の効果は変量効果か
  int<lower=0, upper=1> is_SR_group  ; // SRの傾きはグループごとに異なるか
  int<lower=0, upper=1> is_RH_group  ; // RHの傾きはグループごとに異なるか
  int<lower=0, upper=1> is_mat_ested ; // 登熟期間を推定するか
  int<lower=0, upper=1> is_mat_group ; // 登熟期間はグループごとに異なるか
  
}

transformed data {
  
  // 播種後日数のベクトル
  vector[D] DAF = linspaced_vector(D, 1, D);

}

parameters {
  
  // 変量効果のベクトル
  vector[K] intercept;
  vector[K] eff_MET26;
  vector[K] eff_SR_group;
  vector[K] eff_RH_group;
  vector[M] eff_refsite;
  
  // 変量効果の平均
  real mu_intercept;
  real mu_eff_MET26;
  real mu_eff_SR;
  real mu_eff_RH;
  
  // 変量効果の分散
  real<lower=0> sigma_intercept;
  real<lower=0> sigma_eff_MET26;
  real<lower=0> sigma_eff_SR;
  real<lower=0> sigma_eff_RH;
  real<lower=0> sigma_eff_refsite;
  
  // 固定効果
  real eff_SR_common;
  real eff_RH_common;
  
  // 登熟期間
  vector<lower=0, upper=40>[K] maturity_center_group;
  real mu_maturity_center;
  real<lower=0> sigma_maturity_center;
  real<lower=0, upper=40> maturity_center_common;

  // 気象要因の重み関数の傾き
  real<lower=0, upper=20> maturity_exp;
  
  // 残差の分散
  real<lower=0> sigma;
  
}

transformed parameters {
  
  // グループごとに象要因への重みを計算する
  vector[D] weight[K];
  real sum_weight[K];
  for (k in 1:K) {
    
    real maturity_center_k;
    real maturity_exp_k;
    
    // モデルに応じてパラメータを選択する
    if (is_mat_ested) {
      
      maturity_exp_k    = maturity_exp;
      if (is_mat_group) {
        maturity_center_k = maturity_center_group[k];
      } else {
        maturity_center_k = maturity_center_common;
      }
      
    } else {
      
      maturity_center_k = 20.5;
      maturity_exp_k    = 0.01;
      
    }
    
    // 重みを計算する
    weight[k] = 1 / (1 + exp((DAF - maturity_center_k) / maturity_exp_k));
    sum_weight[k] = sum(weight[k]);
    
  }
  
  // 気象要因の説明変数の計算
  vector[N] MET26_wm;
  vector[N] SR_wm;
  vector[N] RH_wm;
  for (n in 1:N) {
    int k = group[n];
    MET26_wm[n] = AT26[n, ] * weight[k] / sum_weight[k];
    SR_wm   [n] = SR  [n, ] * weight[k] / sum_weight[k];
    RH_wm   [n] = RH  [n, ] * weight[k] / sum_weight[k];
  }
  
  // 気象要因のセンタリング
  vector[N] MET26_wms = centering(MET26_wm);
  vector[N] SR_wms    = centering(SR_wm);
  vector[N] RH_wms    = centering(RH_wm);
  
  // SR, RHの効果の選択
  vector[K] eff_SR;
  vector[K] eff_RH;
  for (k in 1:K) {
    eff_SR[k] = is_SR_group ? eff_SR_group[k] : eff_SR_common;
    eff_RH[k] = is_RH_group ? eff_RH_group[k] : eff_RH_common;
  }
  
  // logitCGの平均
  vector[N] mu;
  for (n in 1:N) {
    
    int k = group[n]; // サンプルnの品種
    int m = refsite [n]; // サンプルnのReference-site
    
    // logitCGに対する環境の効果
    real met_eff_n = 0;
    met_eff_n += eff_MET26[k] * MET26_wms[n];
    met_eff_n += eff_SR   [k] * SR_wms   [n];
    met_eff_n += eff_RH   [k] * RH_wms   [n];
    
    // logitCGの平均
    mu[n] = intercept[k] + met_eff_n + eff_refsite[m];

  }
  
}

model {
  
  // 変量効果 - 切片、MET26、文献場所
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
  
  // 変量効果 - SR
  eff_SR_group ~ normal(mu_eff_SR, sigma_eff_SR);
  if (is_SR_group) {
    eff_SR_common ~ normal(0, 1);
  } else {
    mu_eff_SR    ~ normal(0, 1);
    sigma_eff_SR ~ normal(10, 1);
  }
  
  // 変量効果 - RH
  eff_RH_group ~ normal(mu_eff_RH, sigma_eff_RH);
  if (is_RH_group) {
    eff_RH_common ~ normal(0, 1);
  } else {
    mu_eff_RH    ~ normal(0, 1);
    sigma_eff_RH ~ normal(10, 1);
  }
  
  // 登熟期間
  maturity_center_group ~ normal(mu_maturity_center, sigma_maturity_center);
  if (is_mat_group || (!is_mat_ested)) {
    maturity_center_common ~ normal(20, 1);
  }
  if (!is_mat_group) {
    mu_maturity_center    ~ normal(20, 1);
    sigma_maturity_center ~ normal(10, 1);
  }
  
  // データ
  logitCG ~ normal(mu, sigma);
  
}

generated quantities {
  
  // 気象要因の平均値
  real mean_MET26 = mean(MET26_wm);
  real mean_SR    = mean(SR_wm);
  real mean_RH    = mean(RH_wm);
  
  // データごとの対数尤度
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(logitCG[n] | mu[n], sigma);
  }
  
}
