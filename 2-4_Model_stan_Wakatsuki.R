# パラメータ推定におけるMCMCサンプリング
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# MCMCサンプルを用いたWAICの計算
library(loo)

# RMSEの計算
library(TodaFunc)

# 可視化
library(gridExtra)
library(ggplot2)

# 進捗の表示
library(timechecker)

# ロジット変換・逆ロジット変換
source(file.path('Script', '0_CommonFunctions', '0-1-1_Logit.R'))



#### 1. パスの設定 ####

# 結果の保存先
path_out <- file.path(
  '2_ModelWithPubData', '2-4_Model_stan_Wakatsuki', 
  paste0(Sys.Date(), '_ver.2.010.2_Init'))
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

# 白未熟粒率データ
path_dat <- file.path(
  '2_ModelWithPubData', '2-4_Model_stan', 
  '2025-12-04_ver.2.010.2_AllModInOneStan_LittleModify', 'Data_CG.rds')

# 気象データ
path_met <- file.path(
  '2_ModelWithPubData', '2-4_Model_stan', 
  '2025-12-04_ver.2.010.2_AllModInOneStan_LittleModify', 'Data_Met.rds')

# stanスクリプト
path_stan <- file.path('Script', '2_ModelWithPubData', '2-4_Model_stan.stan')




#### 2. パラメーターの設定 ####

# MCMCの設定
N_CHAIN <- 4
N_ITER <- 4500
N_WARMUP <- 500
N_THIN <- 2
# N_CHAIN <- 2; N_ITER <- 100; N_WARMUP <- 20

# 結果を要約するときの分位点
PROBS <- c(0.025, 0.05, 0.5, 0.95, 0.975)

# パラメータ名（データ間で共通のもの）
PARS_SINGLE <- c(
  'mu_intercept', 'sigma_intercept', 
  'mu_eff_MET26', 'sigma_eff_MET26', 
  'mu_eff_SR'   , 'sigma_eff_SR'   , 'eff_SR_common', 
  'mu_eff_RH'   , 'sigma_eff_RH'   , 'eff_RH_common', 
  'sigma_eff_refsite', 'sigma',
  'mu_maturity_center', 'maturity_center_common', 'sigma_maturity_center', 
  'maturity_exp')

# パラメータ名（データ内で複数あるもの）
PARS_MULTIPLE <- c(
  'intercept', 'eff_MET26', 'eff_RH_group', 'eff_SR_group', 
  'maturity_center_group', 'eff_refsite')



#### 3. データの読み込み ####

# 白未熟粒率データを読み込み
dat <- readRDS(path_dat)

# 気象データを読み込み
met <- readRDS(path_met)

# 高温登熟性を因子化する
dat$HTR <- as.factor(paste0('HTR', dat$HTR))

# stanモデルを作成しておく
stan_mod <- stan_model(path_stan)



#### 4. データの準備 ####

# 基本データを準備する
standat <- with(dat, list(
  
  # 基礎情報
  N = nrow(dat),         # サンプル数
  D = dim(met)[1],       # 気象データの日数
  K = nlevels(HTR),      # 品種数（ここでは高温登熟性ランクの数）
  M = nlevels(RefSite),  # Reference-siteの数
  
  # 目的変数
  logitCG = logitCG,
  
  # 説明変数
  AT26 = t(pmax(met[, , 'tmp'] - 26, 0)),
  SR = t(met[, , 'srd']),
  RH = t(met[, , 'rhu']),
  
  # 各データに付随する情報
  group = as.integer(HTR),
  refsite = as.integer(RefSite),

  # モデル構造を定義
  is_interc_ran = FALSE,
  is_MET26_ran  = FALSE,
  is_SR_group   = FALSE,
  is_RH_group   = FALSE,
  is_mat_ested  = FALSE,
  is_mat_group  = FALSE
  
))



#### 5. 推定 ####

# stanの警告を保存する関数
save_warning <- function(w, file) {
  
  # ログ行を作成（時刻・関数名・メッセージ）
  log_line <- sprintf(
    '%s | WARNING | call=%s\n%s\n',
    format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
    deparse(conditionCall(w)),
    conditionMessage(w))
  
  # 追記モードで書き込み
  cat(log_line, '\n', file = file, append = TRUE)
  
  # 以降の警告出力をコンソールに出したくないので有効化
  invokeRestart('muffleWarning')
  
}

# フォルダをを作成する
filename <- paste0(path_out, .Platform$file.sep)

# stanを用いてパラメータを推定する
withCallingHandlers(
  {
    stan_res <- sampling(
      stan_mod, data = standat, 
      chains = N_CHAIN, warmup = N_WARMUP, iter = N_ITER, thin = N_THIN)
  },
  warning = function(w) save_warning(w, paste0(filename, '01_stan_warn.txt'))
)

# stanの結果を保存しておく
saveRDS(stan_res, paste0(filename, '01_stan.rds'))

# summaryを保存しておく
stan_summary <- summary(stan_res, pars = PARS_SINGLE)
write.csv(stan_summary$summary, paste0(filename, '01_stan_summary.csv'))
for (j in seq_len(dim(stan_summary$c_summary)[3])) {
  filename_j <- sprintf('%s01_stan_summary_chain%02i.csv', filename, j)
  write.csv(stan_summary$c_summary[, , j], filename_j)
}

# 収束したか確認する
pdf(paste0(filename, '02_Fig_convergence.pdf'), width = 10)
for (parname_j in PARS_SINGLE) {
  g_list <- lapply(c(TRUE, FALSE), function(inc_warmup) {
    stan_trace(stan_res, pars = parname_j, inc_warmup = inc_warmup) +
      labs(title = parname_j) +
      theme(legend.position = 'bottom')
  })
  do.call(grid.arrange, c(g_list, list(nrow = 1)))
}
for (parname_j in PARS_MULTIPLE) {
  g <- stan_trace(stan_res, pars = parname_j, inc_warmup = FALSE) +
    labs(title = parname_j)
  print(g)
}
dev.off()

# pairsプロットで変数間相関を確認する
pars_pairs <- PARS_SINGLE
pars_pairs <- pars_pairs[!grepl('_eff_SR$', pars_pairs)]
pars_pairs <- pars_pairs[!grepl('_eff_RH$', pars_pairs)]
pars_pairs <- pars_pairs[!grepl('maturity', pars_pairs)]
png(paste0(filename, '03_Fig_pairs.png'), 
    width = 10, height = 10, unit = 'in', res = 300)
pairs(stan_res, pars = pars_pairs, cex.labels = 0.7)
dev.off()

# MCMCサンプルを抽出する
stan_res_list <- rstan::extract(stan_res)

# パラメータの要約を保存する（全データ共通）
summary_single <- sapply(stan_res_list[PARS_SINGLE], function(vec) {
  quantile(vec, probs = PROBS)
})
write.csv(summary_single, paste0(filename, '03_Res_Summary_Common.csv'))
saveRDS(summary_single, paste0(filename, '03_Res_Summary_Common.rds'))

# パラメータの要約を保存する（複数あるもの）
summary_multiple <- lapply(PARS_MULTIPLE, function(parname) {
  
  mat <- stan_res_list[[parname]]
  res <- apply(mat, 2, quantile, probs = PROBS)
  if (parname == 'eff_refsite') {
    colnames(res) <- levels(dat$RefSite)
  } else {
    colnames(res) <- levels(dat$HTR)
  }
  write.csv(res, paste0(filename, '04_Res_Summary_', parname, '.csv'))
  res
  
})
names(summary_multiple) <- PARS_MULTIPLE
saveRDS(summary_multiple, paste0(filename, '04_Res_Summary_Multiple.rds'))

# RMSEを計算する
est_logitcg <- colMeans(stan_res_list$mu)
rmse_logitcg <- rmse(dat$logitCG, est_logitcg)
rmse_cg <- rmse(inv_logit(dat$logitCG), inv_logit(est_logitcg)) * 100

# WAIC, ELPDを計算する
log_lik <- extract_log_lik(stan_res, merge_chains = FALSE)
withCallingHandlers(
  {waic <- waic(log_lik)},
  warning = function(w) save_warning(w, paste0(filename, '05_WAIC_warn.txt')))
withCallingHandlers(
  {elpd <- loo(log_lik)},
  warning = function(w) save_warning(w, paste0(filename, '05_ELPD_warn.txt')))
stan_eval <- list(WAIC = waic, ELPD = elpd)
saveRDS(stan_eval, paste0(filename, '05_Res_WAIC.rds'))

# WAIC, ELPD, RMSEをcsvファイルに保存する
stan_eval_csv <- rbind(
  stan_eval$WAIC$estimates, stan_eval$ELPD$estimates, 
  cbind(c(RMSE_logitCG = rmse_logitcg, RMSE_CG = rmse_cg), NA))
write.csv(stan_eval_csv, paste0(filename, '05_Res_WAIC.csv'))
