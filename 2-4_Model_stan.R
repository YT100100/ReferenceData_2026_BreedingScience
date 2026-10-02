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
  '2_ModelWithPubData', '2-4_Model_stan', 
  paste0(Sys.Date(), '_ver.2.011.4_DataID'))
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

# 白未熟粒率データ
path_dat <- file.path(
  '..', 'Data', '2-2_Paper', '2024_Wakatsuki', 'CleanData', 
  '2025-12-02_Mutsuhomare', 'Data_Main.rds')

# 気象データ
path_met <- file.path(
  '1_DataCheck', '1-1_IshigookaMeshData', '2025-02-25_ver.2.002.1_Init', 
  'meshdat.rds')

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

# データを読み込み、整形
dat_raw <- readRDS(path_dat)
dat <- dat_raw$df[, c('ID', 'HTR', 'Cultivar', 'Pref', 'RefNo', 'CG')]
dat <- within(dat, {
  
  # ref-siteを作成
  RefSite  <- sprintf('%02i_%s', RefNo, Pref)
  
  # 白未熟粒率をロジット変換
  logitCG <- logit(CG / 100)
  
  # 因子化
  HTR      <- as.factor(HTR)
  Cultivar <- as.factor(Cultivar)
  RefSite  <- as.factor(RefSite)
  
})

# 気象データを読み込み
met <- readRDS(path_met)

# データ行数と品種数を確認する
cat(sprintf('nrow: %i, n. cultivar: %i, n. refsite: %i\n', 
            nrow(dat), nlevels(dat$Cultivar), nlevels(dat$RefSite)))

# NAを削除する
sel <- !is.na(dat$CG)
sel <- sel & apply(met, 2, function(x) all(!is.na(x)))
dat <- dat[sel, ]
met <- met[, sel, ]
dat$Cultivar <- droplevels(dat$Cultivar)
dat$Refsite  <- droplevels(dat$RefSite)

# データ行数と品種数を確認する
cat(sprintf('nrow: %i, n. cultivar: %i, n. refsite: %i\n', 
            nrow(dat), nlevels(dat$Cultivar), nlevels(dat$RefSite)))

# データを保存しておく
write.csv(dat, file.path(path_out, 'Data_CG.csv'))
saveRDS(dat, file.path(path_out, 'Data_CG.rds'))
saveRDS(met, file.path(path_out, 'Data_Met.rds'))

# stanモデルを作成しておく
stan_mod <- stan_model(path_stan)



#### 4. データの準備 ####

# 基本データを準備する
standat <- with(dat, list(
  
  # 基礎情報
  N = nrow(dat),         # サンプル数
  D = dim(met)[1],       # 気象データの日数
  K = nlevels(Cultivar), # 品種数
  M = nlevels(RefSite),  # Reference-siteの数
  
  # 目的変数
  logitCG = logitCG,
  
  # 説明変数
  AT26 = t(pmax(met[, , 'tmp'] - 26, 0)),
  SR = t(met[, , 'srd']),
  RH = t(met[, , 'rhu']),
  
  # 各データに付随する情報
  group = as.integer(Cultivar),
  refsite  = as.integer(RefSite),
  
  # モデル構造を定義
  is_interc_ran = TRUE,
  is_MET26_ran  = TRUE

))

# モデルのパターンを生成する
pat <- expand.grid(
  is_SR_group = c(TRUE, FALSE),
  is_RH_group = c(TRUE, FALSE),
  maturity_mode  = c('20days', 'fixed', 'cultivar'),
  stringsAsFactors = FALSE)
pat <- within(pat, {
  is_mat_ested <- maturity_mode != '20days'
  is_mat_group <- maturity_mode == 'cultivar'
})

# モデルの名前を作成する
pat <- within(pat, {
  modname <- sprintf(
    'SRcult=%s_RHcult=%s_mat=%s',
    is_SR_group, is_RH_group, maturity_mode)
})



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

tc <- set_loop_timechecker(nrow(pat))
out <- lapply(seq_len(nrow(pat)), function(i) {
  # i <- 1
  
  # フォルダをを作成する
  filename_i <- paste0(
    path_out, .Platform$file.sep, pat[i, 'modname'], .Platform$file.sep)
  dir.create(filename_i, showWarnings = FALSE)
  
  # モデルに応じたデータを作成する
  standat_i <- c(standat, list(
    is_SR_group  = pat[i, 'is_SR_group' ],
    is_RH_group  = pat[i, 'is_RH_group' ],
    is_mat_ested = pat[i, 'is_mat_ested'],
    is_mat_group = pat[i, 'is_mat_group']))
  
  # stanを用いてパラメータを推定する
  withCallingHandlers(
    {
      stan_res_i <- sampling(
        stan_mod, data = standat_i, 
        chains = N_CHAIN, warmup = N_WARMUP, iter = N_ITER, thin = N_THIN)
    },
    warning = function(w) save_warning(w, paste0(filename_i, '01_stan_warn.txt'))
  )
  
  # stanの結果を保存しておく
  saveRDS(stan_res_i, paste0(filename_i, '01_stan.rds'))
  
  # summaryを保存しておく
  stan_summary_i <- summary(stan_res_i, pars = PARS_SINGLE)
  write.csv(stan_summary_i$summary, paste0(filename_i, '01_stan_summary.csv'))
  for (j in seq_len(dim(stan_summary_i$c_summary)[3])) {
    filename_j <- sprintf('%s01_stan_summary_chain%02i.csv', filename_i, j)
    write.csv(stan_summary_i$c_summary[, , j], filename_j)
  }
  
  # 収束したか確認する
  pdf(paste0(filename_i, '02_Fig_convergence.pdf'), width = 10)
  for (parname_j in PARS_SINGLE) {
    g_list <- lapply(c(TRUE, FALSE), function(inc_warmup) {
      stan_trace(stan_res_i, pars = parname_j, inc_warmup = inc_warmup) +
        labs(title = parname_j) +
        theme(legend.position = 'bottom')
    })
    do.call(grid.arrange, c(g_list, list(nrow = 1)))
  }
  for (parname_j in PARS_MULTIPLE) {
    g <- stan_trace(stan_res_i, pars = parname_j, inc_warmup = FALSE) +
      labs(title = parname_j)
    print(g)
  }
  dev.off()
  
  # pairsプロットで変数間相関を確認する
  pars_pairs <- PARS_SINGLE
  pars_pairs <- with(standat_i, {
    
    if (is_SR_group) {
      pars_pairs <- pars_pairs[!grepl('^eff_SR_common$', pars_pairs)]
    } else {
      pars_pairs <- pars_pairs[!grepl('_eff_SR$', pars_pairs)]
    }
    
    if (is_RH_group) {
      pars_pairs <- pars_pairs[!grepl('^eff_RH_common$', pars_pairs)]
    } else {
      pars_pairs <- pars_pairs[!grepl('_eff_RH$', pars_pairs)]
    }
    
    if (is_mat_ested) {
      if (is_mat_group) {
        pars_pairs <- pars_pairs[!grepl('^maturity_center_common$', pars_pairs)]
      } else {
        pars_pairs <- pars_pairs[!grepl('_maturity_', pars_pairs)]
      }
    } else {
      pars_pairs <- pars_pairs[!grepl('maturity', pars_pairs)]
    }
    
  })
  png(paste0(filename_i, '03_Fig_pairs.png'), 
      width = 10, height = 10, unit = 'in', res = 300)
  pairs(stan_res_i, pars = pars_pairs, cex.labels = 0.7)
  dev.off()
  
  # MCMCサンプルを抽出する
  stan_res_list_i <- rstan::extract(stan_res_i)
  
  # パラメータの要約を保存する（全データ共通）
  summary_single_i <- sapply(stan_res_list_i[PARS_SINGLE], function(vec) {
    quantile(vec, probs = PROBS)
  })
  write.csv(summary_single_i, paste0(filename_i, '03_Res_Summary_Common.csv'))
  saveRDS(summary_single_i, paste0(filename_i, '03_Res_Summary_Common.rds'))
  
  # パラメータの要約を保存する（複数あるもの）
  summary_multiple_i <- lapply(PARS_MULTIPLE, function(parname) {
    
    mat <- stan_res_list_i[[parname]]
    res <- apply(mat, 2, quantile, probs = PROBS)
    if (parname == 'eff_refsite') {
      colnames(res) <- levels(dat$RefSite)
    } else {
      colnames(res) <- levels(dat$Cultivar)
    }
    write.csv(res, paste0(filename_i, '04_Res_Summary_', parname, '.csv'))
    res
    
  })
  names(summary_multiple_i) <- PARS_MULTIPLE
  saveRDS(summary_multiple_i, paste0(filename_i, '04_Res_Summary_Multiple.rds'))
  
  # RMSEを計算する
  est_logitcg_i <- colMeans(stan_res_list_i$mu)
  rmse_logitcg_i <- rmse(dat$logitCG, est_logitcg_i)
  rmse_cg_i <- rmse(inv_logit(dat$logitCG), inv_logit(est_logitcg_i)) * 100
  
  # WAIC, ELPDを計算する
  log_lik_i <- extract_log_lik(stan_res_i, merge_chains = FALSE)
  withCallingHandlers(
    {waic_i <- waic(log_lik_i)},
    warning = function(w) save_warning(w, paste0(filename_i, '05_WAIC_warn.txt')))
  withCallingHandlers(
    {elpd_i <- loo(log_lik_i)},
    warning = function(w) save_warning(w, paste0(filename_i, '05_ELPD_warn.txt')))
  stan_eval_i <- list(WAIC = waic_i, ELPD = elpd_i)
  saveRDS(stan_eval_i, paste0(filename_i, '05_Res_WAIC.rds'))
  
  # WAIC, ELPD, RMSEをcsvファイルに保存する
  stan_eval_csv_i <- rbind(
    stan_eval_i$WAIC$estimates, stan_eval_i$ELPD$estimates, 
    cbind(c(RMSE_logitCG = rmse_logitcg_i, RMSE_CG = rmse_cg_i), NA))
  write.csv(stan_eval_csv_i, paste0(filename_i, '05_Res_WAIC.csv'))
  
  tc()
  
})
