# 
# Fitting a model proposed by Wakatsuki et al. (2024)
# for chalky grain proportions of rice
# 
# Author: Yusuke Toda (yusuke0320.research@gmail.com)
# 

# MCMC sampling for parameter estimation
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Calculation of WAIC using MCMC samples
library(loo)

# Visualization
library(gridExtra)
library(ggplot2)

# Progress monitoring
library(timechecker)

# Logit and inverse-logit transformations
logit <- function(x) log(x / (1 - x))
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Create a directory for storing the output files
path_out <- file.path('Result', '2-2_ParameterEstimation_Wakatsuki2024')
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)



#### 1. Set analysis parameters ####

# MCMCの設定
N_CHAIN <- 4
N_ITER <- 4500
N_WARMUP <- 500
N_THIN <- 2

# Quantiles used to summarize posterior distributions
PROBS <- c(0.025, 0.05, 0.5, 0.95, 0.975)

# Names of scalar parameters
PARS_SINGLE <- c(
  'mu_intercept', 'sigma_intercept', 
  'mu_eff_MET26', 'sigma_eff_MET26', 
  'mu_eff_SR'   , 'sigma_eff_SR'   , 'eff_SR_common', 
  'mu_eff_RH'   , 'sigma_eff_RH'   , 'eff_RH_common', 
  'sigma_eff_refsite', 'sigma',
  'mu_maturity_center', 'maturity_center_common', 'sigma_maturity_center', 
  'maturity_exp')

# Names of vector parameters
PARS_MULTIPLE <- c(
  'intercept', 'eff_MET26', 'eff_RH_group', 'eff_SR_group', 
  'maturity_center_group', 'eff_refsite')



#### 2. Import data ####

# Import and format the chalky grain proportion data
dat <- readRDS(file.path('Result', '2-1_ParameterEstimation_MainModels', 'Data_CG.rds'))
dat$HTR <- as.factor(paste0('HTR', dat$HTR))

# Import weather data
met <- readRDS(file.path('Result', '2-1_ParameterEstimation_MainModels', 'Data_Met.rds'))

# Compile the Stan model
stan_mod <- stan_model(file.path('Script', '2_MainModels.stan'))



#### 3. Prepare data for Stan ####

# 基本データを準備する
standat <- with(dat, list(
  
  # Basic dimensions
  N = nrow(dat),         # Number of observations
  D = dim(met)[1],       # Number of days in the weather data
  K = nlevels(HTR),      # Number of heat tolerance ranks
  M = nlevels(RefSite),  # Number of reference-site combinations
  
  # Response variable
  logitCG = logitCG,
  
  # Explanatory variables
  AT26 = t(pmax(met[, , 'tmp'] - 26, 0)),
  SR = t(met[, , 'srd']),
  RH = t(met[, , 'rhu']),
  
  # Grouping information associated with each observation
  group = as.integer(HTR),
  refsite = as.integer(RefSite),

  # Define the common model structure
  is_interc_ran = FALSE,
  is_MET26_ran  = FALSE,
  is_SR_group   = FALSE,
  is_RH_group   = FALSE,
  is_mat_ested  = FALSE,
  is_mat_group  = FALSE
  
))



#### 4. Fit models ####

# Function for saving warnings generated during Stan sampling
save_warning <- function(w, file) {
  
  # Create a log entry containing the time, function call, and warning message
  log_line <- sprintf(
    '%s | WARNING | call=%s\n%s\n',
    format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
    deparse(conditionCall(w)),
    conditionMessage(w))
  
  # Append the warning to the log file
  cat(log_line, '\n', file = file, append = TRUE)
  
  # Suppress subsequent display of the warning in the console
  invokeRestart('muffleWarning')
  
}

rmse <- function(obs, fit) {
  
  selector <- (!is.na(obs)) & (!is.na(fit))
  if (sum(selector) == 0) stop('No pair of data point.')
  obs <- obs[selector]
  fit <- fit[selector]
  
  n <- length(obs)
  if (n != length(fit)) stop('Length of two vectors are not equal.')
  
  sse <- sum((obs - fit) ^ 2)
  rmse <- sqrt(sse / n)
  return(rmse)
  
}

# Create a directory for the current model
filename <- paste0(path_out, .Platform$file.sep)

# Estimate model parameters using Stan
withCallingHandlers(
  {
    stan_res <- sampling(
      stan_mod, data = standat, 
      chains = N_CHAIN, warmup = N_WARMUP, iter = N_ITER, thin = N_THIN)
  },
  warning = function(w) save_warning(w, paste0(filename, '01_stan_warn.txt'))
)

# Save the fitted Stan model
saveRDS(stan_res, paste0(filename, '01_stan.rds'))

# Save posterior summaries
stan_summary <- summary(stan_res, pars = PARS_SINGLE)
write.csv(stan_summary$summary, paste0(filename, '01_stan_summary.csv'))
for (j in seq_len(dim(stan_summary$c_summary)[3])) {
  filename_j <- sprintf('%s01_stan_summary_chain%02i.csv', filename, j)
  write.csv(stan_summary$c_summary[, , j], filename_j)
}

# Examine MCMC convergence using trace plots
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

# Examine posterior correlations using a pairs plot
pars_pairs <- PARS_SINGLE
pars_pairs <- pars_pairs[!grepl('_eff_SR$', pars_pairs)]
pars_pairs <- pars_pairs[!grepl('_eff_RH$', pars_pairs)]
pars_pairs <- pars_pairs[!grepl('maturity', pars_pairs)]
png(paste0(filename, '03_Fig_pairs.png'), 
    width = 10, height = 10, unit = 'in', res = 300)
pairs(stan_res, pars = pars_pairs, cex.labels = 0.7)
dev.off()

# Extract MCMC samples
stan_res_list <- rstan::extract(stan_res)

# Save summaries of scalar parameters
summary_single <- sapply(stan_res_list[PARS_SINGLE], function(vec) {
  quantile(vec, probs = PROBS)
})
write.csv(summary_single, paste0(filename, '03_Res_Summary_Common.csv'))
saveRDS(summary_single, paste0(filename, '03_Res_Summary_Common.rds'))

# Save summaries of vector parameters
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

# Calculate RMSE
est_logitcg <- colMeans(stan_res_list$mu)
rmse_logitcg <- rmse(dat$logitCG, est_logitcg)
rmse_cg <- rmse(inv_logit(dat$logitCG), inv_logit(est_logitcg)) * 100

# Calculate WAIC and ELPD
log_lik <- extract_log_lik(stan_res, merge_chains = FALSE)
withCallingHandlers(
  {waic <- waic(log_lik)},
  warning = function(w) save_warning(w, paste0(filename, '05_WAIC_warn.txt')))
withCallingHandlers(
  {elpd <- loo(log_lik)},
  warning = function(w) save_warning(w, paste0(filename, '05_ELPD_warn.txt')))
stan_eval <- list(WAIC = waic, ELPD = elpd)
saveRDS(stan_eval, paste0(filename, '05_Res_WAIC.rds'))

# Save WAIC, ELPD, and RMSE as a CSV file
stan_eval_csv <- rbind(
  stan_eval$WAIC$estimates, stan_eval$ELPD$estimates, 
  cbind(c(RMSE_logitCG = rmse_logitcg, RMSE_CG = rmse_cg), NA))
write.csv(stan_eval_csv, paste0(filename, '05_Res_WAIC.csv'))
