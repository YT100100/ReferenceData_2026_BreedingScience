# 
# Fitting main models for chalky grain proportions of rice
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
path_out <- file.path('Result', '2-1_ParameterEstimation_MainModels')
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)



#### 1. Set analysis parameters ####

# MCMC settings
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
dat_raw <- readRDS(file.path('Result', '1_DataPreparation', 'Data.rds'))
dat <- dat_raw$df[, c('ID', 'HTR', 'Cultivar', 'Pref', 'RefNo', 'CG')]
dat <- within(dat, {
  
  # Create reference-site identifiers
  RefSite  <- sprintf('%02i_%s', RefNo, Pref)
  
  # Apply the logit transformation to chalky grain percentage
  logitCG <- logit(CG / 100)
  
  # Convert categorical variables to factors
  HTR      <- as.factor(HTR)
  Cultivar <- as.factor(Cultivar)
  RefSite  <- as.factor(RefSite)
  
})

# Import weather data
met <- readRDS(file.path('Data', 'meshdat.rds'))

# Print the numbers of observations, cultivars, and reference-site combinations
cat(sprintf('nrow: %i, n. cultivar: %i, n. refsite: %i\n', 
            nrow(dat), nlevels(dat$Cultivar), nlevels(dat$RefSite)))

# Remove observations containing missing values
sel <- !is.na(dat$CG)
sel <- sel & apply(met, 2, function(x) all(!is.na(x)))
dat <- dat[sel, ]
met <- met[, sel, ]
dat$Cultivar <- droplevels(dat$Cultivar)
dat$RefSite  <- droplevels(dat$RefSite)

# Print the numbers of observations, cultivars, and reference-site combinations
cat(sprintf('nrow: %i, n. cultivar: %i, n. refsite: %i\n', 
            nrow(dat), nlevels(dat$Cultivar), nlevels(dat$RefSite)))

# Save the analysis datasets
write.csv(dat, file.path(path_out, 'Data_CG.csv'))
saveRDS(dat, file.path(path_out, 'Data_CG.rds'))
saveRDS(met, file.path(path_out, 'Data_Met.rds'))

# Compile the Stan model
stan_mod <- stan_model(file.path('Script', '2_MainModels.stan'))



#### 3. Prepare data for Stan ####

# Prepare data shared among all model specifications
standat <- with(dat, list(
  
  # Basic dimensions
  N = nrow(dat),         # Number of observations
  D = dim(met)[1],       # Number of days in the weather data
  K = nlevels(Cultivar), # Number of cultivars
  M = nlevels(RefSite),  # Number of reference-site combinations
  
  # Response variable
  logitCG = logitCG,
  
  # Explanatory variables
  AT26 = t(pmax(met[, , 'tmp'] - 26, 0)),
  SR = t(met[, , 'srd']),
  RH = t(met[, , 'rhu']),
  
  # Grouping information associated with each observation
  group   = as.integer(Cultivar),
  refsite = as.integer(RefSite),
  
  # Define the common model structure
  is_interc_ran = TRUE,
  is_MET26_ran  = TRUE

))

# Generate combinations of model specifications
pat <- expand.grid(
  is_SR_group = c(TRUE, FALSE),
  is_RH_group = c(TRUE, FALSE),
  maturity_mode  = c('20days', 'fixed', 'cultivar'),
  stringsAsFactors = FALSE)
pat <- within(pat, {
  is_mat_ested <- maturity_mode != '20days'
  is_mat_group <- maturity_mode == 'cultivar'
})

# Create model names
pat <- within(pat, {
  modname <- sprintf(
    'SRcult=%s_RHcult=%s_mat=%s',
    is_SR_group, is_RH_group, maturity_mode)
})



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

tc <- set_loop_timechecker(nrow(pat))
invisible(lapply(seq_len(nrow(pat)), function(i) {
  
  # Create a directory for the current model
  filename_i <- paste0(
    path_out, .Platform$file.sep, pat[i, 'modname'], .Platform$file.sep)
  dir.create(filename_i, showWarnings = FALSE)
  
  # Add model-specific settings to the Stan data
  standat_i <- c(standat, list(
    is_SR_group  = pat[i, 'is_SR_group' ],
    is_RH_group  = pat[i, 'is_RH_group' ],
    is_mat_ested = pat[i, 'is_mat_ested'],
    is_mat_group = pat[i, 'is_mat_group']))
  
  # Estimate model parameters using Stan
  withCallingHandlers(
    {
      stan_res_i <- sampling(
        stan_mod, data = standat_i, 
        chains = N_CHAIN, warmup = N_WARMUP, iter = N_ITER, thin = N_THIN)
    },
    warning = function(w) save_warning(w, paste0(filename_i, '01_stan_warn.txt'))
  )
  
  # Save the fitted Stan model
  saveRDS(stan_res_i, paste0(filename_i, '01_stan.rds'))
  
  # Save posterior summaries
  stan_summary_i <- summary(stan_res_i, pars = PARS_SINGLE)
  write.csv(stan_summary_i$summary, paste0(filename_i, '01_stan_summary.csv'))
  for (j in seq_len(dim(stan_summary_i$c_summary)[3])) {
    filename_j <- sprintf('%s01_stan_summary_chain%02i.csv', filename_i, j)
    write.csv(stan_summary_i$c_summary[, , j], filename_j)
  }
  
  # Examine MCMC convergence using trace plots
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
  
  # Examine posterior correlations using a pairs plot
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
  
  # Extract MCMC samples
  stan_res_list_i <- rstan::extract(stan_res_i)
  
  # Save summaries of scalar parameters
  summary_single_i <- sapply(stan_res_list_i[PARS_SINGLE], function(vec) {
    quantile(vec, probs = PROBS)
  })
  write.csv(summary_single_i, paste0(filename_i, '03_Res_Summary_Common.csv'))
  saveRDS(summary_single_i, paste0(filename_i, '03_Res_Summary_Common.rds'))
  
  # Save summaries of vector parameters
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
  
  # Calculate RMSE
  est_logitcg_i <- colMeans(stan_res_list_i$mu)
  rmse_logitcg_i <- rmse(dat$logitCG, est_logitcg_i)
  rmse_cg_i <- rmse(inv_logit(dat$logitCG), inv_logit(est_logitcg_i)) * 100
  
  # Calculate WAIC and ELPD
  log_lik_i <- extract_log_lik(stan_res_i, merge_chains = FALSE)
  withCallingHandlers(
    {waic_i <- waic(log_lik_i)},
    warning = function(w) save_warning(w, paste0(filename_i, '05_WAIC_warn.txt')))
  withCallingHandlers(
    {elpd_i <- loo(log_lik_i)},
    warning = function(w) save_warning(w, paste0(filename_i, '05_ELPD_warn.txt')))
  
  stan_eval_i <- list(WAIC = waic_i, ELPD = elpd_i)
  saveRDS(stan_eval_i, paste0(filename_i, '05_Res_WAIC.rds'))
  
  # Save WAIC, ELPD, and RMSE as a CSV file
  stan_eval_csv_i <- rbind(
    stan_eval_i$WAIC$estimates, stan_eval_i$ELPD$estimates, 
    cbind(c(RMSE_logitCG = rmse_logitcg_i, RMSE_CG = rmse_cg_i), NA))
  write.csv(stan_eval_csv_i, paste0(filename_i, '05_Res_WAIC.csv'))
  
  tc()
  
}))
