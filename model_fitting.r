
# Description -------------------------------------------------------------

# Code to reproduce statistical analysis in Read et al. manuscript
# This code reproduces the analysis described under the subheading "Model Fitting" in the Methods section of the manuscript.

# Script created by QDR, 4 October 2018
# Script last modified by QDR, 9 October 2018
# Code last tested (with reduced number of iterations) by QDR, 9 October 2018
# Tested under R version 3.5.1, brms version 2.5.0
# Contact: qread@sesync.org

# Define model fitting function -------------------------------------------

fit_mv_mm <- function(pred_df, resp_df, pred_vars, resp_vars, id_var, region_var, adj_matrix, distribution = 'gaussian', priors = NULL, n_chains = 2, n_iter = 2000, n_warmup = 1000, delta = 0.9, random_effect_type = 'spatial', force_zero_intercept = FALSE) {
  require(dplyr)
  require(brms)
  require(reshape2)
  # Create data frame with site IDs and ecoregions, then bind it to the predictor data frame.
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_vars)]
  resp_df[,-1] <- scale(resp_df[,-1, drop = FALSE])
  resp_df[,-1] <- lapply(resp_df[,-1, drop = FALSE], as.numeric)
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  resp_var_names <- paste0('cbind(', paste(resp_vars, collapse = ','), ')')

  if (length(pred_vars) > 0) {
    # Create formula string containing all fixed and random effects.
    pred_df <- pred_df[, c(id_var, pred_vars)]
    pred_df[,-1] <- scale(pred_df[,-1, drop = FALSE])
    pred_df[,-1] <- lapply(pred_df[,-1, drop = FALSE], as.numeric)
    pred_var_names <- names(pred_df)[-1]
    fixed_effects <- paste(pred_var_names, collapse = '+')
    intercepts <- if (force_zero_intercept) '0' else paste('(1|', region_name, ')', sep = '')
    random_effects <- paste(c(intercepts, paste('(', pred_var_names, ' - 1|', region_name, ')', sep = '')), collapse = '+')
    formula_string <- paste(resp_var_names, '~', fixed_effects, '+', random_effects)
    dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  } else {
    # If fitting a null model (no fixed predictors) create formula string with intercepts only.
    intercepts <- if (force_zero_intercept) '0 +' else ''
    formula_string <- paste(resp_var_names, '~', intercepts, paste('(1|', region_name, ')', sep = ''))
    dat <- Reduce(left_join, list(id_df, resp_df)) %>% filter(complete.cases(.))
  }
  
  dat <- dat %>% group_by(region) %>% filter(n() >= 5) # Get rid of any region with <5 sites.
  
  # If any region no longer has an adjacent neighbor at this point, get rid of it too.
  reduced_adj_matrix <- adj_matrix[rownames(adj_matrix) %in% dat$region, rownames(adj_matrix) %in% dat$region]
  nneighb <- rowSums(reduced_adj_matrix)
  keep_regions <- names(nneighb)[nneighb > 0]
  dat <- filter(dat, region %in% keep_regions)
  
  # Fit model
  if (random_effect_type == 'spatial') {
    mm <- brm(formula = formula_string, data = dat, family = distribution, autocor = cor_car(adj_matrix, formula = ~ 1|region),
              chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  } else {
    mm <- brm(formula = formula_string, data = dat, family = distribution,
              chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  }
  # Extract random effects from fit model object.
  random_effects <- ranef(mm)
  random_effects <- cbind(effect = 'random', melt(random_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
  if (!force_zero_intercept | length(pred_vars) > 0) {
    # Extract fixed effects and coefficients (fixed + random effects) from fit model object
    fixed_effects <- fixef(mm)
    region_effects <- coef(mm)
    fixed_effects <- cbind(effect = 'fixed', region = as.character(NA), melt(fixed_effects, varnames = c('parameter', 'stat')))
    region_effects <- cbind(effect = 'coefficient', melt(region_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
    mm_coef <- fixed_effects %>% full_join(random_effects) %>% full_join(region_effects)
  } else {
    # Skip extraction of fixed effects if fitting a null model with only random effects.
    mm_coef <- random_effects
  }
  return(list(model = mm, coef = mm_coef))
}
##### END definition of fit_mv_mm()


# Define function to do single fold of the k-fold cross-validation --------

onefold <- function(fit, k, ksub, n_chains, n_iter, n_warmup, delta = 0.8, seed = 101) {
  require(brms)
  require(purrr)
  require(dplyr)
  require(reshape2)
  
  # The group 'fold' is specified ahead of time.
  assign_fold <- function(n, k) {
    sample(rep_len(sample(1:k), n))
  }
  
  set.seed(seed)
  folds <- fit$data %>% group_by(region) %>% transmute(fold = assign_fold(n(), k))
  fit$data$fold <- factor(folds$fold)
  
  # Create MCMC seed based on time
  now <- as.numeric(Sys.time())
  (mcmc_seed <- trunc((now-floor(now))*10000))
  
  # Fit only the specified fold.
  kf <- kfold(fit, Ksub = as.array(ksub), chains = n_chains, cores = n_chains, iter = n_iter, warmup = n_warmup, control = list(adapt_delta = delta), save_fits = TRUE, group = 'fold', seed = mcmc_seed)
  
  # Get names of response variables. Note that underscore characters are removed by brms, causing issues.
  resp_idx <- match(fit$formula$responses, gsub('_', '', names(fit$data)))
  resp_names <- names(fit$data)[resp_idx]
  
  # Predicted values for the subset not in the specified fold.
  # Since this is multivariate, we need to rewrite this code to get multiple y obs and y pred columns
  oos_pred <- map2(kf$fits[,'fit'], kf$fits[,'omitted'], function(fit_fold, idx) {
    pred_raw <- predict(fit_fold, newdata = fit$data[idx,], summary = FALSE)
    dimnames(pred_raw)[[3]] <- resp_names
    obs_raw <- fit$data[idx, resp_names]
    sweep(pred_raw, 2:3, as.matrix(obs_raw), FUN = '-') %>% # Subtract predicted - observed
      melt(varnames = c('iter', 'idx', 'response'))
  })

  # Add fold ID
  oos_pred <- data.frame(fold = ksub, bind_rows(oos_pred$fit))
  
  return(list(kfold_estimates = kf$estimates, oos_pred = oos_pred))
  
}
#### END definition of onefold()


# Load data ---------------------------------------------------------------

library(brms) # For instructions on installation please see https://github.com/paul-buerkner/brms
library(purrr)
library(dplyr)
library(tidyr)

# Load predictor and response variables for each site.
bbsbio <- read.csv('bbs_biodiversity.csv', stringsAsFactors = FALSE)
bbsgeo <- read.csv('bbs_geodiversity.csv', stringsAsFactors = FALSE)
fiabio <- read.csv('fia_biodiversity.csv', stringsAsFactors = FALSE)
fiageo <- read.csv('fia_geodiversity.csv', stringsAsFactors = FALSE)

# Load adjacency matrix.
tnc_bin <- read.csv('tnc_adjacencymatrix.csv', row.names = 1)
tnc_bin <- as.matrix(tnc_bin)
dimnames(tnc_bin)[[2]] <- NULL

# Load table of non-default priors, specified to aid model convergence for some of the models.
prior_table <- read.csv('prior_table.csv', stringsAsFactors = FALSE)

# Specify model fitting options for brms ----------------------------------

NC <- 3       # 3 chains
NI <- 5000    # 5000 total steps
NW <- 3000    # 3000 warmup steps (out of 5000)
delta <- 0.9  # target acceptance rate for Hamiltonian Monte Carlo (brms default is 0.8)


# Create table of candidate models to fit ---------------------------------

# There are 12 possible combinations of 2 taxa * 3 diversity levels * 4 candidate models.
model_table <- expand.grid(taxon = c('bbs', 'fia'),
                           rv = c('alpha', 'beta', 'gamma'),
                           model = c('full','climate','space', 'geo'),
                           stringsAsFactors = FALSE)

# Specify predictor and response variable names.
prednames <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')
climate_prednames <- c('bio1_5k_50_mean', 'bio12_5k_50_mean')
geo_prednames <- c('elevation_5k_tri_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'dhi_gpp_5k_tri_50_mean')
alpha_resp <- c('alpha_richness', 'alpha_phy_pa', 'alpha_func_pa')
beta_resp <- c('beta_td_sorensen_pa', 'beta_phy_pa', 'beta_func_pa')
gamma_resp <- c('gamma_richness', 'gamma_phy_pa', 'gamma_func_pa')

# Fit models --------------------------------------------------------------

# Note: this will take a fairly long time to run sequentially.
# The models were originally fit in parallel: 3 parallel chains and 24 models running simultaneously on different cores (72 cores total)

# Here, the model fitting is done by looping through each combination of taxon, response variable, and predictor set sequentially.

n_fits <- nrow(model_table)
model_fits <- list()

# Construct BRMS prior objects to be used for the appropriate models.
model_priors <- prior_table %>%
  group_by(rv, taxon) %>%
  do(prior = pmap(., function(prior, class, resp, ...) set_prior(prior = prior, class = class, resp = resp)))

model_table <- model_table %>% as_tibble %>% left_join(model_priors)

for (i in 1:n_fits) {
  
  model_fits[[i]] <- fit_mv_mm(pred_df = switch(model_table$taxon[i], fia = fiageo, bbs = bbsgeo), 
                               resp_df = switch(model_table$taxon[i], fia = fiabio, bbs = bbsbio), 
                               pred_vars = switch(model_table$model[i], full = prednames, climate = climate_prednames, geo = geo_prednames, space = character(0)), 
                               resp_vars = switch(model_table$rv[i], alpha = alpha_resp, beta = beta_resp, gamma = gamma_resp), 
                               id_var = switch(model_table$taxon[i], fia = 'PLT_CN', bbs = 'rteNo'), 
                               region_var = 'TNC', 
                               distribution = 'gaussian', 
                               adj_matrix = tnc_bin,
                               priors = do.call(c, model_table$prior[[i]]),
                               n_chains = NC,
                               n_iter = NI,
                               n_warmup = NW,
                               delta = delta
  )
}


# Perform cross-validation ------------------------------------------------

# Note: as above, this was originally done in parallel. 24 models * 5 folds per model * 2 chains per fold = 240 cores 

model_kfolds <- replicate(n_fits, list())

for (i in 1:n_fits) {
  for (j in 1:5) {
    model_kfolds[[i]][[j]] <- onefold(model_fits[[i]]$model, k = 5, ksub = j, n_chains = 2, n_iter = NI, n_warmup = NW, delta = delta, seed = i + 303)
  }
}

# Extract summary information from model fits -----------------------------

# Note: as above, this was originally done in parallel. (24 cores, 1 for each model fit)

model_stats <- list()

for (i in 1:n_fits) {
  fit <- model_fits[[i]]
  
  # Get the correct variable names and apply them where needed.
  raw_resp_names <- fit$model$formula$responses
  resp_idx <- match(raw_resp_names, gsub('_', '', names(fit$model$data)))
  resp_names <- names(fit$model$data)[resp_idx]
  
  fit$coef <- fit$coef %>%
    separate(parameter, into = c('response', 'parameter'), extra = 'merge') %>%
    mutate(response = resp_names[match(response, raw_resp_names)])
  
  model_coef <- fit$coef # Includes fixed, random, and coefficient.
  model_pred <- predict(fit$model) # Returns raw array: n data x 4 stats x n response variables.
  # The predict() call takes a long time (~20 min or so in some cases)
  
  dimnames(model_pred)[[3]] <- resp_names
  model_pred <- melt(model_pred, varnames=c('idx','stat','response'))
  
  # Join predicted with observed values
  model_obs <- melt(cbind(idx = 1:nrow(fit$model$data), fit$model$data[, resp_names]), id.vars = 1, value.name = 'observed', variable.name = 'response')
  model_pred <- dcast(model_pred, idx + response ~ stat) %>%
    left_join(model_obs)
  
  # Here, do the RMSE for the model.
  # Prediction raw values. 
  pred_raw <- predict(fit$model, summary = FALSE)
  dimnames(pred_raw)[[3]] <- resp_names
  
  # Observed raw values
  obs_raw <- fit$model$data[, resp_names]
  
  # Get RMSE for each iteration and their quantiles
  rmse_quantiles <- sweep(pred_raw, 2:3, as.matrix(obs_raw), FUN = '-') %>% # Subtract predicted - observed
    melt(varnames = c('iter', 'idx', 'response')) %>%
    group_by(response, iter) %>%
    summarize(RMSE = sqrt(mean(value^2))) %>%
    ungroup %>% group_by(response) %>%
    summarize(RMSE_mean = sqrt(mean(RMSE^2)), 
              RMSE_q025 = quantile(RMSE, probs = 0.025), 
              RMSE_q975 = quantile(RMSE, probs = 0.975))
  
  # Generate ranges of observed data and divide this by the RMSE values to get the relative RMSE values
  model_rmse <- model_pred %>%
    group_by(response) %>%
    summarize(range_obs = diff(range(observed))) %>%
    left_join(rmse_quantiles) %>%
    mutate_at(vars(starts_with('RMSE')), funs(relative = ./range_obs))
  
  # Bayesian R-squared
  model_r2 <- cbind(model_table[i, 1:2], ecoregion = 'TNC', model = model_table[i, 3], response = resp_names, bayes_R2(fit$model))
  
  # Spatial variation in coefficients across regions
  model_coef_var <- summary(fit$model)$random$region[, c('Estimate', 'l-95% CI', 'u-95% CI')]
  dimnames(model_coef_var)[[2]] <-  c('Estimate', 'q025', 'q975')
  
  # Parse names in the spatial variation data frame to match the other names.
  parse_names <- function(ns) {
    ns <- map_chr(strsplit(ns, '\\(|\\)'), 2) # get rid of parentheses around name
    response <- unlist(regmatches(ns, gregexpr('^[^_]*', ns)))
    parameter <- unlist(regmatches(ns, gregexpr('_.*$', ns)))
    parameter <- substr(parameter, 2, nchar(parameter))
    return(data.frame(response = response, parameter = parameter))
  }
  
  model_coef_var <- cbind(parse_names(dimnames(model_coef_var)[[1]]), model_coef_var)
  
  model_stats[[i]] <- list(coef = model_coef, pred = model_pred, rmse = model_rmse, r2 = model_r2, coef_var = model_coef_var)
  
}

model_coef <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = 'TNC', model = model_table$model[y], as.data.frame(x$coef)))
model_pred <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = 'TNC', model = model_table$model[y], as.data.frame(x$pred)))
model_rmse <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = 'TNC', model = model_table$model[y], as.data.frame(x$rmse)))
model_r2 <- map_dfr(model_stats, 'r2')
model_coef_var <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = 'TNC', model = model_table$model[y], as.data.frame(x$coef_var)))

# Extract summary information from cross-validation objects ---------------

# This was not necessary to do in parallel.

model_kfold_pred <- list()
model_kfold_stats <- list()

# For each fit, load each k-fold subset (5) and combine the outputs.

for (i in 1:n_fits) {
  # Combine all the predicted values into a single data frame.
  model_kfold_pred[[i]] <- map_dfr(model_kfolds[[i]], 'oos_pred')
  
  # Calculate RMSE for each fold and each iteration
  RMSE_all <- model_kfold_pred[[i]] %>%
    group_by(fold, iter, response) %>%
    summarize(RMSE = sqrt(mean(value^2)))
  RMSE_quantiles <- RMSE_all %>%
    ungroup %>%
    group_by(response) %>%
    summarize(kfold_RMSE_mean = sqrt(mean(RMSE^2)), 
              kfold_RMSE_q025 = quantile(RMSE, probs = 0.025), 
              kfold_RMSE_q975 = quantile(RMSE, probs = 0.975))
  RMSE_quantiles_byfold <- RMSE_all %>%
    ungroup %>%
    group_by(fold, response) %>%
    summarize(kfold_RMSE_mean = sqrt(mean(RMSE^2)), 
              kfold_RMSE_q025 = quantile(RMSE, probs = 0.025), 
              kfold_RMSE_q975 = quantile(RMSE, probs = 0.975))
  
  
  # Combine all the k-fold ICs and RMSEs into a single object.
  model_kfold_stats[[i]] <- bind_rows(data.frame(fold = NA, RMSE_quantiles),
                                      data.frame(RMSE_quantiles_byfold, 
                                                 map_dfr(model_kfolds[[i]], function(x) data.frame(kfoldic = rep(x$kfold_estimates['kfoldic','Estimate'],3),
                                                                                                   kfoldic_se = rep(x$kfold_estimates['kfoldic','SE'],3)))))
  
}

model_kfold_stats <- map2_dfr(model_kfold_stats, 1:n_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = 'TNC', model = model_table$model[y], as.data.frame(x)))


# Write model fit summaries to CSVs ---------------------------------------

write.csv(model_coef, 'model_coef.csv', row.names = FALSE)                 # Coefficient estimates and credible intervals
write.csv(model_pred, 'model_pred.csv', row.names = FALSE)                 # Fitted values and credible intervals    
write.csv(model_rmse, 'model_rmse.csv', row.names = FALSE)                 # Root mean squared errors for models fit to full dataset
write.csv(model_r2, 'model_r2.csv', row.names = FALSE)                     # Bayesian R-squared values for models
write.csv(model_kfold_stats, 'model_kfold_stats.csv', row.names = FALSE)   # Root mean squared errors for models fit to cross-validation datasets
write.csv(model_coef_var, 'model_coef_var.csv', row.names = FALSE)         # Spatial variation in model coefficients by ecoregion

