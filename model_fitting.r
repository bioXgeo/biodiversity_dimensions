
# Description -------------------------------------------------------------

# Code to reproduce statistical analysis in Read et al. manuscript
# This code reproduces the analysis described under the subheading "Model Fitting" in the Methods section of the manuscript.
# PLEASE NOTE: This code was originally run in parallel but is not written to run in parallel here.

# Script created by QDR, 4 October 2018
# Script last modified by QDR, 20 May 2019 (analysis revised based on initial review from Global Ecology and Biogeography)
# Code last tested (with reduced number of iterations) by QDR, 20 May 2019
# Tested under R version 3.5.1, brms version 2.5.0
# Contact: qread@sesync.org

# Define model fitting function -------------------------------------------

fit_mv_mm <- function(pred_df, resp_df, pred_vars, resp_vars, id_var, region_var, adj_matrix, distribution = 'gaussian', priors = NULL, n_chains = 2, n_iter = 2000, n_warmup = 1000, delta = 0.9, random_effect_type = 'spatial', force_zero_intercept = FALSE, missing_data = FALSE, exclude_locations = character(0)) {
  require(dplyr)
  require(brms)
  require(reshape2)
  # Build formula and data
  if (missing_data) {
	missing_df <- resp_df[,c(id_var, 'missing')]
  }
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_vars)]
  resp_df[,-1] <- scale(resp_df[,-1, drop = FALSE])
  resp_df[,-1] <- lapply(resp_df[,-1, drop = FALSE], as.numeric)
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  
  # Specify missing-data model if there are missing response variables (used for cross-validation step).
  if (missing_data) {
	resp_var_names <- paste(resp_vars, '| mi()')
  } else {
	resp_var_names <- resp_vars
  }
  # below, create full formula string for the model with predictors
  # create much simpler one if there aren't predictors (edited 13 June)
  if (length(pred_vars) > 0) {
	pred_df <- pred_df[, c(id_var, pred_vars)]
	pred_df[,-1] <- scale(pred_df[,-1, drop = FALSE])
	pred_df[,-1] <- lapply(pred_df[,-1, drop = FALSE], as.numeric)
	pred_var_names <- names(pred_df)[-1]
	fixed_effects <- paste(pred_var_names, collapse = '+')
	intercepts <- if (force_zero_intercept) '0' else paste('(1|', region_name, ')', sep = '')
	random_effects <- paste(c(intercepts, paste('(', pred_var_names, ' - 1|', region_name, ')', sep = '')), collapse = '+')
	formula_string <- mvbf(flist = paste(resp_var_names, '~', fixed_effects, '+', random_effects), rescor = FALSE)
	dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  } else {
	intercepts <- if (force_zero_intercept) '0 +' else ''
	formula_string <- mvbf(flist = paste(resp_var_names, '~', intercepts, paste('(1|', region_name, ')', sep = '')), rescor = FALSE)
	dat <- Reduce(left_join, list(id_df, resp_df)) %>% filter(complete.cases(.))
  }
  
  # Added 1 May 2019: Change any missing data rows to NA
  if (missing_data) {
	dat <- dat %>% left_join(missing_df)
	dat[dat$missing, resp_vars] <- NA
  }
  
  # Added 2 May 2018: get rid of any region that has less than 5 sites.
  dat <- dat %>% group_by(region) %>% filter(n() >= 5)
  
  # Added 2 May 2019: exclude entire regions explicitly if they are listed in the exclude_locations argument
  if (length(exclude_locations) > 0) {
	dat <- dat %>% filter(!grepl(paste(exclude_locations, collapse = '|'), region))
  }
  
  # Added 4 May 2018: if any region no longer has a neighbor at this point, get rid of it too.
  reduced_adj_matrix <- adj_matrix[rownames(adj_matrix) %in% dat$region, rownames(adj_matrix) %in% dat$region]
  nneighb <- rowSums(reduced_adj_matrix)
  keep_regions <- names(nneighb)[nneighb > 0]
  dat <- filter(dat, region %in% keep_regions)
  
  # Fit model, extract coefficients, and format them
  if (random_effect_type == 'spatial') {
	  mm <- brm(formula = formula_string, data = dat, family = distribution, autocor = cor_car(adj_matrix, formula = ~ 1|region),
				chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  } else {
	  mm <- brm(formula = formula_string, data = dat, family = distribution,
				chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  }
  # Edit 16 Aug: do not extract fixed effects (and combined fixed+random effects) if it is a null model without fixed effects.
  random_effects <- ranef(mm)
  random_effects <- cbind(effect = 'random', melt(random_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
  if (!force_zero_intercept | length(pred_vars) > 0) {
	fixed_effects <- fixef(mm)
    region_effects <- coef(mm)
	fixed_effects <- cbind(effect = 'fixed', region = as.character(NA), melt(fixed_effects, varnames = c('parameter', 'stat')))
	region_effects <- cbind(effect = 'coefficient', melt(region_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
    mm_coef <- fixed_effects %>% full_join(random_effects) %>% full_join(region_effects)
  } else {
	mm_coef <- random_effects
  }
  return(list(model = mm, coef = mm_coef))
}
##### END definition of fit_mv_mm()

# Define functions to extract output and summary statistics from models -------------------------------------

# get_correct_variable_names(): Function to get correct response variable names and apply where needed.
get_correct_variable_names <- function(fit) {
	raw_resp_names <- fit$model$formula$responses
	resp_idx <- match(raw_resp_names, gsub('_', '', names(fit$model$data)))
	names(fit$model$data)[resp_idx]
}
##### END definition of get_correct_variable_names()

# get_relative_rmse(): Function to get relative RMSE for the model, including its MCMC quantiles
get_relative_rmse <- function(fit, resp_var_names, predicted_values) {
	# Prediction raw values. 
	pred_raw <- predict(fit$model, summary = FALSE)
	dimnames(pred_raw)[[3]] <- resp_var_names

	# Observed raw values
	obs_raw <- fit$model$data[, resp_var_names]

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
	predicted_values %>%
		group_by(response) %>%
		summarize(range_obs = diff(range(observed))) %>%
		left_join(rmse_quantiles) %>%
		mutate_at(vars(starts_with('RMSE')), funs(relative = ./range_obs))
}
##### END definition of get_relative_rmse()

# get_kfold_rmse(): function to get k-fold relative root mean squared error for each model
get_kfold_rmse <- function(fit_ids, K) {
  fit_folds <- model_fits[fit_ids[-1]]

	resp_names <- get_correct_variable_names(fit_folds[[1]])
	
	# Load full fit so we can get the data out
	fit <- model_fits[[fit_ids[1]]]

	pred_folds_raw <- map(fit_folds, function(x) {
		preds <- predict(x$model, summary = FALSE)
		dimnames(preds)[[3]] <- resp_names
		preds
	})
	
	# Extract slices of predicted and observed that correspond to the holdout data points for each fold.
	pred_folds_holdout <- map(1:K, function(i) {
		holdout_idx <- fit_folds[[i]]$model$data$region %in% region_folds[i]
		pred_folds_raw[[i]][, holdout_idx, ]
	})
	
	# also change this so it just reorders the main data df to the same order as the holdout
	obs_folds_holdout <- map(1:K, function(i) {
		holdout_idx <- fit$model$data$region %in% region_folds[i]
		fit$model$data[holdout_idx, c(resp_names)]
	})
	
	# Bind the slices into the proper dimensions 
	pred_all_holdout <- abind(pred_folds_holdout, along = 2)
	obs_all_holdout <- do.call('rbind', obs_folds_holdout)
	
	# sweep out observed from fitted and calculate RMSE
	rmse_quantiles <- sweep(pred_all_holdout, 2:3, as.matrix(obs_all_holdout), FUN = '-') %>% # Subtract predicted - observed
		melt(varnames = c('iter', 'idx', 'response')) %>%
		group_by(response, iter) %>%
		summarize(RMSE = sqrt(mean(value^2))) %>%
		ungroup %>% group_by(response) %>%
		summarize(RMSE_mean = sqrt(mean(RMSE^2)), 
				  RMSE_q025 = quantile(RMSE, probs = 0.025), 
				  RMSE_q975 = quantile(RMSE, probs = 0.975))

	# Generate ranges of observed data and divide this by the RMSE values to get the relative RMSE values
	obs_folds_holdout %>%
		melt(variable.name = 'response') %>%
		group_by(response) %>%
		summarize(range_obs = diff(range(value))) %>%
		left_join(rmse_quantiles) %>%
		mutate_at(vars(starts_with('RMSE')), funs(relative = ./range_obs))
	
}	
##### END definition of get_kfold_rmse()

# Load data ---------------------------------------------------------------

library(brms) # For instructions on installation please see https://github.com/paul-buerkner/brms
library(purrr)
library(reshape2)
library(dplyr)
library(tidyr)
library(abind)

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

# Logit transformation of beta taxonomic diversity.
bbsbio$beta_td_sorensen_pa <- qlogis(bbsbio$beta_td_sorensen_pa)
fiabio$beta_td_sorensen_pa <- qlogis(fiabio$beta_td_sorensen_pa)

# The following six ecoregions should not be used in any model fitting because they have too few data points. 
# They are primarily in Canada or Mexico with only a small portion of area in the USA, once buffer is deducted
exclude_regions <- c('NA0801', 'NA0808', 'NA0417', 'NA0514', 'NA1202', 'NA1301')

# Include the ecoregion folds, less the excluded ones
fold_df <- read.csv('ecoregion_folds.csv', stringsAsFactors = FALSE)
region_folds <- fold_df$TNC
region_folds <- region_folds[!grepl(paste(exclude_regions, collapse = '|'), region_folds)]

# Specify model fitting options for brms ----------------------------------

NC <- 3       # 3 chains
NI <- 5000    # 5000 total steps
NW <- 3000    # 3000 warmup steps (out of 5000)
delta <- 0.9  # target acceptance rate for Hamiltonian Monte Carlo (brms default is 0.8)

K <- 63		  # Number of cross-validation folds

# Create table of candidate models to fit ---------------------------------

# There are 12 possible combinations of 2 taxa * 3 diversity levels * 4 candidate models.
model_table <- expand.grid(taxon = c('fia','bbs'),
                           rv = c('alpha', 'beta', 'gamma'),
                           ecoregion = 'TNC',
                           model = c('full','climate','space', 'geo'),
                           fold = 0:K,
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
# The models were originally fit in parallel: 3 parallel chains and 24 models, each of which is fit to the full dataset and to the 63 cross-validation folds.
# 3 chains * 24 models * (63+1) folds = 4608 tasks.

# Here, the model fitting is done by looping through each combination of taxon, response variable, predictor set, and cross-validation fold sequentially.

n_fits <- nrow(model_table)
model_fits <- list()

# Construct BRMS prior objects to be used for the appropriate models.
model_priors <- prior_table %>%
  group_by(rv, taxon) %>%
  do(prior = pmap(., function(prior, class, resp, ...) set_prior(prior = prior, class = class, resp = resp)))

model_table <- model_table %>% as_tibble %>% left_join(model_priors)

for (i in 1:n_fits) {
  
  biodat <- switch(model_table$taxon[i], fia = fiabio, bbs = bbsbio)
  geodat <- switch(model_table$taxon[i], fia = fiageo, bbs = bbsgeo)
  siteid <- switch(model_table$taxon[i], fia = 'PLT_CN', bbs = 'rteNo')
  
  if (model_table$fold[i] != 0) {
    # Join response variable data with the region ID, then set the appropriate values to NA
    biodat <- biodat %>% left_join(geodat[, c(siteid, 'TNC')])
    biodat$missing <- biodat$TNC == region_folds[model_table$fold[i]]
  }
  
  model_fits[[i]] <- fit_mv_mm(pred_df = geodat, 
                               resp_df = biodat, 
                               pred_vars = switch(model_table$model[i], full = prednames, climate = climate_prednames, geo = geo_prednames, space = character(0)), 
                               resp_vars = switch(model_table$rv[i], alpha = alpha_resp, beta = beta_resp, gamma = gamma_resp), 
                               id_var = siteid, 
                               region_var = 'TNC', 
                               distribution = 'gaussian', 
                               adj_matrix = tnc_bin,
                               priors = do.call(c, model_table$prior[[i]]),
                               n_chains = NC,
                               n_iter = NI,
                               n_warmup = NW,
                               delta = delta,
                               missing_data = model_table$fold[i] > 0,
                               exclude_locations = exclude_regions
  )
}

# Extract summary information from model fits -----------------------------

# Note: as above, this was originally done in parallel. (24 cores, 1 for each model fit)

model_stats <- list()
n_full_fits <- sum(model_table$fold == 0)

for (i in 1:n_full_fits) {
  fit <- model_fits[[i]]
  
  # Get the correct variable names and apply them where needed.
  raw_resp_names <- fit$model$formula$responses
  resp_names <- get_correct_variable_names(fit)
  
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
  model_rmse <- get_relative_rmse(fit, resp_names, model_pred)
  
  # Bayesian R-squared
  model_r2 <- cbind(model_table[i, c('taxon', 'rv', 'ecoregion', 'model')], response = resp_names, bayes_R2(fit$model))
  
  # Information criterion (WAIC)
  model_waic <- waic(fit$model)
  model_waic <- cbind(model_table[i, c('taxon', 'rv', 'ecoregion', 'model')], WAIC = model_waic$estimates['waic','Estimate'], WAIC_SE = model_waic$estimates['waic','SE'])
  
  # RMSE from K-fold cross validation
  fit_ids <- with(model_table, which(taxon == taxon[i] & rv == rv[i] & model == model[i]))
  kfold_rmse <- get_kfold_rmse(fit_ids, K = K)
	
  model_stats[[i]] <- list(coef = model_coef, pred = model_pred, rmse = model_rmse, r2 = model_r2, waic = model_waic, kfold_rmse = kfold_rmse)
  
}

model_coef <- map2_dfr(model_stats, 1:n_full_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = 'TNC', model = model_table$model[y], as.data.frame(x$coef)))
model_pred <- map2_dfr(model_stats, 1:n_full_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = 'TNC', model = model_table$model[y], as.data.frame(x$pred)))
model_rmse <- map2_dfr(model_stats, 1:n_full_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = 'TNC', model = model_table$model[y], as.data.frame(x$rmse)))
model_r2 <- map_dfr(model_stats, 'r2')
model_waic <- map_dfr(model_stats, 'waic')
model_kfold_rmse <- map2_dfr(model_stats, 1:n_full_fits, function(x, y) cbind(taxon = model_table$taxon[y], rv = model_table$rv[y], ecoregion = model_table$ecoregion[y], model = model_table$model[y], as.data.frame(x$kfold_rmse)))

# Extract spatial variability of coefficients in full models only ---------

model_sds <- list()

for (i in 1:6) {
  fit <- model_fits[[i]]
  sds <- summary(fit$model)$random$region[,c(1,3,4)]
  dimnames(sds)[[2]] <-  c('Estimate', 'q025', 'q975')
  model_sds[[i]] <- sds
}

# Parse the response and parameter names out.
parse_names <- function(ns) {
  ns <- map_chr(strsplit(ns, '\\(|\\)'), 2) # get rid of parentheses around name
  response <- unlist(regmatches(ns, gregexpr('^[^_]*', ns)))
  parameter <- unlist(regmatches(ns, gregexpr('_.*$', ns)))
  parameter <- substr(parameter, 2, nchar(parameter))
  return(data.frame(response = response, parameter = parameter))
}

model_sds <- map(model_sds, function(x) cbind(parse_names(dimnames(x)[[1]]), x))

model_coef_var <- cbind(taxon = rep(model_table$taxon[1:6], map_int(model_sds, nrow)),
                        do.call(rbind, model_sds))

# Write model fit summaries to CSVs ---------------------------------------

write.csv(model_coef, 'model_coef.csv', row.names = FALSE)                 # Coefficient estimates and credible intervals
write.csv(model_pred, 'model_pred.csv', row.names = FALSE)                 # Fitted values and credible intervals    
write.csv(model_rmse, 'model_rmse.csv', row.names = FALSE)                 # Root mean squared errors for models fit to full dataset
write.csv(model_r2, 'model_r2.csv', row.names = FALSE)                     # Bayesian R-squared values for models
write.csv(model_kfold_rmse, 'model_kfold_rmse.csv', row.names = FALSE)     # Root mean squared errors for models fit to cross-validation datasets
write.csv(model_waic, 'model_waic.csv', row.names = FALSE)                 # Widely Available Information Criterion values for models 
write.csv(model_coef_var, 'model_coef_var.csv', row.names = FALSE)         # Spatial variation in model coefficients by ecoregion

