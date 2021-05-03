## #######################################################################################
##
## Run space-time modeling estimating proportion of women who have never attended school
##
## Author: Nat Henry
## Created: 22 April 2021
## Purpose: Estimate proportion of women who have never attended school by 5x5km gridded
##  pixel across Uganda.
##
## Executed in R, called on the command line. Example model run script call:
## > Rscript --vanilla run_uga_model.R --config-fp /path/to/config.yaml --model-version v1
##     --use-covs intercept access lights --spatial
##
## #######################################################################################


## 0) Setup ----------------------------------------------------------------------------->

  library(argparse)
  library(yaml)

  ## Get arguments from command line
  ap <- argparse::ArgumentParser(description='Uganda education spatial model')
  ap$add_argument('--config-fp', type='character', help='Filepath for config YAML file')
  ap$add_argument('--model_version', type='character', help='Version ID for model settings')
  ap$add_argument('--holdout', type='integer', default=0, help='Holdout number')
  ap$add_argument(
    '--use-covs', type='character', nargs='+',
    help='Covariates to use, including intercept', default=c('intercept','access','lights')
  )
  ap$add_argument('--spatial', help='Use spatial random effect?', action='store_true')
  args <- ap$parse_args(commandArgs(TRUE))
  # Example command line arguments
  # args <- list(
  #   config_fp = '~/temp_data/hss_coding_test/config.yaml', model_version = 'full_model',
  #   holdout = 0, use_covs = c('intercept','access','lights'), spatial=TRUE
  # )

  ## Load config, all default libraries, and R helper functions
  config <- yaml::read_yaml(args$config_fp)
  invisible(lapply(config$settings$load_libraries, require, character.only=TRUE))
  source(file.path(config$folders$code, config$scripts$r_functions))

  # Set random seed for script execution
  set.seed(config$random_seed)


## 1) Load and format input data -------------------------------------------------------->

  # Load tabular data
  raw_data <- data.table::fread(file.path(config$folders$in_data, config$in_files$data))
  # Randomly assign holdouts to each observation for out-of-sample estimation
  raw_data$idx_holdout <- assign_holdouts(
    dataset = raw_data, num_holdouts = config$settings$num_holdouts
  )

  # Load covariate raster brick and use it to create a "index raster" to index estimates
  # Rescale covariates to N(0, 1) to improve model fit
  cov_template_list <- load_rasters(
    covar_file = file.path(config$folders$in_data, config$in_files$covariates),
    raster_covar_names = c('access', 'lights'),
    rescale_covars = TRUE
  )
  covar_rasters <- cov_template_list$covars
  index_raster <- cov_template_list$index

  # Assign covariate values to each data observation based on location
  data_with_covars <- merge_covariates_with_data(
    in_data = raw_data, covar_rasters = covar_rasters, index_raster = index_raster
  )

  # Subset to only data with non-missing numerators, denominators, and covariates
  full_data <- na.omit(data_with_covars)
  dropped_rows <- nrow(data_with_covars) - nrow(full_data)
  if(dropped_rows > 0) message('Dropped ', dropped_rows, ' rows due to missing pixel ID.')


## 2) Prepare data for model fit, estimated in Template Model Builder (TMB) ------------->

  # Create spatial mesh that will be used to project the spatial random effect, along with
  #  corresponding precision and projection matrices
  spatial_mats <- create_spatial_matrices(
    in_data = full_data, index_raster = index_raster
  )
  prior_mats <- spatial_mats$spde$param.inla

  # TMB: input data list
  tmb_data_stack <- list(
    holdout = args$holdout, no_edu_i = full_data$edu_0, n_i = full_data$N,
    X = as.matrix(full_data[, c(args$use_covs), with=F]),
    idx_holdout = full_data$idx_holdout,
    A_projmat = spatial_mats$A_proj_data, M0 = prior_mats$M0, M1 = prior_mats$M1,
    M2 = prior_mats$M2
  )
  # TMB: input parameters list
  params_list <- list(
    beta_covs = rep(0, length(args$use_covs)), log_kappa = 0, log_tau = 0,
    z_space = rep(0, spatial_mats$mesh_s$n)
  )
  # TMB 'map' to fix select parameters during model fit, as needed
  tmb_map = list()
  if(!args$spatial){
    tmb_map$log_kappa <- as.factor(NA)
    tmb_map$log_tau <- as.factor(NA)
    tmb_map$z_space <- rep(as.factor(NA), length(params_list$z_space))
  }
  # Vector of TMB random effects
  tmb_random <- character(0)
  if(args$spatial) tmb_random <- c(tmb_random, 'z_space')


## 3) Fit model in TMB ------------------------------------------------------------------>

  # Run model fitting
  model_fit <- setup_run_tmb(
    tmb_data_stack = tmb_data_stack, params_list = params_list, tmb_random = tmb_random,
    tmb_map = tmb_map, optimization_method = 'nlminb',
    template_file = file.path(config$folders$code, config$scripts$tmb_template)
  )
  # Get joint precision matrix of all parameters
  sdrep <- TMB::sdreport(model_fit$obj, bias.correct = args$spatial, getJointPrecision = TRUE)


## 4) Predict no-education ratios by pixel ---------------------------------------------->

  # Calculate parameter and prediction draws -> summarize mean & UI bounds
  draws_list <- generate_draws(
    tmb_sdreport = sdrep, num_draws = 250, covar_names = args$use_covs,
    covar_rasters = covar_rasters, index_raster = index_raster,
    A_proj_raster = spatial_mats$A_proj_raster, spatial_re = args$spatial
  )


## 5) Summarize and save outputs -------------------------------------------------------->

  out_dir <- file.path(config$folders$outputs, args$model_version)
  dir.create(out_dir)
  setwd(out_dir)

  if(args$holdout == 0){
    # Generate raster summaries of pixel draws
    pixel_summaries <- summarize_pixel_draws(draws_list$pixel_draws, index_raster)

    # Save model objects, draws, and summaries
    data.table::fwrite(draws_list$param_draws, row.names=TRUE, file='param_draws.csv')
    for(item_name in c('pixel_draws','fe_draws','re_draws')){
      data.table::fwrite(draws_list[[item_name]], file=paste0(item_name,'.csv'))
    }
    for(item_name in names(pixel_summaries)){
      writeRaster(pixel_summaries[[item_name]], file=paste0(item_name,'.tif'))
    }
    # Save some input data objects
    writeRaster(index_raster, file='index_raster.tif')
    saveRDS(spatial_mats$mesh_s, file='mesh_s.RDS')
    fwrite(full_data, file='full_data.csv')

  } else {
    # Save pixel draws only for out-of-sample draws
    data.table::fwrite(
      draws_list$pixel_draws, row.names=TRUE,
      file=paste0('pixel_draws_holdout',holdout,'.csv')
    )
  }
