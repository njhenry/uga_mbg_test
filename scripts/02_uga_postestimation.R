## #######################################################################################
##
## Visualize model results for Ugandan education
##
## Author: Nat Henry
## Created: 22 April 2021
## Purpose: Estimate proportion of women who have never attended school by 5x5km gridded
##  pixel across Uganda.
##
## #######################################################################################


## 0) Setup ----------------------------------------------------------------------------->

  model_version <- 'full_model'
  config_fp <- '~/temp_data/hss_coding_test/config.yaml'

  config <- yaml::read_yaml(config_fp)
  invisible(lapply(config$settings$load_libraries, require, character.only=TRUE))
  source(file.path(config$folders$code, config$scripts$r_functions))

## 1) Load data ------------------------------------------------------------------------->

  out_dir <- file.path(config$folders$outputs, model_version)
  if(!dir.exists(out_dir)) stop("Missing outputs folder")
  setwd(out_dir)

  full_data <- data.table::fread('full_data.csv')
  summ_stats <- c('mean', 'median', 'lower', 'upper')
  rasters <- lapply(summ_stats, function(x) raster::raster(paste0(x,'.tif')))
  rnames <- paste0('r',summ_stats)
  names(rasters) <- rnames
  pixel_draws <- data.table::fread('pixel_draws.csv')
  param_draws <- data.table::fread('param_draws.csv')
  mesh_s <- readRDS("mesh_s.RDS")
  covariates <- raster::brick('../../in_data/covariates_brick.tif')
  access <- covariates[[1]]
  lights <- covariates[[2]]


## 2) Map mean, lower, and upper rasters using the same color scheme -------------------->

  plot_rasters <- brick(rasters[c('rmean','rlower','rupper')])
  names(plot_rasters) <- c('Mean','Lower','Upper')
  png('raster_summaries.png', height=8, width=8, units='in', res=300)
  plot(plot_rasters, col = viridis::magma(20))
  dev.off()

  # Also plot the spatial mesh
  plot_png <- function(fn, fig){
    png(fn, height=6, width=6, units='in', res=300); plot(fig); dev.off()
    invisible()
  }
  plot_png('spatial_mesh.png', mesh_s)

  # Also plot spatial covariates
  plot_png('access.png', covariates[[1]])
  plot_png('lights.png', covariates[[2]])
  png('access_hist.png', height=6, width=6, units='in', res=300)
  hist(values(covariates[[1]]))
  dev.off()
  png('lights_hist.png', height=6, width=6, units='in', res=300)
  hist(values(covariates[[2]]))
  dev.off()

## 3) Get in-sample model fit metrics: RMSE and coverage

  # Get summary statistics (mean, lower, upper) for each pixel associated with a data point
  full_data$est_mean <- na.omit(values(rasters$rmean))[full_data$pixel_idx]
  full_data$est_lower <- na.omit(values(rasters$rlower))[full_data$pixel_idx]
  full_data$est_upper <- na.omit(values(rasters$rupper))[full_data$pixel_idx]

  full_data[, raw_est := edu_0 / N
    ][, rmse_contrib := sqrt((raw_est - est_mean)^2)
    ][, coverage_contrib := as.integer((raw_est > est_lower) & (raw_est < est_upper))]

  error_data <- full_data[, .(rmse=sum(rmse_contrib)/.N, coverage=sum(coverage_contrib)/.N)]
  print(error_data)
  data.table::fwrite(error_data, file='is_metrics.csv')

