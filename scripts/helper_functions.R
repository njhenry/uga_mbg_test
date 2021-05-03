## #######################################################################################
##
## Helper functions for Uganda education modeling code
##
## Author: Nat Henry
## Created: 22 April 2021
##
## #######################################################################################


## FUNCTION: Assign holdout IDs to a new field, "idx_holdout"
assign_holdouts <- function(dataset, num_holdouts){
  sample(rep_len(1:num_holdouts, length.out = nrow(dataset)), replace = FALSE)
}


## FUNCTION: insert a vector of values into a raster layer or brick
vec_to_raster <- function(new_vals, rast){
  # Ensure matrix format
  if(!'matrix' %in% class(new_vals)) new_vals <- matrix(new_vals, ncol=1)
  # Get cell index
  cell_idx <- which(!is.na(getValues(rast[[1]])))
  # check the index makes superficial sense
  stopifnot(length(cell_idx) == nrow(new_vals))
  stopifnot(max(cell_idx) <= ncell(rast))
  # create results raster
  num_cols <- ncol(new_vals)
  raster_new <- raster::brick(replicate(num_cols, rast[[1]], simplify = FALSE))
  names(raster_new) <- colnames(new_vals)
  # update the values
  for (col_idx in 1:num_cols) raster_new[[col_idx]][cell_idx] <- new_vals[, col_idx]
  return(raster_new)
}


## FUNCTION: Given a filepath to a raster .tif, load covariate rasters and use them to
##  create a pixel-based index for modeling
load_rasters <- function(
  covar_file, raster_covar_names, rescale_covars = TRUE
){
  # Load raster brick from file
  raw_brick <- raster::brick(covar_file)
  names(raw_brick) <- raster_covar_names
  # It's clear that the light intensity covariate is flipped and heavily skewed
  # Flip, rescale the values, and log so output values are in (0,1) and brighter -> higher
  if('lights' %in% raster_covar_names){
    val_range <- range(values(raw_brick$lights), na.rm=T)
    max_val <- max(val_range)
    val_diff <- diff(val_range)
    raw_brick$lights <- calc(
      raw_brick$lights,
      function(x) log( (1-exp(1)) * (x - max_val)/val_diff + 1 )
    )
  }

  # If needed, rescale to N(0, 1) for each layer
  if(rescale_covars){
    for(lyr in 1:dim(raw_brick)[3]){
      lyr_mean <- mean(values(raw_brick[[lyr]]), na.rm=T)
      lyr_sd <- sd(values(raw_brick[[lyr]]), na.rm=T)
      raw_brick[[lyr]] <- raster::calc(raw_brick[[lyr]], function(x) (x-lyr_mean)/lyr_sd)
    }
  }

  # Create an index raster, which gives the integer cell location for each non-NA cell
  non_na_rast <- raster::calc(raw_brick, function(x) 0 * sum(x) + 1)
  cropped_brick <- raster::mask(x=raw_brick, mask=non_na_rast)
  names(cropped_brick) <- raster_covar_names
  index_raster <- vec_to_raster(1:sum(getValues(non_na_rast), na.rm=T), non_na_rast)

  return(list(covars = cropped_brick, index = index_raster))
}


## FUNCTION: Merge covariates onto the tabular data
merge_covariates_with_data <- function(in_data, covar_rasters, index_raster){
  # Convert latitude and longitude fields to SpatialPoints
  data_pts <- sp::SpatialPoints(in_data[, .(longitude, latitude)])
  # Get covariate raster values at each point
  covars_dt <- cbind(
    data.table(intercept = rep(1, nrow(in_data))),
    data.table::as.data.table(raster::extract(x=covar_rasters, y=data_pts))
  )
  names(covars_dt)[2:ncol(covars_dt)] <- names(covar_rasters)
  covars_dt$pixel_idx <- raster::extract(x=index_raster, y=data_pts)
  return(cbind(in_data, covars_dt))
}


## FUNCTION: Create the spatial mesh, spatial projection matrix, and SPDE precision
##  objects
create_spatial_matrices <- function(in_data, index_raster){
  simple_raster <- raster::calc(index_raster, function(x) 0*x+1)
  simple_polygon <- raster::rasterToPolygons(simple_raster, dissolve = TRUE)
  simple_polygon <- rgeos::gSimplify(
    raster::buffer(simple_polygon, width=0.2), tol = 0.1
  )
  # Data points -> Spatial mesh
  mesh_s <- INLA::inla.mesh.2d(
    boundary = INLA::inla.sp2segment(simple_polygon),
    loc = cbind(in_data$longitude, in_data$latitude),
    max.edge = c(0.1, 0.5), offset = 0.5, cutoff = 0.1
  )
  # Project from spatial mesh points to data observations
  A_proj_data <- INLA::inla.spde.make.A(
    mesh = mesh_s,
    loc = as.matrix(in_data[, .(longitude, latitude)])
  )
  # Project from spatial mesh points to ALL raster locations
  A_proj_raster <- INLA::inla.spde.make.A(
    mesh = mesh_s,
    loc = raster::rasterToPoints(x=index_raster)[, 1:2]
  )
  # Get range term
  # Get matrices used to make INLA precision matrices
  spde <- INLA::inla.spde2.matern(mesh = mesh_s,  alpha = 2, constr = FALSE)

  # Return all spatial matrix objects
  return(list(
    mesh_s = mesh_s, A_proj_data = A_proj_data, A_proj_raster = A_proj_raster, spde = spde
  ))
}


## FUNCTION: Execute TMB model
setup_run_tmb <- function(
  tmb_data_stack, params_list, tmb_random, tmb_map, optimization_method, template_file
){
  # Compile and load C++ code
  TMB::compile(template_file, flags='-w')
  setwd(dirname(template_file))
  fp_no_ext <- tools::file_path_sans_ext(basename(template_file))
  dyn.load(TMB::dynlib(fp_no_ext))
  # Create objective function
  obj <- TMB::MakeADFun(
    data = tmb_data_stack, parameters = params_list, random = tmb_random,
    map = tmb_map, DLL = fp_no_ext
  )
  # Optimize
  tictoc::tic("  Optimization")
  opt <- optimx::optimx(
    par = obj$par, fn = obj$fn, gr = obj$gr,
    method = optimization_method,
    itnmax = 3000,
    hessian = FALSE,
    control = list(trace = TRUE, dowarn = 0, maxit = 3000, starttests = FALSE, kkt = FALSE)
  )
  conv_code <- opt$convcode
  if(conv_code != 0) stop("Optimization failed with code", conv_code)
  tictoc::toc()
  # TODO
  return(list(obj = obj, opt = opt))
}


## FUNCTION: Take multivariate normal draws given a mean vector and precision matrix
#' @param mu vector of parameter means
#' @param prec joint precision matrix
#' @param n.sims number of draws
#' @return length(mu) by n.sims matrix of parameter draws
rmvnorm_prec <- function(mu, prec, n.sims) {
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv = Matrix::Cholesky(prec, super=TRUE)
  return(mu + solve(as(L_inv, 'pMatrix'), solve(t(as.matrix(as(L_inv, 'Matrix'))), z)))
}


## FUNCTION: Predict draw-level outputs by pixel
generate_draws <- function(
  tmb_sdreport, num_draws, covar_names, covar_rasters, index_raster, A_proj_raster,
  spatial_re = TRUE
){
  # Check that joint precision matrix was created
  if(!"jointPrecision" %in% names(tmb_sdreport)) stop("Missing joint precision matrix")

  # Get names and mean values for all parameters
  mu <- c(tmb_sdreport$par.fixed, tmb_sdreport$par.random)
  parnames <- names(mu)
  # Generate posterior predictive draws for all parameters
  message(sprintf(" - Generating %i parameter draws...", num_draws))
  param_draws <- rmvnorm_prec(
    mu = mu,
    prec = tmb_sdreport$jointPrecision,
    n.sims = num_draws
  )
  rownames(param_draws) <- parnames

  ## Prepare full table of covariates
  covars_dt <- data.table(cbind(1, values(covar_rasters), values(index_raster)))
  colnames(covars_dt) <- c('intercept', names(covar_rasters), 'pixel_idx')
  covars_dt <- na.omit(covars_dt)

  # Fixed effects
  fes <- as.matrix(covars_dt[, ..covar_names]) %*% param_draws[parnames=='beta_covs', ]
  # Spatial random effect
  if(spatial_re){
    res <- A_proj_raster %*% param_draws[parnames == 'z_space', ]
  } else {
    res <- 0
  }
  # Output probability = expit(fes + res)
  logit_prob <- fes + res
  prob <- exp(logit_prob) / (1 + exp(logit_prob))

  return(list(
    param_draws=param_draws, fe_draws=fes, re_draws=as.matrix(res),pixel_draws=as.matrix(prob)
  ))
}


## FUNCTION: Summarize pixel draws into (mean, median, lower, upper)
summarize_pixel_draws <- function(pixel_draws, index_raster){
  summ_stats <- c('mean', 'lower', 'median', 'upper')
  summary_dt <- data.table(cbind(
    rowMeans(pixel_draws),
    matrixStats::rowQuantiles(pixel_draws, probs = c(0.025, 0.5, 0.975))
  ))
  colnames(summary_dt) <- summ_stats
  # Convert to 4 raster layers
  summ_raster_list <- lapply(
    summ_stats,
    function(summ_name) vec_to_raster(summary_dt[[summ_name]], index_raster)
  )
  names(summ_raster_list) <- summ_stats
  return(summ_raster_list)
}
