// -------------------------------------------------------------------------------------->
//
// TMB template for Uganda education spatial modeling
// Author: Nat Henry
// Created: 22 April 2021
// Purpose: Estimate proportion of women who have never attended school by 5x5km gridded
//  pixel across Uganda.
//
// NOTE: This is a TMB model template to be called from R
//
// -------------------------------------------------------------------------------------->

// include libraries
#include <TMB.hpp>
using Eigen::SparseMatrix;


// FUNCTION: Create a spatial precision matrix corresponding with an SPDE model
//
// parameter log_kappa: Log of SPDE spatial range term
// parameter log_tau: Log of SPDE spatial variance term
// parameters M0, M1, M2: Sparse matrix outputs from INLA::inla.spde2.matern()
//
// Output: A sparse precision matrix to evaluate against lat-long data
//
template<class Type>
SparseMatrix<Type> spde_precision(
  Type log_kappa, Type log_tau, SparseMatrix<Type> M0, SparseMatrix<Type> M1,
  SparseMatrix<Type> M2
){
  // Q = tau^2 * (kappa^4 * M0 + 2 * kappa^2 * M1 + M2)
  SparseMatrix<Type> Q = (
    exp(log_tau * 2.) *
    (exp(4. * log_kappa) * M0 + 2. * exp(2. * log_kappa) * M1 + M2)
  );
  return Q;
}


// OBJECTIVE FUNCTION: returning joint negative log-likelihood (jnll)
template<class Type>
Type objective_function<Type>::operator() () {

  // 1) Data inputs --------------------------------------------------------------------->

    // Holdout number for this run
    DATA_INTEGER(holdout);

    // Data observations and metadata
    DATA_VECTOR(no_edu_i); // Numerator: number with no education
    DATA_VECTOR(n_i); // Denominator
    DATA_MATRIX(X); // Covariate design matrix
    DATA_IVECTOR(idx_holdout); // Holdout identifier for each observation

    // Projection matrix: SPDE mesh nodes -> observations
    DATA_SPARSE_MATRIX(A_projmat);

    // Objects used to calculate precision matrix
    DATA_SPARSE_MATRIX(M0);
    DATA_SPARSE_MATRIX(M1);
    DATA_SPARSE_MATRIX(M2);


  // 2) Parameter inputs ---------------------------------------------------------------->

    // Fixed effects
    PARAMETER_VECTOR(beta_covs); // Vector of covariate fixed effects, including intercept
    PARAMETER(log_kappa); // Log of SPDE spatial range parameter
    PARAMETER(log_tau); // Log of SPDE spatial variance parameter

    // Spatial random effect
    PARAMETER_VECTOR(z_space);

  // 3) Transform data and parameters as needed ----------------------------------------->

    // Spatial range and variance objects
    Type kappa = exp(log_kappa);
    Type tau = exp(log_tau);

    // Spatial precision matrix, evaluated against mesh nodes
    SparseMatrix<Type> Q_spatial = spde_precision(log_kappa, log_tau, M0, M1, M2);

    // Spatial random effect projected onto observation points
    vector<Type> spatial_effects = A_projmat * z_space;
    // Fixed effects projected onto observation points
    vector<Type> fixed_effects = X * beta_covs.matrix();
    // Estimated underlying probabilities in logit (-Inf, Inf) and observation (0,1) space
    vector<Type> logit_prob = fixed_effects + spatial_effects;
    vector<Type> prob = exp(logit_prob)/(Type(1.0) + exp(logit_prob));


  // 4) Likelihood evaluation ----------------------------------------------------------->

    Type jnll = 0.0;

    // 4a) Likelihood of parameters given priors

      // N(0,3) priors on covariate fixed effects
      jnll -= dnorm(beta_covs, Type(0.0), Type(3.0), true).sum();
      // Wide gamma priors on tau and kappa
      jnll -= dlgamma(kappa, Type(1.0), Type(1000.0), true);
      jnll -= dlgamma(tau, Type(1.0), Type(1000.0), true);
      // Evaluate the density of the spatial random effects against the spatial precision
      //  matrix
      jnll += density::GMRF(Q_spatial)(z_space);

    // 4b) Likelihood of data given parameters

      for(int obs_i = 0; obs_i < no_edu_i.size(); obs_i++){
        // Only evaluate observations that are NOT in the holdout fold
        if(idx_holdout(obs_i) != holdout){
          jnll -= dbinom(no_edu_i(obs_i), n_i(obs_i), prob(obs_i), true);
        }
      }

    return jnll;
}
