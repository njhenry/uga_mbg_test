# Uganda MBG model example

Purpose: Estimate proportion of adult women with no prior education by 5x5 km pixel across Uganda. Note that this is a continuous spatial model for a single year of covariate and outcome data.

Note that although this repository was tested on education in Uganda, the same estimation framework could be adapted to other outcomes and country settings.


## Setup

This repository is split into the following files:

```
config.yaml: File paths and model settings that may change by user
scripts/
  01_run_uga_model.R: Script, callable from the command line, that runs a geostatistical model with covariates to estimate prior education across space
  02_uga_postestimation.R: Create maps and calculate summary statistics and predictive validity metrics for a given model run
  helper_functions.R: Functions called by the two executable scripts
  tmb_template.cpp: Template Model Builder (TMB) objective function for the geostatistical model

```

#### Dependencies

This code depends on the packages specified under "Settings > Load_Libraries" specified in the `config.yaml` file.


## Input data

The script expects the following inputs:


#### CSV of survey points, numerator, and denominator

A CSV of point-geolocated survey data should contain the following fields:
- 'latitude' (float)
- 'longitude' (float)
- 'edu0' (numerator, float)
- 'N' (denominator, float, strictly >= 'edu0')


#### Raster covariate data

A raster brick, in .tif format, where each layer corresponds to a named covariate that may be associated with the outcome in the same year that data was observed.
