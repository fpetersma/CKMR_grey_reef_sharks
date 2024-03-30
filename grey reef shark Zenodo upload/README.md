## Close-kin mark-recapture simulation study on the effects of incorrect ageing

Contains the code and data used to run the simulation study. 
The results and details of this study are described in the manuscripted "Age is not just a number: how incorrect ageing impacts close-kin mark-recapture estimates of population size", which is currently under review at Ecology and Evolution. The results also appear in Chapter 4 of the PhD thesis "Advancements in methods for estimating the abundance of marine megafauna using novel sampling techniques" by Felix T Petersma. 
We provide some information in support of the files presented here; however, for extensive background we kindly refer you to the PhD chapter or upcoming publication.

Most of the functionality of the model fitting is contained in the 'CKMRcpp'-package. 
This package has been included pre-built and can be installed using the file 'CKMRcpp_1.0.tar.gz'.
Once this package is installed, all scripts included in this location should run fine, given that all other packages that are used will be installed as well.

### General files

- extract_performance_metrics_new_no_growth.R
	Extracts the performance metrics used in the study from the fitted models. It uses data files created by the fitting and simulation scripts.
 	

### Simulation files

- simple_sims.R. 
	Contain the code to simulate 1000 populations with the demographic traits of the simple species.

- complex_sims.R. 
	Contain the code to simulate 1000 populations with the demographic traits of the complex species.

- simulated_data_with_recaptures.R
	Takes the simulated populations and prepares the data for the model fitting by e.g., extracting POPs and reducing the number of comparisons that require evaluation.

### Fitting files


### Result summaries

These are some RData files with the results that appear in the manuscript. 


### General files

There are some general scripts and data files that are used by the case study and the simulation study. These are as follows.

- LICENSE.md
	Information on the license under which these data are shared.

- nllR.R
	This R script derives the negative log likelihood of an acoustic spatial capture-recapture model for a set of parameters given the data. This script is prepares data and runs the C++ functionality from the 'ascrRcpp'-package. 

- alaska_albers_grid_adaptive_levels=2_inner=10k_outer=50k_maxD2C=Inf_area=8450_n=438.csv

	This file contains a spatial grid using the Albers projection. Every grid point contains information on the ocean depth (bathymetry) and the southward distance to the coast in meters. The grid contains of a high resolution inner grid for improved precision near the DASAR array, and a coarser outer grid farther out to improve runtime. 

- create_grid_Albers.R
	
	Create a spatial grid that is used for the fitting, based on the Albers projection.

### Case study

The following case study files are included. The data files contain the relevant information from 438 calls detected on at least two Directional Autonomous Seafloor Acoustic Recorders (DASARs) from DASAR array 5 on 31-08-2010, which consisted of 6 equally spaced DASARs. A DASAR not only records information about the frequency and loudness of calls, but also the direction the call came from. This means that for every call, information was recorded by every DASAR that was involved in the detection (at least 2, at most 6). This information consisted of the received sound level, the bearing estimated by the DASAR, and an overview of which DASAR detected a call and which one missed it. We also included a file containing information on the placement of the DASARs. We also included scripts to fit the 35 candidate density models to the case study data, and the script that used bootstrapping to estimate the uncertainty around the parameter estimates.

#### Data files

- detections_31-08-2010_successful_loc.csv

	A matrix of detections for every included call for every detector that was involved with the detection, where TRUE indicates a positive detection;

- received_levels_31-08-2010_successful_loc.csv

	A matrix of received sound levels for every included call for every detector that was involved with the detection, and NA otherwise;

- bearings_31-08-2010_successful_loc.csv

	A matrix of bearings for every included call for every detector that was involved with the detection, and -1 otherwise;

- DASARs.txt

	Placement information for the DASARs.	

- fits_1_35_nlminb_n=443.RData
	
	An RData-file containing the results of 35 ASCR models with varying density specifications fitted to the case study data.

- results_model_33_1_999.RData

	Relevant data extracted from fitting the best model (model 33) to the 999 bootstrapped data sets.


#### Scripts

- varying_fits_to_real_data.R
	
	This script fits 35 models

- bootstrap_real_data.R

	This scripts creates 999 bootstrapped data sets, and fits the best model to all.

- extracting_data_from_bootstraps.R
	
	This script extracts relevant information from the fits to the bootstrapped data.


### Simulation study

The study also included a simulation study to evaluate the performance of the method in various scenarios. We include the scripts that were used to simulate the data used in the study, as well as the scripts that perform the fitting. Moreover, we include the actual data files that contain the results from the simulation study, i.e., the simulated data and fit results.

- simulate_1000_data_sets.R

	Simulates data that can be used for simulation studies.
	
- fitting_to_simulated_data.R
	
	Fits models to simulated data for the simlation study.
	
- extracting_data_from_simulations.R

	This script extracts relevant information from the fits to the simulated data for various scenarios.


#### Data files

The following files are associated with the scenarios for which the data were simulated with a variable source level.

- 1000_simulations_fixedSL=FALSE.RData

	A thousand simulated data sets with a variable source level.

- fits_with_fixedSL=FALSE_USE_BEARINGS=0_1-100_on_simulations_with_fixedSL=FALSE.RData

	Results from fitting a model without using the bearing information to the first 100 simulated variable source level data sets.

- fits_with_fixedSL=FALSE_USE_BEARINGS=1_1-100_on_simulations_with_fixedSL=FALSE.RData

	Results from fitting a model with a single-precision bearing model to the first 100 simulated variable source level data sets.

- fits_with_fixedSL=FALSE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=FALSE.RData

	Results from fitting a correct model to the first 100 simulated variable source level data sets.

- fits_with_fixedSL=FALSE_USE_BEARINGS=2_constant_density_1-100_on_simulations_with_fixedSL=FALSE.RData

	Results from fitting a model with an incorrect homogeneous density model to the first 100 simulated variable source level data sets.

- fits_with_fixedSL=TRUE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=FALSE.RData

	Results from fitting a fixed source level model to the first 100 simulated variable source level data sets.


The following files are associated with the scenarios for which the data were simulated with a fixed source level.

- 1000_simulations_fixedSL=TRUE.RData

	A thousand simulated data sets with a fixed source level.

- fits_with_fixedSL=TRUE_USE_BEARINGS=0_1-100_on_simulations_with_fixedSL=TRUE.RData

	Results from fitting a model without using the bearing information to the first 100 simulated fixed source level data sets.

- fits_with_fixedSL=TRUE_USE_BEARINGS=1_1-100_on_simulations_with_fixedSL=TRUE.RData

	Results from fitting a model with a single-precision bearing model to the first 100 simulated fixed source level data sets.

- fits_with_fixedSL=TRUE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=TRUE.RData

	Results from fitting a correct model to the first 100 simulated fixed source level data sets.

- fits_with_fixedSL=TRUE_USE_BEARINGS=2_constant_density_1-100_on_simulations_with_fixedSL=TRUE.RData

	Results from fitting a model with an incorrect homogeneous density model to the first 100 simulated fixed source level data sets.

- fits_with_fixedSL=FALSE_USE_BEARINGS=2_1-100_on_simulations_with_fixedSL=TRUE_density_pars=c(-16,57,-68.5).RData

	Results from fitting a variable source level model to the first 100 simulated fixed source level data sets.


