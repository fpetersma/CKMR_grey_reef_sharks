# CKMR_grey_reef_sharks

## Introduction

This folder contains all scripts used for a close-kin mark-recapture simulation study. 
These scripts underpin the analysis presented an article currently under review for Ecology and Evolution, and of Chapter 4 of my PhD thesis. 
The study is collaborative effort with Len Thomas, Danielle Harris, Darcy Bradley and Yannis Papastamatiou.

If anything remains unclear about the files, running of the scripts, and reproducibility of the results, please contact me at felix.petersma@gmail.com or ftp@st-andrews.ac.uk. 

## Overview

This repository contains several folders, only some of which are relevant for the formal analysis. I describe those folders below:
+ images_no_growth: a folder containing images associated with the results from fitting a model without population growth.
+ images_old_years: a folder with images using an old year structure that is no longer relevant
+ images_with_growth: a folder containing images associated with the results from fitting a model with population growth.
+ reports: a folder containing intermediate analyses and reports as part of the research process.
+ source: 
	+ CKMRcpp: an R package that includes functionality for fitting of a CKMR model using R and C++.
	+ fitting: scripts associated with the fitting of models.
	+ functions: old functions that were replaced by CKMRcpp.
	+ general: general functionality not contained in CKMRcpp.
	+ simulating: functionality for running the simulation using CKMRcpp underpinning the analysis.
+ result_summaries: Rdata files that contain data created in steps 4--6 below. 

Some folders also contain a subfolder called "extra". This folder mostly contains scripts that have become redundant but have not been removed yet. The contents of these folders can be ignored, unless indicated otherwise. 

## Simulation procedure

Due to the size of the data files (several gigabites per datafile) these were not shared on github. 
Instead, I explain the simulation process below with the seeds that were used to produce the data
underpinning the analysis. The scripts were used on a machine with 30+ cores for parallel processing. Be careful running them on any personal machine and make sure to first alter scripts to match your machine's specifications.

The simulation procedure was as follows:
1. Data was simulated using "simulating/simple_sims.R" for the simple population and "simulating/complex_sims.R" for the complex population. The seeds were hardcoded in (see l.49 in "simple_sims.R" and l.57 in "complex_sims.R"). 
2. The script "simulating/finding_pairs_bluewhale.R" was used to extract the kin pairs from the data simulated in step 1. This script was written to run in parallel (similar to the previous two scripts) on a server with 30+ cores. 
3. Next, the script "fitting/fit_and_visualise_scenarios.R" is used to fit the CKMR models to the simulated data. 
4. Descriptive statistics were extracted using the functions "general/extract_descriptives.R" and "general/some_simulation_descriptives.R".
5. An overview of maximum likelihood variance estimates was created using "general/estimate_variance_mle.R"
6. The performance metrics were derived and overviews/tables were created in "general/extract_performance_metrics_new.R". 

## Visualisations

The visualisations for the analysis were created using the followin scripts:
+ general/plot_abundance_through_time.R, to simulated the historic abundances for the fitted models
+ general/plot_growth_curves.R, to plot the true growth curve with deviations and true measurement error with deviations