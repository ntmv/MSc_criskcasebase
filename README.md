# 💡 Penalized Competing Risks Analysis using Case-base Sampling

📍 Repository with my (Nirupama Tamvada) code and my final report for my MSc. thesis

🎯 We develop an alternative penalized cause-specific hazards model that extends on the `casebase` package for competing risks survival analysis of high-dimensional biological data

## 📂 Repository organization

- `bash_script_template` - Template for bash scripts for Compute Canada. Creates individual bash scripts to be submitted to the cluster for one simulation run. Requires `runscripts` and `logs` folder to be setup

- `doc` - Thesis write up using `ubcdiss` template

- `Final_results` - Contains all final scripts,figures and results generated for the thesis for variable selection and prediction performance

- `mtool` - `R` package for fitting a penalized multinomial regression model with the K-1 logit parameterization 

- `mtool_fit_diagnostics` - Documentation of known issues with `mtool`

- `Papers` - Papers for background information. To be updated. 

- `simulation_scripts` - Template script for simulation from a two-cause model using the `replicate` function 

- `src` - Helper fitting functions including the cross-validation function, simulation functions, and survival performance metrics functions 

- `survsim_mod` - Contains modified funtions from the `survsim` `R` package to generate competing risks data with normally distributed covariates having 1) AR(1) correlation structure and 2) Block correlation structure

- `updates` - Weekly meeting updates

# 📃 Abstract

Genome-wide transcriptome profiling and advances in experimental technologies have greatly increased the generation of high-dimensional genomic data, particularly microarray data correlated with survival outcomes such as patient survival time or time to cancer relapse. The analysis of such genomic time-to-event data becomes more complicated when there are competing events, i.e., the failure of a patient can occur due to one of multiple distinct causes. We develop an alternate elastic-net penalized competing risks analysis method that is able to produce easily interpretable hazard ratios akin to the Cox regression model. This approach is also able accurately produce smooth-in-time predicted estimates patient risk, in a variety of settings, such as non-proportional hazards as well. We examine the performance of this method in a simulation study, looking at both variable selection, as well as patient risk estimation performance in both the $N > p$ as well as the $p > N$ scenarios.

🎺 The original model is based on the [casebase](https://arxiv.org/abs/2009.10264) paper by Sahir Bhatnagar, Maxime Turgeon, Jesse Islam and James Hanley, which is based on the [sampling methodology](https://www.degruyter.com/document/doi/10.2202/1557-4679.1125/html?lang=en) proposed by James Hanley and Olli Miettinen. 

💻 The optimization for the stochastic gradient descent in `mtool` was written by Dr. Yi Lian.

📌 The relaxed LASSO branch contains a relaxed LASSO implementation for the casebase penalized model (WIP) by Alex Romanus.



## 💡 Relaxed LASSO as a Penalization Method for Competing Risks Analysis using Case-Base Sampling

📍 relaxed_LASSO branch contains the addition of Relaxed LASSO penalization for fitting Competing Risks models using high-dimensional data

🎯 To address issues posed by fitting models to high-dimensional data, such as biased coefficient estimates, we developed a version of Relaxed LASSO regularization for efficient variable selection during model fitting

### 📂 Relaxed_LASSO-specific branch organization

- `practice` - Contains files created during development of relaxed LASSO method
  - `development` - Files used to edit and check results of relaxed LASSO implementations for each data type (linear, multinomial, and Case-base) during development
  - `simulation` - Mock simulation scripts and results

- `renv` - Virtual environment directory containing all packages needed to run relaxed LASSO project files (except mtool - more on that below)

- `simulation_final` - Scripts and results from final simulation used in report for each data type
 - Note: Data in `final_results` sub-directories are results averaged over all iterations of each simulation for each data type
 
- `src`
  - `final_relaxed_implementation.R` - Contains implementation of Relaxed LASSO used in final simulations in multinomial and Case-base data settings and its helper functions
  - `practice_relaxed_implementation`- All versions of Relaxed LASSO implementations made during development and their helper functions
  - simulation_helper_functions.R files - Helper functions intended for use during simulation, but replaced with scripts in final simulation runs (except in linear case)
  

### 💻 Instructions on running files from Relaxed LASSO project and installing repository packages
- To run any file and install the packages properly, please set your working directory to the project directory (casebase_relaxed_LASSO)

- Installing packages:
  - Published packages - To install most of the packages needed for this project, use renv to install the virtual environment.
    
    1) Install renv using the following command in your RStudio console:
    
      install.packages("renv")
    
    2) Use the following command to install all published packages:
    
      renv::restore()
    
  - mtool package - mtool is the only unpublished package required in this project, and hence must be downloaded from source. Use the following lines to install the package:
  
      install_path = paste(getwd(), "/mtool_1.0.tar.gz", sep = "")
      install.packages(install_path, repos = NULL, type = "source")
    
    Finally, restart R to run the files using the required packages with:
    
      .rs.restartR()