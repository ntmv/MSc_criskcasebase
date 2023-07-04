# 💡 Penalized Competing Risks Analysis using Case-base Sampling

📍 Repository with my code and my final report for my MSc. thesis

🎯 We develop an alternative penalized cause-specific hazards model that extends on the `casebase` package for competing risks survival analysis of high-dimensional biological data

## 📂 Repository organization 

- `bash_script_template` - Template for bash scripts for Compute Canada. Creates individual bash scripts to be submitted to the cluster for one simulation run. Requires `runscripts` and `logs` folder to be setup

- `doc` - Thesis write up using `ubcdiss` template

- `Final_results` - Contains all final scripts,figures and results generated for the thesis for variable selection and prediction performance

- `mtool` - `R` package for fitting a penalized multinomial regression model with the K-1 logit parameterization 

- `mtool_fit_diagnostics` - Documentation of known issues with `mtool`

- `Papers` - Papers for background information. To be updated. 

-`simulation_scripts` - Template script for simulation from a two-cause model using the `replicate` function 

- `src` - Fitting functions and survival performance metrics functions 

- `survsim_mod` - Contains modified funtions from the `survsim` `R` package to generate competing risks data with normally distributed covariates having 1) AR(1) correlation structure and 2) Block correlation structure

- `updates` - Weekly meeting updates

# 📃 Abstract

Genome-wide transcriptome profiling and advances in experimental technologies have greatly increased the generation of high-dimensional genomic data, particularly microarray data correlated with survival outcomes such as patient survival time or time to cancer relapse. The analysis of such genomic time-to-event data become more complicated when there are competing events, i.e., the failure of a patient can occur due to one of multiple distinct causes. We develop an alternate elastic-net penalized competing risks analysis method that is able to produce easily interpretable hazard ratios akin to the Cox regression model, and also able accurately produce smooth-in-time predicted estimates patient risk, in a variety of settings, such as non-proportional hazards as well. We examine the performance of this method in a simulation study. 

🎺 The original model is based on the [casebase](https://arxiv.org/abs/2009.10264) paper by Sahir Bhatnagar, Maxime Turgeon, Jesse Islam and James Hanley, which is based on the [sampling methodology](https://www.degruyter.com/document/doi/10.2202/1557-4679.1125/html?lang=en) proposed by James Hanley and Olli Miettinen. 

💻 The optimization for the stochastic gradient descent in `mtool` was written by Dr. Yi Lian.
