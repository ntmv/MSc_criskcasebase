library(casebase)
library(future.apply)
library(glmnet)
#library(mtool)
library(parallelly)
library(timereg)
library(parallel)
library(tictoc)
library(tidyverse)
#library(riskRegression)
library(cmprsk)
library(survsim)
library(caret)
library(Matrix)
library(dplyr)

# Helper functions 
source("src/helper_functions.R")

# Run simulation
sim_results = runSim(20, 400, 1000, TRUE)


# Get coefficient biases
bias_table = formatCoefficientBiasTable(sim_results)

# Get average test prediction MSE
prediction_MSE_table = formatAverageTestMSETable(sim_results)

write.csv(bias_table, file = "results/iid_coefficient_relaxed.csv")
write.csv(prediction_MSE_table, file = "results/iid_MSE_relaxed.csv")
