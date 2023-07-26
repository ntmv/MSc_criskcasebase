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
source("../src/helper_functions.R")

n = 400
p = 120
N = 1000

# Run simulation
sim_results = runSim(p, n, N, TRUE)


# Get coefficient biases
bias_table = formatCoefficientBiasTable(sim_results, p)

# Get average test prediction MSE
prediction_MSE_table = formatAverageTestMSETable(sim_results)

write.csv(bias_table, file = paste("results/", as.character(runif(1)), "iid_coefficient_relaxed.csv"))
write.csv(prediction_MSE_table, file = paste("results/", as.character(runif(1)), "_iid_MSE_relaxed.csv"))
