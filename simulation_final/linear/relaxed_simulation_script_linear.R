library(casebase)
library(future.apply)
library(glmnet)
library(mtool)
library(parallelly)
library(timereg)
library(parallel)
library(tictoc)
library(tidyverse)
library(cmprsk)
library(survsim)
library(caret)
library(Matrix)
library(dplyr)
library(glue)

# Helper functions 
source("src/linear_simulation_helper_functions.R")

n = 400
p = 20
N = 100
# Run simulation
sim_results = runSim(p, n, N, TRUE)


# Get coefficient biases
bias_table = formatCoefficientBiasTable(sim_results, p)

# Get average test prediction MSE
prediction_MSE_table = formatAverageTestMSETable(sim_results)

write.csv(prediction_MSE_table, file = glue("simulation_final/linear/final_simulation_results/summarized_results_MSE_{runif(1)}.csv"))

coefficient_MSEs <- (rowMeans(bias_table[, c(2:(p+1))]^2))

write_csv(as.data.frame(coefficient_MSEs), file = glue("simulation_final/linear/final_simulation_results/summarized_results_coefficients_{runif(1)}.csv"))