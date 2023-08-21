########################### n = 400, p = 120, Tp = 60  ######################
library(casebase)
library(future.apply)
library(glmnet)
library(devtools)
# install.packages("mtool_1.0.tar.gz")
library(mtool)
library(parallel)
library(tictoc)
library(tidyverse)
library(foreach)
library(survival)
library(cmprsk)
library(glue)
library(pec)
library(pamr)

# Fitting functions
source("src/fitting_functions.R")
source("src/helper_functions.R")


seed <- as.integer(Sys.time())
# take the last five digits of the initial seed
the_seed = seed %% 100000

n = 400
p = 120
N = 5

############################## TEST CODE ######################

# set.seed(the_seed)
# # Set seed
# 
# num_true <- 20
# beta1 <- c(rep(0, p))
# beta2 <- c(rep(0, p))
# nu_ind <- seq(num_true)
# # Here out of 20 predictors, 10 should be non-zero
# beta1[nu_ind] <- c(rep(1, p/2), rep(0, p/2))
# beta2[nu_ind] <- c(rep(-1, p/2), rep(0, p/2))
# 
# # Simulate data
# sim.data <- cause_hazards_sim(n = n, p = p, nblocks = 4,
#                               beta1 = beta1, beta2 = beta2, rate_cens = 0.25,
#                               h1 = 0.55, h2 = 0.10, gamma1 = 1.5, gamma2 = 1.5)
# 
# 
# # Censoring proportion
# cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)
# 
# # Training-test split
# # We only do this (instead of generating datasets for train and test like Anthony mentioned because it is faster computationally
# # as casebase resamples) + proportion of censoring can be quite random in each run of the simulation so we want to maintain the same in validation and test set
# train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
# train <- sim.data[train.index,]
# test <- sim.data[-train.index,]
# 
# start_cv <- Sys.time()
# res_cv <- mtool.multinom.cv(train = train, seed = seed, nfold = 5)
# end_cv <- Sys.time()
# start_post <- Sys.time()
# res_post <- multinom.post_enet_old(train = train, test = test, nfold = 5, seed = seed)
# end_post <- Sys.time()
# start_relaxed <- Sys.time()
# res_relaxed = multinom.relaxed_enet(train = train, nfold = 5, seed = seed)
# end_relaxed <- Sys.time()
# 
# print(paste(end_cv - start_cv, " taken to run multinomial cv"))
# print(paste(end_post - start_post, " taken to run post LASSO"))
# print(paste(end_relaxed - start_relaxed, " taken to run relaxed LASSO"))
# 
# # Write to csv
# write.csv(as.data.frame(res_relaxed$deviance_grid),
#           file = paste("simulation_script_example/results/", as.character(runif(1)), "deviance_grid.csv"))
############################## END TEST CODE ######################


############################## FULL SIMULATION CODE ######################

results_table = runCasebaseSim(n, p, N, nfolds = 5, seed = the_seed)
summarizedSimCaseBaseTable = formatCaseBaseTable(results_table)
write.csv(summarizedSimCaseBaseTable,
          file = paste("simulation_script_example/results/", as.character(runif(1)), "test_results_N=", as.character(N), "_p=",
          as.character(p), ".csv", sep = ""))

############################## END SIMULATION CODE ######################

