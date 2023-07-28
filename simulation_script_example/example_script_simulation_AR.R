########################### n = 400, p = 120, Tp = 20  ######################
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


############################### TEST CODE ######################
# res1 <- multinom.post_enet_old(train, test)
# res2 <- multinom.post_enet(train, test)
# res1_cb_post_lasso <- varsel_perc(res1$coefficients, beta1)
# res2_cb_post_lasso <- varsel_perc(res2$coefficients, beta1)


n = 400
p = 20
N = 5
  # Set seed
seed <- as.integer(Sys.time())

# take the last five digits of the initial seed
the_seed= seed %% 100000
set.seed(the_seed)

num_true <- 20
beta1 <- c(rep(0, p))
beta2 <- c(rep(0, p))
nu_ind <- seq(num_true)
# Here out of 20 predictors, 10 should be non-zero 
beta1[nu_ind] <- c(rep(1, 10), rep(0, 10))
beta2[nu_ind] <- c(rep(-1, 10), rep(0, 10))

# Simulate data
sim.data <- cause_hazards_sim(n = n, p = p, nblocks = 4, 
                              beta1 = beta1, beta2 = beta2, rate_cens = 0.25, 
                              h1 = 0.55, h2 = 0.10, gamma1 = 1.5, gamma2 = 1.5)
  
  
  # Censoring proportion
cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)

# Training-test split 
# We only do this (instead of generating datasets for train and test like Anthony mentioned because it is faster computationally 
# as casebase resamples) + proportion of censoring can be quite random in each run of the simulation so we want to maintain the same in validation and test set
train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
train <- sim.data[train.index,]
test <- sim.data[-train.index,]

res_relaxed = multinom.relaxed_enet(train = train, seed = 2023)


# Write to csv
write.csv(as.data.frame(res_relaxed$deviance_grid), file = paste("simulation_script_example/results/", as.character(runif(1)), "deviance_grid.csv"))

# results_table = runCasebaseSim(400, 20, 5)
# summarizedSimCaseBaseTable = formatCaseBaseTable(results_table)


