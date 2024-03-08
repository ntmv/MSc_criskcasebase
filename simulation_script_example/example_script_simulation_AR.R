<<<<<<< HEAD
########################### n = 400, p = 120, Tp = 10  ######################
=======
########################### N = 400, p = 20 bias simulation  ######################
>>>>>>> 1d82361 (new sample script)
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
<<<<<<< HEAD
library(pamr)

# Fitting functions
source("src/fitting_functions.R")
source("src/helper_functions.R")
=======
library(survminer)
>>>>>>> 1d82361 (new sample script)


<<<<<<< HEAD
=======

# Set seed
>>>>>>> 1d82361 (new sample script)
seed <- as.integer(Sys.time())
# take the last five digits of the initial seed
the_seed = seed %% 100000

n = 400
p = 20
N = 1

<<<<<<< HEAD
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
=======
# Simulate data
sim.data <- cause_hazards_sim(n = n, p = p, 
                              beta1 = beta1, beta2 = beta2, rate_cens = 0.25, 
                              h1 = 0.55, h2 = 0.05, gamma1 = 1.5, gamma2 = 1.5, exchangeable = TRUE)
>>>>>>> 1d82361 (new sample script)


############################## FULL SIMULATION CODE ######################

<<<<<<< HEAD
results_table = runCasebaseSim(n, p, N, nfolds = 5, seed = the_seed)
summarizedSimCaseBaseTable = formatCaseBaseTable(results_table)
write.csv(summarizedSimCaseBaseTable,
          file = paste("simulation_script_example/results/", as.character(runif(1)), "test_results_N=", as.character(N), "_p=",
          as.character(p), ".csv", sep = ""))
=======

# Let us plot and visualize the competing risks curves
#cif <- cuminc(ftime = sim.data$ftime, fstatus = sim.data$fstatus)

#ggcompetingrisks(cif)

# Training-test split 
# We only do this (instead of generating datasets for train and test like Anthony mentioned because it is faster computationally 
# as casebase resamples) + proportion of censoring can be quite random in each run of the simulation so we want to maintain the same in validation and test set
train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
train <- sim.data[train.index,]
test <- sim.data[-train.index,]
>>>>>>> 1d82361 (new sample script)

############################## END SIMULATION CODE ######################

<<<<<<< HEAD
=======
x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 

# Censor competing event
y_test <- Surv(time = test$ftime, event = test$fstatus == 1)

x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 

# Fit cause-specific cox model with glmnet on training set 
cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7)

# Fit on validation set 
cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                      lambda = cox_mod$lambda.min)

cc_min <- coef(cox_val_min)

res_cox_min1 <- varsel_perc(cc_min, beta1)

# let's calculate the bias for all the competitors as well (a task could be turning this one line into a function as well)
# Only for the true non-zero variables
cox_mse1 <- mse_bias(cc_min, beta1)
########################## Cause 2 #####################################
# Censor competing event
y_train <- Surv(time = train$ftime, event = train$fstatus == 2)

x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 

# Censor competing event
y_test <- Surv(time = test$ftime, event = test$fstatus == 2)

x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 

# Fit cause-specific cox model with glmnet on training set 
cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7)

# Fit on validation set 
cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                      lambda = cox_mod$lambda.min)

cc_min <- coef(cox_val_min)

# Function to output variable selection performance metrics
res_cox_min2 <- varsel_perc(cc_min, beta2)

########################## Fit PenCR model ##################################
penCR = cv.glmnet.CR(data = train, family="cox", alpha= 0.7, standardize= TRUE,
                     nlambda = 20, t.BS = median(train$ftime), seed = 115, causeOfInt = 1,
                     nfold = 5)

cc_min_penCR1 <- penCR$glmnet.fits$models$`Cause 1`$glmnet.res$lambda[penCR$min.index[1]]
cc_min_penCR2 <- penCR$glmnet.fits$models$`Cause 2`$glmnet.res$lambda[penCR$min.index[2]]

# Fit on validation set 
penCR_val_min1 <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                         lambda = cc_min_penCR1)

cc_min_penCR1 <- coef(penCR_val_min1)

penCR_val_min2 <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                         lambda = cc_min_penCR2)


cc_min_penCR2 <- coef(penCR_val_min2)

res_pencr_min1 <- varsel_perc(cc_min_penCR1, beta1)
res_pencr_min2 <- varsel_perc(cc_min_penCR2, beta2)

# Calculate MSE here as well (try and fill it out!)
penCR_mse1 <-mse_bias(cc_min_penCR1, beta1)
########################## Fit casebase model #################################
# Train case-base model through cross-validation 
tic()
cv.lambda <- mtool.multinom.cv_cluster(train, seed = 1, nfold = 5, alpha = 0.7)
toc()

# Test set 
surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)

# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = cv.lambda$lambda.min , alpha = 0.7, unpen_cov = 2)


res_cb_min1 <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)

res_cb_min2 <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 2], beta2)

# Calculate MSE here as well
casebase_mse1 <-mse_bias(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)
############################ Case-base with post-"elastic net" ###########################
cb_postenet <- multinom.post_enet(fit_val_min, cause = 1)

# Calculate MSE here as well
casebase_mse_enet1 <-mse_bias(coef = cb_postenet$coef_selected[,1][1:6], true_coefs = beta1[cb_postenet$non_zero_coefs])
############################# Format and export results ###############################
Res_bias <- rbind(cox_mse1, penCR_mse1, casebase_mse1, casebase_mse_enet1)

coef_all <- rbind(as.vector(cc_min), as.vector(cc_min_penCR1), 
        as.vector(fit_val_min$coefficients[1:eval(parse(text="p")), 1]), cb_postenet$coefs_all)

colnames(coef_all) <- paste("X", 1:eval(parse(text="p")), sep = "")

rownames(Res_bias) <- c("iCR.bias", "penCRbias", "CB-bias", "postCB-bias")

colnames(Res_bias) <- "Coefficient Bias"

write.csv(Res, file = glue("bias_setting1{runif(1)}.csv"))

write.csv(coef_all, file = glue("coefficients_setting1{runif(1)}.csv"))
>>>>>>> 1d82361 (new sample script)
