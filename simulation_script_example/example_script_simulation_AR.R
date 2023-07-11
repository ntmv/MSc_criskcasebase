########################### n = 400, p = 120, Tp = 20  ######################
library(casebase)
library(future.apply)
library(glmnet)
library(mtool)
library(parallel)
library(tictoc)
library(tidyverse)
library(foreach)
library(survival)
library(cmprsk)
library(glue)
library(pec)

# Fitting functions 
source("../src/fitting_functions.R")

# Control number of simulations here
# For each simulation you want to output the MSE
Results <- replicate(5, {
# Set seed
seed <- as.integer(Sys.time())

# take the last five digits of the initial seed
the_seed= seed %% 100000
set.seed(the_seed)

# Setup
n <- 400
# Let's start yours with p = 20
p <- 20
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

##############################################################
# We have two competitor models for variable selection:
# 1) Independent cox-regression model 
# 2) penCR cox regression model - where the lambda penalties are trained together 
######################## Fit indepedent cox-regression model ###############################
######################### Cause-1 #########################################
# Censor competing event
y_train <- Surv(time = train$ftime, event = train$fstatus == 1)

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
mean((cc_min[nu_ind] - beta1[nu_ind])^2)
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

# let's calculate the bias for all the competitors as well (a task could be turning this one line into a function as well)
cc_min_bias <- which(coef(cox_val_min) != 0)

# MSE (bias)
mean((cc_min[cc_min_bias] - beta2[cc_min_bias])^2)
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

# Free up memory for case-base
rm(penCR)

# Calculate MSE here as well (try and fill it out!)
#...
########################## Fit casebase model #################################
# Train case-base model 
cv.lambda <- mtool.multinom.cv(train, seed = 1, nfold = 5)

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

# Calculate MSE here as well!
#....
###################################################################################
Res <- rbind(res_cb_min1, res_cb_min2, res_cox_min1, res_cox_min2, res_pencr_min1, res_pencr_min2, cen.prop)

rownames(Res) <- c("casebase.lambda.min_cause1", "casebase.lambda.min_cause2", "cox.lambda.min_cause1",
                   "cox.lambda.min_cause2", "pencr.lambda.mincause1", "pencr.lambda.mincause2", "cens.prop")

Res

}, simplify = FALSE)

Results <- do.call(rbind, Results)


############ Sketch of function for post-LASSO (or post elastic net in this case) #########
# Look into ... argument to pass parameters from other functions because you want to pass cross-validation parameters
multinom.post_enet <- function(train, test) {
  # Train case-base model 
  cv.lambda <- mtool.multinom.cv(train, seed = 1, nfold = 5)
  # This fit (with lambda.min) needs to be de-biased
  # Fit on test set 
  # Covariance matrix
  cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
  
  # Case-base dataset
  cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)
  
  # Case-base fits 
  # Lambda.min
  fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                             lambda = cv.lambda$lambda.min , alpha = 0.7, unpen_cov = 2)
  
  # Obtain all non-zero selected covariates across both classes
  non_zero_coefs_cause1 <- which(fit_val_min$coefficients[1:eval(parse(text="p")), 1] != 0)
  non_zero_coefs_cause2 <- which(fit_val_min$coefficients[1:eval(parse(text="p")), 2] != 0)
  # Combine them 
  non_zero_coefs <- union(non_zero_coefs_cause1, non_zero_coefs_cause2)
  non_zero_coefs <- paste("X", non_zero_coefs , sep = "")
  # Create new subsetted dataset 
  testnew <- cbind(test[, c(colnames(test) %in% non_zero_coefs)], ftime = (test$ftime), fstatus = test$fstatus)
  # Fit "OLS" (unparameterized multinomial model)
  # For working of this function see: http://sahirbhatnagar.com/casebase/articles/competingRisk.html
  model_cb <- fitSmoothHazard(fstatus ~. +log(ftime) -fstatus,
                              data = testnew,
                              ratio = 100,
                              time = "ftime")
  # Return coefficients
}

# Calculate MSE for this as well (fill in here)
#...