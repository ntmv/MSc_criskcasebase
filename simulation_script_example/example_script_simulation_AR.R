########################### N = 400, p = 20 bias simulation  ######################
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
library(survminer)

# Fitting functions 
source("../src/fitting_functions.R")


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
sim.data <- cause_hazards_sim(n = n, p = p, 
                              beta1 = beta1, beta2 = beta2, rate_cens = 0.25, 
                              h1 = 0.55, h2 = 0.05, gamma1 = 1.5, gamma2 = 1.5, exchangeable = TRUE)


# Censoring proportion
cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)


# Let us plot and visualize the competing risks curves
#cif <- cuminc(ftime = sim.data$ftime, fstatus = sim.data$fstatus)

#ggcompetingrisks(cif)

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
cox_mse1 <- mse_bias(cc_min, beta1)

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

res_pencr_min1 <- varsel_perc(cc_min_penCR1, beta1)


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
