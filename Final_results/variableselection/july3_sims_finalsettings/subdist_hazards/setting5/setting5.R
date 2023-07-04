########################### n = 400, p = 120, Tp = 20 b1 = b2 simulation ######################
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

set.seed(110)
# Fitting functions 
source("../../src/fitting_functions.R")

# Setup
n <- 400
p <- 120
num_true <- 50
beta1 <- c(rep(0, p))
beta2 <- c(rep(0, p))
nu_ind <- seq(num_true)
beta1[nu_ind] <- 1
beta2[nu_ind] <- 0.8

# Simulate data
sim.data <- cause_subdist_sim(n = n, p = p, nblocks = 5, 
                                  beta1 = beta1, mix_p = 0.6, 
                                  beta2 = beta2, u.max = 1.5, lambda1 = 1, 
                              lambda2 = 100, rho1 = 1.5, rho2 = 1.5, 
                              cor_vals = c(0.7, 0.4, 0.6, 0.5, 0.8))


# Censoring proportion
cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)

# Training-test split 
train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
train <- sim.data[train.index,]
test <- sim.data[-train.index,]
######################## Fit cox-regression model ###############################
# Censor competing event
y_train <- Surv(time = train$ftime, event = train$fstatus == 1)

x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 

# Censor competing event
y_test <- Surv(time = test$ftime, event = test$fstatus == 1)

x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 

# Fit cause-specific cox model with glmnet on training set 
cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.5)

# Fit on validation set 
cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 1, 
                      lambda = cox_mod$lambda.min)

cc_min <- coef(cox_val_min)

res_cox_min <- varsel_perc(cc_min, beta1)

########################## Fit casebase model #################################
# Split training set into train and validation
train.index <- caret::createDataPartition(train$fstatus, p = 0.70, list = FALSE)
train_new <- train[train.index,]
validation <- train[-train.index,]

# Fit case-base model 
# Convert to case-base dataset
surv_obj_train <- with(train_new, Surv(ftime, as.numeric(fstatus), type = "mstate"))

cov_train <- as.matrix(cbind(train_new[, c(grepl("X", colnames(train)))], time = log(train_new$ftime)))

# Create case-base dataset
cb_data_train <- create_cbDataset(surv_obj_train, cov_train, ratio = 20)

# Create case base dataset for validation
surv_obj_validation <- with(validation, Surv(ftime, as.numeric(fstatus), type = "mstate"))

cov_validation <- as.matrix(cbind(validation[, c(grepl("X", colnames(validation)))], time = log(validation$ftime)))

# Create case-base dataset
cb_data_validation <- create_cbDataset(surv_obj_validation, cov_validation, ratio = 20)


# Fit lambda grid to training set 
epsilon <- 0.001
lambda_max <- 0.1
grid_size <- 100
lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))

# Set the number of cores to be used for parallel processing
num_cores <- parallelly::availableCores()  # Adjust the number of cores as per your system's capacity
# Create a parallel cluster using the specified number of cores
cl <- parallel::makeCluster(num_cores, setup_strategy = "sequential")

# Load necessary packages on each cluster
clusterEvalQ(cl, {
  library(casebase)
  library(future.apply)
  library(mtool)
  library(parallel)
  library(dplyr)
})

# Export necessary objects to all clusters
objects_to_export <- c("lambdagrid", "fit_cbmodel", "multi_deviance", "cb_data_train")
clusterExport(cl, objects_to_export)

cv_res <- foreach(lambda_val = lambdagrid, .packages = "mtool") %dopar% {
  fit_cbmodel(cb_data_train, regularization = 'elastic-net', lambda = lambda_val, alpha = 0.5)
}

stopCluster(cl)

# Fit to validation test to choose lambda min
mult_deviance <- unlist(lapply(cv_res, multi_deviance, cb_data = cb_data_validation))

# Lambda min 
lambda.min <- lambdagrid[which.min(mult_deviance)]
#})

# Test set 
surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)

# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = lambda.min , alpha = 0.5, unpen_cov = 2)

res_cb_min <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)

###################################################################################
Res <- rbind(res_cb_min, res_cox_min, cen.prop)

rownames(Res) <- c("casebase.lambda.min", "cox.lambda.min", "cens.prop")

Res

write.csv(Res, file = glue("setting5{runif(1)}.csv"))
