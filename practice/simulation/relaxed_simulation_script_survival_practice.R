########################## Practice survival simulation for relaxed LASSO #######################
# Libraries 
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

source("src/practice_relaxed_implementations.R")
source("src/fitting_functions.R")

# Setup
n <- 400
p <- 20
seed = 2023
beta1 <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, rep(0, p/2))
beta2 <-  c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, rep(0, p/2))

# Simulate data
sim.data <- cause_hazards_sim(n = n, p = p,  
                              beta1 = beta1, beta2 = beta2, rate_cens = 0.15, 
                              h1 = 0.55, h2 = 0.15, gamma1 = 1.5, gamma2 = 1.5)

# Training-test split 
train.index <- caret::createDataPartition(as.factor(sim.data$fstatus), p = 0.75, list = FALSE)
train <- sim.data[train.index,]
test <- sim.data[-train.index,]

test_list = list("time" = test$ftime,
                 "event_ind" = test$fstatus,
                 "covariates" = test[, paste("X", as.character(seq(1:p)), sep = "")],
                 "offset" = test$offset)  
  
  
############################ Coxnet model cause 1 #############################
# Censor competing event
y_train1 <- Surv(time = train$ftime, event = train$fstatus == 1)

x_train1 <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 

# Censor competing event
y_test1 <- Surv(time = test$ftime, event = test$fstatus == 1)

x_test1 <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 

set.seed(seed)
# Fit cause-specific cox model with glmnet on training set 
cox_mod1 <- cv.glmnet(x = x_train1, y = y_train1, family = "cox", alpha = 1)

# Fit on validation set 
cox_val_min1 <- glmnet(x = x_test1, y = y_test1, family = "cox", alpha = 1, 
                       lambda = cox_mod1$lambda.min)

cc_min1 <- coef(cox_val_min1)

res_cox_min1 <- varsel_perc(cc_min1, beta1)

############################ Coxnet model cause 2 #############################
# Censor competing event
y_train2 <- Surv(time = train$ftime, event = train$fstatus == 2)

x_train2 <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 

# Censor competing event
y_test2 <- Surv(time = test$ftime, event = test$fstatus == 2)

x_test2 <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 

set.seed(seed)
# Fit cause-specific cox model with glmnet on training set 
cox_mod2 <- cv.glmnet(x = x_train2, y = y_train2, family = "cox", alpha = 1)

# Fit on validation set 
cox_val_min2 <- glmnet(x = x_test2, y = y_test2, family = "cox", alpha = 1, 
                       lambda = cox_mod1$lambda.min)

cc_min2 <- coef(cox_val_min2)

res_cox_min2 <- varsel_perc(cc_min2, beta2)

glmnet_cox_coefs = list("coefficients" = as.matrix(cbind(cc_min1, cc_min2)))


# Test multi-deviance
dev_cox <- multi_deviance_cox(cb_data = test_list, fit_object = glmnet_cox_coefs)

######################### casebase model (without relaxed fit) ####################
# Train case-base model 
cv.lambda <- mtool.multinom.cv(train, seed = seed, alpha = 1, nfold = 5, train_ratio = 20, lambda_max = 0.05)

# Test set 
surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)

# Case-base fits 
# Lambda.min
fit_val_min_cv_cb <- fit_cbmodel(cb_data_val, regularization = 'l1',
                           lambda = cv.lambda$lambda.min , alpha = 1, unpen_cov = 2)

coef_val1 <- fit_val_min_cv_cb$coefficients[1:eval(parse(text="p")), 1]
coef_val2 <- fit_val_min_cv_cb$coefficients[1:eval(parse(text="p")), 2]

res_cb_min1_cv <- varsel_perc(coef_val1, beta1)
res_cb_min2_cv <- varsel_perc(coef_val2, beta2)


# Test multi-deviance
dev_cb_mtool_cv <- multi_deviance(cb_data = cb_data_val, fit_object = fit_val_min_cv_cb)



######################### casebase model (without relaxed fit) ####################
# Train case-base model 
relaxed.lambda <- multinom.relaxed_lasso_cb(train, seed = seed, alpha = 1, gamma = 0,
                                       nfold = 5, train_ratio = 20, lambda_max = 0.05)


# Case-base relaxed fit 
# Lambda.min
fit_val_min_relaxed_cb <- fit_cbmodel(cb_data_val, regularization = 'l1',
                              lambda = relaxed.lambda$lambda.min , alpha = 1, unpen_cov = 2)

coef_val1 <- fit_val_min_relaxed_cb$coefficients[1:eval(parse(text="p")), 1]
coef_val2 <- fit_val_min_relaxed_cb$coefficients[1:eval(parse(text="p")), 2]

res_cb_min1_relaxed <- varsel_perc(coef_val1, beta1)
res_cb_min2_relaxed <- varsel_perc(coef_val2, beta2)


# Test multi-deviance
dev_cb_mtool_relaxed <- multi_deviance(cb_data = cb_data_val, fit_object = fit_val_min_cv_cb)




#############################################################
Res <- rbind(res_cb_min1, res_cox_min1)

rownames(Res) <- c("casebase.lambda.min_cause1", "cox.lambda.min_cause1")

write.csv(Res, file = paste("practice/simulation/results/", as.character(runif(1)), "_survival_devs_n=", as.character(n), "_p=",
      as.character(p), ".csv", sep = ""))

