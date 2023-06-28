########################### P > n simulation ######################
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

# Fitting functions 
source("../src/fitting_functions_nonparallel.R")


Results <- replicate(10, {
# Simulate data from cause-specific subdistribution hazards 
n <- 400
p <- 1000
num_true <- 20
beta1 <- c(rep(0, p))
beta2 <- c(rep(0, p))
nu_ind <- seq(num_true)
beta1[nu_ind] <- 1
beta2[nu_ind] <- 0.8

sim.data <- simulateTwoCauseModel(n = n, p = p, nblocks = 4, 
                                                beta1 = beta1, mix_p = 0.6, 
                                                beta2 = beta2, u.max = 1)

cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)

# Training-test split 
train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
train <- sim.data[train.index,]
test <- sim.data[-train.index,]

cif <- cuminc(fstatus = train$fstatus, ftime = train$ftime)
plot(cif)
######################## Fit cox-regression model ###############################
# Censor competing event
y_train <- Surv(time = train$ftime, event = train$fstatus == 1)

x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 

# Censor competing event
y_test <- Surv(time = test$ftime, event = test$fstatus == 1)

x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 

# Fit cause-specific cox model with glmnet on training set 
cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 1)

# Fit on validation set 
cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 1, 
                      lambda = cox_mod$lambda.min)

cc_min <- coef(cox_val_min)

res_cox_min <- varsel_perc(cc_min, beta1)
########################## Fit casebase model #################################
# Fit case-base model 
# Convert to case-base dataset
surv_obj_train <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))

cov_train <- cbind(train[, c(grepl("X", colnames(train)))], time =  log(train$ftime))

# Create case-base dataset
cb_data_train <- create_cbDataset(surv_obj_train, cov_train, ratio = 10)

cv.alpha <- mtool.multinom.cv(cb_data_train, alpha = 1, nfold = 5, lambda_max = 0.1, constant_covariates = 2)

cv.alpha

# Cross-validation plot 
p1 <- plot_cv.multinom(cv.alpha$deviance_grid, cv.alpha$lambdagrid, cv.alpha$lambda.min, cv.alpha$lambda.1se, nfold = 5)

# validation set
surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = bs(test$ftime))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))

# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = cv.alpha$lambda.min , alpha = 1, unpen_cov = 4)

res_cb_min <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)

###################################################################################
Res <- rbind(res_cb_min, res_cox_min, cen.prop)

rownames(Res) <- c("casebase.lambda.min", "cox.lambda.min", "cens.prop")

Res

}, simplify = FALSE)


Result3 <- do.call(rbind, Results)

write.csv(Res, file = paste0(runif(1), "iid_sparse.csv"))

png(filename = paste0(runif(1), "cv.png"), height = 15, width = 20, res = 300, units = "cm")
p1
dev.off()

