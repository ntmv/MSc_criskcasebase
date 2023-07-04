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
library(survminer)
# Fitting functions 
source("../src/fitting_functions_nonparallel.R")


#Results <- replicate(5, {
# Simulate data from cause-specific sub-distribution hazards 
n <-400
p <- 120
num.true <- 20
beta1 <- c(rep(1, 10), rep(-1, 10),rep(0, p-num.true))
beta2 <- c(rep(1, 10), rep(-1, 10), rep(0, p-num.true))


beta1 <- c(1, 0.8, 0.5, 1,1,  1, 1, 1, 0.6, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1,1, rep(0, p-num.true))
beta2 <- c(1, 1, 0.2, 1, 1,  -0.6, -1, -1, -1, -1,  1, 1, 1, 1, 1, 0.8, 0.8, 0.8, 0.8, 0.8, rep(0, p-num.true))

sim.data <- cause_subdist_sim(n = n, p = p, nblocks = 4, 
                                                beta1 = beta1, mix_p = 0.6, 
                                                beta2 = beta2, u.max = 1)

cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)

# Training-test split 
train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
train <- sim.data[train.index,]
test <- sim.data[-train.index,]

cif <- cuminc(fstatus = train$fstatus, ftime = train$ftime)
plot(cif)

png(filename = "~/Desktop/cif.png", width = 15, height = 10, units = "cm", res = 300)
ggcompetingrisks(cif, palette = "Dark2",
                 legend = "top",
                 ggtheme = theme_bw())
dev.off()
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
cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.5, 
                      lambda = cox_mod$lambda.min)

cc_min <- coef(cox_val_min)

res_cox_min <- varsel_perc(cc_min, beta1)
########################## Fit casebase model #################################
#   S
# Fit case-base model 
# Convert to case-base dataset
surv_obj_train <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))

cov_train <- cbind(train[, c(grepl("X", colnames(train)))], time = log(train$ftime))

# Create case-base dataset
cb_data_train <- create_cbDataset(surv_obj_train, cov_train, ratio = 1)

cv.alpha <- mtool.multinom.cv(cb_data_train, alpha = 1, nfold = 5,
                              constant_covariates = 2, seed = 115)

#cv.alpha

# Cross-validation plot 
p1 <- plot_cv.multinom(cv.alpha)


train.index <- caret::createDataPartition(train$fstatus, p = 0.70, list = FALSE)
train_new <- train[train.index,]
validation <- train[-train.index,]

# Create case base dataset for validation
surv_obj_validation <- with(validation, Surv(ftime, as.numeric(fstatus), type = "mstate"))

cov_validation <- as.matrix(cbind(validation[, c(grepl("X", colnames(validation)))],
                                  time = log(validation$ftime)))

# Create case-base dataset
cb_data_validation <- create_cbDataset(surv_obj_validation, cov_validation, ratio = 3)


# Fit lambda grid to training set 
epsilon <- 0.001
lambda_max <- 0.0001
lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
cv_res <- mclapply(lambdagrid, 
                 function(lambda_val) {
                   fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                               lambda = lambda_val, alpha = 0.5)}, mc.cores = 4, mc.set.seed = 112)


# Fit to validation test to choose lambda min
mult_deviance <- unlist(lapply(cv_res, multi_deviance, cb_data = cb_data_validation))

# Lambda min 
lambda.min <- lambdagrid[which.min(mult_deviance)]
cv_se <- sqrt(var(mult_deviance))
dev.1se <- mult_deviance[which.min(mult_deviance)] + cv_se
range.1se <- lambdagrid[which(mult_deviance <= dev.1se)]
lambda.1se <- max(range.1se)
#})

# validation set
surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)

# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = 1e-16 , alpha = 1e-16, unpen_cov = 2)

res_cb_min <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)

###################################################################################
Res <- rbind(res_cb_min, res_cox_min, cen.prop)

rownames(Res) <- c("casebase.lambda.min", "cox.lambda.min", "cens.prop")

Res
toc()

#}, simplify = FALSE)


#Result3 <- do.call(rbind, Results)

write.csv(Res, file = paste0(runif(1), "iid_sparse.csv"))

#png(filename = paste0(runif(1), "cv.png"), height = 15, width = 20, res = 300, units = "cm")
#p1
#dev.off()
#####################################################################################

