################ Variable selection simulation  ##############
library(casebase)
library(future.apply)
library(glmnet)
library(mtool)
library(timereg)
library(parallel)
library(tictoc)
library(tidyverse)
library(riskRegression)
library(cmprsk)

# Fitting functions 
source("../src/fitting_functions_nonparallel.R")

n <- 400
p <- 20
beta <- list(c(0.5, rep(0, 18), 0.5),c(0.2, rep(0, 18), 0.2))
dist.ev <- c("weibull", "weibull")
anc.ev <- c(0.8, 0.3)
beta0.ev <- c(0.1, 0.1)


# Generating survival data 
sim.data <- crisk.sim_mvn(n = n, p = p, rho = 0.5, foltime = 4, dist.ev = dist.ev, 
                      anc.ev = anc.ev, beta0.ev = beta0.ev, beta0.cens = 0.05, anc.cens = 4, nsit = 2, 
                      beta = beta)

# fix status variable
sim.data$cause <- with(sim.data, ifelse(is.na(sim.data$cause), 0, sim.data$cause))
colnames(sim.data)[grepl("x", colnames(sim.data))]   <- paste0("X", seq_len(p))

# Format data
sim.data <- sim.data %>%
  select(-nid, -status, -start, -stop, -z) %>%
  rename(status = cause)

# Average estimates of incidence and censoring rate
prop.table(table(sim.data$status))

# True cumulative incidence 
cif <- cuminc(ftime = sim.data$time, fstatus = sim.data$status)
cif

#################################################################
# Split into training and validation sets (stratified)
train.index <- caret::createDataPartition(sim.data$status, p = 0.80, list = FALSE)
train <- sim.data[train.index,]
test <- sim.data[-train.index,]
######################### Cause-specific proportional hazards model ###############
# Censor competing event
y_train <- Surv(time = train$time, event = train$status == 1)

x_train <- model.matrix(~ . -time -status, data = train)[, -1] 

# Censor competing event
y_test <- Surv(time = test$time, event = test$status == 1)

x_test <- model.matrix(~ . -time -status, data = test)[, -1] 

# Fit cause-specific cox model with glmnet on training set 
cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 1)

# Fit on validation set 
cox_val <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 1, 
                  lambda = cox_mod$lambda.min)

cc <- coef(cox_val)[1:20]

# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(cc != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(cc == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_cox <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)

################ Comparing case-base fit ###########################
# Convert to case-base dataset
surv_obj_train <- with(train, Surv(time, as.numeric(status), type = "mstate"))

cov_train <- as.matrix(cbind(train[, c(3:22)], time = log(train$time)))

# Create case-base dataset
cb_data_train <- create_cbDataset(surv_obj_train, cov_train, ratio = 10)

tic()
res <- replicate(5,{
cv.alpha <- mtool.multinom.cv(cb_data_train, lambda_max = 0.3, alpha = 1, nfold = 10)
})
toc()

# Check one run output 
cv.alpha <- mtool.multinom.cv(cb_data_train, lambda_max = 0.3, alpha = 1, nfold = 10)

cv.alpha

# Cross-validation plot 
p1 <- plot_cv.multinom(cv.alpha$deviance_grid, cv.alpha$lambdagrid, cv.alpha$lambda.min, cv.alpha$lambda.1se, nfold = 10)

res

# Average over lambda 
res1 <- data.frame(res[c(1, 3, 4, 5, 6),])

res1 <- data.frame(lapply(res1, as.numeric))

res1 <- rowMeans(res1)

res1[1]; res1[2]; res1[3]; res1[4]; res1[5]


# validation set
surv_obj_val <- with(test, Surv(time, as.numeric(status), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(3:22)], time = log(test$time))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))


# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                           lambda = res1[1] , alpha = 1)

# Lambda min 1se
fit_val_min1se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = res1[2], alpha = 1)

# Lambda min 0.5se
fit_val_min0.5se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                             lambda = res1[3] , alpha = 1)

# Lambda 1 se
fit_val_1se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                         lambda = res1[4], alpha = 1)

# Lambda 0.5 se
fit_val_0.5se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                             lambda = res1[5] , alpha = 1)

# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(fit_val_min$coefficients[1:20, 1] != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(fit_val_min$coefficients[1:20, 1] == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_cb_min <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)

# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(fit_val_min1se$coefficients[1:20, 1] != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(fit_val_min1se$coefficients[1:20, 1] == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_cb_min1se <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)


# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(fit_val_min0.5se$coefficients[1:20, 1] != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(fit_val_min0.5se$coefficients[1:20, 1] == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_cb_min0.5se <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)


# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(fit_val_1se$coefficients[1:20, 1] != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(fit_val_1se$coefficients[1:20, 1] == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_cb_0.5se <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)

# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(fit_val_0.5se$coefficients[1:20, 1] != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(fit_val_0.5se$coefficients[1:20, 1] == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_cb_1se <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)


##########################################################################
Res <- rbind(res_cb_min, res_cb_1se, res_cb_0.5se, res_cb_min0.5se, res_cb_min1se, res_cox)

write.csv(Res, file = paste0(runif(1), "iid_sparse.csv"))

png(filename = paste0(runif(1), "cv.png"), height = 15, width = 20, res = 300, units = "cm")
p1
dev.off()
