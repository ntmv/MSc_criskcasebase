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
library(survsim)

# Fitting functions 
source("../src/fitting_functions_nonparallel.R")

n <- 400
p <- 20
beta <- list(c(0.5, rep(0, 8), 0.5, rep(0, 10)), c(0.2, rep(0, 8), 0.2, rep(0, 10)))
dist.ev <- c("weibull", "weibull")
anc.ev <- c(0.8, 0.3)
beta0.ev <- c(0.1, 0.1)

# Generating survival data 
sim.data <- crisk.sim_block(n =n, p = p, nblocks = 4, foltime = 4, dist.ev = dist.ev, 
                          anc.ev = anc.ev, beta0.ev = beta0.ev, beta0.cens = 0.05, anc.cens = 4, nsit = 2, 
                          beta = beta)

# fix status variable
sim.data$cause <- with(sim.data, ifelse(is.na(sim.data$cause), 0, sim.data$cause))
colnames(sim.data)[grepl("x", colnames(sim.data))]   <- paste0("X", seq_len(p))

# Format data
sim.data <- sim.data %>%
  dplyr::select(-nid, -status, -start, -stop, -z) %>%
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

# 0.5 se values 
cvse <- sqrt(var(cox_mod$cvm))
dev.0.5se <- min(cox_mod$cvm) + cvse/2
dev.1se <- min(cox_mod$cvm) + cvse
range.1se <- cox_mod$lambda[which((cox_mod$cvm <= dev.1se))]
range.0.5se <- cox_mod$lambda[which((cox_mod$cvm <= dev.0.5se))]
lambda.0.5se <- max(range.0.5se)
lambda.min0.5se <- min(range.0.5se)
lambda.min1se <- min(range.1se)

# Fit on validation set 
cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 1, 
                  lambda = cox_mod$lambda.min)

cox_val_1se <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 1, 
                      lambda = cox_mod$lambda.1se)

cox_val_0.5se <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 1, 
                      lambda = lambda.0.5se)

cox_val_min1se <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 1, 
                      lambda = lambda.min1se)

cox_val_min0.5se <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 1, 
                         lambda = lambda.min0.5se)

cc_min <- coef(cox_val_min)

res_cox_min <- varsel_perc(cc_min, beta[[1]])

cc_1se <- coef(cox_val_1se)

res_cox_1se <- varsel_perc(cc_1se, beta[[1]])

cc_0.5se <- coef(cox_val_0.5se)

res_cox_0.5se <- varsel_perc(cc_0.5se, beta[[1]])

cc_min1se <- coef(cox_val_min1se)

res_cox_min1se <- varsel_perc(cc_min1se, beta[[1]])

cc_min0.5se <- coef(cox_val_min0.5se)

res_cox_min0.5se <- varsel_perc(cc_min0.5se, beta[[1]])


################ Comparing case-base fit ###########################
# Convert to case-base dataset
surv_obj_train <- with(train, Surv(time, as.numeric(status), type = "mstate"))

cov_train <- as.matrix(cbind(train[, c(grepl("X", colnames(train)))], time = log(train$time)))

# Create case-base dataset
cb_data_train <- create_cbDataset(surv_obj_train, cov_train, ratio = 10)

tic()
cv.alpha <- mtool.multinom.cv(cb_data_train, alpha = 1, nfold = 5)
toc()

cv.alpha

# Cross-validation plot 
p1 <- plot_cv.multinom(cv.alpha$deviance_grid, cv.alpha$lambdagrid, cv.alpha$lambda.min, cv.alpha$lambda.1se, nfold = 10)

# validation set
surv_obj_val <- with(test, Surv(time, as.numeric(status), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$time))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))

# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = cv.alpha$lambda.min , alpha = 1)

# Lambda min 1se
fit_val_min1se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                              lambda = cv.alpha$lambda.min1se, alpha = 1)

# Lambda min 0.5se
fit_val_min0.5se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                                lambda = cv.alpha$lambda.min0.5se , alpha = 1)

# Lambda 1 se
fit_val_1se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = cv.alpha$lambda.1se, alpha = 1)

# Lambda 0.5 se
fit_val_0.5se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                             lambda = cv.alpha$lambda.0.5se , alpha = 1)

res_cb_min <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta[[1]])

res_cb_min1se <-  varsel_perc(fit_val_min1se$coefficients[1:eval(parse(text="p")), 1], beta[[1]])

res_cb_min0.5se <-  varsel_perc(fit_val_min0.5se$coefficients[1:eval(parse(text="p")), 1], beta[[1]])

res_cb_0.5se <- varsel_perc(fit_val_0.5se$coefficients[1:eval(parse(text="p")), 1], beta[[1]])

res_cb_1se <- varsel_perc(fit_val_1se$coefficients[1:eval(parse(text="p")), 1], beta[[1]])


##########################################################################
Res <- rbind(res_cb_min, res_cb_1se, res_cb_0.5se, res_cb_min0.5se, res_cb_min1se, res_cox, 
             res_cox_1se, res_cox_0.5se, res_cox_min1se, res_cox_min0.5se)

rownames(Res) <- c("casebase.lambda.min", "casebase.lambda.1se", "casebase.lambda.0.5se",
                   "casebase.lambda.min0.5se", "casebase.lambda.min1se", "cox.lambda.min", 
                   "cox.lambda.1se",  "cox.lambda.0.5se", "cox.lambda.min1se",
                   "cox.lambda.min0.5se")

write.csv(Res, file = paste0(runif(1), "block_sparse.csv"))

png(filename = paste0(runif(1), "cv.png"), height = 15, width = 20, res = 300, units = "cm")
p1
dev.off()
