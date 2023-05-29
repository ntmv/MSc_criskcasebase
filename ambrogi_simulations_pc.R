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
source("../src/fitting_functions.R")

p <- 20
n <- 400

# Very sparse case 
beta <- c(0.5, rep(0,18), 0.5)
beta1 <- -beta

zero_ind1 <- which(beta == 0)
nonzero_ind1 <- which(beta != 0)

# Generate X (iid case)
X <- matrix(rnorm(n*p), nrow = n, ncol = p)

# XB matrix
suma1 <- X %*% beta
suma2 <- X %&% beta1

# Function to generate survival times 
create.times <- function(n, ch, sup.int = 100) { 
  times <- numeric(n) 
  i <- 1 
  while (i <= n) 
  { u <- runif(1) 
  if (ch(0, -log(u)) * ch(sup.int, -log(u)) < 0) 
  { times[i] <- uniroot(ch, c(0, sup.int), tol = 0.0001, y= -log(u))$root 
  i <- i + 1 
  } 
  else { 
    cat("pos")
  }} 
  times
}

# binomial probability of cause 1 
binom.status <- function(ftime, n, a01, a02, size = 1) 
{ prob <- a01(ftime) / (a01(ftime) + a02(ftime))
out <- rbinom(n, size, prob) 
out }


# Cause-specific proportional hazards 
times <- vector()
f.status <- vector()
for (i in seq_len(n)) {
  alpha.1 <- function(t) { ((0.5*t)*exp(suma1[i])) }
  alpha.2 <- function(t) { t*exp(suma2[i]) }
  
  cum.haz <- function(t, y) { stats::integrate(alpha.1, lower=0.001, upper=t, 
                                               subdivisions=1000)$value + 
      stats::integrate(alpha.2, lower=0.001, upper=t, 
                       subdivisions=1000)$value - y } 
  times[i] <- create.times(1, cum.haz)
  f.status[i]<- binom.status(times, 1, alpha.1, alpha.2) + 1
}


# Censoring 
# Set seed to be somwehat consistent between simulations 
set.seed(121)
cens.times <- runif(n, 0, 6)

# Censoring in status variable 
f.status <- as.numeric(times <= cens.times) * f.status

# times with censoring 
times <- pmin(times, cens.times) 

# Dataset 
sim.dat <- data.frame(time = times, status = f.status)

sim.dat <- cbind(sim.dat, X)

# True cumulative incidence 
cif <- cuminc(ftime = sim.dat$time, fstatus = sim.dat$status)
cif

colnames(sim.dat)[3:22] <- paste0("X", seq_len(p))
#################################################################
# Split into training and validation sets (stratified)
train.index <- caret::createDataPartition(sim.dat$status, p = 0.70, list = FALSE)
train <- sim.dat[train.index,]
test <- sim.dat[-train.index,]
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
