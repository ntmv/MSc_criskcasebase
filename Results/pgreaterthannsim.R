#!/usr/bin/env Rscript

#################### Simulating competing risks data ##########################################
##### Packages ######
library(casebase)
library(survsim)
library(tidyverse)
library(mvtnorm)
library(tidymodels)
library(mtool)
library(future.apply)
library(fastcmprsk)
library(timereg)
library(CoxBoost)
library(devtools)
library(CoxBoost)
library(timereg)
library(riskRegression)
library(penalized)
library(glmnet)
####################### Setting simulation parameters ##############################
source("~/Desktop/src/fitting_functions.R")

results <- replicate(3, {
beta1 <- c(rep(-1.5, 5), rep(1, 30), rep(1.5, 5), rep(-1.5, 5), rep(0, 30), rep(1.5, 5), rep(1.5, 20))
beta2 <- c(rep(-1.5, 5), rep(0, 10), rep(1.5, 5), rep(-1, 5), rep(1.5, 5), rep(1.5, 5), rep(-1.5, 5), rep(1, 30),
           rep(1, 10), rep(1.5, 5),
           rep(1, 15))
rho <- 0.6
p <- 100
n <- 80
zero_ind1 <- which(beta1 == 0)
nonzero_ind1 <- which(beta1 != 0)
zero_ind2 <- which(beta2 == 0)
nonzero_ind2 <- which(beta2 != 0)
alpha = 0.5

################################################################################
# Generating competing risks data (without outliers)
sim.dat <- crisk.sim_mvn(n = n, p = p, rho = rho, foltime = 100, dist.ev = c("weibull", "weibull"), 
                         anc.ev = c(1.9, 1.2), beta0.ev = c(1, 1), dist.cens = "unif", anc.cens= 4,
                         beta0.cens= 1, nsit = 2, beta=list(c(beta1),c(beta2)))

# Check if you can recover parameters with survsim 

# Fix column names 
colnames(sim.dat)[8:107] <- paste0("X", seq_len(p))

# Remove columns we don't need
sim.dat <- sim.dat %>%
  select(-start, -stop, -z)

# Fix cause column for censored observations 
sim.dat <- sim.dat %>%
  mutate(cause = ifelse(is.na(cause), 0, cause))

# Remove status variable
sim.dat$status <- NULL
###############################################################################
# Split into training and validation sets (stratified)
train.index <- caret::createDataPartition(sim.dat$cause, p = 0.50, list = FALSE)
train <- sim.dat[train.index,]
validation <- sim.dat[-train.index,]
############################################################################
# Training dataset 

# Create survival object
surv_obj_train <- with(train, Surv(time, as.numeric(cause), type = "mstate")) 

# Create covariate matrix
cov_train <- cbind(train[4:103])

# Create case-base sampled dataset
cb_data_train <- create_cbDataset(surv_obj_train, as.matrix(cov_train))


# Cross-validation for lambda 
lambda.val <- mtool.multinom.cv(cb_data_train, regularization = 'elastic-net', lambdagrid = seq(0.003, 0.007, 0.0001), 
                                alpha = 0.5, nfold = 10)

# Create survival object
surv_obj_val <- with(validation, Surv(time, as.numeric(cause), type = "mstate")) 

# Create covariate matrix
cov_val <- cbind(validation[4:103])

# Create case-base sampled dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))


# Apply to validation set
fit_val <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                        lambda = lambda.val$lambda.min, alpha = 0.5)


# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(fit_val$coefficients[2:100, 1] != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(fit_val$coefficients[2:101, 1] == 0), zero_ind1))/length(zero_ind1)

1-num_zeroes_1
1-wrong_zeroes1

res_cb <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)

# Cox boost model
# Penalty parameter
optim.res <- optimCoxBoostPenalty(time = train$time, status = train$cause, x = as.matrix(train[, c(4:103)]),
                                  trace = TRUE, start.penalty = 500)
cbfit <- iCoxBoost(Hist(time,cause) ~ ., data= validation, stepno = optim.res$cv.res$optimal.step)

coef(cbfit)
summary(cbfit)


# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(coef(cbfit) != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(coef(cbfit) == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_bfg <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)

# Fitting penalized binomial model
times1 <- sort(unique(sim.dat$time[sim.dat$cause == 1]))
times2 <- quantile(times1, probs=seq(0.01,1,length.out=10))

# Generate stacked data for competing risks
stacked_dat <- prep.glm.comprisk(train,time="time",cause="cause",
                                 times=times2,censmod=0,cens.code=0)

p <- length(times2)
nn <- nrow(stacked_dat)/p
fmla <- reformulate(paste0("X",  seq_len(100)))
mm <- model.matrix(fmla, data = stacked_dat)
mm <- mm[,-1] ### removes intercept
pp=100
timep <- ncol(mm)-pp ### number of unpenalized variables

penalty <- c(rep(1,pp))

# Adaptive LASSO weights for the unpenalized variables
ridge.cv <- cv.glmnet(mm, stacked_dat$Nt,weights= stacked_dat$weights,
                      standardize=TRUE,
                      type.measure="deviance",
                      penalty.factor=penalty,
                      family="binomial",alpha=0)

bhat<-as.matrix(coef(ridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
if(all(bhat==0)){
  ## if bhat is all zero then assign very close to zero weight to all.
  ## Amounts to penalizing all of the second stage to zero.
  bhat<-rep(.Machine$double.eps*2,length(bhat))
}
adpen<-(penalty/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight

cvudglmnet_train <- cv.glmnet(mm,stacked_dat$Nt,weights= stacked_dat$weights,
                        ###        intercept=FALSE,
                        standardize=TRUE,
                        type.measure="deviance",
                        penalty.factor=adpen,
                        family="binomial",alpha=1)

# Fit to validation set 
# Generate stacked data for competing risks
stacked_dat <- prep.glm.comprisk(validation,time="time",cause="cause",
                                 times=times2,censmod=0,cens.code=0)

p <- length(times2)
nn <- nrow(stacked_dat)/p
fmla <- reformulate(paste0("X",  seq_len(20)))
mm <- model.matrix(fmla, data = stacked_dat)
mm <- mm[,-1] ### removes intercept
pp=100
timep <- ncol(mm)-pp ### number of unpenalized variables

penalty <- c(rep(1,pp))

# Adaptive LASSO weights for the unpenalized variables
ridge.cv <- cv.glmnet(mm, stacked_dat$Nt,weights= stacked_dat$weights,
                      standardize=TRUE,
                      type.measure="deviance",
                      penalty.factor=penalty,
                      family="binomial",alpha=0)

bhat<-as.matrix(coef(ridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
if(all(bhat==0)){
  ## if bhat is all zero then assign very close to zero weight to all.
  ## Amounts to penalizing all of the second stage to zero.
  bhat<-rep(.Machine$double.eps*2,length(bhat))
}
adpen<-(penalty/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight

cvudglmnet <- glmnet(mm,stacked_dat$Nt,weights= stacked_dat$weights,
                        ###        intercept=FALSE,
                        standardize=TRUE,
                        type.measure="deviance",
                        penalty.factor=adpen,
                        family="binomial",alpha=1, lambda = cvudglmnet_train$lambda.min)

cc  <- as.numeric(coef(cvudglmnet, s="lambda.min"))

# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(cc != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(cc == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_bm <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)
Res <- rbind(res_cb, res_bfg, res_bm)

}, simplify = FALSE) 

res <- do.call(rbind.data.frame, results)

write.csv(res, file = "~/Desktop/results.csv")
