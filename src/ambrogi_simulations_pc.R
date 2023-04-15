################ Replicating Ambrogi et al. 2016 simulation ##############
# Libraries 
library(casebase)
library(future.apply)
library(glmnet)
library(mtool)
library(timereg)
library(parallel)
library(foreach)
library(doParallel)

# Fitting functions 
source("../src/fitting_functions.R")

# Cluster set up cores
cl <- parallel::makeCluster(2, timeout = 60)
plan(cluster, workers = cl)

# Set a different seed on each member of the cluster (just in case)
results <-  future_replicate(1000, {
# independent normal variables 
p <- 20
n <- 400
rho <- 0.5

Sigma <- covAR1(p = p, rho = rho)
X <- mvtnorm::rmvnorm(n, sigma = Sigma)
# Very sparse case
# Common betas for both competing risks
beta <- c(0.5, rep(0, 18), 0.5)
zero_ind1 <- which(beta == 0)
nonzero_ind1 <- which(beta != 0)


# XB matrix
suma <- X %*% beta

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
alpha.1 <- function(t) { ((0.5*t)*exp(suma[i])) }
alpha.2 <- function(t) { t*exp(suma[i]) }

cum.haz <- function(t, y) { stats::integrate(alpha.1, lower=0.001, upper=t, 
                                             subdivisions=1000)$value + 
 stats::integrate(alpha.2, lower=0.001, upper=t, 
                 subdivisions=1000)$value - y } 
times[i] <- create.times(1, cum.haz)
f.status[i]<- binom.status(times, 1, alpha.1, alpha.2) + 1
}

# Censoring 
cens.times <- runif(n, 0, 6)

# Censoring in status variable 
f.status <- as.numeric(times <= cens.times) * f.status

prop.table(table(f.status))

# times with censoring 
times <- pmin(times, cens.times) 

# Dataset 
sim.dat <- data.frame(time = times, status = f.status)

sim.dat <- cbind(sim.dat, X)

colnames(sim.dat)[3:22] <- paste0("X", seq_len(p))

# Let's try fitting binomial regression to this dataset
# Split into training and validation sets (stratified)
train.index <- caret::createDataPartition(sim.dat$status, p = 0.70, list = FALSE)
train <- sim.dat[train.index,]
validation <- sim.dat[-train.index,]


# Fitting penalized binomial model
times1 <- sort(unique(train$time[train$status == 1]))
times2 <- quantile(times1, probs=seq(0.01,1,length.out=10))

# Generate stacked data for competing risks
stacked_dat <- prep.glm.comprisk(train,time="time",cause= "status",
                                 times= times2,censmod=0,cens.code=0)

stacked_dat <- stacked_dat[order(stacked_dat$id),]

p <- length(times2)
nn <- nrow(stacked_dat)/p

fmla <- fmla <- as.formula(paste(" ~ factor(h) +", paste(colnames(train)[3:22], collapse= "+")))
mm <- model.matrix(fmla, data = stacked_dat)
mm <- mm[,-1] ### removes intercept
pp=20
timep <- ncol(mm)-pp ### number of unpenalized variables

penalty <- c(rep(0, timep), rep(1,pp))

# Keep id's together
uno  <- rep(1:9, each=round(nn/10)*p)
qui <- rep(10, each=(nrow(mm)-length(uno)))
foldID = c(uno, qui)


# Adaptive LASSO weights for the unpenalized variables
ridge.cv <- cv.glmnet(mm, stacked_dat$Nt, weights= stacked_dat$weights,
                      standardize=TRUE,
                      type.measure="deviance",
                      penalty.factor=penalty,
                      family="binomial",alpha= 0, foldid=foldID, parallel = TRUE)

bhat<-as.matrix(coef(ridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
if(all(bhat==0)){
  ## if bhat is all zero then assign very close to zero weight to all.
  ## Amounts to penalizing all of the second stage to zero.
  bhat<-rep(.Machine$double.eps*2,length(bhat))
}
adpen<-(penalty/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight

cvudglmnet_train <- cv.glmnet(mm,stacked_dat$Nt, weights= stacked_dat$weights,
                              intercept=FALSE,
                              standardize=TRUE,
                              type.measure="deviance",
                              penalty.factor=adpen,
                              family="binomial",alpha=1, parallel = TRUE)


# Fitting penalized binomial model on validation dataset 
times1 <- sort(unique(validation$time[validation$status == 1]))
times2 <- quantile(times1, probs=seq(0.01,1,length.out=10))

# Generate stacked data for competing risks
stacked_dat <- prep.glm.comprisk(validation,time="time",cause= "status",
                                 times= times2,censmod=0,cens.code=0)

stacked_dat <- stacked_dat[order(stacked_dat$id),]

p <- length(times2)
nn <- nrow(stacked_dat)/p

fmla <- fmla <- as.formula(paste(" ~ factor(h) +", paste(colnames(validation)[3:22], collapse= "+")))
mm <- model.matrix(fmla, data = stacked_dat)
mm <- mm[,-1] ### removes intercept
pp=20
timep <- ncol(mm)-pp ### number of unpenalized variables

penalty <- c(rep(0, timep), rep(1,pp))

# Validation fit 
val.fit <- glmnet(mm,stacked_dat$Nt, weights= stacked_dat$weights,
                              intercept=FALSE,
                              standardize=TRUE,
                              type.measure="deviance",
                              family="binomial",alpha=1, lambda = cvudglmnet_train$lambda.1se)

cc <- coef(val.fit)[11:30]

# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(cc != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(cc == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_bm <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)

################ Comparing case-base fit ###########################
surv_obj_train <- with(train, Surv(time, as.numeric(status), type = "mstate"))

cov_train <- cbind(train[3:22])

# Create case-base dataset
cb_data_train <- create_cbDataset(surv_obj_train, as.matrix(cov_train))

# Cross-validation for lambda
lambda.val_0.5 <- mtool.multinom.cv(cb_data_train, regularization = 'elastic-net', lambdagrid = c(0.001, 0.1),
                         alpha = 1, nfold = 10)
# validation set
surv_obj_val <- with(validation, Surv(time, as.numeric(status), type = "mstate"))

# Covariance matrix
cov_val <- cbind(validation[3:22])

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))

# Apply to validation set
fit_val_0.5 <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                                          lambda = 0.009, alpha = 1)


# Check number of right zeroes (sensitivity)
num_zeroes_1 <- length(intersect(which(fit_val_0.5$coefficients[1:20, 1] != 0), nonzero_ind1))/length(nonzero_ind1)

# Check number of wrong zeroes (specificity)
wrong_zeroes1 <- length(intersect(which(fit_val_0.5$coefficients[1:20, 1] == 0), zero_ind1))/length(zero_ind1)

1-wrong_zeroes1
1-num_zeroes_1

res_cb <- c(wrong_zeroes1, num_zeroes_1, 1-wrong_zeroes1, 1-num_zeroes_1)
##########################################################################

Res <- rbind(res_bm, res_cb)

}, simplify = FALSE)

Res <- do.call(rbind.data.frame, results)

## CLEANUP
parallel::stopCluster(cl)

write.csv(../results/"iidsparseresults_ambrogi.csv")
