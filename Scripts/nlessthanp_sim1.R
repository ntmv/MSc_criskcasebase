#################### Simulating competing risks data ##########################################
##### Packages ######
library(casebase)
library(survsim)
library(tidyverse)

####################### Setting simulation parameters ##############################
beta1 <- c(rep(-1, 4), rep(0, 5), rep(1, 4), rep(0, 16), rep(1, 5), rep(0, 16))
beta2 <- c(rep(0, 4), rep(-1, 5), rep(1, 4), rep(0, 16), rep(1, 5), rep(0, 16))
rho <- 0.6
p <- 50
n <- 250
zero_ind1 <- which(beta == 0)
nonzero_ind1 <- which(beta != 0)
zero_ind2 <- which(beta == 0)
nonzero_ind2 <- which(beta != 0)
alpha = 0.5
################### Generating MVN covariates with an AR1 correlation (N > p) #################
# AR1 correlation function 
covAR1 <- function(p, rho) {
  stopifnot(p >= 2 && rho >= 0)
  # https://statisticaloddsandends.wordpress.com/2020/02/07/generating-correlation-matrix-for-ar1-model/
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  
  return(rho^exponent)
}

Sigma <- covAR1(p = p,
                rho = rho)

X <- rmvnorm(n = n, sigma = Sigma)
colnames(X) <- paste0("X", seq_len(p))
################## Simple experimental generation of competing risks data #############################
set.seed(100)
# Event 1
# Generate failure times according to 
# exponential model
fail_times1 <- rexp(n = n,
                   rate = exp(X %*% beta1))

# Competing event failure times 
fail_times2 <- rexp(n = n,
                    rate = exp(X %*% beta2))

# Binomial experiment for which event 
z <- rbinom(n, 1, 0.5)

dat_times <- cbind(fail_times1, fail_times2, z)

# Generate censoring times 
# U * exp(x^T beta), where U ~ Unif(1, 3)
# which should lead to ~30% censoring
cens_times1 <- runif(n = n, min = 1, max = 3) *
  rexp(n = n, rate = exp(X %*% beta1))
cens_times2 <- runif(n = n, min = 1, max = 3) *
  rexp(n = n, rate = exp(X %*% beta2))

times <- pmin(fail_times1, fail_times2)

dat_times <- as.data.frame(cbind(fail_times1, fail_times2, times, z))

dat_times <- dat_times %>%
  mutate(status = ifelse(fail_times2 > fail_times1, "1", "2"))

status_cens <- as.numeric(fail_times1 < cens_times1)
status_cens <- which(status_cens == 0)

#Induce censoring 
dat_times$status[status_cens] <- "0"

dat_times <- dat_times %>%
  select(-fail_times1, -fail_times2)

dat_times <- cbind(dat_times, X)

# Split into train and test sets
spec = c(train = .6, validation = .2, test = .2)

cut_dat = sample(cut(
  seq(nrow(dat_times)), 
  nrow(dat_times)*cumsum(c(0,spec)),
  labels = names(spec)
))

res = split(dat_times, cut_dat)
######################################################################
# training dataset 

train <- res$train

# Create survival object
surv_obj_train <- with(train, Surv(times, as.numeric(status), type = "mstate")) 

# Create covariate matrix
cov_train<- cbind(train[4:53], train$z)

# Create case-base sampled dataset
cb_data_train <- create_cbDataset(surv_obj_train, as.matrix(cov_train))

# testing dataset

test <- res$test

# Create survival object
surv_obj_test <- with(test, Surv(times, as.numeric(status), type = "mstate")) 

# Create covariate matrix
cov_test<- cbind(test[4:53], test$z)

# Create case-base sampled dataset
cb_data_test <- create_cbDataset(surv_obj_test, as.matrix(cov_test))


# validation dataset


validation <- res$validation

# Create survival object
surv_obj_val <- with(validation, Surv(times, as.numeric(status), type = "mstate")) 

# Create covariate matrix
cov_val <- cbind(validation[4:53], validation$z)

# Create case-base sampled dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))



# Fit casebase penalized model for different values of lambda and calculate multinomial deviance
fit_object <- list()
lambda_vals <- seq(0, 1, 0.01)
for (i in 1:length(lambda_vals)) {
  fit_object[[i]] <- fit_cbmodel(cb_data_train, lambda = lambda_vals[i], alpha = alpha)
}
# Compute deviance
dev <- c()
for (i in 1:101) {
  dev[i] <- multi_deviance(cb_data_test, fit_object[[i]])
}

lambda_vals[which.min(dev)]

