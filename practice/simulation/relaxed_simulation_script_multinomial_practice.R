########################## Practice multinomial simulation for relaxed LASSO #######################
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
n = 400; p = 20; N = 1; nfolds = 5

tryCatch({
# For each simulation you want to output the deviance + coefficient table
  Results <- replicate(N, {
    
    # Simulate data
    X <- matrix(rnorm(n * p), n, p)
    
    # coefficients for each choice
    beta0 <- rep(0, p)
    beta1 <- c(rep(3, p/2), rep(0, p/2))
    zero_beta1 <- which(beta1 == 0)
    
    beta2 <- c(rep(3, p/2), rep(0, p/2))
    zero_beta2 <- which(beta2 == 0)
    
    
    # vector of probabilities
    vProb = cbind(exp(X%*%beta0), exp(X%*%beta1), exp(X%*%beta2))
    
    # multinomial draws
    mChoices <- t(apply(vProb, 1, rmultinom, n = 1, size = 1))
    dfM <- cbind.data.frame(y = apply(mChoices, 1, function(x) which(x == 1)), X)
    # Rename covariates 
    colnames(dfM)[2:(p+1)] <- paste0('x', 1:p)
    
    train_index_y = caret::createDataPartition(dfM$y, p = 0.80, list = FALSE)
    
    dfM_train = dfM[train_index_y, ]
    dfM_test = dfM[-train_index_y, ]
    
    x_train = as.matrix(dfM_train[, 2:(p+1)])
    y_train = dfM_train$y - 1
    
    x_test = as.matrix(dfM_test[, 2:(p+1)])
    y_test = dfM_test$y - 1
    
    dfM_train_list = list("covariates" = as.matrix(cbind(dfM_train[, c(2:ncol(dfM_train))])),
                         "y" = dfM_train$y - 1)
    
    dfM_test_list = list("covariates" = as.matrix(cbind(dfM_test[, c(2:ncol(dfM_test))])),
                         "y" = dfM_test$y - 1)
    
    ############################ glmnet model #############################
    fit.glmnet.relaxed <- glmnet::cv.glmnet(
      x = x_train, y = y_train,
      family = "multinomial",
      intercept = TRUE,
      type.multinomial = "grouped",  # same sparsity pattern for all outcome classes
      alpha = 1, nfolds = 10,
      relax = TRUE, gamma = 0)
    
    glmnet_lambda_max = fit.glmnet.relaxed$lambda[1]
    glmnet_lambda_min = fit.glmnet.relaxed$relaxed$lambda.min
    
    fit_val_min_glmnet <- glmnet(x = x_test, y = y_test,
                                 lambda = glmnet_lambda_min,
                                 family = "multinomial",
                                 intercept = TRUE,
                                 type.multinomial = "grouped")
    
    
    glmnet_min1 <- coef(fit_val_min_glmnet)$`1` - coef(fit_val_min_glmnet)$`0`
    glmnet_min2 <- coef(fit_val_min_glmnet)$`2` - coef(fit_val_min_glmnet)$`0`
    
    res_glmnet_min1 <- varsel_perc(glmnet_min1[-1], beta1)
    res_glmnet_min2 <- varsel_perc(glmnet_min2[-1], beta2)
    
    reformat_glmnet_min1 = c(glmnet_min1[-1], glmnet_min1[1])
    reformat_glmnet_min2 = c(glmnet_min2[-1], glmnet_min2[1])
    glmnet_coefs = list("coefficients" = as.matrix(cbind(reformat_glmnet_min1, reformat_glmnet_min2)))
    
    # Test multi-deviance
    dev_glmnet <- multi_deviance_multinomial(data = dfM_test_list, fit_object = glmnet_coefs)
    
    
    ######################### mtool model no penalization ####################
    # Train mtool no penalization
    fit.mtool.relaxed_no_penalization <- 
      multinom.relaxed_lasso_multinomial(train = dfM_train, 
                                         lambda_max = glmnet_lambda_max, gamma = 0, seed = seed, epsilon = 0.001, nfolds = nfolds)
    
    
    # Refit on lambda.min 
    fit_val_min_no_pen <- fit_model_lasso(dfM_test_list, regularization = 'l1',
                               lambda = fit.mtool.relaxed_no_penalization$lambda.min, alpha = 1, unpen_cov = 1)
    
    coef_val1 <- fit_val_min_no_pen$coefficients[1:p, 1]
    coef_val2 <- fit_val_min_no_pen$coefficients[1:p, 2]
    
    # Res
    res_mtool_relaxed_min1_no_pen <- varsel_perc(coef_val1, beta1)
    res_mtool_relaxed_min2_no_pen <- varsel_perc(coef_val2, beta2)
    
    
    # Test multi-deviance
    dev_mtool_relaxed_no_pen <- multi_deviance_multinomial(data = dfM_test_list, fit_object = fit_val_min_no_pen)
    
    
    ######################### mtool model small penalization ####################
    # Train mtool small penalization
    fit.mtool.relaxed_small_penalization <- 
      multinom.relaxed_lasso_multinomial(train = dfM_train, 
                                         lambda_max = glmnet_lambda_max, gamma = 0.00001, seed = seed, epsilon = 0.001, nfolds = nfolds)
    
    
    # Refit on lambda.min 
    fit_val_min_small_pen <- fit_model_lasso(dfM_test_list, regularization = 'l1',
                                   lambda = fit.mtool.relaxed_small_penalization$lambda.min, alpha = 1, unpen_cov = 1)
    
    coef_val1 <- fit_val_min_small_pen$coefficients[1:p, 1]
    coef_val2 <- fit_val_min_small_pen$coefficients[1:p, 2]
    
    # Res
    res_mtool_relaxed_min1_small_pen <- varsel_perc(coef_val1, beta1)
    res_mtool_relaxed_min2_small_pen <- varsel_perc(coef_val2, beta2)
    
    
    # Test multi-deviance
    dev_mtool_relaxed_small_pen <- multi_deviance_multinomial(data = dfM_test_list, fit_object = fit_val_min_small_pen)
    
    
    ######################### mtool model significant penalization ####################
    # Train mtool significant penalization
    fit.mtool.relaxed_significant_penalization <- 
      multinom.relaxed_lasso_multinomial(train = dfM_train, 
                                         lambda_max = glmnet_lambda_max, gamma = 0.01, seed = seed, epsilon = 0.001, nfolds = nfolds)
    
    
    # Refit on lambda.min 
    fit_val_min_sig_pen <- fit_model_lasso(dfM_test_list, regularization = 'l1',
                                   lambda = fit.mtool.relaxed_significant_penalization$lambda.min, alpha = 1, unpen_cov = 1)
    
    coef_val1 <- fit_val_min_sig_pen$coefficients[1:p, 1]
    coef_val2 <- fit_val_min_sig_pen$coefficients[1:p, 2]
    
    # Res
    res_mtool_relaxed_min1_significant_pen <- varsel_perc(coef_val1, beta1)
    res_mtool_relaxed_min2_significant_pen <- varsel_perc(coef_val2, beta2)
    
    
    # Test multi-deviance
    dev_mtool_relaxed_significant_pen <- multi_deviance_multinomial(data = dfM_test_list, fit_object = fit_val_min_sig_pen)
    
    
    
    #############################################################
    Res <- rbind(res_glmnet_min1, res_glmnet_min2, res_mtool_relaxed_min1_no_pen, res_mtool_relaxed_min2_no_pen,
                 res_mtool_relaxed_min1_small_pen, res_mtool_relaxed_min2_small_pen, res_mtool_relaxed_min1_significant_pen,
                 res_mtool_relaxed_min2_significant_pen)
    
    rownames(Res) <- c("glmnet.lambda.min_reference1", "glmnet.lambda.min_reference2", "mtool.relaxed.no.pen.lambda.min_reference1", 
                       "mtool.relaxed.no.pen.lambda.min_reference2", "mtool.relaxed.small.pen.lambda.min_reference1", 
                       "mtool.relaxed.small.pen.lambda.min_reference2", "mtool.relaxed.significant.pen.lambda.min_reference1",
                       "mtool.relaxed.significant.pen.lambda.min_reference2")
    
    devs <- rbind(dev_glmnet, dev_mtool_relaxed_no_pen, dev_mtool_relaxed_small_pen,
                  dev_mtool_relaxed_significant_pen)
    
    rownames(devs) <- c("glmnet.lambda.min",  "mtool.relaxed.no.pen.lambda.min1", 
                       "mtool.relaxed.small.pen.lambda.min", 
                       "mtool.relaxed.significant.pen.lambda.min")
    
    list(Res, devs)
  
  }, simplify = FALSE)
  
    write.csv(Results$Res,
              file = paste("results/", as.character(runif(1)), "_multinomial_results_n=", as.character(n), "_p=",
                           as.character(p), ".csv", sep = ""))
    
    write.csv(Results$devs,
              file = paste("results/", as.character(runif(1)), "_multinomial_devs_n=", as.character(n), "_p=",
                           as.character(p), ".csv", sep = ""))

    # Results <- do.call(rbind, Res)
    # Results
  },
  error = function(e) {
    print(e)
})


