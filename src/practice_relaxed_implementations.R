library(tictoc)
library(caret)
library(dplyr)
library(survival)
library(pec)
library(casebase)
library(future.apply)
library(glmnet)
library(parallelly)
library(timereg)
library(parallel)
library(tictoc)
library(tidyverse)
library(cmprsk)
library(survsim)
library(caret)
library(Matrix)
library(dplyr)




#' ################################################# Practice Relaxed LASSO implementations ##################################################


###########################################################
#' relaxed LASSO function for mtool 
multinom.relaxed_lasso_nnet_old_mtool <- function(train, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
                                        constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                        ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20) {
  # FOR TESTING
  # regularization = 'l1'; lambda_max = NULL; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.009; i = 10; train = train1;
  # p = ncol(train)
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  
  cb_data_train = train
  cb_data_train <- as.data.frame(cb_data_train)
  cb_data_train <- cb_data_train %>%
    select(-time)
  
  
  # Create folds 
  folds <- caret::createFolds(y = cb_data_train$event_ind, k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs_matrix <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  lowest_deviance = .Machine$double.xmax
  
  current_deviance = .Machine$double.xmax
  best_fit = NULL
  min_lambda_index = 0
  
  
  cb_data_all = list("time" = cb_data_train$time,
                     "event_ind" = cb_data_train$event_ind,
                     "covariates" = cb_data_train[, grepl("covariates", names(cb_data_train))],
                     "offset" = cb_data_train$offset)
  
  #Perform 10 fold cross validation
  for(i in 1:nfold){
    
    train_cv <- cb_data_train[which(folds != i), ] #Set the training set
    test_cv <- cb_data_train[which(folds == i), ] #Set the validation set
    # Create X and Y
    train_cv <- list("time" = train_cv$time,
                     "event_ind" = train_cv$event_ind,
                     "covariates" = train_cv[, grepl("covariates", names(train_cv))],
                     "offset" = train_cv$offset)
    
    test_cv <- list("time" = test_cv$time,
                    "event_ind" = test_cv$event_ind,
                    "covariates" = test_cv[, grepl("covariates", names(test_cv))],
                    "offset" = test_cv$offset)
    
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    
    
    res <- unlist(
      mclapply(lambdagrid,
               function(lambda_v) {
                 out = 0
                 # lambda1 <- lambda_v*alpha
                 # lambda2 <- 0.5*lambda_v*(1 - alpha)
                 # mtool
                 fit.mtool.parameterized <- fit_cbmodel(train_cv, regularization = 'l1',
                                                        lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                 coefs_cause1 = fit.mtool.parameterized$coefficients[, 1]
                 coefs_cause2 = fit.mtool.parameterized$coefficients[, 2]
                 selected_cause1 <- which(coefs_cause1 != 0)
                 selected_cause2 <- which(coefs_cause2 != 0)
                 
                 # Maybe try just using first cause 1
                 non_zero_coef_names <- paste("covariates.x", 
                                              as.character(sort(union(selected_cause1, selected_cause2))), 
                                              sep = "")
                 non_zero_coef_names_no_int = non_zero_coef_names[1:length(non_zero_coef_names) - 1]
                 
                 all_zeros = length(non_zero_coef_names) == 1
                 
                 if (all_zeros){
                   # If no predictors selected, skip OLS with cross validation and continue with next lambda
                   
                   test_cv$covariates <- as.matrix(rep(1, nrow(test_cv$covariates)))
                   
                   coef_subset = t(as.matrix(fit.mtool.parameterized$coefficients[nrow(fit.mtool.parameterized$coefficients), ]))
                   out = multi_deviance_fsh(cb_data = test_cv, coefs = coef_subset, add_intercept = FALSE)
                   
                 } else {
                   train_cv$covariates <- train_cv$covariates[, colnames(train_cv$covariates) %in% non_zero_coef_names]
                   test_cv$covariates <- test_cv$covariates[, colnames(test_cv$covariates) %in% non_zero_coef_names]                    
                   
                   p_new = length(non_zero_coef_names_no_int)
                   
                   dfM_subset = cbind.data.frame(train_cv$covariates, train_cv$event_ind)
                   
                   colnames(dfM_subset) = c(non_zero_coef_names_no_int, "y")
                   
                   formula = ""
                   for (j in non_zero_coef_names_no_int) {
                     formula = paste(formula, j, "+")
                   }
                   formula = as.formula(paste("y ~", substr(formula, 2, str_length(formula) - 2)))
                   
                   
                   fit.subset <- nnet::multinom(formula,
                                                data = dfM_subset)
                   
                   coef_subset = rbind(t(coef(fit.subset))[2:(p_new + 1), ],
                                       t(coef(fit.subset))[1, ])
                   
                   out = multi_deviance_fsh(cb_data = test_cv, coefs = coef_subset, add_intercept = TRUE)
                 }
                 
                 if(i == nfold) {
                   fit.all.data = fit_cbmodel(cb_data_all, regularization = 'l1',
                                              lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                   non_zero_coefs = length(union(which(fit.all.data$coefficients[, 1] != 0), 
                                                 which(fit.all.data$coefficients[, 2] != 0))) 
                   out = list(non_zero_coefs, out)
                 }
                 
                 out
               },
               mc.cores = ncores, mc.set.seed = seed))
    if(i == nfold) {
      all_deviances[, i] = res[seq(2, 200, by = 2)]
      all_non_coefs = res[seq(1, 200, by = 2)]
    } else {
      all_deviances[, i] = res
    }
    cat("Completed Fold", i, "\n")
  }
  
  mean_dev <- rowMeans(all_deviances)
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, coef_grid = all_non_coefs,
              cv_se = cv_se, deviance_grid = all_deviances))
}





###########################################################
#' relaxed LASSO function for mtool 
multinom.relaxed_lasso_mtool_old_mtool <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10, 
                                         constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                         ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = NULL) {
  # FOR TESTING
  # regularization = 'l1'; lambda_max = NULL; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.009; i = 10; train = train1;
  
  
  # Create lambda grid
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  
  
  # Create dataset to fit full penalized model to for each value of lambda
  data_all <- train
  
  # Create dataframe object to fit entire dataset to for each value of lambda
  data_df <- as.data.frame(train)
  
  
  # Create folds 
  folds <- caret::createFolds(y = train$event_ind, k = nfold, list = FALSE)
  
  # Initialize matrices of results of cross validation
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs_matrix <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  
  #Perform 10 fold cross validation
  for(i in 1:nfold){
    
    #Segment data by fold
    train_cv <- data_df[which(folds != i), ] #Set the training set
    test_cv <- data_df[which(folds == i), ] #Set the validation set
    
    # Create cb-styled object to pass training + test data to
    train_cv <- list("event_ind" = train_cv$event_ind,
                     "covariates" = train_cv[, grepl("covariates", names(train_cv))],
                     "offset" = train_cv$offset)
    test_cv <- list("event_ind" = test_cv$event_ind,
                    "covariates" = test_cv[, grepl("covariates", names(test_cv))],
                    "offset" = test_cv$offset)
    
    # Scale and center covariates before fitting
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    
    
    res <- unlist(
      mclapply(lambdagrid,
               function(lambda_v) {
                 out = 0
                 
                 # Fit penalized model to get covariate subset for each lambda
                 #TRY fit_cbmodel_lasso here
                 fit.mtool.penalized <- fit_cbmodel(train_cv, regularization = 'l1',
                                                    lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                 
                 # Get list of names of selected covariates
                 coefs_cause1 = fit.mtool.penalized$coefficients[, 1]
                 coefs_cause2 = fit.mtool.penalized$coefficients[, 2]
                 selected_cause1 <- which(coefs_cause1 != 0)
                 selected_cause2 <- which(coefs_cause2 != 0)
                 non_zero_cov_names <- paste("covariates.x", 
                                             as.character(sort(union(selected_cause1, selected_cause2))), 
                                             sep = "")
                 
                 # Subset columns of covariates in training and test set by those that were selected by penalized model
                 train_cv$covariates <- train_cv$covariates[, colnames(train_cv$covariates) %in% non_zero_cov_names]
                 test_cv$covariates <- test_cv$covariates[, colnames(test_cv$covariates) %in% non_zero_cov_names]                    
                 
                 # Fit a multinomial model on each subset (where gamma is penalization value)
                 #TRY fit_cbmodel_lasso here
                 fit.mtool.subset <- fit_cbmodel(train_cv, regularization = 'l1',
                                                 lambda = gamma, alpha = alpha, unpen_cov = 1)
                 
                 # Get selected coefficient estimates from fit on selected covariates
                 coef_subset = fit.mtool.subset$coefficients
                 
                 # Calculate the resulting deviance on the test set
                 out = multi_deviance_fsh(cb_data = test_cv, coefs = coef_subset, add_intercept = TRUE)
                 
                 
                 # Fit a penalized model on entire dataset (all covariates, all observations) during last fold of cv
                 if(i == nfold) {
                   
                   train_col_names = paste("covariates.", colnames(data_all$covariates), sep = "")
                   
                   # Subset columns of covariates 
                   data_all$covariates <- data_all$covariates[, train_col_names %in% non_zero_cov_names]
                   
                   # Fit a multinomial model on subsetted data with a gamma penalization
                   fit.all.data = fit_cbmodel(data_all, regularization = 'l1',
                                              lambda = gamma, alpha = alpha, unpen_cov = 1)
                   num_non_zero_cov = length(union(which(fit.all.data$coefficients[, 1] != 0), 
                                                   which(fit.all.data$coefficients[, 2] != 0)))
                   
                   # return both the number of selected covariates and the test deviance
                   out = list(num_non_zero_cov, out)
                 }
                 
                 out
               },
               mc.cores = ncores, mc.set.seed = seed))
    
    if(i == nfold) {
      all_deviances[, i] = res[seq(2, 200, by = 2)]
      all_non_coefs = res[seq(1, 200, by = 2)]
    } else {
      all_deviances[, i] = res
    }
    cat("Completed Fold", i, "\n")
  }
  
  mean_dev <- rowMeans(all_deviances)
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, coef_grid = all_non_coefs,
              cv_se = cv_se, deviance_grid = all_deviances))
}


# Subsets coefficients correctly
multinom.relaxed_lasso_correct_coefs_old_mtool <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10,
                                                           constant_covariates = 1, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE,
                                                           ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0.009) {

  # FOR TESTING
  # regularization = 'l1'; lambda_max = NULL; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.009; i = 10; train = dfM;
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  
  # Create folds
  folds <- caret::createFolds(y = train$y, k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  
  # Initialize matrices of results of cross validation
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs_matrix <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  
  train_list = list("covariates" = as.matrix(cbind(train[, c(2:ncol(train))])),
                    "event_ind" = train$y,
                    "offset" = rep(0, length(train$y)))
  
  
  
  #Perform 10 fold cross validation
  for(i in 1:nfold){
    
    #Segment your data by fold using the which() function
    train_cv <- train[which(folds != i), ] #Set the training set
    test_cv <- train[which(folds == i), ] #Set the validation set
    
    train_cv = list("covariates" = as.matrix(cbind(train_cv[, c(2:ncol(train_cv))])),
                    "event_ind" = train_cv$y,
                    "offset" = rep(0, length(train_cv$y)))
    test_cv = list("covariates" = as.matrix(cbind(test_cv[, c(2:ncol(test_cv))])),
                   "event_ind" = test_cv$y,
                   "offset" = rep(0, length(test_cv$y)))
    
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    
    # Fit the relaxed LASSO on each value in lambda grid for fold i
    res <- unlist(
      mclapply(lambdagrid,
               function(lambda_v) {
                 dev = 0
                 
                 # Fit penalized model to get covariate subset for each lambda
                 fit.mtool.penalized <- fit_cbmodel(train_cv, regularization = 'l1',
                                                               lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                 
                 # Extract non-zero coefficients
                 coefs_cause1 = fit.mtool.penalized$coefficients[, 1]
                 coefs_cause2 = fit.mtool.penalized$coefficients[, 2]
                 selected_cause1 <- which(coefs_cause1 != 0)
                 selected_cause2 <- which(coefs_cause2 != 0)
                 
                 # Get list of names of selected covariates
                 non_zero_cov_names <- c(paste("x",
                                               as.character(sort(union(selected_cause1, selected_cause2))),
                                               sep = ""))
                 
                 # Subset columns of covariates in training and test set by those that were selected by penalized model
                 train_cv$covariates <- train_cv$covariates[, colnames(train_cv$covariates) %in% non_zero_cov_names]
                 test_cv$covariates <- test_cv$covariates[, colnames(test_cv$covariates) %in% non_zero_cov_names]     
                 
                 # Fit a multinomial model on each subset (where gamma is penalization value)
                 fit.mtool.subset <- fit_cbmodel(train_cv, regularization = 'l1',
                                                            lambda = gamma, alpha = alpha, unpen_cov = 1)
                 
                 
                 # Calculate the resulting deviance on the test set
                 dev = multi_deviance_fsh(cb_data = test_cv, coefs = fit.mtool.subset$coefficients, add_intercept = TRUE)
                 
                 # Fit a penalized model on entire dataset (all covariates, all observations) during last fold of cv
                 if(i == nfold) {
                   # Subset columns of covariates
                   train_list$covariates <- train_list$covariates[, colnames(train_list$covariates) %in% non_zero_cov_names]
                   
                   # Fit a multinomial model on subsetted data with a gamma penalization
                   fit.all.data = fit_cbmodel(train_list, regularization = 'l1',
                                                         lambda = gamma, alpha = alpha, unpen_cov = 1)
                   num_non_zero_cov = length(union(which(fit.all.data$coefficients[, 1] != 0),
                                                   which(fit.all.data$coefficients[, 2] != 0)))
                   
                   # return both the number of selected covariates and the test deviance
                   dev = list(num_non_zero_cov, dev)
                 }
                 
                 # Return the deviance
                 dev
               }, mc.cores = ncores, mc.set.seed = seed))
    if(i == nfold) {
      all_deviances[, i] = res[seq(2, 200, by = 2)]
      all_non_coefs = res[seq(1, 200, by = 2)]
    } else {
      all_deviances[, i] = res
    }
    cat("Completed Fold", i, "\n")
  }
  
  # Return mean deviance across all folds for each lambda
  mean_dev <- rowMeans(all_deviances)
  
  # Find lambda min
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, coef_grid = all_non_coefs,
              cv_se = cv_se, deviance_grid = all_deviances))
}



# Subsets coefficients correctly
multinom.relaxed_lasso_correct_coefs_new_mtool <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10,
                                                           constant_covariates = 1, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE,
                                                           ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0.009) {
  
  # FOR TESTING
  # regularization = 'l1'; lambda_max = NULL; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.009; i = 10; train = dfM;
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  
  # Create folds
  folds <- caret::createFolds(y = train$y, k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  
  # Initialize matrices of results of cross validation
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs_matrix <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  
  train_list = list("covariates" = as.matrix(cbind(train[, c(2:ncol(train))])),
                    "y" = train$y)
  
  
  
  #Perform 10 fold cross validation
  for(i in 1:nfold){
    
    #Segment your data by fold using the which() function
    train_cv <- train[which(folds != i), ] #Set the training set
    test_cv <- train[which(folds == i), ] #Set the validation set
    
    train_cv = list("covariates" = as.matrix(cbind(train_cv[, c(2:ncol(train_cv))])),
                    "y" = train_cv$y)
    test_cv = list("covariates" = as.matrix(cbind(test_cv[, c(2:ncol(test_cv))])),
                   "y" = test_cv$y)
    
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    
    # Fit the relaxed LASSO on each value in lambda grid for fold i
    res <- unlist(
      mclapply(lambdagrid,
               function(lambda_v) {
                 dev = 0
                 
                 # Fit penalized model to get covariate subset for each lambda
                 fit.mtool.penalized <- fit_cbmodel_lasso_mine(train_cv, regularization = 'l1',
                                                               lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                 
                 # Extract non-zero coefficients
                 coefs_cause1 = fit.mtool.penalized$coefficients[, 1]
                 coefs_cause2 = fit.mtool.penalized$coefficients[, 2]
                 selected_cause1 <- which(coefs_cause1 != 0)
                 selected_cause2 <- which(coefs_cause2 != 0)
                 
                 # Get list of names of selected covariates
                 non_zero_cov_names <- c(paste("x",
                                               as.character(sort(union(selected_cause1, selected_cause2))),
                                               sep = ""))
                 
                 # Subset columns of covariates in training and test set by those that were selected by penalized model
                 train_cv$covariates <- train_cv$covariates[, colnames(train_cv$covariates) %in% non_zero_cov_names]
                 test_cv$covariates <- test_cv$covariates[, colnames(test_cv$covariates) %in% non_zero_cov_names]     
                 
                 # Fit a multinomial model on each subset (where gamma is penalization value)
                 fit.mtool.subset <- fit_cbmodel_lasso_mine(train_cv, regularization = 'l1',
                                                            lambda = gamma, alpha = alpha, unpen_cov = 1)
                 
                 
                 # Calculate the resulting deviance on the test set
                 dev = multi_deviance_mine(data = test_cv, fit.mtool.subset)
                 
                 # Fit a penalized model on entire dataset (all covariates, all observations) during last fold of cv
                 if(i == nfold) {
                   # Subset columns of covariates
                   train_list$covariates <- train_list$covariates[, colnames(train_list$covariates) %in% non_zero_cov_names]
                   
                   # Fit a multinomial model on subsetted data with a gamma penalization
                   fit.all.data = fit_cbmodel_lasso_mine(train_list, regularization = 'l1',
                                                         lambda = gamma, alpha = alpha, unpen_cov = 1)
                   num_non_zero_cov = length(union(which(fit.all.data$coefficients[, 1] != 0),
                                                   which(fit.all.data$coefficients[, 2] != 0)))
                   
                   # return both the number of selected covariates and the test deviance
                   dev = list(num_non_zero_cov, dev)
                 }
                 
                 # Return the deviance
                 dev
               }, mc.cores = ncores, mc.set.seed = seed))
    if(i == nfold) {
      all_deviances[, i] = res[seq(2, 200, by = 2)]
      all_non_coefs = res[seq(1, 200, by = 2)]
    } else {
      all_deviances[, i] = res
    }
    cat("Completed Fold", i, "\n")
  }
  
  # Return mean deviance across all folds for each lambda
  mean_dev <- rowMeans(all_deviances)
  
  # Find lambda min
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, coef_grid = all_non_coefs,
              cv_se = cv_se, deviance_grid = all_deviances))
}

###########################################################
#' relaxed LASSO function for mtool 
multinom.relaxed_enet <- function(X, Y, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
                                  constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                  ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20) {
  # FOR TESTING
  # regularization = 'elastic-net'; lambda_max = NULL; alpha = 0.7; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .0001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; i = 1; l = 1;
  p = ncol(train) - 2
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  
  
  # Create folds 
  folds <- caret::createFolds(y = Y, k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs_matrix <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  lowest_deviance = .Machine$double.xmax
  
  all_non_coefs <- c()
  
  
  current_deviance = .Machine$double.xmax
  best_fit = NULL
  min_lambda_index = 0
  
  cb_data_all = list("event_ind" = Y,
                     "covariates" = X,
                     "offset" = rep(0, length(Y)))
  
  #Perform 10 fold cross validation
  for(i in 1:nfold){
    #Segment your data by fold using the which() function 
    X_train <- as.matrix(X[which(folds != i), ]) #Set the training set
    X_val <- as.matrix(X[which(folds == i), ]) #Set the validation set
    Y_train <- as.numeric(Y[which(folds != i)]) - 1 #Set the training set
    Y_val <- as.numeric(Y[which(folds == i)]) - 1 #Set the validation set
    
    # Standardize
    X_train <- scale(X_train, center = T, scale = T)
    X_val <- scale(X_val, center = T, scale = T)
    
    cb_data_train = list("event_ind" = Y_train,
                         "covariates" = X_train,
                         "offset" = rep(0, length(Y_train)))
    
    
    all_deviances[, i] <- unlist(mclapply(lambdagrid,
                                          function(lambda_v) {
                                            fit.mtool.parameterized <- fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                                                                                   lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                                            coefs_cause1 = fit.mtool.parameterized$coefficients[, 1]
                                            coefs_cause2 = fit.mtool.parameterized$coefficients[, 2]
                                            
                                            
                                            selected_cause1 <- which(coefs_cause1 != 0)
                                            selected_cause2 <- which(coefs_cause2 != 0)
                                            non_zero_coef_namesunion(selected_cause1, selected_cause2)
                                            
                                            # Maybe try just using first cause 1
                                            non_zero_coef_names <- paste("x", as.character(sort(union(selected_cause1, selected_cause2))), 
                                                                         sep = "")
                                            
                                            
                                            X_train_new <- X_train[, colnames(X_train) %in% non_zero_coef_names]
                                            X_val_new <- X_val[, colnames(X_train) %in% non_zero_coef_names]
                                            
                                            
                                            cb_data_train_new <- list("event_ind" = Y_train,
                                                                      "covariates" = X_train_new,
                                                                      "offset" = rep(0, length(Y_train)))
                                            
                                            
                                            fit.mtool.subset <- fit_cbmodel(cb_data_train_new, regularization = 'elastic-net',
                                                                            lambda = 0.009, alpha = alpha, unpen_cov = 1)
                                            
                                            multi_deviance_fsh(covs = X_val_new, response = Y_val, coefs = fit.mtool.subset$coefficients, add_intercept = TRUE)
                                          },
                                          mc.cores = ncores, mc.set.seed = seed))
    cat("Completed Fold", i, "\n")
  }
  mean_dev <- rowMeans(all_deviances)
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, non_zero_coefs = non_zero_coefs, lambdagrid = lambdagrid, 
              cv_se = cv_se, deviance_grid = all_deviances, selected_coefs_grid = all_non_coefs))
}




#' Fits linear regression using relaxed LASSO regularization to given data
#' @param: train_data - dataframe containing values of features
#' @param: response - list of response values
#' @param: cv - boolean indicating whether to fit lasso with cross validation or not
#' @param: print_time - boolean indicating whether to print time function takes to run
myRelaxedFinal = function(train_data, response, folds = 5, print_time = FALSE) {
  tryCatch({ 
    if(print_time)
      tic()
    
    coefficient_names = colnames(as.data.frame(train_data))
    response_name = names(as.data.frame(response))
    
    # Data frame and list to keep track of coefficients and CV-MSEs of all lambdas
    coefs_all_lambdas = data.frame(matrix(nrow = 1, ncol = length(coefficient_names) + 1))
    colnames(coefs_all_lambdas) = c("(Intercept)", coefficient_names)
    MSEs_all_lambda = c()
    
    fit_lambdas = glmnet(train_data, response, family="gaussian", keep=TRUE, alpha=1)
    # all_coef_fit_lambdas = as.data.frame.matrix(coef(fit_lambdas))[-1, ]
    all_lambdas = fit_lambdas$lambda
    
    
    # Local variables to track lowest CV-MSE, best fit, and index of lambda of lowest MSE while iterating over all lambdas
    current_MSE = .Machine$double.xmax;
    best_fit = lm(1~1)
    best_lambda_index = 0
    mse_matrix = matrix(NA_real_, nrow = length(all_lambdas), ncol = folds)
    
    # Iterate over all lambdas
    indices = c(1:length(all_lambdas))
    
    for(i in indices) {
      rows = nrow(coefs_all_lambdas)
      # TODO: Look into using predict() with type = "nonzero" as a way to get nonzero coefficients from a glmnet fit for particular lambda values
      # (might be more efficent that using which and subsetting the whole training set, then refitting with lm())
      
      # TODO: replace lm() with: solve(t(X) %*% X) %*% t(X) %*% y
      
      # Perform cross-validation
      # TODO: Check if there's a difference between creating folds outside lambda loop or inside lambda loop
      folds_list <- createFolds(response, k = folds, list = TRUE, returnTrain = TRUE)
      res = list()
      mse_values = c()
      
      for(m in 1:folds) {
        # print(paste("Current lambda and fold indices: ", as.character(i), ", ", as.character(m)))
        train_indices <- folds_list[[m]]
        if (is.matrix(train_data)) {
          k_minus_1_train_data = train_data[train_indices, ]
          validation_fold_data = train_data[-train_indices, ]
        } else {
          k_minus_1_train_data = train_data[train_indices]
          validation_fold_data = train_data[-train_indices]
        }
        y_train = response[train_indices]
        y_valid = response[-train_indices]
        
        fit_en = glmnet(k_minus_1_train_data, y_train, family="gaussian", alpha=1, lambda = all_lambdas[i])
        
        
        # TODO: Replace coef(fit_en) with list of coefficients found using lambda min
        current_coef = coef(fit_en)[-1, ]
        non_zero_coef = current_coef[current_coef != 0]
        all_zeros = length(non_zero_coef) == 0
        
        # If no predictors selected, skip OLS with cross validation and continue with next lambda
        if (all_zeros){
          coefs_all_lambdas[rows,] = 0
          mse_values[m] = mean((y_valid)^2)
          mse_matrix[i, m] = mse_values[m]
          next
        }
        
        non_zero_coef_indices = which(current_coef != 0)
        non_zero_coef_names = coefficient_names[non_zero_coef_indices]
        # print(head(k_minus_1_train_data))
        # print(non_zero_coef_indices)
        new_x_train = k_minus_1_train_data[, non_zero_coef_indices]
        
        if(is.null(ncol(validation_fold_data))) {
          validation_fold_data = as.list(validation_fold_data)
        } else {
          validation_fold_data = as.data.frame.matrix(validation_fold_data)
        }
        
        # lm() requires single dataframe of covariates and response for training
        data = data.frame(new_x_train, y_train)
        # Probably can get rid of this check in case of 0 selected covariates
        if(length(data) < 3) {
          names(data) = c(non_zero_coef_names, response_name)
          names(validation_fold_data) = non_zero_coef_names
        } else {
          colnames(data) = c(non_zero_coef_names, response_name)
          colnames(validation_fold_data) = non_zero_coef_names
        }
        
        form = as.formula(paste(response_name, "~."))
        fit_OLS_on_LASSO_subset = lm(formula = form, data = data, x = TRUE, y = TRUE)
        pred = predict(fit_OLS_on_LASSO_subset, newdata = validation_fold_data)
        mse_values[m] = mean((y_valid - pred)^2)
        mse_matrix[i, m] = mse_values[m]
      }
      
      MSE = mean(mse_values)
      
      # Format new row of coefficient table for current lambda
      coefs_all_lambdas[rows+1,] = 0
      j = 1
      
      for (l in c(1:length(coefs_all_lambdas))) {
        if(is.na(fit_OLS_on_LASSO_subset$coefficients[j]))
          break
        if(colnames(coefs_all_lambdas)[l] == names(fit_OLS_on_LASSO_subset$coefficients[j])) {
          coefs_all_lambdas[i, l] = fit_OLS_on_LASSO_subset$coefficients[j]
          j = j + 1
        }
      } 
      
      MSEs_all_lambda = rbind(MSEs_all_lambda, MSE)
      
      if(MSE < current_MSE) {
        current_MSE = MSE
        best_fit = fit_OLS_on_LASSO_subset
        best_lambda_index = i
      }
    }
  }
  ,
  error = function(e) {
    print(e)
  }
  ,
  warning = function(w) {
  }
  )
  if(print_time)
    toc()
  
  result = list(coefficients = coefs_all_lambdas, CV_MSEs = MSEs_all_lambda, all_folds_MSEs = mse_matrix, all_lambdas = all_lambdas, 
                min_lambda = all_lambdas[best_lambda_index], min_lambda_index = best_lambda_index, 
                best_fit = best_fit)
  return (result)
}



#' Fits linear regression using relaxed LASSO regularization to given data
#' @param: train_data - dataframe containing values of features
#' @param: response - list of response values
#' @param: cv - boolean indicating whether to fit lasso with cross validation or not
#' @param: print_time - boolean indicating whether to print time function takes to run
myRelaxedOld = function(train_data, response, cv, folds = 5, print_time) {
  tryCatch({ 
    if(print_time)
      tic()
    
    coefficient_names = colnames(as.data.frame(train_data))
    response_name = names(as.data.frame(response))
    
    # Data frame and list to keep track of coefficients and CV-MSEs of all lambdas
    coefs_all_lambdas = data.frame(matrix(nrow = 1, ncol = length(coefficient_names) + 1))
    colnames(coefs_all_lambdas) = c("(Intercept)", coefficient_names)
    MSEs_all_lambda = c()
    
    if (cv) {
      fit_lasso = cv.glmnet(train_data, response, family="gaussian", keep=TRUE, alpha=1)
      all_coef_fit_lasso = as.data.frame.matrix(fit_lasso$glmnet.fit$beta)
    } else {
      fit_lasso = glmnet(train_data, response, family="gaussian", alpha=1)
      all_coef_fit_lasso = as.data.frame.matrix(coef(fit_lasso))[-1, ]
    }
    
    # Local variables to track lowest CV-MSE, best fit, and index of lambda of lowest MSE while iterating over all lambdas
    current_MSE = .Machine$double.xmax;
    best_fit = lm(1~1)
    best_lambda_index = 0
    
    # Iterate over all lambdas
    indices = c(1:length(all_coef_fit_lasso))
    
    for(i in indices) {
      rows = nrow(coefs_all_lambdas)
      current_coef = all_coef_fit_lasso[i]
      non_zero_coef = current_coef[all_coef_fit_lasso[i] != 0]
      all_zeros = length(non_zero_coef) == 0
      
      # If no predictors selected, skip OLS with cross validation and continue with next lambda
      if (all_zeros){
        coefs_all_lambdas[rows,] = 0
        next
      }
      
      
      # TODO: Look into using predict() with type = "nonzero" as a way to get nonzero coefficients from a glmnet fit for particular lambda values
      # (might be more efficent that using which and subsetting the whole training set, then refitting with lm())
      
      # Create new training dataframe to pass to lm() for OLS fit on selected predictors
      non_zero_coef_indices = which(current_coef != 0)
      non_zero_coef_names = coefficient_names[non_zero_coef_indices]
      new_x_train = train_data[, non_zero_coef_indices]
      
      # TODO: replace lm() with: solve(t(X) %*% X) %*% t(X) %*% y
      
      #Perform cross-validation
      folds_list <- createFolds(response, k = folds, list = TRUE, returnTrain = TRUE)
      res = list()
      mse_values = c()
      for(m in 1:folds) {
        train_indices <- folds_list[[m]]
        # TODO: Fix code so column names
        # if(is.numeric(new_x_train)) {
        #   names(new_x_train) = c(non_zero_coef_names)
        # } else {
        #   colnames(new_x_train) = c(non_zero_coef_names)
        # }
        if (is.matrix(new_x_train)) {
          traindata = new_x_train[train_indices, ]
          testdata = new_x_train[-train_indices, ]
        } else {
          traindata = new_x_train[train_indices]
          testdata = new_x_train[-train_indices]
        }
        ytrain = response[train_indices]
        ytest = response[-train_indices]
        
        
        if(is.null(ncol(testdata))) {
          testdata = as.list(testdata)
        } else {
          testdata = as.data.frame.matrix(testdata)
        }
        
        
        data = data.frame(traindata, ytrain)
        if(length(data) < 3) {
          names(data) = c(non_zero_coef_names, response_name)
          names(testdata) = non_zero_coef_names
        } else {
          colnames(data) = c(non_zero_coef_names, response_name)
          colnames(testdata) = non_zero_coef_names
        }
        
        form = as.formula(paste(response_name, "~."))
        fit_OLS_on_LASSO_subset = lm(formula = form, data = data, x = TRUE, y = TRUE)
        pred = predict(fit_OLS_on_LASSO_subset, newdata = testdata)
        
        mse_values[m] = mean((ytest - pred)^2)
      }
      
      MSE = mean(mse_values)
      
      # Format new row of coefficient table for current lambda
      coefs_all_lambdas[rows+1,] = 0
      j = 1
      
      for (l in c(1:length(coefs_all_lambdas))) {
        if(is.na(fit_OLS_on_LASSO_subset$coefficients[j]))
          break
        if(colnames(coefs_all_lambdas)[l] == names(fit_OLS_on_LASSO_subset$coefficients[j])) {
          coefs_all_lambdas[i, l] = fit_OLS_on_LASSO_subset$coefficients[j]
          j = j + 1
        }
      } 
      
      MSEs_all_lambda = rbind(MSEs_all_lambda, MSE)
      
      if(MSE < current_MSE) {
        current_MSE = MSE
        best_fit = fit_OLS_on_LASSO_subset
        best_lambda_index = i
      }
    }
  }
  ,
  error = function(e) {
    print(e)
  }
  ,
  warning = function(w) {
  }
  )
  if(print_time)
    toc()
  
  result = list(coefficients = coefs_all_lambdas, CV_MSEs = MSEs_all_lambda, all_lambdas = fit_lasso$lambda, 
                min_lambda = fit_lasso$lambda[best_lambda_index], min_lambda_index = best_lambda_index, 
                best_fit = best_fit)
  
  return (result)
}


#' ################################################# Nirupama's implementations functions ################################################' 


multinom.relaxed_lasso_Nirupama <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10, 
                                            constant_covariates = 1, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                            ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0.009) {
  # # FOR TESTING
  # regularization = 'l1'; lambda_max = NULL; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.009; i = 10; train = dfM;
  
  print("RELAXED")
  p = ncol(train)
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  
  #cb_data_train = train
  #cb_data_train <- as.data.frame(cb_data_train)
  #cb_data_train <- cb_data_train %>%
  #  select(-time)
  # Create folds 
  folds <- caret::createFolds(y = train$y, k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs_matrix <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  lowest_deviance = .Machine$double.xmax
  
  current_deviance = .Machine$double.xmax
  best_fit = NULL
  min_lambda_index = 0
  
  #Perform 10 fold cross validation
  for(i in 1:nfold){
    
    #Segment your data by fold using the which() function 
    train_cv <- train[which(folds != i), ] #Set the training set
    test_cv <- train[which(folds == i), ] #Set the validation set
    
    train_cv = list("covariates" = as.matrix(cbind(train_cv[, c(2:ncol(train_cv))])),
                    "y" = train_cv$y)
    test_cv = list("covariates" = as.matrix(cbind(test_cv[, c(2:ncol(test_cv))])),
                   "y" = test_cv$y)
    
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    
    res <- unlist(
      mclapply(lambdagrid,
               function(lambda_v) {
                 dev = 0
                 fit.mtool.penalized <- fit_cbmodel_lasso(train_cv, regularization = 'l1',
                                                          lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                 coefs_cause1 = fit.mtool.penalized$coefficients[, 1]
                 coefs_cause2 = fit.mtool.penalized$coefficients[, 2]
                 selected_cause1 <- which(coefs_cause1 != 0)
                 selected_cause2 <- which(coefs_cause2 != 0)
                 
                 # Maybe try just using first cause 1
                 non_zero_coef_names <- paste("x", 
                                              as.character(sort(union(selected_cause1, selected_cause2))), 
                                              sep = "")
                 non_zero_coef_names_no_int = non_zero_coef_names[1:length(non_zero_coef_names) - 1]
                 
                 all_zeros = length(non_zero_coef_names) == 1
                 
                 train_cv$covariates <- train_cv[, colnames(train_cv) %in% non_zero_coef_names]
                 test_cv$covariates <- test_cv[, colnames(test_cv) %in% non_zero_coef_names]                    
                 
                 p_new = length(non_zero_coef_names_no_int)
                 
                 dfM_subset = cbind.data.frame(train_cv$covariates, train_cv$y)
                 
                 colnames(dfM_subset) = c(non_zero_coef_names_no_int, "y")
                 
                 fit.mtool.unpenalized <- fit_cbmodel_lasso(train_cv, regularization = 'l1',
                                                            lambda = gamma, alpha = alpha, unpen_cov = 1) 
                 dev = multi_deviance(data = test_cv, fit.mtool.unpenalized)
                 
                 # Fit a penalized model on entire dataset (all covariates, all observations) during last fold of cv
                 if(i == nfold) {
                   # Subset columns of covariates
                   train$covariates <- train[, colnames(train) %in% non_zero_coef_names]
                   
                   # Fit a multinomial model on subsetted data with a gamma penalization
                   fit.all.data = fit_cbmodel_lasso(train, regularization = 'l1',
                                                    lambda = gamma, alpha = alpha, unpen_cov = 1)
                   num_non_zero_cov = length(union(which(fit.all.data$coefficients[, 1] != 0),
                                                   which(fit.all.data$coefficients[, 2] != 0)))
                   
                   # return both the number of selected covariates and the test deviance
                   dev = list(num_non_zero_cov, dev)
                 }
                 
                 
                 dev
               }, mc.cores = ncores, mc.set.seed = seed))
    if(i == nfold) {
      all_deviances[, i] = res[seq(2, 200, by = 2)]
      all_non_coefs = res[seq(1, 200, by = 2)]
    } else {
      all_deviances[, i] = res
    }
    
    cat("Completed Fold", i, "\n")
  }
  
  print(all_deviances)
  mean_dev <- rowMeans(all_deviances)
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, coef_grid = all_non_coefs,
              cv_se = cv_se, deviance_grid = all_deviances))
}



#' ################################################# Final implementations functions ##################################################



#' @param train: dataframe containing training data to fit model to, with first column being the response variable and the rest being the covariate values
# Subsets coefficients correctly
multinom.relaxed_lasso_multinomial <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10,
                                               constant_covariates = 1, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE,
                                               ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0) {
  
  # FOR TESTING
  # regularization = 'l1'; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.009; i = 10; train = dfM;
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  
  # Create folds
  folds <- caret::createFolds(y = train$y, k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  
  # Initialize matrices of results of cross validation
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs_matrix <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  
  train_list = list("covariates" = as.matrix(cbind(train[, c(2:ncol(train))])),
                    "y" = train$y)
  
  
  
  #Perform 10 fold cross validation
  for(i in 1:nfold){
    
    #Segment your data by fold using the which() function
    train_cv <- train[which(folds != i), ] #Set the training set
    test_cv <- train[which(folds == i), ] #Set the validation set
    
    train_cv = list("covariates" = as.matrix(cbind(train_cv[, c(2:ncol(train_cv))])),
                    "y" = train_cv$y)
    test_cv = list("covariates" = as.matrix(cbind(test_cv[, c(2:ncol(test_cv))])),
                   "y" = test_cv$y)
    
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    
    # Fit the relaxed LASSO on each value in lambda grid for fold i
    res <- unlist(
      mclapply(lambdagrid,
               function(lambda_v) {
                 dev = 0
                 
                 # Fit penalized model to get covariate subset for each lambda
                 fit.mtool.penalized <- fit_model_lasso(train_cv, regularization = 'l1',
                                                        lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                 
                 # Extract non-zero coefficients
                 coefs_cause1 = fit.mtool.penalized$coefficients[, 1]
                 coefs_cause2 = fit.mtool.penalized$coefficients[, 2]
                 selected_cause1 <- which(coefs_cause1 != 0)
                 selected_cause2 <- which(coefs_cause2 != 0)
                 
                 # Get list of names of selected covariates
                 non_zero_cov_names <- c(paste("x",
                                               as.character(sort(union(selected_cause1, selected_cause2))),
                                               sep = ""))
                 
                 # Subset columns of covariates in training and test set by those that were selected by penalized model
                 train_cv$covariates <- train_cv$covariates[, colnames(train_cv$covariates) %in% non_zero_cov_names]
                 test_cv$covariates <- test_cv$covariates[, colnames(test_cv$covariates) %in% non_zero_cov_names]     
                 
                 # Fit a multinomial model on each subset (where gamma is penalization value)
                 fit.mtool.subset <- fit_model_lasso(train_cv, regularization = 'l1',
                                                     lambda = gamma, alpha = alpha, unpen_cov = 1)
                 
                 
                 # Calculate the resulting deviance on the test set
                 dev = multi_deviance_multinomial(data = test_cv, fit.mtool.subset)
                 
                 # Fit a penalized model on entire dataset (all covariates, all observations) during last fold of cv
                 if(i == nfold) {
                   # Subset columns of covariates
                   train_list$covariates <- train_list$covariates[, colnames(train_list$covariates) %in% non_zero_cov_names]
                   
                   # Fit a multinomial model on subsetted data with a gamma penalization
                   fit.all.data = fit_model_lasso(train_list, regularization = 'l1',
                                                  lambda = gamma, alpha = alpha, unpen_cov = 1)
                   num_non_zero_cov = length(union(which(fit.all.data$coefficients[, 1] != 0),
                                                   which(fit.all.data$coefficients[, 2] != 0)))
                   
                   # return both the number of selected covariates and the test deviance
                   dev = list(num_non_zero_cov, dev)
                 }
                 
                 # Return the deviance
                 dev
               }, mc.cores = ncores, mc.set.seed = seed))
    if(i == nfold) {
      all_deviances[, i] = res[seq(2, 200, by = 2)]
      all_non_coefs = res[seq(1, 200, by = 2)]
    } else {
      all_deviances[, i] = res
    }
    cat("Completed Fold", i, "\n")
  }
  
  # Return mean deviance across all folds for each lambda
  mean_dev <- rowMeans(all_deviances)
  
  # Find lambda min
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, coef_grid = all_non_coefs,
              cv_se = cv_se, deviance_grid = all_deviances))
}




#' @param train: dataframe containing training data to fit model to, with first column being the response variable and the rest being the covariate values
# Subsets coefficients correctly
multinom.relaxed_lasso_cb <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10,
                                      constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE,
                                      ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0) {
  
  # FOR TESTING
  # regularization = 'l1'; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max_survival; gamma = 0; i = 10; train = train_surv;
  
  
  surv_obj_train <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))
  cov_train <- as.matrix(cbind(train[, c(grepl("X", colnames(train)))], time = log(train$ftime)))
  # Create case-base dataset
  cb_data_train <- create_cbDataset(surv_obj_train, cov_train, ratio =  train_ratio)
  # Default lambda grid
  if(is.null(lambda_max)) {
    # Lambda max grid for bisection search
    if(is.null(initial_max_grid)) {
      initial_max_grid <-  c(0.9, 0.5, 0.1, 0.07, 0.05, 0.01, 0.009, 0.005)
      fit_val_max <- mclapply(initial_max_grid,
                              function(lambda_val) {
                                fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                                            lambda = lambda_val, alpha = alpha, unpen_cov = constant_covariates)}, mc.cores = ncores)
      non_zero_coefs <-  unlist(lapply(fit_val_max, function(x) {return(x$no_non_zero)}))
      if(!isTRUE(any(non_zero_coefs == (constant_covariates*2)))){
        warning("Non-zero coef value not found in default grid. Re-run function and specify initial grid")
      }
      upper <- initial_max_grid[which(non_zero_coefs > (constant_covariates*2 + 1))[1]-1]
      lower <- initial_max_grid[which(non_zero_coefs > (constant_covariates*2 + 1))[1]]
      new_max_searchgrid <- seq(lower, upper, precision)
      fit_val_max <-  mclapply(new_max_searchgrid,
                               function(lambda_val) {
                                 fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                                             lambda = lambda_val, alpha = alpha,
                                             unpen_cov = constant_covariates)}, mc.cores = ncores, mc.set.seed = seed)
      non_zero_coefs <-  unlist(mclapply(fit_val_max, function(x) {return(x$no_non_zero)}, mc.cores = ncores, mc.set.seed = seed))
      lambda_max <- new_max_searchgrid[which.min(non_zero_coefs)]
    }
    # else {
    #   lambda_max = max(initial_max_grid)
    # }
  }
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  # lambdagrid <- lambdagrid[1:60]
  
  cb_data_train <- as.data.frame(cb_data_train)
  cb_data_train <- cb_data_train %>%
    select(-time)
  cb_data_train_list = train_cv <- list("time" = cb_data_train$time,
                                        "event_ind" = cb_data_train$event_ind,
                                        "covariates" = cb_data_train[, grepl("covariates", names(cb_data_train))],
                                        "offset" = cb_data_train$offset)
  # Create folds 
  folds <- caret::createFolds(factor(cb_data_train$event_ind), k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  
  
  
  #Perform 10 fold cross validation
  for(i in 1:nfold){
    
    #Segment your data by fold using the which() function 
    train_cv <- cb_data_train[which(folds != i), ] #Set the training set
    test_cv <- cb_data_train[which(folds == i), ] #Set the validation set
    # Create X and Y
    train_cv <- list("time" = train_cv$time,
                     "event_ind" = train_cv$event_ind,
                     "covariates" = train_cv[, grepl("covariates", names(train_cv))],
                     "offset" = train_cv$offset)
    test_cv <- list("time" = test_cv$covariates.time,
                    "event_ind" = test_cv$event_ind,
                    "covariates" = as.matrix(test_cv[, grepl("covariates", names(test_cv))]),
                    "offset" = test_cv$offset)
    
    
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    
    # Fit the relaxed LASSO on each value in lambda grid for fold i
    res <- unlist(
      mclapply(lambdagrid,
               function(lambda_v) {
                 dev = 0
                 
                 # Fit penalized model to get covariate subset for each lambda
                 fit.mtool.penalized <- fit_cbmodel_lasso(train_cv, regularization = 'l1',
                                                          lambda = lambda_v, alpha = alpha, unpen_cov = constant_covariates)
                 
                 # Extract non-zero coefficients
                 coefs_cause1 = fit.mtool.penalized$coefficients[, 1]
                 coefs_cause2 = fit.mtool.penalized$coefficients[, 2]
                 selected_cause1 <- which(coefs_cause1 != 0)
                 selected_cause2 <- which(coefs_cause2 != 0)
                 
                 # Get list of names of selected covariates
                 non_zero_cov_names <- c(paste("covariates.X",
                                               as.character(sort(union(selected_cause1, selected_cause2))),
                                               sep = ""))
                 
                 # Subset columns of covariates in training and test set by those that were selected by penalized model
                 train_cv$covariates <- train_cv$covariates[, colnames(train_cv$covariates) %in% non_zero_cov_names]
                 test_cv$covariates <- test_cv$covariates[, colnames(test_cv$covariates) %in% non_zero_cov_names]     
                 
                 # Fit a multinomial model on each subset (where gamma is penalization value)
                 fit.mtool.subset <- fit_cbmodel_lasso(train_cv, regularization = 'l1',
                                                       lambda = gamma, alpha = alpha, unpen_cov = constant_covariates)
                 
                 
                 # Calculate the resulting deviance on the test set
                 dev = multi_deviance_cb(cb_data = test_cv, fit.mtool.subset)
                 
                 # Fit a penalized model on entire dataset (all covariates, all observations) during last fold of cv
                 if(i == nfold) {
                   # Subset columns of covariates
                   cb_data_train_list$covariates <- cb_data_train_list$covariates[, colnames(cb_data_train_list$covariates) %in% non_zero_cov_names]
                   
                   # Fit a multinomial model on subsetted data with a gamma penalization
                   fit.all.data = fit_cbmodel_lasso(cb_data_train_list, regularization = 'l1',
                                                    lambda = gamma, alpha = alpha, unpen_cov = constant_covariates)
                   num_non_zero_cov = length(union(which(fit.all.data$coefficients[, 1] != 0),
                                                   which(fit.all.data$coefficients[, 2] != 0)))
                   
                   # return both the number of selected covariates and the test deviance
                   dev = list(num_non_zero_cov, dev)
                 }
                 
                 # Return the deviance
                 dev
               }, mc.cores = ncores, mc.set.seed = seed))
    if(i == nfold) {
      all_deviances[, i] = res[seq(2, 2 * length(lambdagrid), by = 2)]
      all_non_coefs = res[seq(1, 2 * length(lambdagrid), by = 2)]
    } else {
      all_deviances[, i] = res
    }
    cat("Completed Fold", i, "\n")
  }
  
  # Return mean deviance across all folds for each lambda
  mean_dev <- rowMeans(all_deviances)
  
  # Find lambda min
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, coef_grid = all_non_coefs,
              cv_se = cv_se, deviance_grid = all_deviances))
}



#' ################################################# Relaxed implementatson helper functions ##################################################



#' Multinomial deviance
#' 
#' @param cb_data Output of \code{create_cbDataset}
#' @param fit_object Output of \code{fit_cbmodel}
#' @return Multinomial deviance
multi_deviance <- function(data, fit_object) {
  X <- as.matrix(cbind(data[, c(2:ncol(data))], 1))
  print(head(X))
  print(dim(X))
  fitted_vals <- as.matrix(X %*% fit_object$coefficients)
  pred_mat <- VGAM::multilogitlink(fitted_vals, 
                                   inverse = TRUE)
  # Turn event_ind into Y_mat
  Y_fct <- factor(data$y)
  Y_levels <- levels(Y_fct)
  Y_mat <- matrix(NA_integer_, ncol = length(Y_levels),
                  nrow = nrow(X))
  for (i in seq_len(length(Y_levels))) {
    Y_mat[,i] <- (Y_fct == Y_levels[i])
  }
  
  dev <- VGAM::multinomial()@deviance(pred_mat, Y_mat, 
                                      w = rep(1, nrow(X)))
  
  return(dev)
}




############################################
#' Fit case-base sampling model
#'
#' @param cb_data Output of \code{create_cbDataset}
fit_cbmodel_lasso <- function(data, regularization = 'l1',
                              lambda, alpha = 1, unpen_cov = 1) {
  stopifnot(alpha >= 0 && alpha <= 1)
  
  # Prepare covariate matrix with intercept
  X <- as.matrix(cbind(data[, c(2:ncol(data))], 1))
  print(head(X))
  print(colnames(X))
  print(class(data$y))
  out <- fit.mtool <- mtool::mtool.MNlogistic(
    X = as.matrix(X),
    Y = data$y,
    offset =  rep(0, length(data$y)),
    N_covariates = unpen_cov,
    regularization = 'l1',
    transpose = FALSE,
    lambda1 = lambda, lambda2 = 0,
    lambda3 = 0
  )
  
  return(out)
}


############################################
#' Fit case-base sampling model
#'
#' @param cb_data Output of \code{create_cbDataset}
fit_cbmodel_lasso_mine <- function(data, regularization = 'l1',
                                   lambda, alpha = 1, unpen_cov = 1) {
  stopifnot(alpha >= 0 && alpha <= 1)
  
  # Prepare covariate matrix with intercept
  X <- as.matrix(cbind(data$covariates, 1))
  print(head(X))
  print(dim(X))
  print(length(data$y))
  out <- fit.mtool <- mtool::mtool.MNlogistic(
    X = as.matrix(X),
    Y = data$y,
    offset =  rep(0, length(data$y)),
    N_covariates = unpen_cov,
    regularization = 'l1',
    # learning_rate = 0.001,
    transpose = FALSE,
    lambda1 = lambda, lambda2 = 0,
    lambda3 = 0
  )
  
  return(out)
}


#' Multinomial deviance
#' 
#' @param cb_data Output of \code{create_cbDataset}
#' @param fit_object Output of \code{fit_cbmodel}
#' @return Multinomial deviance
multi_deviance_mine <- function(data, fit_object) {
  X <- as.matrix(cbind(data$covariates, 1))
  print(head(X))
  print(dim(X))
  print(length(data$y))
  fitted_vals <- as.matrix(X %*% fit_object$coefficients)
  pred_mat <- VGAM::multilogitlink(fitted_vals, 
                                   inverse = TRUE)
  # Turn event_ind into Y_mat
  Y_fct <- factor(data$y)
  Y_levels <- levels(Y_fct)
  Y_mat <- matrix(NA_integer_, ncol = length(Y_levels),
                  nrow = nrow(X))
  for (i in seq_len(length(Y_levels))) {
    Y_mat[,i] <- (Y_fct == Y_levels[i])
  }
  
  dev <- VGAM::multinomial()@deviance(pred_mat, Y_mat, 
                                      w = rep(1, nrow(X)))
  
  return(dev)
}


################################################
#' Fit case-base sampling model
#' 
#' @param data Output of \code{create_cbDataset}
fit_cbmodel_mine <- function(data, regularization = 'elastic-net',
                             lambda, alpha = 0.5, unpen_cov = 1) {
  stopifnot(alpha >= 0 && alpha <= 1)
  # Elastic-net reparametrization
  lambda1 <- lambda*alpha
  lambda2 <- 0.5*lambda*(1 - alpha)
  # Prepare covariate matrix with intercept
  X <- as.matrix(cbind(data$covariates, 1))
  print(head(X))
  print(head(data$y))
  out <- fit.mtool <- mtool::mtool.MNlogistic(
    X = as.matrix(X),
    Y = data$y,
    offset = rep(0, length(data$y)),
    N_covariates = unpen_cov,
    regularization = 'elastic-net',
    transpose = FALSE,
    lambda1 = lambda1, lambda2 = lambda2, 
    lambda3 = 0
  )
  
  return(out)
}


################################################
#' Fit case-base sampling model
#' 
#' @param cb_data Output of \code{create_cbDataset}
fit_cbmodel <- function(cb_data, regularization = 'elastic-net',
                        lambda, alpha = 0.5, unpen_cov = 1) {
  stopifnot(alpha >= 0 && alpha <= 1)
  # Elastic-net reparametrization
  lambda1 <- lambda*alpha
  lambda2 <- 0.5*lambda*(1 - alpha)
  # Prepare covariate matrix with intercept
  X <- as.matrix(cbind(cb_data$covariates, 1))
  print(head(X))
  print(head(cb_data$event_ind))
  out <- fit.mtool <- mtool::mtool.MNlogistic(
    X = as.matrix(X),
    Y = cb_data$event_ind,
    offset = cb_data$offset,
    N_covariates = unpen_cov,
    regularization = 'elastic-net',
    transpose = FALSE,
    lambda1 = lambda1, lambda2 = lambda2, 
    lambda3 = 0
  )
  
  return(out)
}



###########################################################
#' Calculate multinomial deviance on fitSmoothHazard object
multi_deviance_fsh <- function(cb_data, coefs, add_intercept) {
  if(add_intercept) {
    X <- as.matrix(cbind(cb_data$covariates, 1))
  } else {
    X = cb_data$covariates
  }
  
  print(head(X))
  print(dim(X))
  print(dim(coefs))
  print(class(cb_data$event_ind))
  
  fitted_vals <- as.matrix(X %*% coefs)
  pred_mat <- VGAM::multilogitlink(fitted_vals, 
                                   inverse = TRUE)
  # Turn event_ind into Y_mat
  Y_fct <- factor(cb_data$event_ind)
  Y_levels <- levels(Y_fct)
  Y_mat <- matrix(NA_integer_, ncol = length(Y_levels),
                  nrow = nrow(X))
  for (i in seq_len(length(Y_levels))) {
    Y_mat[,i] <- (Y_fct == Y_levels[i])
  }
  
  dev <- VGAM::multinomial()@deviance(pred_mat, Y_mat, 
                                      w = rep(1, nrow(X)))
  
  return(dev)
} 


#' Multinomial deviance
#' 
#' @param cb_data Output of \code{create_cbDataset}
#' @param fit_object Output of \code{fit_cbmodel}
#' @return Multinomial deviance
multi_deviance_cox <- function(cb_data, fit_object) {
  X <- as.matrix(cb_data$covariates)
  fitted_vals <- as.matrix(X %*% fit_object$coefficients)
  pred_mat <- VGAM::multilogitlink(fitted_vals, 
                                   inverse = TRUE)
  # Turn event_ind into Y_mat
  Y_fct <- factor(cb_data$event_ind)
  Y_levels <- levels(Y_fct)
  Y_mat <- matrix(NA_integer_, ncol = length(Y_levels),
                  nrow = nrow(X))
  for (i in seq_len(length(Y_levels))) {
    Y_mat[,i] <- (Y_fct == Y_levels[i])
  }
  
  dev <- VGAM::multinomial()@deviance(pred_mat, Y_mat, 
                                      w = rep(1, nrow(X)))
  
  return(dev)
} 



############################################
#' Fit multinomial logistic regression
#'
#' @param data: list of data with elements:
#'    1) covariates: a dataframe of covariate values
#'    2) y: a list of response values

fit_model_lasso <- function(data, regularization = 'l1',
                            lambda, alpha = 1, unpen_cov = 1) {
  stopifnot(alpha >= 0 && alpha <= 1)
  
  # Prepare covariate matrix with intercept
  X <- as.matrix(cbind(data$covariates, 1))
  print(head(X))
  print(dim(X))
  print(length(data$y))
  out <- fit.mtool <- mtool::mtool.MNlogistic(
    X = as.matrix(X),
    Y = data$y,
    offset =  rep(0, length(data$y)),
    N_covariates = unpen_cov,
    regularization = 'l1',
    # learning_rate = 0.001,
    transpose = FALSE,
    lambda1 = lambda, lambda2 = 0,
    lambda3 = 0
  )
  
  return(out)
}

#' Multinomial deviance for multinomial data
#' 
#' @param data: list of data with elements:
#'    1) covariates: a dataframe of covariate values
#'    2) y: a list of response values
#' @param fit_object Output of \code{fit_model_lasso}
#' @return Multinomial deviance
multi_deviance_multinomial <- function(data, fit_object) {
  X <- as.matrix(cbind(data$covariates, 1))
  print(head(X))
  print(dim(X))
  print(length(data$y))
  fitted_vals <- as.matrix(X %*% fit_object$coefficients)
  pred_mat <- VGAM::multilogitlink(fitted_vals, 
                                   inverse = TRUE)
  # Turn event_ind into Y_mat
  Y_fct <- factor(data$y)
  Y_levels <- levels(Y_fct)
  Y_mat <- matrix(NA_integer_, ncol = length(Y_levels),
                  nrow = nrow(X))
  for (i in seq_len(length(Y_levels))) {
    Y_mat[,i] <- (Y_fct == Y_levels[i])
  }
  
  dev <- VGAM::multinomial()@deviance(pred_mat, Y_mat, 
                                      w = rep(1, nrow(X)))
  
  return(dev)
}



############################################
#' Fit case-base sampling model
#'
#' @param cb_data Output of \code{create_cbDataset}
fit_cbmodel_lasso <- function(cb_data, regularization = 'l1',
                              lambda, alpha = 1, unpen_cov = 1) {
  stopifnot(alpha >= 0 && alpha <= 1)
  
  # Prepare covariate matrix with intercept
  X <- as.matrix(cbind(cb_data$covariates, 1))
  # print(head(X))
  # print(dim(X))
  # print(length(cb_data$event_ind))
  # print(head(cb_data$offset))
  out <- fit.mtool <- mtool::mtool.MNlogistic(
    X = as.matrix(X),
    Y = cb_data$event_ind,
    offset = cb_data$offset,
    N_covariates = unpen_cov,
    regularization = regularization,
    # learning_rate = 0.001,
    transpose = FALSE,
    lambda1 = lambda, lambda2 = 0,
    lambda3 = 0
  )
  
  return(out)
}

#' Multinomial deviance
#' 
#' @param cb_cb_data Output of \code{create_cbcb_dataset}
#' @param fit_object Output of \code{fit_cbmodel}
#' @return Multinomial deviance
multi_deviance_cb <- function(cb_data, fit_object) {
  X <- as.matrix(cbind(cb_data$covariates, 1))
  print(head(X))
  print(dim(X))
  print(length(cb_data$event_ind))
  fitted_vals <- as.matrix(X %*% fit_object$coefficients)
  pred_mat <- VGAM::multilogitlink(fitted_vals, 
                                   inverse = TRUE)
  # Turn event_ind into Y_mat
  Y_fct <- factor(cb_data$event_ind)
  Y_levels <- levels(Y_fct)
  Y_mat <- matrix(NA_integer_, ncol = length(Y_levels),
                  nrow = nrow(X))
  for (i in seq_len(length(Y_levels))) {
    Y_mat[,i] <- (Y_fct == Y_levels[i])
  }
  
  dev <- VGAM::multinomial()@deviance(pred_mat, Y_mat, 
                                      w = rep(1, nrow(X)))
  
  return(dev)
}




plot_post_relaxed.multinom_test <- function(cv_object) {
  num_lambdas = length(cv_object$lambdagrid)
  no_coef = cv_object$coef_grid[seq(1, num_lambdas, by = 5)]
  nfold <- ncol(cv_object$deviance_grid)
  mean_dev <- rowMeans(cv_object$deviance_grid)
  row_stdev <- apply(cv_object$deviance_grid, 1, function(x) {sd(x)/sqrt(nfold)})
  plot.dat.p <- data.frame(lambdagrid = cv_object$lambdagrid, mean.dev = mean_dev, 
                           upper = mean_dev +row_stdev, lower = mean_dev - row_stdev)
  p <- ggplot(plot.dat.p, aes(log(lambdagrid), mean.dev)) + geom_point(colour = "red", size = 3) + theme_bw() + 
    geom_errorbar(aes(ymin= lower, ymax=upper), width=.2, colour = "grey") + 
    labs(x = "log(lambda)", y = "Multinomial Deviance")  + 
    geom_vline(xintercept = log(cv_object$lambda.min), linetype = "dotted", colour = "blue")
  
  positions = log(cv_object$lambdagrid)[seq(1, num_lambdas, by = 5)]
  labels = cv_object$coef_grid[seq(1, num_lambdas, by = 5)]
  
  p1 <- add_secondary_axis(p, positions, labels)
  
  return(p1)
}

# Define a function to add secondary axis
add_secondary_axis <- function(p, positions, labels) {
  print(positions)
  print(labels)
  n <- length(positions)
  
  if (n != length(labels)) {
    stop("positions and labels must have the same length")
  }
  
  # Get the current x-axis limits
  x_limits <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
  
  # Filter positions and labels that are within x-axis limits
  valid_indices <- which(positions >= x_limits[1] & positions <= x_limits[2])
  
  if (length(valid_indices) < n) {
    warning("Some positions are outside the x-axis limits and will be ignored.")
    positions <- positions[valid_indices]
    labels <- labels[valid_indices]
  }
  
  # Dynamically get the y-limit based on the data in the plot
  max_error_bar = max(ggplot_build(p)$data[[2]]$ymax)
  y_data <- rep(max_error_bar, length(ggplot_build(p)$data[[2]]$ymax))
  y_limit <- max(y_data) + 0.01 * diff(range(y_data))
  
  # Add text labels and tick marks
  for (i in seq_along(positions)) {
    # Text labels
    p <- p + 
      annotation_custom(
        grob = grid::textGrob(label = labels[i], hjust = 0.5, gp = gpar(cex = 0.9)),
        xmin = positions[i], xmax = positions[i],
        ymin = y_limit, ymax = y_limit
      )
  }
  
  return(p)
}