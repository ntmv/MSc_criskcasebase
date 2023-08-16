################################################# My implementation functions ##################################################



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
  regularization = 'l1'; lambda_max = NULL; alpha = 1; nfold = 10;
  constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.009; i = 10; train = train1;
  
  
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
  regularization = 'l1'; lambda_max = NULL; alpha = 1; nfold = 10;
  constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.009; i = 10; train = dfM;
  
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







plot_post_relaxed.multinom_test <- function(cv_object) {
  no_coef = cv_object$coef_grid[seq(1, 100, by = 5)]
  nfold <- ncol(cv_object$deviance_grid)
  mean_dev <- rowMeans(cv_object$deviance_grid)
  row_stdev <- apply(cv_object$deviance_grid, 1, function(x) {sd(x)/sqrt(nfold)})
  plot.dat.p <- data.frame(lambdagrid = cv_object$lambdagrid, mean.dev = mean_dev, 
                           upper = mean_dev +row_stdev, lower = mean_dev - row_stdev)
  p <- ggplot(plot.dat.p, aes(log(lambdagrid), mean.dev)) + geom_point(colour = "red", size = 3) + theme_bw() + 
    geom_errorbar(aes(ymin= lower, ymax=upper), width=.2, colour = "grey") + 
    labs(x = "log(lambda)", y = "Multinomial Deviance")  + 
    geom_vline(xintercept = log(cv_object$lambda.min), linetype = "dotted", colour = "blue")
  
  positions = log(cv_object$lambdagrid)[seq(1, 100, by = 5)]
  labels = cv_object$coef_grid[seq(1, 100, by = 5)]
  
  p1 <- add_secondary_axis(p, positions, labels)
  
  return(p1)
}

# Define a function to add secondary axis
add_secondary_axis <- function(p, positions, labels) {
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


#' ################################################# Nirupama's implementations functions ##################################################
#' 





multinom.relaxed_lasso_Nirupama <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10, 
                                   constant_covariates = 1, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                   ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0.009) {
  # FOR TESTING
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

  mean_dev <- rowMeans(all_deviances)
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, coef_grid = all_non_coefs,
              cv_se = cv_se, deviance_grid = all_deviances))
}




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

