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



generateDataset = function(p, n) {
  beta = c(rep(5, p/2 + 1), rep(0, p/2))
  beta = as.matrix(beta)
  beta_neg = -beta
  
  
  zero_ind1 <- which(beta == 0)
  nonzero_ind1 <- which(beta != 0)
  
  
  # Generate X (iid case)
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  X = cbind(rep(1, times = n), X)
  
  
  #Generate Y
  epsilon <- rnorm(n, sd = 0.5)
  Y1 = X %*% beta + epsilon
  Y2 = X %*% beta_neg + epsilon
  
  result = list(X = X, Y = Y1)
  return(result)
}


partitionData = function(x, y) {
  train_index_y = caret::createDataPartition(y, p = 0.80, list = FALSE)
  
  x_train = x[train_index_y, -1]
  y_train = y[train_index_y]
  
  x_test = x[-train_index_y, -1]
  y_test = y[-train_index_y]
  
  result = list(x_train = x_train, y_train = y_train,
                x_test = x_test, y_test = y_test)
  
  return(result)
}



fitAll = function(train_data, response, folds = 5, print_time) {
  # cv.glmnet with LASSO, cv.glmnet relaxed LASSO, and my relaxed implementation fits
  fit = cv.glmnet(train_data, response, gamma = 0, nfolds = folds, relax = FALSE)
  relaxed_fit = cv.glmnet(train_data, response, gamma = 0, nfolds = folds, relax = TRUE)
  my_relaxed_fit = myRelaxedOld(train_data, response, FALSE, folds = folds, print_time)
  my_relaxed_fit_cv = myRelaxedOld(train_data, response, TRUE, folds = folds, print_time)
  my_relaxed_new_fit_cv = myRelaxedFinal(train_data, response, folds = folds, print_time)
  
  
  # Save lambda values of interest from cv.glmnet LASSO to use for post LASSO fit
  lambda_min = fit$lambda.min

  # Post LASSO fits
  fit_post_lasso_lambda_min = glmnet(train_data, response, lambda = lambda_min)

  result = list(cv_glmnet_relaxed = relaxed_fit, my_relaxed = my_relaxed_fit,
                my_relaxed_cv = my_relaxed_fit_cv, my_relaxed_new = my_relaxed_new_fit_cv, 
                post_lasso_lambda_min = fit_post_lasso_lambda_min)
  
  return(result)
}


computeMSEs = function(resulting_fits, test_data) {
  pred_fit_relaxed_min = predict(resulting_fits$cv_glmnet_relaxed, newx = test_data$x_test, s = resulting_fits$cv_glmnet_relaxed$relaxed$lambda.min)

  pred_myRelaxed = predict.lm(resulting_fits$my_relaxed$best_fit, newdata = as.data.frame(test_data$x_test))
  pred_myRelaxed_cv = predict.lm(resulting_fits$my_relaxed_cv$best_fit, newdata = as.data.frame(test_data$x_test))
  pred_myRelaxed_new_cv = predict.lm(resulting_fits$my_relaxed_new$best_fit, newdata = as.data.frame(test_data$x_test))
  
  coefficient_names = rownames(resulting_fits$post_lasso_lambda_min$beta)
  indices = which(colnames(as.data.frame.matrix(test_data$x_test)) %in% coefficient_names)
  new_x_test = test_data$x_test[, indices]
  colnames(new_x_test) = coefficient_names
  pred_fit_post_lasso_min = predict(resulting_fits$post_lasso_lambda_min, newx = new_x_test)

  
  # Compute MSE
  mse_relaxed_min <- mean((pred_fit_relaxed_min - test_data$y_test)^2)
  mse_myRelaxed <- mean((pred_myRelaxed - test_data$y_test)^2)
  mse_myRelaxed_cv <- mean((pred_myRelaxed_cv - test_data$y_test)^2)
  mse_myRelaxed_new_cv <- mean((pred_myRelaxed_new_cv - test_data$y_test)^2)
  mse_post_lasso_min <- mean((pred_fit_post_lasso_min - test_data$y_test)^2)

  MSEs = data.frame(matrix(nrow = 5, ncol = 2))
  colnames(MSEs) = c("Model", "Test MSE")
  
  model_fit_labels = c("cv.glmnet relax = TRUE, lambda min",
                       "Relaxed LASSO implementation, glmnet(), start CV AFTER predictors selected", 
                       "Relaxed LASSO implementation, cv.glmnet(), start CV AFTER predictors selected",
                       "Relaxed LASSO implementation, glmnet(), start CV BEFORE predictors selected",
                       "Post LASSO lambda min")
  
  MSEs[, 1] = model_fit_labels
  MSEs[, 2] = c(mse_relaxed_min, mse_myRelaxed, mse_myRelaxed_cv,
                mse_myRelaxed_new_cv, mse_post_lasso_min)
  
  return(MSEs)
}



formatCoefficientTable = function(p, n, coefficient_names, resulting_fits, print_time) {
  model_fit_labels = c("cv.glmnet relax = TRUE, lambda min",
                       "Relaxed LASSO implementation, glmnet(), start CV AFTER predictors selected", 
                       "Relaxed LASSO implementation, cv.glmnet(), start CV AFTER predictors selected",
                       "Relaxed LASSO implementation, glmnet(), start CV BEFORE predictors selected",
                       "Post LASSO lambda min")
  
  # Create dataframes to keep track of coefficient values and MSEs per each iteration of simulation
  coefficient_values = data.frame(matrix(nrow = 0, ncol = length(coefficient_names)))
  
  col_names = rownames(as.data.frame.matrix(coef(resulting_fits$cv_glmnet_relaxed)))
  coef_values = as.numeric(coef(resulting_fits$post_lasso_lambda_min))
  coef_names = rownames(coef(resulting_fits$post_lasso_lambda_min))
  names(coef_values) = coef_names
  y_post_lasso_lambda_min_coefficient_row = formatCoefficientTableRow(coef_values, col_names)
  
  
  
  coefficient_values = rbind(coefficient_values, 
                             c(0, resulting_fits$cv_glmnet_relaxed$glmnet.fit$relaxed$beta[, resulting_fits$cv_glmnet_relaxed$relaxed$index[1]]),
                             c(0, resulting_fits$my_relaxed$coefficients[resulting_fits$my_relaxed$min_lambda_index, -1]),
                             c(0, resulting_fits$my_relaxed_cv$coefficients[resulting_fits$my_relaxed_cv$min_lambda_index, -1]),
                             c(0, resulting_fits$my_relaxed_new$coefficients[resulting_fits$my_relaxed_new$min_lambda_index, -1]),
                             c(0, y_post_lasso_lambda_min_coefficient_row))
  options(digits = 5)
  colnames(coefficient_values) = coefficient_names
  coefficient_values[, 1] = model_fit_labels
  
  return(coefficient_values)
}


oneIteration = function(p, n, training_data, test_data, print_time) {
  tryCatch({
    # Fit training data with different models
    resulting_fits = fitAll(training_data$x_train, training_data$y_train, folds = 5, FALSE)
    
    # Get table of coefficient values for each model
    coefficient_names = colnames(as.data.frame(training_data$x_train))
    coefficient_names = c("Model", coefficient_names)
    coefs = formatCoefficientTable(p, n, coefficient_names, resulting_fits, print_time)
    
    # Compute biases
    MSEs = computeMSEs(resulting_fits, test_data)
    result = list(coefficient_table = coefs, MSE_table = MSEs)
    
    return(result)
  },
  error = function(e) {
    print(e)
  })
}



runSim = function(p, n, N, print_time) {
  # Generate/split data
  dataset = generateDataset(p, n)
  partitioned_data = partitionData(dataset$X, dataset$Y)
  test_data = partitioned_data[c("x_test", "y_test")]
  
  # Run simulation N times
  sim_replicates = replicate(N, {
    dataset = generateDataset(p, n)
    partitioned_data = partitionData(dataset$X, dataset$Y)
    training_data = partitioned_data[c("x_train", "y_train")]
    
    result_one_iter = oneIteration(p, n, training_data, test_data, print_time)
    
    result = list(coefficient_table = result_one_iter[[1]], MSE_table = result_one_iter[[2]])
    
    return(result)
  }, simplify = FALSE)
  
  
  #TODO: Fix convoluted method of binding resulting tables from simulation
  sim_coefficient_results = data.frame(matrix(nrow = 0, ncol = ncol(sim_replicates[[1]][[1]])))
  sim_MSE_results = data.frame(matrix(nrow = 0, ncol = ncol(sim_replicates[[1]][[2]])))
  for(i in c(1:N)) {
    sim_coefficient_results = rbind(sim_coefficient_results, sim_replicates[[i]][[1]])
    sim_MSE_results = rbind(sim_MSE_results, sim_replicates[[i]][[2]])
  }
  
  return(list(sim_coefficient_table = sim_coefficient_results, sim_MSE_table = sim_MSE_results))
}


formatCoefficientTableRow = function(non_zero_coefficients, column_names) {
  non_zero_coefficient_names = names(non_zero_coefficients)
  result = c()
  i = 1
  
  if(non_zero_coefficient_names[i] == "(Intercept)") {
    i = i + 1
  }
  
  for (l in c(1:length(column_names))) {
    if(column_names[l] == "(Intercept)") {
      next
    }
    if (column_names[l] %in% non_zero_coefficient_names[i]) {
      result = c(result, as.numeric(non_zero_coefficients[i]))
      i = i + 1
    } else {
      result = c(result, 0)
    }
  }
  return(result)
}


formatCoefficientBiasTable = function(sim_results, p) {
  betas = c(rep(5, p/2), rep(0, p/2))
  
  model_fit_labels = c("cv.glmnet relax = TRUE, lambda min",
                       "Relaxed LASSO implementation, glmnet(), start CV AFTER predictors selected", 
                       "Relaxed LASSO implementation, cv.glmnet(), start CV AFTER predictors selected",
                       "Relaxed LASSO implementation, glmnet(), start CV BEFORE predictors selected",
                       "Post LASSO lambda min")
  
  bias_table = data.frame(matrix(nrow = 0, ncol = length(sim_results$sim_coefficient_table)))
  coefficient_names = colnames(as.data.frame(sim_results$sim_coefficient_table))
  colnames(bias_table) = coefficient_names
  
  for (i in model_fit_labels) {
    average_coef_row = colMeans((sim_results$sim_coefficient_table %>% filter(Model == i))[, -1])
    average_coef_row = average_coef_row - betas
    bias_table = rbind(bias_table, c(i, average_coef_row))
  }
  colnames(bias_table) = coefficient_names
  bias_table[, -1] = lapply(bias_table[, -1], as.numeric)
  
  options(digits = 5)
  
  bias_table[, -1] = round_df(bias_table[, -1], 3)
  
  return(bias_table)
}



formatAverageTestMSETable = function(sim_results) {
  model_fit_labels = c("cv.glmnet relax = TRUE, lambda min",
                       "Relaxed LASSO implementation, glmnet(), start CV AFTER predictors selected", 
                       "Relaxed LASSO implementation, cv.glmnet(), start CV AFTER predictors selected",
                       "Relaxed LASSO implementation, glmnet(), start CV BEFORE predictors selected",
                       "Post LASSO lambda min")
  
  average_MSE_table = data.frame(matrix(nrow = 0, ncol = 2))
  colnames(average_MSE_table) = c("Model", "Test MSE")
  
  for (i in model_fit_labels) {
    average_MSE_row = mean((sim_results$sim_MSE_table %>% filter(Model == i))[, -1])
    average_MSE_table = rbind(average_MSE_table, c(i, average_MSE_row))
  }
  colnames(average_MSE_table) = c("Model", "Test MSE")
  average_MSE_table[, -1] = as.numeric(average_MSE_table[, -1])
  options(digits = 5)
  
  return(average_MSE_table)
}




round_df = function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


formatEstimatedCoefs = function(coefficient_estimates, all_coefficient_names) {
  formated_estimates = c()
  j = 1
  
  for (l in c(1:length(all_coefficient_names))) {
    if(!(is.na(names(coefficient_estimates)[j])) & all_coefficient_names[l] == names(coefficient_estimates)[j]) {
      formated_estimates = c(formated_estimates, coefficient_estimates[j])
      j = j + 1
    } else {
      formated_estimates = c(formated_estimates, 0)
    }
  }
  names(formated_estimates) = all_coefficient_names
  return(formated_estimates)
}


#' ###########################################################
#' Plotting function for relaxed LASSO (MSEs calculated on each fold using unparameterized fit from each set of selected predictors for each lambda)
plot_myRelaxed = function(relaxed_object) {
  nfold <- ncol(relaxed_object$all_folds_MSEs)
  mean_MSE <- rowMeans(relaxed_object$all_folds_MSEs)
  row_stdev <- apply(relaxed_object$all_folds_MSEs, 1, function(x) {sd(x)/sqrt(nfold)})
  plot.dat <- data.frame(lambdagrid = relaxed_object$all_lambdas, mean.MSE = mean_MSE, 
                         upper = mean_MSE +row_stdev, lower = mean_MSE - row_stdev)
  p <- ggplot(plot.dat, aes(log(lambdagrid), mean.MSE)) + geom_point(colour = "red", size = 3) + theme_bw() + 
    geom_errorbar(aes(ymin= lower, ymax=upper), width=.2, position=position_dodge(0.05), colour = "grey") + 
    labs(x = "log(lambda)", y = "Mean-Squared Error")  + 
    geom_vline(xintercept = log(relaxed_object$min_lambda), linetype = "dotted", colour = "blue")  
  return(p)
}



################################################# Multinomial LR helper functions ##################################################


###########################################################
#' relaxed LASSO function for mtool 
multinom.relaxed_enet <- function(X, Y, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
                                  constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                  ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20) {
  # FOR TESTING
  # regularization = 'elastic-net'; lambda_max = NULL; alpha = 0.7; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .0001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; i = 1; l = 1;
  print("RELAXED")
  p = ncol(train) - 2
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  print(lambdagrid)


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
                                                       lambda = 0.000001, alpha = alpha, unpen_cov = 1)

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



###########################################################
#' relaxed LASSO function for mtool 
multinom.relaxed_enet_nnet <- function(train, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
                                  constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                  ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20) {
  # FOR TESTING
  # regularization = 'elastic-net'; lambda_max = NULL; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .0001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; i = 1; l = 1; lambda_max = 0.2
  print("RELAXED")
  p = ncol(train)
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  print(lambdagrid)
  
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
    
    train_cv <- train[which(folds != i), ] #Set the training set
    test_cv <- train[which(folds == i), ] #Set the validation set
    # Create X and Y
    # Standardize
    train_cv <- scale(train_cv, center = T, scale = T)
    test_cv <- scale(test_cv, center = T, scale = T)
    
    #Segment your data by fold using the which() function 
    X_train <- as.matrix(X[which(folds != i), ]) #Set the training set
    X_val <- as.matrix(X[which(folds == i), ]) #Set the validation set
    Y_train <- as.numeric(Y[which(folds != i)]) - 1 #Set the training set
    Y_val <- as.numeric(Y[which(folds == i)]) - 1 #Set the validation set
    
    # Standardize
    X_train <- scale(X_train, center = T, scale = T)
    X_val <- scale(X_val, center = T, scale = T)
    # #Segment your data by fold using the which() function 
    # X_train <- as.matrix(X[which(folds != i), ]) #Set the training set
    # X_val <- as.matrix(X[which(folds == i), ]) #Set the validation set
    # Y_train <- as.numeric(Y[which(folds != i)]) - 1 #Set the training set
    # Y_val <- as.numeric(Y[which(folds == i)]) - 1 #Set the validation set
    
    # # Standardize
    # X_train <- scale(X_train, center = T, scale = T)
    # X_val <- scale(X_val, center = T, scale = T)
    
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    

    
    all_deviances[, i] <- unlist(
      mclapply(lambdagrid,
                function(lambda_v) {
                  dev = 0
                  # lambda1 <- lambda_v*alpha
                  # lambda2 <- 0.5*lambda_v*(1 - alpha)
                  # mtool
                  fit.mtool.parameterized <- fit_cbmodel(train_cv, regularization = 'elastic-net',
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
                    dev = multi_deviance_fsh(cb_data = test_cv, coefs = coef_subset, add_intercept = FALSE)
                    
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
                    
                    dev = multi_deviance_fsh(cb_data = test_cv, coefs = coef_subset, add_intercept = TRUE)
                  }
                  
                  dev
                },
                mc.cores = ncores, mc.set.seed = seed))
    cat("Completed Fold", i, "\n")
  }
  print(all_deviances)
  
  mean_dev <- rowMeans(all_deviances)
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, lambdagrid = lambdagrid, 
              cv_se = cv_se, deviance_grid = all_deviances))
}



###########################################################
#' relaxed LASSO function for mtool 
multinom.relaxed_enet_coefs_each_lambda <- function(X, Y, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
                                                    constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                                    ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20) {
  # FOR TESTING
  # regularization = 'elastic-net'; lambda_max = NULL; alpha = 0.7; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .0001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; i = 1; l = 1;
  print("RELAXED")
  p = ncol(train) - 2
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  print(lambdagrid)
  
  
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
    
    
    res <- unlist(mclapply(lambdagrid,
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
                                                             lambda = 0.000001, alpha = alpha, unpen_cov = 1)
                             
                             
                             if(i == nfold) {
                               fit.all.data = fit_cbmodel(cb_data_all, regularization = 'elastic-net',
                                                          lambda = lambda_v, alpha = alpha, unpen_cov = 1)
                               non_zero_coefs = length(union(which(fit.all.data$coefficients[, 1] != 0), 
                                                             which(fit.all.data$coefficients[, 2] != 0)))
                               list(length(non_zero_coefs), 
                                    multi_deviance_fsh(covs = X_val_new, response = Y_val, coefs = fit.mtool.subset$coefficients), add_intercept = TRUE)
                               
                             } else {
                               multi_deviance_fsh(covs = X_val_new, response = Y_val, coefs = fit.mtool.subset$coefficients, add_intercept = TRUE)
                             }
                           },
                           mc.cores = ncores, mc.set.seed = seed))
    if(i == nfold) {
      all_deviances[, i] = res[seq(2, 200, by = 2)]
    } else {
      all_deviances[, i] = res
    }
    cat("Completed Fold", i, "\n")
  }
  all_non_coefs[, nfold] = res[seq(1, 200, by = 2)]
  
  
  mean_dev <- rowMeans(all_deviances)
  lambda.min <- lambdagrid[which.min(mean_dev)]
  
  
  cv_se <- sqrt(var(mean_dev))
  rownames(all_deviances) <- lambdagrid
  return(list(lambda.min = lambda.min, non_zero_coefs = non_zero_coefs, lambdagrid = lambdagrid, 
              cv_se = cv_se, deviance_grid = all_deviances, selected_coefs_grid = all_non_coefs))
}




###########################################################
#' Calculate multinomial deviance on fitSmoothHazard object
multi_deviance_fsh <- function(cb_data, coefs, add_intercept) {
  if(add_intercept) {
    X <- as.matrix(cbind(cb_data$covariates, 1))
  } else {
    X = cb_data$covariates
  }
  
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
#' Runs simulation of hazard function fits

runCasebaseSim = function(n = 400, p = 20, N = 5, nfolds = 5, seed) {
  # n = 400
  # p = 20
  # N = 5
  # nfolds = 5
  # seed = 2023
  
  # tryCatch({
    # For each simulation you want to output the MSE
    Results <- replicate(N, {
      # Set seed
      # take the last five digits of the initial seed
      # seed = as.integer(Sys.time())
      # the_seed= seed %% 100000
      # set.seed(the_seed)
      
      set.seed(seed)
      
      num_true <- p/2
      beta1 <- c(rep(0, p))
      beta2 <- c(rep(0, p))
      nu_ind <- seq(num_true)
      # Here out of 20 predictors, 10 should be non-zero 
      beta1[nu_ind] <- c(rep(1, p/2), rep(0, p/2))
      beta2[nu_ind] <- c(rep(-1, p/2), rep(0, p/2))
      
      # Simulate data
      sim.data <- cause_hazards_sim(n = n, p = p, nblocks = 4, 
                                    beta1 = beta1, beta2 = beta2, rate_cens = 0.25, 
                                    h1 = 0.55, h2 = 0.10, gamma1 = 1.5, gamma2 = 1.5)
      
      
      # Censoring proportion
      cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)
      
      # Training-test split 
      # We only do this (instead of generating datasets for train and test like Anthony mentioned because it is faster computationally 
      # as casebase resamples) + proportion of censoring can be quite random in each run of the simulation so we want to maintain the same in validation and test set
      train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
      train <- sim.data[train.index,]
      test <- sim.data[-train.index,]

      
      ##############################################################
      # We have two competitor models for variable selection:
      # 1) Independent cox-regression model 
      # 2) penCR cox regression model - where the lambda penalties are trained together 
      ######################## Fit indepedent cox-regression model ###############################
      ######################### Cause-1 #########################################
      # Censor competing event
      y_train <- Surv(time = train$ftime, event = train$fstatus == 1)
      
      x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
      
      # Censor competing event
      y_test <- Surv(time = test$ftime, event = test$fstatus == 1)
      
      x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
      
      # Fit cause-specific cox model with glmnet on training set 
      cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7, folds = nfolds)
      
      # Fit on validation set 
      cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                            lambda = cox_mod$lambda.min)
      
      cc_min <- coef(cox_val_min)
      
      res_cox_min1 <- varsel_perc(cc_min, beta1)
      
      # let's calculate the bias for all the competitors as well (a task could be turning this one line into a function as well)
      # Only for the true non-zero variables
      mean((cc_min[nu_ind] - beta1[nu_ind])^2)
      ########################## Cause 2 #####################################
      # Censor competing event
      y_train <- Surv(time = train$ftime, event = train$fstatus == 2)
      
      x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
      
      # Censor competing event
      y_test <- Surv(time = test$ftime, event = test$fstatus == 2)
      
      x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
      
      # Fit cause-specific cox model with glmnet on training set
      cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7, folds = nfolds)
      
      # Fit on validation set 
      cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                            lambda = cox_mod$lambda.min)
      
      cc_min <- coef(cox_val_min)
      
      # Function to output variable selection performance metrics
      res_cox_min2 <- varsel_perc(cc_min, beta2)
      
      # let's calculate the bias for all the competitors as well (a task could be turning this one line into a function as well)
      cc_min_bias <- which(coef(cox_val_min) != 0)
      
      # MSE (bias)
      mean((cc_min[nu_ind] - beta2[nu_ind])^2)
      ########################## Fit PenCR model ##################################
      penCR = cv.glmnet.CR(data = train, family="cox", alpha= 0.7, standardize= TRUE,
                           nlambda = 20, t.BS = median(train$ftime), seed = seed, causeOfInt = 1,
                           nfold = nfolds)
      
      cc_min_penCR1 <- penCR$glmnet.fits$models$`Cause 1`$glmnet.res$lambda[penCR$min.index[1]]
      cc_min_penCR2 <- penCR$glmnet.fits$models$`Cause 2`$glmnet.res$lambda[penCR$min.index[2]]
      
      # Fit on validation set 
      penCR_val_min1 <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                               lambda = cc_min_penCR1)
      
      cc_min_penCR1 <- coef(penCR_val_min1)
      
      penCR_val_min2 <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                               lambda = cc_min_penCR2)
      
      
      cc_min_penCR2 <- coef(penCR_val_min2)
      
      res_pencr_min1 <- varsel_perc(cc_min_penCR1, beta1)
      res_pencr_min2 <- varsel_perc(cc_min_penCR2, beta2)
      
      # Free up memory for case-base
      rm(penCR)
      
      # Calculate MSE here as well (try and fill it out!)
      # MSE Cox model for cause of interest
      mean((cc_min_penCR1[nu_ind] - beta1[nu_ind])^2)
      # MSE Cox model for competing risk
      cc_min_bias_pen <- which(cc_min_penCR2 != 0)
      mean((cc_min_penCR2[nu_ind] - beta2[nu_ind])^2)
      
      ########################## Fit casebase model #################################
      # Test set 
      surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
      
      # Covariance matrix
      cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
      
      # Case-base dataset
      cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)
      
      # Train case-base model
      start_cv <- Sys.time()
      cv.lambda <- mtool.multinom.cv(train, seed = seed, nfold = nfolds, lambda_max =  0.9, alpha = 0.7)
      end_cv <- Sys.time()
      
      
      # Case-base fits 
      # Lambda.min
      fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                                 lambda = cv.lambda$lambda.min , alpha = 0.7, unpen_cov = 2)
      
      
      res_cb_min1 <- varsel_perc(fit_val_min$coefficients[1:p, 1], beta1)
      
      res_cb_min2 <- varsel_perc(fit_val_min$coefficients[1:p, 2], beta2)
      
      # Calculate MSE here as well!
      # MSE for casebase model for cause of interest
      fit_val_coef_1 <- fit_val_min$coefficients[1:p, 1]
      mean((beta1[nu_ind] - fit_val_coef_1[nu_ind])^2)
      
      # MSE for casebase model for competing risk
      fit_val_coef_2 <- fit_val_min$coefficients[1:p, 2]
      cb_min_bias <- which(fit_val_coef_2 != 0)
      mean((fit_val_coef_2[nu_ind] - beta2[nu_ind])^2)
      
      
      ########################## Fit casebase model with post LASSO#################################
      
      start_post <- Sys.time()
      res <- multinom.post_enet_old(train, test, nfold = nfolds, seed = seed, lambda_max =  0.9, alpha = 0.7)
      end_post <- Sys.time()
      
      # Calculate MSE for this as well (fill in here)
      mean((res$coefficients[nu_ind, 1]- beta1[nu_ind])^2)
      mean((res$coefficients[nu_ind, 2] - beta2[nu_ind])^2)
      
      res_cb_post_lasso1 <- varsel_perc(res$coefficients[nu_ind, 1], beta1)
      res_cb_post_lasso2 <- varsel_perc(res$coefficients[nu_ind, 2], beta2)
      
      
      ########################## Fit casebase model with relaxed LASSO#################################
      
      start_relaxed <- Sys.time()
      res <- multinom.relaxed_enet(train, nfold = nfolds, seed = seed, lambda_max =  0.9, alpha = 0.7)
      end_relaxed <- Sys.time()
      print(res$lambda.min)
      
      model_final <- fitSmoothHazard(fstatus ~. +log(ftime) -fstatus,
                                      data = test,
                                      time = "ftime",
                                      lambda = res$lambda.min)
      
      all_coef_names_cause1 =  paste("X", seq(1:p), ":1", sep ="")
      all_coef_names_cause2 =  paste("X", seq(1:p), ":2", sep ="")
      
      exclude_coefs = c("(Intercept):1", "ftime:1",
                        "log(ftime):1", "(Intercept):2", "ftime:2",
                        "log(ftime):2")
      
      est_betas = coef(model_final)[!(names(coef(model_final)) %in% exclude_coefs)]
      est_betas_cause1 = est_betas[names(est_betas) %in% all_coef_names_cause1]
      est_betas_cause2 = est_betas[names(est_betas) %in% all_coef_names_cause2]
      
      
      # Calculate MSE for this as well (fill in here)
      mean((est_betas_cause1 - beta1[nu_ind])^2)
      mean((est_betas_cause2 - beta2[nu_ind])^2)
      
      res_cb_relaxed_lasso1 <- varsel_perc(est_betas_cause1, beta1)
      res_cb_relaxed_lasso2 <- varsel_perc(est_betas_cause2, beta2)
      
      ###################################################################################
      
      print(paste(round(end_cv - start_cv, 3), " taken to run multinomial cv"))
      print(paste(round(end_post - start_post, 3), " taken to run post LASSO"))
      print(paste(round(end_relaxed - start_relaxed, 3), " taken to run relaxed LASSO"))
      
      Res <- rbind(res_cb_relaxed_lasso1, res_cb_relaxed_lasso2, res_cb_post_lasso1, 
                   res_cb_post_lasso2, res_cb_min1, res_cb_min2, 
                   res_cox_min1, res_cox_min2, res_pencr_min1, res_pencr_min2, cen.prop)
      
      rownames(Res) <- c("casebase.relaxed.lasso.lambda.min_cause1", "casebase.relaxed.lasso.lambda.min_cause2", 
                         "casebase.post.lasso.lambda.min_cause1", "casebase.post.lasso.lambda.min_cause2", 
                         "casebase.lambda.min_cause1", 
                         "casebase.lambda.min_cause2", "cox.lambda.min_cause1",
                         "cox.lambda.min_cause2", "pencr.lambda.mincause1", "pencr.lambda.mincause2", "cens.prop")
      
      Res
      
    }, simplify = FALSE)
    
  Results <- do.call(rbind, Results)
  Results
  # },
  # 
  # error = function(e) {
  #   print(e)
  # })
}



plot_post_relaxed.multinom <- function(cv_object) {
  nfold <- ncol(cv_object$deviance_grid)
  mean_dev <- rowMeans(cv_object$deviance_grid)
  row_stdev <- apply(cv_object$deviance_grid, 1, function(x) {sd(x)/sqrt(nfold)})
  plot.dat.p <- data.frame(lambdagrid = cv_object$lambdagrid, mean.dev = mean_dev, 
                         upper = mean_dev +row_stdev, lower = mean_dev - row_stdev)
  p <- ggplot(plot.dat.p, aes(log(lambdagrid), mean.dev)) + geom_point(colour = "red", size = 3) + theme_bw() + 
    geom_errorbar(aes(ymin= lower, ymax=upper), width=.2, position=position_dodge(0.05), colour = "grey") + 
    labs(x = "log(lambda)", y = "Multinomial Deviance")  + 
    geom_vline(xintercept = log(cv_object$lambda.min), linetype = "dotted", colour = "blue")
  
  # plot.dat.g <- data.frame(no_coefs = cv_object$lambdagrid, mean.dev = mean_dev, 
  #                          upper = mean_dev +row_stdev, lower = mean_dev - row_stdev)
  # 
  # g <-
  # 
  return(p)
}




###########################################################
#' Formats resulting table from simulation of hazard function fits

formatCaseBaseTable = function(sim_results) {
  model_fit_labels = c("casebase.relaxed.lasso.lambda.min_cause1", "casebase.relaxed.lasso.lambda.min_cause2", 
                        "casebase.post.lasso.lambda.min_cause1", "casebase.post.lasso.lambda.min_cause2", 
                        "casebase.lambda.min_cause1", 
                        "casebase.lambda.min_cause2", "cox.lambda.min_cause1",
                        "cox.lambda.min_cause2", "pencr.lambda.mincause1", "pencr.lambda.mincause2", "cens.prop")
  
  # stat_names = colnames(as.data.frame(sim_results))
  
  num_models = length(model_fit_labels)
  rows = nrow(sim_results)
  # Add model name to model column
  model_list = c()
  for (i in c(1:rows)) {
    index = i %% num_models
    if (i %% num_models == 0)
      index = num_models
    curr_model = model_fit_labels[index]
    model_list = c(model_list, curr_model)
  }
  
  sim_results["Model"] = model_list
  
  stat_names = c("Model", "Sensitivity",	"Specificity",	"MCC",	"Coefficient_Bias")
  
  formatted_table = data.frame(matrix(nrow = 0, ncol = length(stat_names)))
  colnames(formatted_table) = stat_names
  stat_table = sim_results[stat_names]
  
  for (i in model_fit_labels) {
    average_stat_row = colMeans((stat_table %>% filter(Model == i))[2:ncol(stat_table)])
    formatted_table = rbind(formatted_table, c(i, average_stat_row))
  }
  colnames(formatted_table) = stat_names
  formatted_table = formatted_table %>% replace(is.na(.), NaN)

  options(digits = 5)
  
  return(formatted_table)
}


