library(tictoc)
library(caret)

#' Fits linear regression using relaxed LASSO regularization to given data
#' @param: train_data - dataframe containing values of features
#' @param: response - list of response values
#' @param: cv - boolean indicating whether to fit lasso with cross validation or not
#' @param: print_time - boolean indicating whether to print time function takes to run
myRelaxed = function(train_data, response, cv, print_time) {
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
    }
    else {
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
      k = 5
      folds <- createFolds(response, k = k, list = TRUE, returnTrain = TRUE)
      res = list()
      mse_values = c()
      for(m in 1:k) {
        train_indices <- folds[[m]]
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


generateDataset = function(p, n) {
  beta = c()
  for (i in seq(1:(p+1))) {
    if (i == 1) {
      beta = c(beta, 1)
    }
    else if (i %% 2 == 0) {
      beta = c(beta, 0)
    } else {
      beta = c(beta, 0.5)
    }
  }
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


fitAll = function(train_data, response, print_time) {
  # cv.glmnet with LASSO, cv.glmnet relaxed LASSO, and my relaxed implementation fits
  fit_Y1 = cv.glmnet(train_data, response, gamma = 0, relax = FALSE)
  relaxed_fit_Y1 = cv.glmnet(train_data, response, gamma = 0, relax = TRUE)
  my_relaxed_fit_Y1 = myRelaxed(train_data, response, FALSE, print_time)
  my_relaxed_fit_cv_Y1 = myRelaxed(train_data, response, TRUE, print_time)
  
  # Save lambda values of interest from cv.glmnet LASSO to use for post LASSO fit
  lambda_min = fit_Y1$lambda.min
  lambda_1se = fit_Y1$lambda.1se
  
  # Post LASSO fits
  fit_Y1_post_lasso_lambda_min = glmnet(train_data, response, lambda = lambda_min)
  fit_Y1_post_lasso_lambda_1se = glmnet(train_data, response, lambda = lambda_1se)
  
  result = list(cv_glmnet = fit_Y1, cv_glmnet_relaxed = relaxed_fit_Y1, my_relaxed = my_relaxed_fit_Y1,
                my_relaxed_cv = my_relaxed_fit_cv_Y1, post_lasso_lambda_min = fit_Y1_post_lasso_lambda_min,
                post_lasso_lambda_1se = fit_Y1_post_lasso_lambda_1se)
  
  return(result)
}


computeMSEs = function(resulting_fits, test_data) {
  
  pred_fit_min = predict(resulting_fits$cv_glmnet, newx = test_data$x_test, s = resulting_fits$cv_glmnet$lambda.min)
  pred_fit_1se = predict(resulting_fits$cv_glmnet, newx = test_data$x_test, s = resulting_fits$cv_glmnet$lambda.1se)
  
  pred_fit_relaxed_min = predict(resulting_fits$cv_glmnet_relaxed, newx = test_data$x_test, s = resulting_fits$cv_glmnet_relaxed$relaxed$lambda.min)
  pred_fit_relaxed_1se = predict(resulting_fits$cv_glmnet_relaxed, newx = test_data$x_test, s = resulting_fits$cv_glmnet_relaxed$relaxed$lambda.1se)
  
  pred_myRelaxed = predict.lm(resulting_fits$my_relaxed$best_fit, newdata = as.data.frame(test_data$x_test))
  pred_myRelaxed_cv = predict.lm(resulting_fits$my_relaxed_cv$best_fit, newdata = as.data.frame(test_data$x_test))
  
  coefficient_names = rownames(resulting_fits$post_lasso_lambda_min$beta)
  indices = which(colnames(as.data.frame.matrix(test_data$x_test)) %in% coefficient_names)
  new_x_test = test_data$x_test[, indices]
  colnames(new_x_test) = coefficient_names
  pred_fit_post_lasso_min = predict(resulting_fits$post_lasso_lambda_min, newx = new_x_test)
  
  coefficient_names = rownames(resulting_fits$post_lasso_lambda_1se$beta)
  indices = which(colnames(as.data.frame.matrix(test_data$x_test)) %in% coefficient_names)
  new_x_test = test_data$x_test[, indices]
  colnames(new_x_test) = coefficient_names
  pred_fit_post_lasso_1se = predict(resulting_fits$post_lasso_lambda_1se, newx = new_x_test)
  
  
  
  # Compute MSE
  mse_min <- sum((pred_fit_min - test_data$y_test)^2)
  mse_1se <- sum((pred_fit_1se - test_data$y_test)^2)
  mse_relaxed_min <- sum((pred_fit_relaxed_min - test_data$y_test)^2)
  mse_relaxed_1se <- sum((pred_fit_relaxed_1se - test_data$y_test)^2)
  mse_myRelaxed <- sum((pred_myRelaxed - test_data$y_test)^2)
  mse_myRelaxed_cv <- sum((pred_myRelaxed_cv - test_data$y_test)^2)
  mse_post_lasso_min <- sum((pred_fit_post_lasso_min - test_data$y_test)^2)
  mse_post_lasso_1se <- sum((pred_fit_post_lasso_1se - test_data$y_test)^2)
  
  MSEs = data.frame(matrix(nrow = 8, ncol = 2))
  colnames(MSEs) = c("Model", "Test MSE")
  
  model_fit_labels = c("cv.glmnet relax = FALSE, LASSO, lambda min", "cv.glmnet relax = FALSE, LASSO, lambda 1se", 
                       "cv.glmnet relax = TRUE, lambda min", "cv.glmnet relax = TRUE, lambda 1se", 
                       "Relaxed LASSO implementation no CV", "Relaxed LASSO implementation CV", 
                       "Post LASSO lambda min", "Post LASSO lambda 1se")
  
  MSEs[, 1] = model_fit_labels
  MSEs[, 2] = c(mse_min, mse_1se, mse_relaxed_min, mse_relaxed_1se, mse_myRelaxed, mse_myRelaxed_cv, mse_post_lasso_min, mse_post_lasso_1se)
  
  return(MSEs)
}



formatCoefficientTable = function(p, n, coefficient_names, resulting_fits, print_time) {
  model_fit_labels = c("cv.glmnet relax = FALSE, LASSO, lambda min", "cv.glmnet relax = FALSE, LASSO, lambda 1se", 
                       "cv.glmnet relax = TRUE, lambda min", "cv.glmnet relax = TRUE, lambda 1se", 
                       "Relaxed LASSO implementation no CV", "Relaxed LASSO implementation CV", 
                       "Post LASSO lambda min", "Post LASSO lambda 1se")
  
  # Create dataframes to keep track of coefficient values and MSEs per each iteration of simulation
  coefficient_values = data.frame(matrix(nrow = 0, ncol = length(coefficient_names)))
  
  col_names = rownames(as.data.frame.matrix(coef(resulting_fits$cv_glmnet_relaxed)))
  coef_values = as.numeric(coef(resulting_fits$post_lasso_lambda_min))
  coef_names = rownames(coef(resulting_fits$post_lasso_lambda_min))
  names(coef_values) = coef_names
  y_post_lasso_lambda_min_coefficient_row = formatCoefficientTableRow(coef_values, col_names)
  
  coef_values = as.numeric(coef(resulting_fits$post_lasso_lambda_1se))
  coef_names = rownames(coef(resulting_fits$post_lasso_lambda_1se))
  names(coef_values) = coef_names
  y_post_lasso_lambda_1se_coefficient_row = formatCoefficientTableRow(coef_values, col_names)
  
  
  coefficient_values = rbind(coefficient_values, 
                             c(0, resulting_fits$cv_glmnet$glmnet.fit$beta[, resulting_fits$cv_glmnet$index[1]]),
                             c(0, coef(resulting_fits$cv_glmnet)[-1]),
                             c(0, resulting_fits$cv_glmnet_relaxed$glmnet.fit$relaxed$beta[, resulting_fits$cv_glmnet_relaxed$relaxed$index[1]]),
                             c(0, coef(resulting_fits$cv_glmnet_relaxed)[-1]),
                             c(0, resulting_fits$my_relaxed$coefficients[resulting_fits$my_relaxed$min_lambda_index, -1]),
                             c(0, resulting_fits$my_relaxed_cv$coefficients[resulting_fits$my_relaxed_cv$min_lambda_index, -1]),
                             c(0, y_post_lasso_lambda_min_coefficient_row),
                             c(0, y_post_lasso_lambda_1se_coefficient_row))
  
  options(digits = 5)
  colnames(coefficient_values) = coefficient_names
  coefficient_values[, 1] = model_fit_labels
  
  return(coefficient_values)
}


oneIteration = function(p, n, training_data, test_data, print_time) {
  tryCatch({
    # Fit training data with different models
    resulting_fits = fitAll(training_data$x_train, training_data$y_train, FALSE)
    
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



round_df = function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


formatCoefficientBiasTable = function(sim_results) {
  betas = c()
  p = 20
  for (i in seq(1:p)) {
    if (i %% 2 != 0) {
      betas = c(betas, 0)
    } else {
      betas = c(betas, 0.5)
    }
  }
  
  model_fit_labels = c("cv.glmnet relax = FALSE, LASSO, lambda min", "cv.glmnet relax = FALSE, LASSO, lambda 1se", 
                       "cv.glmnet relax = TRUE, lambda min", "cv.glmnet relax = TRUE, lambda 1se", 
                       "Relaxed LASSO implementation no CV", "Relaxed LASSO implementation CV", 
                       "Post LASSO lambda min", "Post LASSO lambda 1se")
  
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
  model_fit_labels = c("cv.glmnet relax = FALSE, LASSO, lambda min", "cv.glmnet relax = FALSE, LASSO, lambda 1se", 
                       "cv.glmnet relax = TRUE, lambda min", "cv.glmnet relax = TRUE, lambda 1se", 
                       "Relaxed LASSO implementation no CV", "Relaxed LASSO implementation CV", 
                       "Post LASSO lambda min", "Post LASSO lambda 1se")
  
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

