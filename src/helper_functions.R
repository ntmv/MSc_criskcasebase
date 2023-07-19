library(tictoc)
library(caret)

#' Fits linear regression using relaxed LASSO regularization to given data
#' @param: train_data - dataframe containing values of features
#' @param: response - list of response values
#' @param: folds - number of folds used to compute cross validation MSE (default is 5)
#' @param: print_time - boolean indicating whether to print time function takes to run
myRelaxedGlmnetFolds = function(train_data, response, folds = 5, print_time) {
  tryCatch({ 
    if(print_time)
      tic()
  
    coefficient_names = colnames(as.data.frame(train_data))
    response_name = names(as.data.frame(response))
    
    # Data frame and list to keep track of coefficients and CV-MSEs of all lambdas
    coefs_all_lambdas = data.frame(matrix(nrow = 1, ncol = length(coefficient_names) + 1))
    colnames(coefs_all_lambdas) = c("(Intercept)", coefficient_names)
    MSEs_all_lambda = c()
    
    set.seed(2023)
    fit_lasso = cv.glmnet(train_data, response, family="gaussian", keep=TRUE, nfolds = folds, alpha=1)
    all_coef_fit_lasso = as.data.frame.matrix(fit_lasso$glmnet.fit$beta)

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
      
      fold_ids = fit_lasso$foldid
      mse_values = c()
      
      for(m in c(1:folds)) {
        train_indices = which(fold_ids != m)

        # Check if new training data includes one or more covariates to make sure it is subsetted as a matrix or a list correctly
        if (is.matrix(new_x_train)) {
          traindata = new_x_train[train_indices, ]
          testdata = new_x_train[-train_indices, ]
        } else {
          traindata = new_x_train[train_indices]
          testdata = new_x_train[-train_indices]
        }
        ytrain = response[train_indices]
        ytest = response[-train_indices]
        
        # Check what type to treat validation data as
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



#' Fits linear regression using relaxed LASSO regularization to given data
#' @param: train_data - dataframe containing values of features
#' @param: response - list of response values
#' @param: folds - number of folds used to compute cross validation MSE (default is 5)
#' @param: print_time - boolean indicating whether to print time function takes to run
myRelaxedCustomFolds = function(train_data, response, cv, folds = 5, print_time) {
  tryCatch({ 
    if(print_time)
      tic()
    
    # Get lambda grid from glmnet fit on entire training set
    
    
    # TODO: replace lm() with: solve(t(X) %*% X) %*% t(X) %*% y
    
    #Perform cross-validation
    set.seed(2023)
    folds_list <- createFolds(response, k = folds, list = TRUE, returnTrain = TRUE)
    res = list()
    mse_values = c()
    for(m in 1:folds) {
      train_indices <- folds_list[[m]]

      if (is.matrix(new_x_train)) {
        traindata = new_x_train[train_indices, ]
        testdata = new_x_train[-train_indices, ]
      } else {
        traindata = new_x_train[train_indices]
        testdata = new_x_train[-train_indices]
      }
      ytrain = response[train_indices]
      ytest = response[-train_indices]
      
      fit_lasso = glmnet(traindata, ytrain, family="gaussian", alpha=1)
      all_coef_fit_lasso = as.data.frame.matrix(coef(fit_lasso))[-1, ]
      
      
      
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
    
    coefficient_names = colnames(as.data.frame(train_data))
    response_name = names(as.data.frame(response))
    
    # Data frame and list to keep track of coefficients and CV-MSEs of all lambdas
    coefs_all_lambdas = data.frame(matrix(nrow = 1, ncol = length(coefficient_names) + 1))
    colnames(coefs_all_lambdas) = c("(Intercept)", coefficient_names)
    MSEs_all_lambda = c()
    
    if (cv) {
      set.seed(2023)
      fit_lasso = cv.glmnet(train_data, response, family="gaussian", keep=TRUE, alpha=1)
      all_coef_fit_lasso = as.data.frame.matrix(fit_lasso$glmnet.fit$beta)
    }
    else {
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
    
    # Iterate over all lambdas
    indices = c(1:length(all_lambdas))

    for(i in indices) {
      rows = nrow(coefs_all_lambdas)
      # TODO: Look into using predict() with type = "nonzero" as a way to get nonzero coefficients from a glmnet fit for particular lambda values
      # (might be more efficent that using which and subsetting the whole training set, then refitting with lm())
      
      # TODO: replace lm() with: solve(t(X) %*% X) %*% t(X) %*% y
      
      # Perform cross-validation
      # TODO: Check if there's a difference between creating folds outside lambda loop or inside lambda loop
      set.seed(2023)
      folds_list <- createFolds(response, k = folds, list = TRUE, returnTrain = TRUE)
      res = list()
      mse_values = c()
      for(m in 1:folds) {
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
        
        # TODO: Check if setting seed v.s. not setting seed here gives different results
        fit_en = glmnet(k_minus_1_train_data, y_train, family="gaussian", alpha=1, lambda = all_lambdas[i])
        
        
        # TODO: Replace coef(fit_en) with list of coefficients found using lambda min
        current_coef = coef(fit_en)[-1, ]
        non_zero_coef = current_coef[current_coef != 0]
        all_zeros = length(non_zero_coef) == 0
        if(m == 1) {
          print(paste("current lambda index: ", as.character(i)))
          print(current_coef)
        }

        
        # If no predictors selected, skip OLS with cross validation and continue with next lambda
        if (all_zeros){
          coefs_all_lambdas[rows,] = 0
          mse_values[m] = mean((y_valid)^2)
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
      }
      
      print(mse_values)
      
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
  
  result = list(coefficients = coefs_all_lambdas, CV_MSEs = MSEs_all_lambda, all_lambdas = all_lambdas, 
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
      set.seed(2023)
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
  mse_relaxed_min <- sum((pred_fit_relaxed_min - test_data$y_test)^2)
  mse_myRelaxed <- sum((pred_myRelaxed - test_data$y_test)^2)
  mse_myRelaxed_cv <- sum((pred_myRelaxed_cv - test_data$y_test)^2)
  mse_myRelaxed_new_cv <- sum((pred_myRelaxed_new_cv - test_data$y_test)^2)
  mse_post_lasso_min <- sum((pred_fit_post_lasso_min - test_data$y_test)^2)

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

