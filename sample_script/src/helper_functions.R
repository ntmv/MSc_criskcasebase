install.packages( "/Users/alex/Downloads/lmvar_1.5.2.tar.gz", repos=NULL, type="source")
library(lmvar)
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
    indeces = c(1:length(all_coef_fit_lasso))
  
    for(i in indeces) {
      rows = nrow(coefs_all_lambdas)
      current_coef = all_coef_fit_lasso[i]
      non_zero_coef = current_coef[all_coef_fit_lasso[i] != 0]
      all_zeros = length(non_zero_coef) == 0
      
      # If no predictors selected, skip OLS with cross validation and continue with next lambda
      if (all_zeros){
        coefs_all_lambdas[rows,] = 0
        next
      }
  
      # Create new training dataframe to pass to lm() for OLS fit on selected predictors
      non_zero_coef_indeces = which(current_coef != 0)
      non_zero_coef_names = coefficient_names[non_zero_coef_indeces]
      new_x_train = train_data[, non_zero_coef_indeces]
      # print(head(new_x_train))
      # writeLines("")
      # print(non_zero_coef_names)
      # writeLines("")
      # print(i)
      # writeLines("")
      
      # Subset X and y for cross validation
      # TODO: replace lm() with: solve(t(X) %*% X) %*% t(X) %*% y
      
      #Perform cross-validation
      k = 5
      folds <- createFolds(response, k = k, list = TRUE, returnTrain = TRUE)
      res = list()
      mse_values = c()
      for(m in 1:k) {
        train_indices <- folds[[m]]
        # TODO: Fix code so column names
        print(head())
        if(is.numeric(new_x_train)) {
          names(new_x_train) = c(non_zero_coef_names)
        } else {
          colnames(new_x_train) = c(non_zero_coef_names)
        }
        traindata <- new_x_train[train_indices, ]
        testdata  <- new_x_train[-train_indices, ]
        ytrain = response[train_indices]
        ytest = response[-train_indices]
  
        data = data.frame(traindata, ytrain)
        colnames(data) = c(non_zero_coef_names, response_name)
  
        fit_OLS_on_LASSO_subset = lm(response~., data = data, x = TRUE, y = TRUE)
        pred = predict(fit_OLS_on_LASSO_subset, newdata = as.data.frame.matrix(testdata))
  
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
    print(class(new_x_train))
    print(head(new_x_train))
    print(paste("lambda index: ", i))
    # print(paste("length new_x_train: ", length(new_x_train)))
    # print(non_zero_coef_names)
  }
  # ,
  # warning = function(w) {
  #   print(head(new_x_train))
  #   print(non_zero_coef_names)
  #   return(NA)
  # }
  )
  if(print_time)
    toc()
  
  result = list(coefficients = coefs_all_lambdas, CV_MSEs = MSEs_all_lambda, all_lambdas = fit_lasso$lambda, 
                min_lambda = fit_lasso$lambda[best_lambda_index], min_lambda_index = best_lambda_index, 
                best_fit = best_fit)
  return (result)
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