install.packages( "/Users/alex/Downloads/lmvar_1.5.2.tar.gz", repos=NULL, type="source")
library(lmvar)
library(tictoc)
library(caret)


myRelaxed = function(train_data, response, cv) {
  tic()
  coefficient_names = colnames(as.data.frame(train_data))
  response_name = names(as.data.frame(response))
  
  coefs_all_lambdas = data.frame(matrix(nrow = 1, ncol = length(coefficient_names) + 1))
  colnames(coefs_all_lambdas) = c("(Intercept)", coefficient_names)
  MSEs_all_lambda = c()
  
  if (cv) {
    fit_lasso = cv.glmnet(train_data, response, family="gaussian", alpha=1)
    all_coef_fit_lasso = as.data.frame.matrix(fit_lasso$glmnet.fit$beta)
  }
  else {
    fit_lasso = glmnet(train_data, response, family="gaussian", alpha=1)
    all_coef_fit_lasso = as.data.frame.matrix(coef(fit_lasso))[-1, ]
  }

  current_MSE = .Machine$double.xmax;
  best_fit = lm(1~1)
  best_lambda_index = 0

  
  indeces = c(1:length(all_coef_fit_lasso))
  # indeces = c(30:33)
  
  for(i in indeces) {
    rows = nrow(coefs_all_lambdas)
    current_coef = all_coef_fit_lasso[i]
    non_zero_coef_lasso_then_OLS_all_predictors = current_coef[all_coef_fit_lasso[i] != 0]
    all_zeros = length(non_zero_coef_lasso_then_OLS_all_predictors) == 0
    
    if (all_zeros){
      coefs_all_lambdas[rows,] = 0
      next
    }

    non_zero_coef_lasso_then_OLS_all_predictors_indeces = which(current_coef != 0)
    selected_coefficient_names = coefficient_names[non_zero_coef_lasso_then_OLS_all_predictors_indeces]
    new_x_train = train_data[, non_zero_coef_lasso_then_OLS_all_predictors_indeces]
    # data = data.frame(new_x_train, response)
    # colnames(data) = c(selected_coefficient_names, response_name)
    
    # Subset X and y for cross validation
    # TODO: replace lm() with: solve(t(X) %*% X) %*% t(X) %*% y
    # fit_OLS_on_LASSO_subset = lm(response~., data = data, x = TRUE, y = TRUE)
    # fit_OLS_cv_on_LASSO_subset = cv.lm(fit_OLS_on_LASSO_subset, k = 5)
    # MSE = fit_OLS_cv_on_LASSO_subset$MSE$mean
    
    k = 5
    folds <- createFolds(response, k = k, list = TRUE, returnTrain = TRUE)
    res = list()
    mse_values = c()

    for(m in 1:k) {
      train_indices <- folds[[m]]
      colnames(new_x_train) = c(selected_coefficient_names)
      traindata <- new_x_train[train_indices, ]
      testdata  <- new_x_train[-train_indices, ]
      ytrain = response[train_indices]
      ytest = response[-train_indices]

      data = data.frame(traindata, ytrain)
      colnames(data) = c(selected_coefficient_names, response_name)

      fit_OLS_on_LASSO_subset = lm(response~., data = data, x = TRUE, y = TRUE)
      # print(colnames(data))
      # cat("/n")
      # print(fit_OLS_on_LASSO_subset$coefficients)
      # cat("/n")

      pred = predict(fit_OLS_on_LASSO_subset, newdata = as.data.frame.matrix(testdata))

      mse_values[m] = mean((ytest - pred)^2)
    }

    MSE = mean(mse_values)
    cat(paste("Num rows: ", rows))
    writeLines("")
    cat(paste("i: ", i))
    writeLines("")
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
    # if(all_zeros) {
    #   cat(coefs_all_lambdas[rows+1,])
    # }
    # 
    MSEs_all_lambda = rbind(MSEs_all_lambda, MSE)

    if(MSE < current_MSE) {
      current_MSE = MSE
      best_fit = fit_OLS_on_LASSO_subset
      best_lambda_index = i
    }
  }
  
  # rownames(coefs_all_lambdas) = paste("lambda: ", fit_lasso$lambda)
  cat(paste("\nbest fit is for lambda at index: ", best_lambda_index, "\n"))
  cat("\nEND OF FITTING\n")
  toc()
  
  result = list(coefficients = coefs_all_lambdas, CV_MSEs = MSEs_all_lambda, all_lambdas = fit_lasso$lambda, 
                min_lambda = fit_lasso$lambda[best_lambda_index], min_lambda_index = best_lambda_index, 
                best_fit = best_fit)
  return (result)
}

# cross_validation <- function(new_x_train, response, seed_val, model, k = 10) {
#   set.seed(2023)
#   folds <- createFolds(response, k = k, list = TRUE, returnTrain = TRUE)
#   response_name = names(as.data.frame(response))
#     
#   res = list()
# 
#   for(i in 1:k) {
#     train_indices <- folds[[i]]
#     train_data <- new_x_train[train_indices, ]
#     test_data  <- new_x_train[-train_indices, ]
#     y_train = response[train_indices]
#     y_test = response[-train_indices]
#     
#     data = data.frame(train_data, y_train)
#     colnames(data) = c(selected_coefficient_names, response_name)
#     
#     fit_OLS_on_LASSO_subset = lm(response~., data = data, x = TRUE, y = TRUE)
#     
#     
#     pred <- predict(fit_OLS_on_LASSO_subset, newdata = test_data)
#     
#     mse_values[i] <- mean((y_test - pred)^2)
#   }
#   
#   return(mean(mse_values))
# } 


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