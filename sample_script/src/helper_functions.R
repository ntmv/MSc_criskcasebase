install.packages( "/Users/alex/Downloads/lmvar_1.5.2.tar.gz", repos=NULL, type="source")
library(lmvar)
library(tictoc)


myRelaxed = function(train_data, response, cv) {
  tic()
  set.seed(2023)
  coefficient_names = colnames(as.data.frame(train_data))
  response_name = names(as.data.frame(response))
  
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

  for(i in (1: length(all_coef_fit_lasso))) {
    current_coef = all_coef_fit_lasso[i]
    non_zero_coef_lasso_then_OLS_all_predictors = current_coef[all_coef_fit_lasso[i] != 0]
    if (length(non_zero_coef_lasso_then_OLS_all_predictors) == 0){
      next
    }

    non_zero_coef_lasso_then_OLS_all_predictors_indeces = which(current_coef != 0)
    selected_coefficient_names = coefficient_names[non_zero_coef_lasso_then_OLS_all_predictors_indeces]
    new_x_train = train_data[, non_zero_coef_lasso_then_OLS_all_predictors_indeces]
    
    
    data = data.frame(new_x_train, response)
    colnames(data) = c(selected_coefficient_names, response_name)
    fit_OLS_on_LASSO_subset = lm(response~., data = data, x = TRUE, y = TRUE)
    
    fit_OLS_cv_on_LASSO_subset = cv.lm(fit_OLS_on_LASSO_subset, k = 5)
    MSE = fit_OLS_cv_on_LASSO_subset$MSE$mean
    
    if(MSE < current_MSE) {
      current_MSE = MSE
      best_fit = fit_OLS_on_LASSO_subset
    }
  }
  toc()
  return (best_fit)
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