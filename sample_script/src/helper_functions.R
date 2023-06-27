
myRelaxed = function(train_data, response) {
  set.seed(2023)
  fit_lasso = glmnet(train_data, response, family="gaussian", alpha=1)
  
  all_coef_fit_lasso = as.data.frame.matrix(coef(fit_lasso))[-1, ]
  
  current_MSE = .Machine$double.xmax;
  best_fit = lm(1~1)
  non_zero_coef_lasso_then_OLS_all_predictors_indeces_best_fit = c()
  
  for(i in (1: length(all_coef_fit_lasso))) {
    current_coef = all_coef_fit_lasso[i]
    non_zero_coef_lasso_then_OLS_all_predictors = all_coef_fit_lasso[i][all_coef_fit_lasso[i] != 0]
    if (length(non_zero_coef_lasso_then_OLS_all_predictors) == 0){
      next
    }
    
    non_zero_coef_lasso_then_OLS_all_predictors_indeces = which(current_coef != 0)
    
    new_x_train = train_data[, non_zero_coef_lasso_then_OLS_all_predictors_indeces]
    
    data = data.frame(X = new_x_train, y = response)
    
    fit_OLS_on_LASSO_subset = lm(y~., data = data, x = TRUE, y = TRUE)
    
    coef_fit_OLS_on_LASSO_subset = coef(fit_OLS_on_LASSO_subset)[-1]
    
    fit_OLS_cv_on_LASSO_subset = cv.lm(fit_OLS_on_LASSO_subset, k = 5)
    MSE = fit_OLS_cv_on_LASSO_subset$MSE$mean
    
    if(MSE < current_MSE) {
      current_MSE = MSE
      best_fit = fit_OLS_on_LASSO_subset
      
      non_zero_coef_lasso_then_OLS_all_predictors_indeces_best_fit = non_zero_coef_lasso_then_OLS_all_predictors_indeces
      
      
      coef_lasso_then_OLS_all_predictors_final = rep(0, 20)
      j = 1
      for (l in c(1:20)) {
        if (l %in% non_zero_coef_lasso_then_OLS_all_predictors_indeces){
          coef_lasso_then_OLS_all_predictors_final[l] = coef_fit_OLS_on_LASSO_subset[j]
          j = j + 1
        }
      }
    }
  }
  
  fit_lasso_then_OLS_all_predictors = best_fit
  c(fit_lasso_then_OLS_all_predictors, coef_lasso_then_OLS_all_predictors_final)
}