#' ################################################# Helper functions for multinomial simulation ##################################################

generateDatasetMultinomial = function(p, n) {
  # Generate covariates 
  X <- matrix(rnorm(n * p), n, p)
  
  # coefficients for each choice
  beta1 <- rep(0, p)
  beta2 <- c(rep(3, p/2), rep(0, p/2))
  zero_beta2 <- which(beta2 == 0)
  
  beta3 <- c(rep(3, p/2), rep(0, p/2))
  zero_beta3 <- which(beta3 == 0)
  
  
  # vector of probabilities
  vProb = cbind(exp(X%*%beta1), exp(X%*%beta2), exp(X%*%beta3))
  
  # multinomial draws
  mChoices <- t(apply(vProb, 1, rmultinom, n = 1, size = 1))
  dfM <- cbind.data.frame(y = apply(mChoices, 1, function(x) which(x == 1)), X)
  # Rename covariates 
  colnames(dfM)[2:(p+1)] <- paste0('x', 1:p)
  
  # 0, 1, 2 for levels 
  Y <- factor(dfM$y-1)
  
  # Covariate matrix 
  X <- as.matrix(dfM[, c(2:(p+1))])
  
  # Rename covariates 
  colnames(X) <- paste0('x', 1:p)
  
  result = list(X = X, Y = Y)
  
  return(result)
}


partitionDataMultinomial = function(x, y) {
  train_index_y = caret::createDataPartition(y, p = 0.80, list = FALSE)
  
  x_train = x[train_index_y, ]
  y_train = y[train_index_y]
  
  x_test = x[-train_index_y, ]
  y_test = y[-train_index_y]
  
  result = list(x_train = x_train, y_train = y_train,
                x_test = x_test, y_test = y_test)
  
  return(result)
}


fitAllMultinomial = function(train_data, response, folds = 5, print_time, seed) {
  # cv.glmnet with LASSO, cv.glmnet relaxed LASSO, and my relaxed implementation fits
  
  fit.glmnet.relaxed <- cv.glmnet(
    x = train_data, y = response,
    family = "multinomial",
    intercept = TRUE,
    type.multinomial = "grouped",  # same sparsity pattern for all outcome classes
    alpha = 1, nfolds = 10,
    relax = TRUE, gamma = 0)
  
  # mtool
  glmnet_lambda_max = fit.glmnet.relaxed$lambda[1]
  
  data = data.frame(y = as.numeric(response) - 1, train_data)
  
  fit.mtool.relaxed_no_penalization <- multinom.relaxed_lasso_multinomial(train = data, lambda_max = glmnet_lambda_max,  
                                                                          gamma = 0, seed = seed, epsilon = 0.001)
  
  fit.mtool.recolknalaxed_small_penalization <- multinom.relaxed_lasso_multinomial(train = data, lambda_max = glmnet_lambda_max, 
                                                                                   gamma = 0.00001, seed = seed, epsilon = 0.001)
  
  fit.mtool.relaxed_significant_penalization <- multinom.relaxed_lasso_multinomial(train = data, lambda_max = glmnet_lambda_max, 
                                                                                   gamma = 0.01, seed = seed, epsilon = 0.001)
  
  result = list(cv_glmnet_relaxed = relaxed_fit, my_relaxed_no_penalization = fit.mtool.relaxed_no_penalization,
                my_relaxed_small_penalization = fit.mtool.relaxed_small_penalization, 
                my_relaxed_significant_penalization = fit.mtool.relaxed_significant_penalization)
  
  return(result)
}


computeDeviances = function(resulting_fits, test_data) {
  test_data_list = list("covariates" = as.matrix(test_data$x_test),
                        "y" = test_cv$y_test)
  
  pred_fit_relaxed_min = multi_deviance_multinomial(data = test_data_list, resulting_fits$cv_glmnet_relaxed)
  
  pred_myRelaxed_no_pen = multi_deviance_multinomial(data = test_data_list, resulting_fits$my_relaxed_no_penalization)
  pred_myRelaxed_small_pen = multi_deviance_multinomial(data = test_data_list, resulting_fits$my_relaxed_small_penalization)
  pred_myRelaxed_significant_pen = multi_deviance_multinomial(data = test_data_list, resulting_fits$my_relaxed_significant_penalization)
  
  
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




oneIteration = function(p, n, training_data, test_data, print_time, seed) {
  tryCatch({
    # Fit training data with different models
    resulting_fits = fitAllMultinomial(training_data$x_train, training_data$y_train, folds = 5, seed = seed, print_time = FALSE)
    
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



runSimMultinomial = function(p, n, N, print_time, seed) {
  # Generate/split data
  dataset = generateDatasetMultinomial(p, n)
  partitioned_data = partitionDataMultinomial(dataset$X, dataset$Y)
  test_data = partitioned_data[c("x_test", "y_test")]
  
  # Run simulation N times
  sim_replicates = replicate(N, {
    dataset = generateDatasetMultinomial(p, n)
    partitioned_data = partitionDataMultinomial(dataset$X, dataset$Y)
    training_data = partitioned_data[c("x_train", "y_train")]
    
    result_one_iter = oneIterationMultinomial(p, n, training_data, test_data, print_time, seed)
    
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





loadCoefficientBiasTableMultinomial = function(p, N) {
  res_data <- read.csv("results_multinomial_local_sim.csv", header=TRUE, stringsAsFactors=FALSE)
  dev_data <- read.csv("deviances_multinomial_local_sim.csv", header=TRUE, stringsAsFactors=FALSE)
  
  header <- colnames(res_data)
  
  model_fit_labels = c("glmnet.lambda.min_reference1", "glmnet.lambda.min_reference2", "mtool.relaxed.no.pen.lambda.min_reference1",
                       "mtool.relaxed.no.pen.lambda.min_reference2", "mtool.relaxed.small.pen.lambda.min_reference1",
                       "mtool.relaxed.small.pen.lambda.min_reference2", "mtool.relaxed.significant.pen.lambda.min_reference1",
                       "mtool.relaxed.significant.pen.lambda.min_reference2")
  
  
  for (l in model_fit_labels) {
    assign(gsub(" ", "", paste(l, "_df", sep = "")), 
           data.frame(matrix(nrow = 0, ncol = 8)), envir = parent.frame())
  }
  
  colnames(glmnet.lambda.min_reference1_df) <- header
  colnames(glmnet.lambda.min_reference2_df) <- header
  colnames(mtool.relaxed.no.pen.lambda.min_reference1_df) <- header
  colnames(mtool.relaxed.no.pen.lambda.min_reference2_df) <- header
  colnames(mtool.relaxed.small.pen.lambda.min_reference1_df) <- header
  colnames(mtool.relaxed.small.pen.lambda.min_reference2_df) <- header
  colnames(mtool.relaxed.significant.pen.lambda.min_reference1_df) <- header
  colnames(mtool.relaxed.significant.pen.lambda.min_reference2_df) <- header
  
  glmnet.lambda.min_reference1_rows = seq(1, 400, by = 8)
  glmnet.lambda.min_reference2_rows = seq(2, 400, by = 8)
  mtool.relaxed.no.pen.lambda.min_reference1_rows = seq(3, 400, by = 8)
  mtool.relaxed.no.pen.lambda.min_reference2_rows = seq(4, 400, by = 8)
  mtool.relaxed.small.pen.lambda.min_reference1_rows = seq(5, 400, by = 8)
  mtool.relaxed.small.pen.lambda.min_reference2_rows = seq(6, 400, by = 8)
  mtool.relaxed.significant.pen.lambda.min_reference1_rows = seq(7, 400, by = 8)
  mtool.relaxed.significant.pen.lambda.min_reference2_rows = seq(8, 400, by = 8)
  
  glmnet.lambda.min_reference1_df = res_data[glmnet.lambda.min_reference1_rows, ]
  glmnet.lambda.min_reference2_df = res_data[glmnet.lambda.min_reference2_rows, ]
  mtool.relaxed.no.pen.lambda.min_reference1_df = res_data[mtool.relaxed.no.pen.lambda.min_reference1_rows, ]
  mtool.relaxed.no.pen.lambda.min_reference2_df = res_data[mtool.relaxed.no.pen.lambda.min_reference2_rows, ]
  mtool.relaxed.small.pen.lambda.min_reference1_df = res_data[mtool.relaxed.small.pen.lambda.min_reference1_rows, ]
  mtool.relaxed.small.pen.lambda.min_reference2_df = res_data[mtool.relaxed.small.pen.lambda.min_reference2_rows, ]
  mtool.relaxed.significant.pen.lambda.min_reference1_df = res_data[mtool.relaxed.significant.pen.lambda.min_reference1_rows, ]
  mtool.relaxed.significant.pen.lambda.min_reference2_df = res_data[mtool.relaxed.significant.pen.lambda.min_reference2_rows, ]
  
  
  glmnet.lambda.min_reference1_results = colMeans(glmnet.lambda.min_reference1_df[, c(5:8)])
  glmnet.lambda.min_reference2_results = colMeans(glmnet.lambda.min_reference2_df[, c(5:8)])
  mtool.relaxed.no.pen.lambda.min_reference1_results = colMeans(mtool.relaxed.no.pen.lambda.min_reference1_df[, c(5:8)])
  mtool.relaxed.no.pen.lambda.min_reference2_results = colMeans(mtool.relaxed.no.pen.lambda.min_reference2_df[, c(5:8)])
  mtool.relaxed.small.pen.lambda.min_reference1_results = colMeans(mtool.relaxed.small.pen.lambda.min_reference1_df[, c(5:8)])
  mtool.relaxed.small.pen.lambda.min_reference2_results = colMeans(mtool.relaxed.small.pen.lambda.min_reference2_df[, c(5:8)])
  mtool.relaxed.significant.pen.lambda.min_reference1_results = colMeans(mtool.relaxed.significant.pen.lambda.min_reference1_df[, c(5:8)])
  mtool.relaxed.significant.pen.lambda.min_reference2_results = colMeans(mtool.relaxed.significant.pen.lambda.min_reference2_df[, c(5:8)])
  
  final_results_multinomial = rbind(glmnet.lambda.min_reference1_results, glmnet.lambda.min_reference2_results,
                                    mtool.relaxed.no.pen.lambda.min_reference1_results, mtool.relaxed.no.pen.lambda.min_reference2_results,
                                    mtool.relaxed.small.pen.lambda.min_reference1_results, mtool.relaxed.small.pen.lambda.min_reference2_results,
                                    mtool.relaxed.significant.pen.lambda.min_reference1_results, mtool.relaxed.significant.pen.lambda.min_reference2_results)
  
  colnames(final_results_multinomial) <- names(glmnet.lambda.min_reference1_results)
  
  
  # Calculate summarized coefficient statistics
  
  model_deviance_labels = c("glmnet.lambda.min","mtool.relaxed.no.pen.lambda.min",
                            "mtool.relaxed.small.pen.lambda.min", "mtool.relaxed.significant.pen.lambda.min")
  
  glmnet.lambda.min_rows = seq(1, 80, by = 4)
  mtool.relaxed.no.pen.lambda.min_rows = seq(2, 80, by = 4)
  mtool.relaxed.small.pen.lambda.min_rows = seq(3, 80, by = 4)
  mtool.relaxed.significant.pen.lambda.min_rows = seq(4, 80, by = 4)
  
  
  glmnet.lambda.min_df = dev_data[glmnet.lambda.min_rows, ]
  mtool.relaxed.no.pen.lambda.min_df = dev_data[mtool.relaxed.no.pen.lambda.min_rows, ]
  mtool.relaxed.small.pen.lambda.min_df = dev_data[mtool.relaxed.small.pen.lambda.min_rows, ]
  mtool.relaxed.significant.pen.lambda.min_df = dev_data[mtool.relaxed.significant.pen.lambda.min_rows, ]
  
  
  
  glmnet.lambda.min_results = mean(glmnet.lambda.min_df)
  mtool.relaxed.no.pen.lambda.min_results = mean(mtool.relaxed.no.pen.lambda.min_df)
  mtool.relaxed.small.pen.lambda.min_results = mean(mtool.relaxed.small.pen.lambda.min_df)
  mtool.relaxed.significant.pen.lambda.min_results = mean(mtool.relaxed.significant.pen.lambda.min_df)
  
  
  final_deviances_multinomial = rbind(glmnet.lambda.min_results, mtool.relaxed.no.pen.lambda.min_results, 
                                      mtool.relaxed.small.pen.lambda.min_results, mtool.relaxed.significant.pen.lambda.min_results)
  
  rownames(final_deviances_multinomial) <- model_deviance_labels
  
  results_multinomial = list(final_results_multinomial, final_deviances_multinomial)
  
  return(results_multinomial)
}


#' ################################################# Helper functions for survival simulation ##################################################


loadCoefficientBiasTableSurvival = function(p, N) {
  res_data1 <- read.csv("results_survival_local_sim.csv", header=TRUE, stringsAsFactors=FALSE)
  res_data2 <- read.csv("results_survival_local_sim2.csv", header=TRUE, stringsAsFactors=FALSE)
  dev_data <- read.csv("deviances_survival_local_sim.csv", header=TRUE, stringsAsFactors=FALSE)
  
  res_data <- rbind(res_data1, res_data2)
  
  header <- colnames(res_data)
  
  model_fit_labels = c("cox.lambda.min_cause1", "cox.lambda.min_cause2",
                       "casebase.cv.lambda.min_cause1", "casebase.cv.lambda.min_cause2",
                       "casebase.relaxed.lambda.min.no.pen_cause1", "casebase.relaxed.lambda.min.no.pen_cause2",
                       "casebase.relaxed.lambda.min.small.pen_cause1", "casebase.relaxed.lambda.min.small.pen_cause2",
                       "casebase.relaxed.lambda.min.sig.pen_cause1", "casebase.relaxed.lambda.min.sig.pen_cause2")
  
  
  for (l in model_fit_labels) {
    assign(gsub(" ", "", paste(l, "_df", sep = "")), 
           data.frame(matrix(nrow = 0, ncol = 8)), envir = parent.frame())
  }
  
  # Calculate summarized coefficient statistics
  
  colnames(cox.lambda.min_cause1_df) <- header
  colnames(cox.lambda.min_cause2_df) <- header
  colnames(casebase.cv.lambda.min_cause1_df) <- header
  colnames(casebase.cv.lambda.min_cause2_df) <- header
  colnames(casebase.relaxed.lambda.min.no.pen_cause1_df) <- header
  colnames(casebase.relaxed.lambda.min.no.pen_cause2_df) <- header
  colnames(casebase.relaxed.lambda.min.small.pen_cause1_df) <- header
  colnames(casebase.relaxed.lambda.min.small.pen_cause2_df) <- header
  colnames(casebase.relaxed.lambda.min.small.pen_cause1_df) <- header
  colnames(casebase.relaxed.lambda.min.small.pen_cause2_df) <- header
  
  cox.lambda.min_cause1_rows = seq(1, 400, by = 10)
  cox.lambda.min_cause2_rows = seq(2, 400, by = 10)
  casebase.cv.lambda.min_cause1_rows = seq(3, 400, by = 10)
  casebase.cv.lambda.min_cause2_rows = seq(4, 400, by = 10)
  casebase.relaxed.lambda.min.no.pen_cause1_rows = seq(5, 400, by = 10)
  casebase.relaxed.lambda.min.no.pen_cause2_rows = seq(6, 400, by = 10)
  casebase.relaxed.lambda.min.small.pen_cause1_rows = seq(7, 400, by = 10)
  casebase.relaxed.lambda.min.small.pen_cause2_rows = seq(8, 400, by = 10)
  casebase.relaxed.lambda.min.sig.pen_cause1_rows = seq(9, 400, by = 10)
  casebase.relaxed.lambda.min.sig.pen_cause2_rows = seq(10, 400, by = 10)
  
  cox.lambda.min_cause1_df = res_data[cox.lambda.min_cause1_rows, ]
  cox.lambda.min_cause2_df = res_data[cox.lambda.min_cause2_rows, ]
  casebase.cv.lambda.min_cause1_df = res_data[casebase.cv.lambda.min_cause1_rows, ]
  casebase.cv.lambda.min_cause2_df = res_data[casebase.cv.lambda.min_cause2_rows, ]
  casebase.relaxed.lambda.min.no.pen_cause1_df = res_data[casebase.relaxed.lambda.min.no.pen_cause1_rows, ]
  casebase.relaxed.lambda.min.no.pen_cause2_df = res_data[casebase.relaxed.lambda.min.no.pen_cause2_rows, ]
  casebase.relaxed.lambda.min.small.pen_cause1_df = res_data[casebase.relaxed.lambda.min.small.pen_cause1_rows, ]
  casebase.relaxed.lambda.min.small.pen_cause2_df = res_data[casebase.relaxed.lambda.min.small.pen_cause2_rows, ]
  casebase.relaxed.lambda.min.sig.pen_cause1_df = res_data[casebase.relaxed.lambda.min.sig.pen_cause1_rows, ]
  casebase.relaxed.lambda.min.sig.pen_cause2_df = res_data[casebase.relaxed.lambda.min.sig.pen_cause2_rows, ]
  
  cox.lambda.min_cause1_results = colMeans(cox.lambda.min_cause1_df[, c(5:8)])
  cox.lambda.min_cause2_results = colMeans(cox.lambda.min_cause2_df[, c(5:8)])
  casebase.cv.lambda.min_cause1_results = colMeans(casebase.cv.lambda.min_cause1_df[, c(5:8)])
  casebase.cv.lambda.min_cause2_results = colMeans(casebase.cv.lambda.min_cause2_df[, c(5:8)])
  casebase.relaxed.lambda.min.no.pen_cause1_results = colMeans(casebase.relaxed.lambda.min.no.pen_cause1_df[, c(5:8)])
  casebase.relaxed.lambda.min.no.pen_cause2_results = colMeans(casebase.relaxed.lambda.min.no.pen_cause2_df[, c(5:8)])
  casebase.relaxed.lambda.min.small.pen_cause1_results = colMeans(casebase.relaxed.lambda.min.small.pen_cause1_df[, c(5:8)])
  casebase.relaxed.lambda.min.small.pen_cause2_results = colMeans(casebase.relaxed.lambda.min.small.pen_cause2_df[, c(5:8)])
  casebase.relaxed.lambda.min.sig.pen_cause1_results = colMeans(casebase.relaxed.lambda.min.sig.pen_cause1_df[, c(5:8)])
  casebase.relaxed.lambda.min.sig.pen_cause2_results = colMeans(casebase.relaxed.lambda.min.sig.pen_cause2_df[, c(5:8)])
  
  final_results_survival = rbind(cox.lambda.min_cause1_results, cox.lambda.min_cause2_results,
                                 casebase.cv.lambda.min_cause1_results, casebase.cv.lambda.min_cause2_results,
                                 casebase.relaxed.lambda.min.no.pen_cause1_results, casebase.relaxed.lambda.min.no.pen_cause2_results,
                                 casebase.relaxed.lambda.min.small.pen_cause1_results, casebase.relaxed.lambda.min.small.pen_cause2_results,
                                 casebase.relaxed.lambda.min.sig.pen_cause1_results, casebase.relaxed.lambda.min.sig.pen_cause2_results)
  
  colnames(final_results_survival) <- names(cox.lambda.min_cause1_results)
  
  
  # Calculate summarized coefficient statistics
  
  model_deviance_labels = c("casebase.cv.lambda.min","casebase.relaxed.lambda.min.no.pen",
                            "casebase.relaxed.lambda.min.small.pen", "casebase.relaxed.lambda.min.sig.pen")
  
  casebase.cv.lambda.min_rows = seq(1, 80, by = 4)
  casebase.relaxed.lambda.min.no.pen_rows = seq(2, 80, by = 4)
  casebase.relaxed.lambda.min.small.pen_rows = seq(3, 80, by = 4)
  casebase.relaxed.lambda.min.sig.pen_rows = seq(4, 80, by = 4)
  
  
  casebase.cv.lambda.min_df = dev_data[casebase.cv.lambda.min_rows, ]
  casebase.relaxed.lambda.min.no.pen_df = dev_data[casebase.relaxed.lambda.min.no.pen_rows, ]
  casebase.relaxed.lambda.min.small.pen_df = dev_data[casebase.relaxed.lambda.min.small.pen_rows, ]
  casebase.relaxed.lambda.min.sig.pen_df = dev_data[casebase.relaxed.lambda.min.sig.pen_rows, ]
  
  
  
  casebase.cv.lambda.min_results = mean(casebase.cv.lambda.min_df)
  casebase.relaxed.lambda.min.no.pen_results = mean(casebase.relaxed.lambda.min.no.pen_df)
  casebase.relaxed.lambda.min.small.pen_results = mean(casebase.relaxed.lambda.min.small.pen_df)
  casebase.relaxed.lambda.min.sig.pen_results = mean(casebase.relaxed.lambda.min.sig.pen_df)
  
  
  final_deviances_survival = rbind(casebase.cv.lambda.min_results, casebase.relaxed.lambda.min.no.pen_results, 
                                   casebase.relaxed.lambda.min.small.pen_results, casebase.relaxed.lambda.min.sig.pen_results)
  
  rownames(final_deviances_survival) <- model_deviance_labels
  
  results_survival = list(final_results_survival, final_deviances_survival)
  
  return(results_survival)
}










#' ################################################# Test implementation CV ##################################################


###########################################################
#' Cross-validation function for mtool 
mtool.multinom.cv_test <- function(train, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
                                   constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                                   ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20) {
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
  
  cb_data_train_list = cb_data_train
  cb_data_train <- as.data.frame(cb_data_train)
  cb_data_train <- cb_data_train %>%
    select(-time)
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
    # Standardize
    train_cv$covariates <- as.data.frame(scale(train_cv$covariates, center = T, scale = T))
    cvs_res <- mclapply(lambdagrid, 
                        function(lambda_val) {
                          fit_cbmodel(train_cv, regularization = 'elastic-net',
                                      lambda = lambda_val, alpha = alpha, unpen_cov = constant_covariates)
                        }, mc.cores = ncores, mc.set.seed = seed)
    test_cv <- list("time" = test_cv$covariates.time,
                    "event_ind" = test_cv$event_ind,
                    "covariates" = as.matrix(test_cv[, grepl("covariates", names(test_cv))]),
                    "offset" = test_cv$offset)
    # Standardize
    test_cv$covariates <- as.data.frame(scale(test_cv$covariates, center = T, scale = T))
    mult_deviance <- unlist(lapply(cvs_res, multi_deviance_cb, cb_data = test_cv))
    all_deviances[, i] <- mult_deviance
    non_zero_coefs[, i] <-  unlist(lapply(cvs_res, function(x) {return(x$no_non_zero)}))
    
    if(i == nfold) {
      # Return the number of selected covariates when fit on all data
      all_non_coefs <- mclapply(lambdagrid, 
                                function(lambda_val) {
                                  fit.all.data = fit_cbmodel(cb_data_train_list, regularization = 'elastic-net',
                                                             lambda = lambda_val, alpha = alpha, unpen_cov = constant_covariates)
                                  num_non_zero_cov = length(union(which(fit.all.data$coefficients[, 1] != 0),
                                                                  which(fit.all.data$coefficients[, 2] != 0)))
                                  num_non_zero_cov
                                }, mc.cores = ncores, mc.set.seed = seed)
    }
    
    cat("Completed Fold", i, "\n")
  }
  mean_dev <- rowMeans(all_deviances)
  lambda.min <- lambdagrid[which.min(mean_dev)]
  sel_lambda_min <- non_zero_coefs[which.min(mult_deviance)]
  if (sel_lambda_min  == 2*constant_covariates) {
    cat("Null model chosen: choosing first non-null model lambda")
    lambda.min <- lambdagrid[which.min(non_zero_coefs != 2*constant_covariates)-2]
  }
  cv_se <- sqrt(var(mean_dev))
  dev.1se <- mean_dev[which.min(mean_dev)] + cv_se
  dev.0.5se <- mean_dev[which.min(mean_dev)] + cv_se/2
  range.1se <- lambdagrid[which(mean_dev <= dev.1se)]
  lambda.1se <- max(range.1se)
  lambda.min1se <- min(range.1se)
  range.0.5se <- lambdagrid[which((mean_dev <= dev.0.5se))]
  lambda.0.5se <- max(range.0.5se)
  lambda.min0.5se <- min(range.0.5se)
  rownames(all_deviances) <- lambdagrid
  rownames(non_zero_coefs) <- lambdagrid
  return(list(lambda.min = lambda.min,  non_zero_coefs = non_zero_coefs, lambda.min1se = lambda.min1se, lambda.min0.5se = lambda.min0.5se, 
              lambda.1se = lambda.1se, lambda.0.5se = lambda.0.5se, cv.se = cv_se, lambdagrid = lambdagrid, deviance_grid = all_deviances, 
              coef_grid = all_non_coefs))
}