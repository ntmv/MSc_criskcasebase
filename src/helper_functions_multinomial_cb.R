#' ################################################# Relaxed LASSO implementations ##################################################


#' @param train: dataframe containing training data to fit model to, with first column being the response variable and the rest being the covariate values
# Subsets coefficients correctly
multinom.relaxed_lasso_multinomial <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10,
                                                           constant_covariates = 1, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE,
                                                           ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0) {
  
  # FOR TESTING
  # regularization = 'l1'; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; lambda_max = glmnet_lambda_max; gamma = 0.01; i = 10; train = dfM;
  
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




#' @param train: dataframe containing training data to fit model to, with first column being the response variable and the rest being the covariate values
# Subsets coefficients correctly
multinom.relaxed_lasso_cb <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = 10,
                                               constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE,
                                               ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0) {
  
  # FOR TESTING
  # regularization = 'l1'; alpha = 1; nfold = 10;
  # constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .001; grid_size = 100; plot = FALSE;
  # ncores = parallelly::availableCores(); seed = 2023; train_ratio = 20; gamma = 0; i = 10; lambda_max = NULL; train = train_surv;
  
  
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
  cb_data_train_list =  list("time" = cb_data_train$time,
                            "event_ind" = cb_data_train$event_ind,
                            "covariates" = cb_data_train[, grepl("covariates", names(cb_data_train))],
                            "offset" = cb_data_train$offset)
  # Create folds 
  folds <- caret::createFolds(factor(cb_data_train$event_ind), k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  
  p = length(cb_data_train_list$covariates) - 1
  

  
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
                 coefs_cause1 = fit.mtool.penalized$coefficients[1:(p+1), 1]
                 coefs_cause2 = fit.mtool.penalized$coefficients[1:(p+1), 2]
                 selected_cause1 <- which(coefs_cause1 != 0)
                 selected_cause2 <- which(coefs_cause2 != 0)
                 
                 # Get list of names of selected covariates
                 non_zero_cov_names <- c(paste("covariates.X",
                                               as.character(sort(union(selected_cause1[1:length(selected_cause1) - 1], 
                                                                       selected_cause2[1:length(selected_cause2) - 1]))),
                                               sep = ""))
                
                if((p + 1) %in% selected_cause1 | (p + 1) %in% selected_cause2) {
                  non_zero_cov_names <- c(non_zero_cov_names, "covariates.time")
                }
                 
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




#' ################################################# Helper functions for relaxed implementations ##################################################



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
#' Fit case-base sampling model with LASSO
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

#' Multinomial deviance for a casebase fit
#' 
#' @param cb_cb_data Output of \code{create_cb_dataset}
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





#' ################################################# Helper functions for multinomial simulation ##################################################

generateDataset = function(p, n) {
  p = 120
  n = 400
  nfolds = 5
  seed = 2023
  
  # Generate covariates 
  X <- matrix(rnorm(n * p), n, p)
  
  # coefficients for each choice
  X1 <- rep(0, p)
  X2 <- c(rep(3, p/2), rep(0, p/2))
  zero_X2 <- which(X2 == 0)
  
  X3 <- c(rep(3, p/2), rep(0, p/2))
  zero_X3 <- which(X3 == 0)
  
  
  # vector of probabilities
  vProb = cbind(exp(X%*%X1), exp(X%*%X2), exp(X%*%X3))
  
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
  
  train1 = list("covariates" = X, "event_ind" = as.numeric(Y) - 1, "offset" = rep(0, length(Y)), "time" = rep(0, length(Y)))
  
  
  
  
  
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






#' ################################################# Test implementation CV ##################################################


###########################################################
#' Cross-validation function for mtool 
mtool.multinom.cv_test <- function(train, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
                              constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
                              ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20) {
  print("CV")
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