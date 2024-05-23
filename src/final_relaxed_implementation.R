#' ################################################# Relaxed LASSO implementation for multinomial and survival data ##################################################

multinomial_casebase_relaxed_lasso <- function(train, regularization = 'l1', lambda_max = NULL, alpha = 1, nfold = NULL, 
                                            initial_max_grid = NULL, precision = 0.001, epsilon = NULL, grid_size = 100, plot = FALSE,
                                            ncores = parallelly::availableCores(), seed = NULL, train_ratio = 20, gamma = 0.001) {
  
  if(is.null(train$fstatus)) 
    model_type <- "multinomial"
  else
    model_type <- "casebase"
  
  data_and_parameters <- prepare_data_and_parameters(train = train, epsilon = epsilon, nfold = nfold, model_type = model_type,
                                                     train_ratio = train_ratio)
  
  data <- data_and_parameters$data
  epsilon <- data_and_parameters$epsilon
  nfold <- data_and_parameters$nfold
  constant_covariates <- data_and_parameters$constant_covariates
  folds <- data_and_parameters$folds
  
  lambdagrid <- find_lambda_grid(train_data = data, lambda_max = lambda_max, alpha = alpha, 
                                 constant_covariates = constant_covariates,
                                 initial_max_grid = initial_max_grid, precision = precision, 
                                 epsilon = epsilon, ncores = ncores, seed = seed, model_type = model_type)
  
  p = length(data$covariates)
  lambda.min <- rep(NA_real_, nfold)
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs_matrix <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  
  for (i in 1:nfold) {
    split <- split_data(data, folds, i, model_type)
    train_cv <- split$train_data
    validation_cv <- split$validation_data

    train_cv$covariates <- scale(train_cv$covariates, center = TRUE, scale = TRUE)
    validation_cv$covariates <- scale(validation_cv$covariates, center = TRUE, scale = TRUE)
    
    
    res <- unlist(
      mclapply(lambdagrid,
               function(lambda_v) {
                 dev <- 0
                 
                 fit.mtool.penalized <- fit_model(train_cv, regularization = 'l1',
                                                  lambda = lambda_v, alpha = alpha, unpen_cov = constant_covariates)
                 
                 non_zero_cov_names <- find_non_zero_coefficients(p, fit.mtool.penalized, model_type)
                 
                 train_cv$covariates <- train_cv$covariates[, colnames(train_cv$covariates) %in% non_zero_cov_names]
                 validation_cv$covariates <- validation_cv$covariates[, colnames(validation_cv$covariates) %in% non_zero_cov_names]     
                 
                 fit.mtool.subset <- fit_model(train_cv, regularization = 'l1',
                                               lambda = gamma, alpha = alpha, unpen_cov = constant_covariates)
                 
                 dev <- multi_deviance_final(data = validation_cv, fit_object = fit.mtool.subset)
                 
                 # During last fold, fit penalized model with gamma to entire dataset
                 if(i == nfold) {
                   data$covariates <- data$covariates[, colnames(data$covariates) %in% non_zero_cov_names]
                   fit.all.data <- fit_model(data, regularization = 'l1',
                                             lambda = gamma, alpha = alpha, unpen_cov = constant_covariates)
                   num_non_zero_cov <- length(union(which(fit.all.data$coefficients[, 1] != 0),
                                                    which(fit.all.data$coefficients[, 2] != 0)))
                   dev <- list(num_non_zero_cov, dev)
                 }
                 dev
               }, mc.cores = ncores, mc.set.seed = seed))
    
    if(i == nfold) {
      all_deviances[, i] <- res[seq(2, 2 * length(lambdagrid), by = 2)]
      all_non_coefs <- res[seq(1, 2 * length(lambdagrid), by = 2)]
    } else {
      all_deviances[, i] <- res
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


#' ################################################# Helper functions for relaxed implementations ##################################################
prepare_data_and_parameters <- function(train, epsilon, nfold, model_type, train_ratio) {
  if (model_type == 'multinomial') {
    data <- prepare_multinomial_data(train)
    
    if(is.null(epsilon))
      epsilon = 0.001
    
    if(is.null(nfold))
      nfold = 10
    
    constant_covariates = 1
    folds <- caret::createFolds(y = data$y, k = nfold, list = FALSE)
  } else if (model_type == 'casebase') {
    data <- prepare_casebase_data(train, train_ratio)
    
    if(is.null(epsilon))
      epsilon = 0.0001
    
    if(is.null(nfold))
      nfold = 5
    
    constant_covariates = 2
    folds <- caret::createFolds(factor(data$event_ind), k = nfold, list = FALSE)
  }
  list(data = data, epsilon = epsilon, nfold = nfold, constant_covariates = constant_covariates, folds = folds)
}


prepare_multinomial_data <- function(train) {
  covariates <- as.matrix(cbind(train[, -1]))
  y <- train$y
  list(covariates = covariates, y = y)
}

prepare_casebase_data <- function(train, train_ratio) {
  surv_obj_train <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))
  cov_train <- as.matrix(cbind(train[, grepl("X", colnames(train))], time = log(train$ftime)))
  cb_data_train <- create_cbDataset(surv_obj_train, cov_train, ratio = train_ratio)
  
  cb_data_train <- as.data.frame(cb_data_train)
  covariates <- cb_data_train[, grepl("covariates", names(cb_data_train))]
  event_ind <- cb_data_train$event_ind
  offset <- cb_data_train$offset
  
  list(time = NULL, covariates = covariates, event_ind = event_ind, offset = offset)
}


find_lambda_grid <- function(train_data, lambda_max, alpha = 1, constant_covariates, initial_max_grid, precision = 0.001, epsilon, grid_size = 100,
                             ncores = parallelly::availableCores(), seed = NULL, model_type) {
  if(is.null(lambda_max)) {
    # Lambda max grid for bisection search
    if(is.null(initial_max_grid)) {
      initial_max_grid <-  c(0.9, 0.5, 0.1, 0.07, 0.05, 0.01, 0.009, 0.005)
      fit_val_max <- mclapply(initial_max_grid,
                              function(lambda_val) {
                                fit_model(train_data, regularization = 'l1',
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
                                 fit_model(train_data, regularization = 'l1',
                                           lambda = lambda_val, alpha = alpha,
                                           unpen_cov = constant_covariates)}, mc.cores = ncores, mc.set.seed = seed)
      non_zero_coefs <-  unlist(mclapply(fit_val_max, function(x) {return(x$no_non_zero)}, mc.cores = ncores, mc.set.seed = seed))
      lambda_max <- new_max_searchgrid[which.min(non_zero_coefs)]
    }
  }
  
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  
  lambdagrid
}


split_data <- function(data, folds, fold_idx, model_type) {
  train_idx <- which(folds != fold_idx)
  validation_idx <- which(folds == fold_idx)
  
  if (model_type == 'multinomial') {
    train_data <- list(covariates = data$covariates[train_idx, ], y = data$y[train_idx])
    validation_data <- list(covariates = data$covariates[validation_idx, ], y = data$y[validation_idx])
  } else {
    train_data <- list(time = NULL, covariates = data$covariates[train_idx, ], event_ind = data$event_ind[train_idx], offset = data$offset[train_idx])
    validation_data <- list(time = NULL, covariates = data$covariates[validation_idx, ], event_ind = data$event_ind[validation_idx], offset = data$offset[validation_idx])
  }
  
  list(train_data = train_data, validation_data = validation_data)
}


find_non_zero_coefficients <- function(p, mtool_fit, model_type) {
  if(model_type == "multinomial") {
    coefs_cause1 = mtool_fit$coefficients[, 1]
    coefs_cause2 = mtool_fit$coefficients[, 2]
    selected_cause1 <- which(coefs_cause1 != 0)
    selected_cause2 <- which(coefs_cause2 != 0)
    
    non_zero_cov_names <- c(paste("x",
                                  as.character(sort(union(selected_cause1, selected_cause2))),
                                  sep = ""))
  } else {
    coefs_cause1 = mtool_fit$coefficients[1:p, 1]
    coefs_cause2 = mtool_fit$coefficients[1:p, 2]
    selected_cause1 <- which(coefs_cause1 != 0)
    selected_cause2 <- which(coefs_cause2 != 0)
    
    non_zero_cov_names <- c(paste("covariates.X",
                                  as.character(sort(union(selected_cause1[1:length(selected_cause1) - 1], 
                                                          selected_cause2[1:length(selected_cause2) - 1]))),
                                  sep = ""))
    
    if(p %in% selected_cause1 | p %in% selected_cause2) {
      non_zero_cov_names <- c(non_zero_cov_names, "covariates.time")
    }
  }
  
  non_zero_cov_names
}


################################################
#' Fit multinomial logistic regression or casebase model
#' 
#' @param data Output of \code{create_cbDataset}
fit_model <- function(data, regularization = 'l1',
                        lambda, alpha = 1, unpen_cov = 1) {
  stopifnot(alpha >= 0 && alpha <= 1)
  # Elastic-net reparametrization
  lambda1 <- lambda*alpha
  lambda2 <- 0.5*lambda*(1 - alpha)
  # Prepare covariate matrix with intercept
  X <- as.matrix(cbind(data$covariates, 1))
  
  if(is.null(data$event_ind)) {
    print("multinomial data identified")
    Y = data$y
    offset = rep(0, length(data$y))
  } else {
    print("survival data identified")
    Y = data$event_ind
    offset = data$offset
  }
  
  out <- fit.mtool <- mtool::mtool.MNlogistic(
    X = as.matrix(X),
    Y = Y,
    offset = offset,
    N_covariates = unpen_cov,
    regularization = regularization,
    # learning_rate = 0.001,
    transpose = FALSE,
    lambda1 = lambda1, lambda2 = lambda2,
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
multi_deviance_final <- function(data, fit_object, cox = FALSE) {
  # If cox model, don't add bias column to covariate matrix
  if(cox) {
    X <- as.matrix(data$covariates) 
  } else {
    X <- as.matrix(cbind(data$covariates, 1)) 
  }
  fitted_vals <- as.matrix(X %*% fit_object$coefficients)
  pred_mat <- VGAM::multilogitlink(fitted_vals, 
                                   inverse = TRUE)
  # Turn event_ind into Y_mat
  if(is.null(data$event_ind)) {
    print("multinomial data identified")
    Y_fct <- factor(data$y)
  } else {
    print("survival data identified")
    Y_fct <- factor(data$event_ind)
  }
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


plot_post_relaxed <- function(cv_object) {
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
        grob = grid::textGrob(label = labels[i], hjust = 0.5, gp = grid::gpar(cex = 0.9)),
        xmin = positions[i], xmax = positions[i],
        ymin = y_limit, ymax = y_limit
      )
  }
  
  return(p)
}
