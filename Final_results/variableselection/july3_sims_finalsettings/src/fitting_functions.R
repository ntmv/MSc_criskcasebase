############### All fitting functions for multinomial casebase simulations ############################
#' Multinomial deviance
#' 
#' @param cb_data Output of \code{create_cbDataset}
#' @param fit_object Output of \code{fit_cbmodel}
#' @return Multinomial deviance
multi_deviance <- function(cb_data, fit_object) {
  X <- as.matrix(cbind(cb_data$covariates, 1))
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

###############################################
#' Create the case-base sampled dataset
#'
#' @param surv_obj
#' @param cov_matrix
#' @param ratio
#' @return List of 4 elements, containing the necessary data to fit a case-base
#'   regression model.
create_cbDataset <- function(surv_obj, cov_matrix, ratio = 5) {
    n <- nrow(surv_obj)
    B <- sum(surv_obj[, "time"])
    c1 <- sum(surv_obj[, "status"] == 1)
    c2 <- sum(surv_obj[, "status"] == 2)
    b <- ratio * (c1)
    offset <- log(B/b)
    prob_select <-  surv_obj[, "time"]/B
    # Create base series
    which_pm <- sample(n, b, replace = TRUE, prob = prob_select)
    bSeries <- as.matrix(surv_obj[which_pm, ])
    time_bseries <- runif(b) * bSeries[, "time"]
    cov_bseries <- cov_matrix[which_pm, , drop = FALSE]
    event_bseries <- rep(0L, nrow(bSeries))
    # Extract case series
    cSeries <- as.matrix(surv_obj[surv_obj[,"status"] != 0L, ])
    time_cseries <- cSeries[,"time"]
    cov_cseries <- cov_matrix[surv_obj[,"status"] != 0L, , drop = FALSE]
    event_cseries <- cSeries[,"status"]
    # Combine and return
    output <- list("time" = c(time_bseries, time_cseries),
                   "event_ind" = c(event_bseries, event_cseries),
                   "covariates" = rbind(cov_bseries, cov_cseries),
                   "offset" = rep(offset, nrow(bSeries) + nrow(cSeries)))
    
    return(output)
}

################################################
#' Fit case-base sampling model
#' 
#' @param cb_data Output of \code{create_cbDataset}
fit_cbmodel <- function(cb_data, regularization = 'elastic-net',
                        lambda, alpha = 0.5, unpen_cov = 2) {
    stopifnot(alpha >= 0 && alpha <= 1)
    # Elastic-net reparametrization
    lambda1 <- lambda*alpha
    lambda2 <- 0.5*lambda*(1 - alpha)
    # Prepare covariate matrix with intercept
    X <- as.matrix(cbind(cb_data$covariates, 1))
    out <- fit.mtool <- mtool::mtool.MNlogistic(
        X = as.matrix(X),
        Y = cb_data$event_ind,
        offset = cb_data$offset,
        N_covariates = unpen_cov,
        regularization = 'elastic-net',
        transpose = FALSE,
        lambda1 = lambda1, lambda2 = lambda2, 
        lambda3 = 0
        )
    
    return(out)
}


################################################
#' Generate an AR(1) covariance matrix for p variables with correlation rho
#' 
#' @param p Number of covariates
#' @param rho Correlation parameter
#' @return pxp matrix whose (i,j) entry is rho^abs(i - j).
covAR1 <- function(p, rho) {
    stopifnot(p >= 2 && rho >= 0)
    # https://statisticaloddsandends.wordpress.com/2020/02/07/generating-correlation-matrix-for-ar1-model/
    exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                        (1:p - 1))
    
    return(rho^exponent)
}

###########################################################
#' Cross-validation function for mtool 
mtool.multinom.cv <- function(cb_data_train, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
                              constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE, 
    ncores = parallelly::availableCores(), seed = NULL) {
  repeat {
   # Default lambda grid
  if(is.null(lambda_max)) {
    # Lambda max grid for bisection search 
    if(is.null(initial_max_grid)) {
      initial_max_grid <- 
        c(0.9, 0.5, 0.1, 0.07, 0.05, 0.01, 0.009, 0.005)
      fit_val_max <- mclapply(initial_max_grid, 
                               function(lambda_val) {
                                 fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                                             lambda = lambda_val, alpha = alpha, unpen_cov = constant_covariates)}, mc.cores = 4)
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
  }
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
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
    cvs_res <- mclapply(lambdagrid, 
      function(lambda_val) {
      fit_cbmodel(train_cv, regularization = 'elastic-net',
                  lambda = lambda_val, alpha = alpha, unpen_cov = constant_covariates)
    }, mc.cores = ncores, mc.set.seed = seed)
    test_cv <- list("time" = test_cv$covariates.time,
                    "event_ind" = test_cv$event_ind,
                    "covariates" = as.matrix(test_cv[, grepl("covariates", names(test_cv))]),
                    "offset" = test_cv$offset)
    mult_deviance <- unlist(lapply(cvs_res, multi_deviance, cb_data = test_cv))
    all_deviances[, i] <- mult_deviance
    mean_dev <- rowMeans(all_deviances)
    lambda.min <- lambdagrid[which.min(mean_dev)]
    cv_se <- sqrt(var(mean_dev))
    dev.1se <- mean_dev[which.min(mean_dev)] + cv_se
    dev.0.5se <- mean_dev[which.min(mean_dev)] + cv_se/2
    range.1se <- lambdagrid[which(mean_dev <= dev.1se)]
    lambda.1se <- max(range.1se)
    lambda.min1se <- min(range.1se)
    range.0.5se <- lambdagrid[which((mean_dev <= dev.0.5se))]
    lambda.0.5se <- max(range.0.5se)
    lambda.min0.5se <- min(range.0.5se)
    non_zero_coefs[, i] <-  unlist(lapply(cvs_res, function(x) {return(x$no_non_zero)}))
    rownames(all_deviances) <- lambdagrid
    rownames(non_zero_coefs) <- lambdagrid
  }
  if (lambda.min != lambda_max) {
    break
  }
  message("lambda.min is equal to lambda.max. Re-running the function.")
  }
  return(list(lambda.min = lambda.min,  non_zero_coefs = non_zero_coefs, lambda.min1se = lambda.min1se, lambda.min0.5se = lambda.min0.5se, 
              lambda.1se = lambda.1se, lambda.0.5se = lambda.0.5se, cv.se = cv_se, lambdagrid = lambdagrid, deviance_grid = all_deviances))
}

#' ###########################################################
#' Plotting function for cross-validation 
plot_cv.multinom <- function(cv_object) {
  nfold <- ncol(cv_object$deviance_grid)
  mean_dev <- rowMeans(cv_object$deviance_grid)
  row_stdev <- apply(cv_object$deviance_grid, 1, function(x) {sd(x)/sqrt(nfold)})
  plot.dat <- data.frame(lambdagrid = cv_object$lambdagrid, mean.dev = mean_dev, 
                         upper = mean_dev +row_stdev, lower = mean_dev - row_stdev)
  p <- ggplot(plot.dat, aes(log(lambdagrid), mean.dev)) + geom_point(colour = "red", size = 3) + theme_bw() + 
    geom_errorbar(aes(ymin= lower, ymax=upper), width=.2, position=position_dodge(0.05), colour = "grey") + 
    labs(x = "log(lambda)", y = "Multinomial Deviance")  + 
    geom_vline(xintercept = log(cv_object$lambda.min), linetype = "dotted", colour = "blue") + 
    geom_vline(xintercept = log(cv_object$lambda.1se), linetype = "dotted", 
               colour = "purple")
  return(p)
}


################### Function to calculate specificity, sensitivity and MCC ###############
varsel_perc <- function(model_coef, true_coef) {
  zero_ind1 <- which(true_coef == 0)
  nonzero_ind1 <- which(true_coef != 0)
  TP <- length(intersect(which(model_coef != 0), nonzero_ind1))
  TN <- length(intersect(which(model_coef == 0), zero_ind1))
  FP <- length(intersect(which(model_coef != 0), zero_ind1))
  FN <- length(intersect(which(model_coef == 0), nonzero_ind1))
  sens <- TP/(TP+FN)
  spec <- TN/(TN+FP)
  mcc_num <- (TP*TN)-(FP*FN)
  # To avoid numerical overflow issues 
  mcc_den1 <- as.numeric((TP+FP)*(TP+FN))
  mcc_den2 <- as.numeric((TN+FP)*(TN+FN))
  mcc <- mcc_num/sqrt(mcc_den1*mcc_den2)
  dat <- as.data.frame(cbind(TP = TP, TN = TN, FP = FP, FN = FN, Sensitivity = sens, Specificity = spec, MCC = mcc))
  return(dat)
}


###############################################################################
#' # Simulate from sub-distribution hazards
cause_subdist_sim <- function(n, p, beta1, beta2, num.true = 20, mix_p = 0.5
                                     , cor_vals = c(0.7, 0.4, 0.6, 0.5), noise_cor = 0.1, 
                              nblocks = 4,  rho1 = 4, u.max = 1, lambda1 = 1, lambda2 = 0.8, rho2 = 10, max_time = 2) {
  # Warnings
  if(length(beta1) != length(beta2)) stop("Dimension of beta1 and beta2 should be the same")
  # Set the number of variables per block
  vpb <- num.true/nblocks
  # Set the correlation values for each covariate block
  correlation_values <- cor_vals
  # Initialize empty matrix
  correlation_matrix <- matrix(noise_cor, nrow = p, ncol = p)
  # Generate the covariance matrix with block correlations
  for (i in 1:nblocks) {
    start_index <- (i - 1) * vpb + 1
    end_index <- i * vpb
    correlation_matrix[start_index:end_index, start_index:end_index] <- correlation_values[i]
  }
  # Diagonal elements should be 1
  diag(correlation_matrix) <- rep(1, length(diag(correlation_matrix)))
  
  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = correlation_matrix)
  # General indicator for cause
  c.ind <- 1 + rbinom(n, 1, prob = (1 - mix_p)^exp(X %*% beta1))
  
  #Conditional on cause indicators, we simulate the model.
  ftime <- numeric(n)
  eta1 <- X[c.ind == 1, ] %*% beta1 #linear predictor for cause on interest
  eta2 <- X[c.ind == 2, ] %*% beta2 #linear predictor for competing risk
  u1 <- runif(length(eta1))
  u2 <- runif(length(eta2))
  t1 <-  (- log(u1) / (lambda1 * exp(eta1)))^(1 / rho1)
  # Winsorising tiny values for t (smaller than one day on a yearly-scale, e.g. 1 / 365.242), and adding a tiny amount of white noise not to have too many concurrent values
  t1 <- ifelse(t1 < 1 / 365.242, 1 / 365.242, t1)
  t1[t1 == 1 / 365.242] <- t1[t1 == 1 / 365.242] + rnorm(length(t1[t1 == 1 / 365.242]), mean = 0, sd = 1e-4)
  # ...and make sure that the resulting value is positive
  t1 <- abs(t1)
  #t1 <- -log(1 - (1 - (1 - u1 * (1 - (1 - p)^exp(eta1)))^(1 / exp(eta1))) / p)
  #rexp(length(eta2), rate = exp(eta2))
  t2 <-  (- log(u2) / (lambda2 * exp(eta2)))^(1 / rho2)
  # Winsorising tiny values for t (smaller than one day on a yearly-scale, e.g. 1 / 365.242), and adding a tiny amount of white noise not to have too many concurrent values
  t2 <- ifelse(t2 < 1 / 365.242, 1 / 365.242, t2)
  t2[t2 == 1 / 365.242] <- t2[t2 == 1 / 365.242] + rnorm(length(t2[t2 == 1 / 365.242]), mean = 0, sd = 1e-4)
  # Fixing large values of t (to make sure it is within 1-year study range)
  t2 <- ifelse(t2 > max_time, max_time, t2)
  t2[t2 ==  max_time] <- t2[t2 == max_time] + rnorm(length(t2[t2 == max_time]), mean = 0, sd = 1e-4)
  t2 <- abs(t2)
  ftime[c.ind == 1] <- t1
  ftime[c.ind == 2] <- t2
 ci <- runif(n, min = 0, max = u.max)
ftime <- pmin(ftime, ci)
fstatus <- ifelse(ftime == ci, 0, 1)
fstatus <- fstatus * c.ind
# Fixing large values of t (to make sure it is within 1-year study range)
fstatus <- ifelse(ftime >= max_time, 0, fstatus)
ftime <- ifelse(ftime >= max_time, max_time, ftime)
ftime[ftime==   max_time] <- ftime[ftime == max_time] + 
  rnorm(length(max_time[max_time == 1]), mean = 0, sd = 1e-4)
  sim.data <- data.frame(fstatus = fstatus, ftime = ftime)
  X <- as.data.frame(X)
  colnames(X) <- paste0("X", seq_len(p))
  sim.data <- as.data.frame(cbind(sim.data, X))
  return(sim.data)
}
################## Prediction function for coxBoost ###########################
predictRisk.iCoxBoost <- function(object, newdata, times, cause,...){
  p <- predict(object, newdata= newdata,type="CIF",times= newdata$ftime)
 # if (is.null(dim(p))) {
#    if (length(p)!=length(times))
 #     stop("Prediction failed (wrong number of times)")
 # }
#  else{
#    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
#      stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
#  }
 p
}
##################### Prediction function for casebase penalized model ################
predict_CompRisk <- function(object, newdata = NULL) {
  ttob <- terms(object)
  X <- model.matrix(delete.response(ttob), newdata,
                    contrasts = if (length(object@contrasts)) {
                      object@contrasts
                    } else NULL,
                    xlev = object@xlevels)
  coeffs <- matrix(coef(object), nrow = ncol(X),
                   byrow = TRUE)
  preds <- X %*% coeffs
  colnames(preds) <- paste0("log(mu[,",
                            seq(2, length(object@typeEvents)),
                            "]/mu[,1])")
  return(preds)
}


########################## Test train validation function for mtool.multinom ########


################## Brier score prediction function for casebase #######################
predictRisk.CompRisk <- function(object, newdata, times, cause, ...) {
  #get all covariates excluding intercept and time
  coVars=colnames(object@originalData[, c(grepl("X", colnames(object@originalData)))])
  #coVars is used in lines 44 and 50
  newdata=data.matrix(drop(subset(newdata, select=coVars)))
  
  # if (missing(cause)) stop("Argument cause should be the event type for which we predict the absolute risk.")
  # the output of absoluteRisk is an array with dimension depending on the length of the requested times:
  # case 1: the number of time points is 1
  #         dim(array) =  (length(time), NROW(newdata), number of causes in the data)
  if (length(times) == 1) {
    a <- casebase::absoluteRisk.CompRisk(object, newdata = newdata, time = times)
    p <- matrix(a, ncol = 1)
  } else {
    # case 2 a) zero is included in the number of time points
    if (0 %in% times) {
      # dim(array) =  (length(time)+1, NROW(newdata)+1, number of causes in the data)
      a <- casebase::absoluteRisk.CompRisk(object, newdata = newdata, time = times)
      p <- t(a)
    } else {
      # case 2 b) zero is not included in the number of time points (but the absoluteRisk function adds it)
      a <- casebase::absoluteRisk.CompRisk(object, newdata = newdata, time = times)
      ### we need to invert the plot because, by default, we get cumulative incidence
      #a[, -c(1)] <- 1 - a[, -c(1)]
      ### we remove time 0 for everyone, and remove the time column
      a <- a[-c(1), -c(1)] ### a[-c(1), ] to keep times column, but remove time 0 probabilities
      # now we transpose the matrix because in riskRegression we work with number of
      # observations in rows and time points in columns
      p <- t(a)
    }
  }
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
    stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  }
  p
}

########### Function to compute weibull hazard #######################################
weibull_hazard <- Vectorize(function(gamma, lambda, t) {
  return(gamma * lambda * t^(gamma - 1))
})

############ Function to simulate from cause-specific hazards ############################
cause_hazards_sim <- function(p, n, beta1, beta2, hazard1 = 0.15, hazard2 = 0.10, 
                              nblocks = 4, cor_vals = c(0.7, 0.4, 0.6, 0.5), num.true = 20, h1 = 0.55, h2 = 0.10, 
                              gamma1 = 100, gamma2 = 100, max_time = 1.5, noise_cor = 0.1, 
                              rate_cens = 0.05, min_time = 0.002) {
  # Warnings
  if(length(beta1) != length(beta2)) stop("Dimension of beta1 and beta2 should be the same")
  # Set the number of variables per block
  vpb <- num.true/nblocks
  if(nblocks != length(cor_vals)) stop("Dimension of nblocks and correlations for blocks should match")
  # Set the correlation values for each covariate block
  correlation_values <- cor_vals
  # Initialize empty matrix
  correlation_matrix <- matrix(noise_cor, nrow = p, ncol = p)
  # Generate the covariance matrix with block correlations
  for (i in 1:nblocks) {
    start_index <- (i - 1) * vpb + 1
    end_index <- i * vpb
    correlation_matrix[start_index:end_index, start_index:end_index] <- correlation_values[i]
  }
  # Diagonal elements should be 1
  diag(correlation_matrix) <- rep(1, length(diag(correlation_matrix)))
  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = correlation_matrix)
  X <- as.matrix(X)
  # Specify rate parameters
  lambda1 <- h1 * exp(X %*% beta1)
  lambda2 <- h2 * exp(X %*% beta2)
  # Define cdf - U 
  cdf_U <- function(t, gamma1, lambda1, gamma2, lambda2, U) {
    F_min_U <- 1 - exp(-(lambda1 * t^gamma1 + lambda2 * t^gamma2)) - U
    return(F_min_U)
  }
  # Generate uniform values and store in dataframe
  u <- stats::runif(n)
  # Inverse transform sampling 
  dat_roots <- cbind.data.frame(u, gamma1, lambda1, gamma2, lambda2)
  times <- dat_roots %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      t_tilde = stats::uniroot(
        cdf_U,
        interval = c(.Machine$double.eps,
                     max_time),
        extendInt = "yes",
        U = u,
        gamma1 = gamma1,
        lambda1 = lambda1,
        gamma2 = gamma2,
        lambda2 = lambda2
      )$`root`
    ) %>%
    dplyr::pull(t_tilde)
  # Generate event indicators
  hazard1 <- weibull_hazard(gamma= gamma1, lambda = lambda1, t = times)
  hazard2 <- weibull_hazard(gamma = gamma2, lambda = lambda2, t = times)
  event <- stats::rbinom(n = n, size = 1, prob = hazard1 / (hazard1 + hazard2))
  c.ind <- ifelse(event == 1, 1, 2)
  # Add censoring
  cens <- stats::rexp(n = n, rate = rate_cens)
  c.ind <- ifelse(cens < times, 0, c.ind)
  times <- pmin(cens, times)
 # Winsorize time ranges to desired ones to make them more realistic and add some white noise 
   c.ind <- ifelse(times >= max_time, 0, c.ind)
   times <- ifelse(times >= max_time, max_time, times)
   times[times == max_time] <- times[times == max_time] + rnorm(length(times[times == max_time]), mean = 0, sd = 1e-4)
   times <- ifelse(times < min_time, min_time, times)
   times[times == min_time] <- times[times == min_time] + abs(rnorm(length(times[times == min_time]), mean = 0, sd = 1e-4))
   sim.data <- data.frame(fstatus = c.ind, ftime = times)
   X <- as.data.frame(X)
   colnames(X) <- paste0("X", seq_len(p))
   sim.data <- as.data.frame(cbind(sim.data, X))
   return(sim.data)
}



#' ###########################################################
#' #' Bootstrapping function for standard error of multinomial deviance
#' bootstrap_cvse <- function(x, y, alpha = alpha, B = B, lambdagrid, test_cv){
#' cvse <- replicate(B, {
#'   train_index <- sample(1:nrow(x), replace = TRUE)
#'   xboot <-  X[index, ]
#'   yboot <- Y[index]
#'   cvs_res <- mclapply(lambdagrid, function(lambda_val) {
#'     lambda1 <- lambda_val*alpha
#'     lambda2 <- 0.5*lambda_val*(1 - alpha)
#'     # mtool.MNlogistic is too verbose...
#'     mtool::mtool.MNlogistic(
#'       X = as.matrix(xboot),
#'       Y = yboot,
#'       offset = train_cv$offset,
#'       N_covariates = constant_covariates,
#'       regularization = 'elastic-net',
#'       transpose = FALSE,
#'       lambda1 = lambda1, lambda2 = lambda2, 
#'       lambda3 = 0)
#'   }, mc.cores = 4)
#' mult_deviance <- unlist(lapply(cvs_res, multi_deviance, cb_data = test_cv))
#' all_deviances[, i] <- mult_deviance
#' mean_dev <- rowMeans(all_deviances)
#' })
#' }
#' 
#' 
#' ############################################################################
#############################################################################
################  functions to fit and evaluate penCR model       ###########
#############################################################################
# from CSlassofunctionsV3.R

oneCSlasso = function(data, cause, lambdavec, ...){ 
  data = data.frame(data)
  vars = colnames(data)[(!colnames(data) %in% c('ftime', 'fstatus'))]
  X = as.matrix(data[, vars])
  y = Surv(data$ftime, data$fstatus == cause)
  glmnet.res = glmnet(x = X, y =y, alpha= 0.7, standardize=FALSE, 
                      lambda = lambdavec, ...)
  lp = lapply(lambdavec, function(s)
    as.numeric(predict(glmnet.res, newx=X, s=s, type="link")))
  out <- list('glmnet.res' = glmnet.res, 'vars' = vars, 'linear.predictor'=lp, 
              'response' = y)
  out$call <- match.call()
  class(out) <- "oneCSlasso"
  out
}


predictSurvProb.oneCSlasso <- function (object, newdata, times, lambdavec, index, ...){
  newx = data.frame(newdata)
  newx = as.matrix(newx[, object$vars])
  lp <- as.numeric(predict(object$glmnet.res, newx=newx, s=lambdavec[index], type="link"))
  bsurv <- basesurv(object$response, object$linear.predictor[[index]] %>% unlist, sort(unique(times)))$cumBaseHaz
  p  <- exp(exp(lp) %*% -t(bsurv))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))  stop("Prediction failed")
  p
}

twoCSlassos <- function(data, lambdavecs, ...){
  m1 = oneCSlasso(data = data, cause = 1, lambdavec = lambdavecs[[1]], ...)
  m2 = oneCSlasso(data = data, cause = 2, lambdavec = lambdavecs[[2]], ...)
  out <-  list('models' = list('Cause 1' = m1, 'Cause 2' = m2),
               'eventTimes' = m1$response[,'time'] %>% sort %>% unique,
               'causes' = c(1,2),
               'lambdas' = lambdavecs)
  out$call <- match.call()
  class(out) <- "twoCSlassos"
  out
}


predictEventProb.twoCSlassos = function (object, newdata, times, cause, lambdavecs, indices, ...){
  eTimes <- object$eventTimes
  pred = predictSurvProb(object$models[[paste("Cause", cause)]], times = eTimes,
                         newdata = newdata, lambdavec = lambdavecs[[cause]],
                         index = indices[cause])
  pred[pred<.000001] = .000001
  cumHaz1 <- -log(pred)
  Haz1 <- t(apply(cbind(0, cumHaz1), 1, diff))
  causes = object$causes
  cumHazOther <- lapply(causes[-match(cause, causes)], 
                        function(c) {
                          cumHaz.c <- -log(predictSurvProb(object$models[[paste("Cause", c)]],
                                                           times = eTimes, newdata = newdata,
                                                           lambdavec = lambdavecs[[c]], 
                                                           index = indices[c]))
                        })
  lagsurv <- exp(-cumHaz1 - Reduce("+", cumHazOther))
  cuminc1 <- t(apply(lagsurv * Haz1, 1, cumsum))
  pos <- prodlim::sindex(jump.times = eTimes, eval.times = times)
  p <- cbind(0, cuminc1)[, pos + 1, drop = FALSE]
  p
}


two.selpen.CSlassos <- function(object, ind){
  out <-  list('alllambdas'=object, 'indices' = ind)
  out$call <- match.call()
  class(out) <- "two.selpen.CSlassos"
  out
}

predictEventProb.two.selpen.CSlassos = function (object, newdata, times, cause, lambdavecs, ...){
  eTimes <- object$alllambdas$eventTimes
  pred = predictSurvProb(object$alllambdas$models[[paste("Cause", cause)]], times = eTimes,
                         newdata = newdata, lambdavec = lambdavecs[[cause]], 
                         index = object$indices[cause])
  pred[pred<.000001] = .000001
  cumHaz1 <- -log(pred)
  Haz1 <- t(apply(cbind(0, cumHaz1), 1, diff))
  causes = object$alllambdas$causes
  cumHazOther <- lapply(causes[-match(cause, causes)], 
                        function(c) {
                          cumHaz.c <- -log(predictSurvProb(object$alllambdas$models[[paste("Cause", c)]], times = eTimes, newdata = newdata, lambdavec = lambdavecs[[c]], 
                                                           index = object$indices[c]))
                        })
  lagsurv <- exp(-cumHaz1 - Reduce("+", cumHazOther))
  cuminc1 <- t(apply(lagsurv * Haz1, 1, cumsum))
  pos <- prodlim::sindex(jump.times = eTimes, eval.times = times)
  p <- cbind(0, cuminc1)[, pos + 1, drop = FALSE]
  p
}


####from cv.glmnetCR_v5_2.R

cv.glmnet.CR = function (data, causeOfInt, nfolds = 10,alpha=1, standardize=FALSE, 
                         t.BS = NULL, nlambda = 20, seed = 1, draw.heatmap = FALSE, ...) {
  set.seed(seed)
  data = data.frame(data)
  vars = colnames(data)[(!colnames(data) %in% c('ftime', 'fstatus'))]
  x = as.matrix(data[, vars])
  N = nrow(x)
  ID = seq_len(nrow(x))
  
  #penalized cause-specific model for each cause
  glmnet.object1 = glmnet(x = x, y=Surv(data$ftime, data$fstatus == 1), alpha=alpha,
                          standardize=standardize, nlambda = nlambda, ...)
  glmnet.object2 = glmnet(x, y=Surv(data$ftime, data$fstatus == 2), alpha=alpha,
                          standardize=standardize, nlambda = nlambda, ...)
  
  lambdavec = list(glmnet.object1$lambda, glmnet.object2$lambda)
  
  {sink("/dev/null"); foldid <- my.balanced.folds(data$fstatus, nfolds); sink(); }
  outlist = as.list(seq(nfolds))
  
  for (i in seq(nfolds)) {
    which = foldid == i
    data_sub = data[!which, ]
    
    outlist[[i]] =twoCSlassos(data = data_sub, lambdavec = lambdavec, ...)
  }

# helper function to make prediction vector for all patients at fixed 
# index(lambda1, lambda2)
makepred = function(outlist, nfolds, foldid, index){
  p = vector('list', nfolds)
  # make prediction on each fold
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    test =  data[which,]
    
    p[[i]] = predictEventProb.twoCSlassos(fitobj, newdata=test, times = t.BS, 
                                          cause = causeOfInt,lambdavecs= lambdavec,
                                          indices = index) %>% as.numeric 
    names(p[[i]]) = ID[which]
  }    
  #reorder predictions to fit with original data
  p = unlist(p) 
  p = p[order(as.numeric(names(p)))]  %>% matrix
  return(p)
}

index = expand.grid('ind.lambda1'=1:length(lambdavec[[1]]), 
                    'ind.lambda2'=1:length(lambdavec[[2]]))
index = as.list(data.frame(t(index)))

# make list of prediction over list of all indices(lambda1, lambda2)
listofpred = lapply(index, function(i) makepred(outlist, nfolds, foldid, ind =i))


# Brier score for all combinations of lambdas
temp = pec(listofpred, formula=Hist(ftime,fstatus)~1, data=data, cause = causeOfInt,
           times=t.BS, start=NULL,exact=FALSE, cens.model = 'marginal')
BS = array(data = unlist(temp$AppErr[-1]),
           dim = c(length(lambdavec[[1]]), length(lambdavec[[2]])),
           dimnames = list(paste('l1=', round(lambdavec[[1]],5)),
                           paste('l2=', round(lambdavec[[2]],5))))

if(draw.heatmap) pheatmap(BS, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
# pheatmap(R2, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)

glmnet.fits = twoCSlassos(data = data, lambdavec = lambdavec, ...)
min.index = which(BS == min(BS), arr.ind = TRUE) %>% as.numeric

best.fit = two.selpen.CSlassos(glmnet.fits, ind = min.index)
out <- list('glmnet.fits' = glmnet.fits, 'min.index'=min.index, 'best.fit'=best.fit,
            'BS' = BS, 'lambdavec' = lambdavec)

out$call <- match.call()
class(out) <- "cv.glmnet.CR"
out
}


# predictEventProb-function needed for pec
predictEventProb.cv.glmnet.CR  = function (object, newdata, times, cause, ...){
  eTimes <- object$best.fit$alllambdas$eventTimes
  pred = predictSurvProb(object$best.fit$alllambdas$models[[paste("Cause", cause)]], times = eTimes, 
                         newdata = newdata, lambdavec = object$lambdavec[[cause]], 
                         index = object$min.index[cause])
  pred[pred<.000001] = .000001
  cumHaz1 <- -log(pred)
  Haz1 <- t(apply(cbind(0, cumHaz1), 1, diff))
  causes = object$best.fit$alllambdas$causes
  cumHazOther <- lapply(causes[-match(cause, causes)], 
                        function(c) {
                          cumHaz.c <- -log(predictSurvProb(object$best.fit$alllambdas$models[[paste("Cause", c)]], times = eTimes, newdata = newdata, lambdavec = object$lambdavec[[c]], index = object$min.index[c]))
                        })
  lagsurv <- exp(-cumHaz1 - Reduce("+", cumHazOther))
  cuminc1 <- t(apply(lagsurv * Haz1, 1, cumsum))
  pos <- prodlim::sindex(jump.times = eTimes, eval.times = times)
  p <- cbind(0, cuminc1)[, pos + 1, drop = FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))  stop("Prediction failed")
  p
}


basesurv <- function (response, lp, times.eval = NULL, centered = FALSE)
{
  if (is.null(times.eval)) times.eval <- sort(unique(response[,1]))
  
  t.unique <- sort(unique(response[,1][response[,2] == 1]))
  alpha    <- length(t.unique)
  
  for (i in 1:length(t.unique)) {
    alpha[i] <- sum(response[,1][response[,2] == 1] == t.unique[i])/sum(exp(lp[response[,1] >=  t.unique[i]]))
  }
  
  obj   <- approx(t.unique, cumsum(alpha), yleft=0, xout = times.eval, rule=2)
  
  if (centered) obj$y <- obj$y * exp(mean(lp))
  obj$z <- exp(-obj$y)
  
  names(obj) <- c("times","cumBaseHaz","BaseSurv")
  return(obj)
}

#from my.balanced.folds.R

my.balanced.folds <- function(class.column.factor, cross.outer){
  # get balanced folds from pamr
  sampleOfFolds  <- get("balanced.folds",envir=asNamespace("pamr"))(class.column.factor,
                                                                    nfolds=cross.outer)
  permutated.cut <- rep(0,length(class.column.factor))
  for (sample in 1:cross.outer){
    #cat(sample,"\n")
    permutated.cut[sampleOfFolds[[sample]]] <- sample
  }
  invisible(permutated.cut)
}

###################### Generate data from Proportional Sub-distribution hazard models ######
# from gen.Prop.SDH.dataV3.R

gen.Prop.SDH.data <- function(n, p, gamma, binary.or.cont = 'cont', 
                              seed, cens,  q = 0.5, cens.param = c(0, 5)) {
  # n = number of subjects to simulate (numeric)
  # p = number of variables (numeric)
  # gamma = matrix of true regression coefficients 
  # cont.or.bin = generate continuous or binary variables ('cont' or 'bin')
  # seed = random number generator seed (numeric)
  # cens = include censoring (logical)
  # q = between 0 and 1, probability for patient with all covariates = 0 to
  #     experience event of type 1
  # cens.param = parameteres used in uniform censoring distribution
  
  set.seed(seed)
  
  ## Draw the covariate
  if(binary.or.cont == 'cont') X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  if(binary.or.cont == 'binary') X <- matrix(rbinom(n*p,1,q), nrow = n, ncol = p)
  
  colnames(X) = c(paste('r', 1:sum(gamma[,1] !=0), sep = ''), 
                  paste('u', 1:sum(gamma[,1] ==0), sep = '')) 
  #call all variables that have effect on CIF 1 (i.e. first column of gamma !=0) 'related'. 
  
  
  ## Binomial experiment to decide which event happens
  prob.event1 <- 1 - (1 - q)^exp(X %*% gamma[,1])
  event <- rbinom(n, 1, prob.event1)
  event <- ifelse(event == 0, 2, event)
  
  ## More or less the CDF conditional on event types. See e.g. Fine and Gray paper.
  ## -u is for numeric inversion to get the event times
  ## incr is the increment in the loop below
  cdf1 <- function(t, u, incr) {
    ((1 - (1 - q * (1 - exp(-t)))^exp(X[incr,]%*% gamma[,1])) /
       prob.event1[incr,]) - u
  }
  
  cdf2 <- function(t, u, incr) {
    (((1 - q)^exp(X[incr,] %*% gamma[,1]) * (1 - exp(-t * exp(X[incr,]%*% gamma[,2])))) /
       (1 - prob.event1[incr,])) - u
  }
  
  time.wo.c <- numeric(n)
  ## Us is the uniform drawings, for inverse transform sampling
  Us <- runif(n)
  for (i in seq_len(n)) {
    time.wo.c[i] <- switch(as.character(event[i]),
                           "1" = {
                             uniroot(cdf1, c(0, 1000000), tol = 0.0001,
                                     u = Us[i], incr = i)$root
                           },
                           "2" = {
                             uniroot(cdf2, c(0, 1000000), tol = 0.0001,
                                     u = Us[i], incr = i)$root
                           })
  }
  
  time.wo.c[time.wo.c<=0] = runif(sum(time.wo.c==0),0,.05)
  if(cens) {
    cens <- runif(n, cens.param[1], cens.param[2])
    time <- pmin(cens, time.wo.c)
    status <- ifelse(time == cens, 0, event)
    out = data.frame(time, status, X)
  } else out = data.frame('time' = time.wo.c, 'status' = event, X)
  out
}

