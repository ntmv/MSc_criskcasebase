############### All fitting functions for multinomial casebase simulations ############################
#' Survsim's crisk.sim function modified for generating competing risks data with covariates from a multivariate normal distribution
#' #' 
#' @param p number of covariates 
#' @param rho correlation in AR1 structure for covariates
#' Other params in survsim::crisk.sim
#' @return simulated competing risks dataset
crisk.sim_mvn <-
  function (n, p, rho, foltime, dist.ev, anc.ev, beta0.ev, dist.cens = "weibull", 
            anc.cens, beta0.cens, z = NULL, beta = NA, nsit) 
  {
    # Arguments check
    if (length(anc.ev) != length(beta0.ev)) stop("Wrong number of parameters")
    if (length(anc.cens) != length(beta0.cens) && length(anc.cens) != 1) stop("Wrong number of parameters")
    if (length(anc.ev) != length(dist.ev)) stop("Wrong number of parameters")
    #if ( all(is.na(beta))) stop("Wrong specification of covariables")
    #if (any(!is.na(beta))) stop("Wrong specification of covariables")
    if (length(beta0.ev) != nsit) stop("Wrong number of distributions or events")
    if (!is.null(z) && length(z) != nsit && length(z) != 1) stop("Wrong numbers of elements in z")
    if (!is.null(z) && !all(lapply(z, function(x) x[1]) %in% c("unif","weibull","invgauss", "gamma","exp"))) 
      stop("Wrong specification of z")
    if (!is.null(z) && any(lapply(z, function(x) length(x)) != 3))
    {
      if(any(lapply(z[lapply(z, function(x) length(x)) != 3], function(x) length(x)) != 2)) stop("Wrong specification of z")
      if(any(lapply(z[lapply(z, function(x) length(x)) != 3], function(x) length(x)) == 2))
      {
        for (i in 1:length(z[lapply(z, function(x) length(x)) == 2]))
        {
          if (z[lapply(z, function(x) length(x)) != 3][[i]][1] != "exp") stop("Wrong specification of z")
        }
      }
    }
    if(!is.null(z) && any(lapply(z, function(x) x[1]) == "unif"))
    {
      for (i in 1:length(z[lapply(z, function(x) x[1]) == "unif"]))
      {
        if (as.numeric(z[lapply(z, function(x) x[1]) == "unif"][[i]][2])-as.numeric(z[lapply(z, function(x) x[1]) == "unif"][[i]][3]) >= 0) 
          stop("Wrong specification of z")
        if (as.numeric(z[lapply(z, function(x) x[1]) == "unif"][[i]][2]) < 0) stop("Wrong specification of z")
      }
    }
    
    sim.data <- list()
    eff <- vector()
    eff[1] <- 0
    Sigma <- covAR1(p = p, rho = rho)
    for (i in 1:n) {
      eff <- mvtnorm::rmvnorm(1, sigma = Sigma)
      sim.data[[i]] <- survsim::crisk.ncens.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z, beta, eff, 
                                                dist.ev, dist.cens, i, nsit)
    }
    sim.data <- do.call(rbind, sim.data)
    sim.data$cause[sim.data$status==0] <- NA
    class(sim.data) <- c("crisk.data.sim", "data.frame")
    attr(sim.data, "n") <- n
    attr(sim.data, "foltime") <- foltime
    attr(sim.data, "nsit") <- nsit
    return(sim.data)
  }

#########################################################
#' Multinomial deviance
#' 
#' @param cb_data Output of \code{create_cbDataset}
#' @param fit_object Output of \code{fit_cbmodel}
#' @return Multinomial deviance
multi_deviance <- function(cb_data, fit_object) {
  X <- as.matrix(cbind(cb_data$time, cb_data$covariates))
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
create_cbDataset <- function(surv_obj, cov_matrix, ratio = 10) {
  n <- nrow(surv_obj)
  B <- sum(surv_obj[, "time"])
  c <- sum(surv_obj[, "status"] != 0)
  b <- ratio * c
  offset <- log(B/b)
  prob_select <- surv_obj[, "time"]/B
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
  # Prepare covariate matrix
  X <- as.matrix(cb_data$covariates)
  # mtool.MNlogistic is too verbose...
  out <- fit.mtool <- mtool::mtool.MNlogistic(
    X = as.matrix(cbind(X, 1)),
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
                              constant_covariates = 2, initial_max_grid = NULL, precision = 0.001, epsilon = .0001, grid_size = 100, plot = FALSE) {
  # Default lambda grid
  if(is.null(lambda_max)) {
    # Lambda max grid for bisection search 
    if(is.null(initial_max_grid)) {
      initial_max_grid <- 
        c(0.9, 0.5, 0.1, 0.07, 0.05, 0.01, 0.009, 0.005)
      fit_val_max <- mclapply(initial_max_grid, 
                              function(lambda_val) {
                                fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                                            lambda = lambda_val, alpha = alpha)}, mc.cores = 4)
      non_zero_coefs <-  unlist(mclapply(fit_val_max, function(x) {return(x$no_non_zero)}), mc.cores = 4)
      if(!isTRUE(any(non_zero_coefs == (constant_covariates*2)))){
        warning("Non-zero coef value not found in default grid. Re-run function and specify initial grid")
      }
      upper <- initial_max_grid[which(non_zero_coefs > (constant_covariates*2 + 1))[1]-1]
      lower <- initial_max_grid[which(non_zero_coefs > (constant_covariates*2 + 1))[1]]
      new_max_searchgrid <- seq(lower, upper, precision)
      fit_val_max <-  lapply(new_max_searchgrid, 
                               function(lambda_val) {
                                 fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                                             lambda = lambda_val, alpha = alpha)})
      non_zero_coefs <-  unlist(mclapply(fit_val_max, function(x) {return(x$no_non_zero)}, mc.cores = 4))
      lambda_max <- new_max_searchgrid[which.min(non_zero_coefs)]
    }
  }
  lambdagrid <- rev(round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = grid_size)), digits = 10))
  cb_data_train <- as.data.frame(cb_data_train)
  # Create folds 
  folds <- caret::createFolds(factor(cb_data_train$event_ind), k = nfold, list = FALSE)
  lambda.min <- rep(NA_real_, nfold)
  all_deviances <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  non_zero_coefs <- matrix(NA_real_, nrow = length(lambdagrid), ncol = nfold)
  
  # Set the number of cores to be used for parallel processing
 # num_cores <- parallelly::availableCores()  # Adjust the number of cores as per your system's capacity
  # Create a parallel cluster using the specified number of cores
  #cl <-  parallel::makeCluster(num_cores, setup_strategy = "sequential")
  
  # Load necessary packages on each cluster
 #clusterEvalQ(cl, {
 #  library(casebase)
 #  library(future.apply)
 #  library(mtool)
 #  library(parallel)
 #  library(dplyr)
 #})
 # 
  #clusterEvalQ(cl, source("../src/fitting_functions_nonparallel.R"))
  
  # Export multiple objects to all clusters
  #objects_to_export <- c("train_cv", "lambdagrid", "test_cv", "alpha")
 # clusterExport(cl, objects_to_export)
  
  #Perform 10 fold cross validation
  for (i in 1:nfold) {
    #Segment your data by fold using the which() function 
    train_cv <- cb_data_train[which(folds != i), ] #Set the training set
    test_cv <- cb_data_train[which(folds == i), ] #Set the validation set
    # Create X and Y
    train_cv <- list("time" = train_cv$time,
                     "event_ind" = train_cv$event_ind,
                     "covariates" = train_cv[, grepl("covariates", names(train_cv))],
                     "offset" = train_cv$offset)
    
    # Define the function to be applied in parallel
    cv_res <- mclapply(lambdagrid, function(lambda_val) {
      fit_cbmodel(train_cv, regularization = 'elastic-net', lambda = lambda_val, alpha = alpha)
    }, mc.cores = 4)
    
    test_cv <- list("time" = test_cv$time,
                    "event_ind" = test_cv$event_ind,
                    "covariates" = test_cv[, grepl("covariates", names(test_cv))],
                    "offset" = test_cv$offset)
    mult_deviance <- unlist(mclapply(cv_res, multi_deviance, cb_data = test_cv, mc.cores = 4))
    all_deviances[, i] <- mult_deviance
    mean_dev <- rowMeans(all_deviances)
    lambda.min <- lambdagrid[which.min(mean_dev)]
    cv_se <- sqrt(var(mean_dev)/nfold)
    dev.1se <- mean_dev[which.min(mean_dev)] + cv_se
    dev.0.5se <- mean_dev[which.min(mean_dev)] + cv_se/2
    range.1se <- lambdagrid[which(mean_dev <= dev.1se)]
    lambda.1se <- tail(range.1se, n = 1)
    lambda.min1se <- range.1se[1]
    range.0.5se <- lambdagrid[which((mean_dev <= dev.0.5se))]
    lambda.0.5se <- tail(range.0.5se, n = 1)
    lambda.min0.5se <- range.0.5se[1]
    non_zero_coefs[, i] <-  unlist(lapply(cv_res, function(x) {return(x$no_non_zero)}))
    rownames(all_deviances) <- lambdagrid
    rownames(non_zero_coefs) <- lambdagrid
  }
  return(list(lambda.min = lambda.min,  non_zero_coefs = non_zero_coefs, lambda.min1se = lambda.min1se, lambda.min0.5se = lambda.min0.5se, 
              lambda.1se = lambda.1se, lambda.0.5se = lambda.0.5se, cv.se = cv_se, lambdagrid = lambdagrid, deviance_grid = all_deviances))
}

#' ###########################################################
#' Plotting function for cross-validation 
plot_cv.multinom <- function(all_deviances, lambdagrid, lambda.min, lambda.1se, nfold = 10) {
  mean_dev <- rowMeans(all_deviances)
  row_stdev <- apply(all_deviances, 1, function(x) {sd(x)/sqrt(nfold)})
  plot.dat <- data.frame(lambdagrid = lambdagrid, mean.dev = mean_dev, 
                         upper = mean_dev +row_stdev, lower = mean_dev - row_stdev)
  p <- ggplot(plot.dat, aes(log(lambdagrid), mean.dev)) + geom_point(colour = "red", size = 3) + theme_bw() +  geom_errorbar(aes(ymin= lower, ymax=upper), width=.2, position=position_dodge(0.05), colour = "grey") + labs(x = "log(lambda)", y = "Multinomial Deviance")  + geom_vline(xintercept = log(lambda.min), linetype = "dotted", colour = "blue") + 
    geom_vline(xintercept = log(lambda.1se), linetype = "dotted", colour = "purple")
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
  mcc_den <- (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  mcc <- mcc_num/sqrt(mcc_den)
  dat <- as.data.frame(cbind(TP = TP, TN = TN, FP = FP, FN = FN, Sensitivity = sens, Specificity = spec, MCC = mcc))
  return(dat)
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