########################### n = 400, p = 120, Tp = 20  ######################
library(casebase)
library(future.apply)
library(glmnet)
library(devtools)
# install.packages("mtool_1.0.tar.gz")
library(mtool)
library(parallel)
library(tictoc)
library(tidyverse)
library(foreach)
library(survival)
library(cmprsk)
library(glue)
library(pec)
library(pamr)

# Fitting functions 
source("src/fitting_functions.R")

runCasebaseSim = function(n = 400, p = 20, N = 5) {
  # For each simulation you want to output the MSE
  n = 400
  p = 20
  N = 5
  Results <- replicate(N, {
  # Set seed
  seed <- as.integer(Sys.time())
  
  # take the last five digits of the initial seed
  the_seed= seed %% 100000
  set.seed(the_seed)
  
  num_true <- 20
  beta1 <- c(rep(0, p))
  beta2 <- c(rep(0, p))
  nu_ind <- seq(num_true)
  # Here out of 20 predictors, 10 should be non-zero 
  beta1[nu_ind] <- c(rep(1, 10), rep(0, 10))
  beta2[nu_ind] <- c(rep(-1, 10), rep(0, 10))
  
  # Simulate data
  sim.data <- cause_hazards_sim(n = n, p = p, nblocks = 4, 
                                beta1 = beta1, beta2 = beta2, rate_cens = 0.25, 
                                h1 = 0.55, h2 = 0.10, gamma1 = 1.5, gamma2 = 1.5)
  
  
  # Censoring proportion
  cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)
  
  # Training-test split 
  # We only do this (instead of generating datasets for train and test like Anthony mentioned because it is faster computationally 
  # as casebase resamples) + proportion of censoring can be quite random in each run of the simulation so we want to maintain the same in validation and test set
  train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
  train <- sim.data[train.index,]
  test <- sim.data[-train.index,]
  
  ##############################################################
  # We have two competitor models for variable selection:
  # 1) Independent cox-regression model 
  # 2) penCR cox regression model - where the lambda penalties are trained together 
  ######################## Fit indepedent cox-regression model ###############################
  ######################### Cause-1 #########################################
  # Censor competing event
  y_train <- Surv(time = train$ftime, event = train$fstatus == 1)
  
  x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
  
  # Censor competing event
  y_test <- Surv(time = test$ftime, event = test$fstatus == 1)
  
  x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
  
  # Fit cause-specific cox model with glmnet on training set 
  cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7)
  
  # Fit on validation set 
  cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                        lambda = cox_mod$lambda.min)
  
  cc_min <- coef(cox_val_min)
  
  res_cox_min1 <- varsel_perc(cc_min, beta1)
  
  # let's calculate the bias for all the competitors as well (a task could be turning this one line into a function as well)
  # Only for the true non-zero variables
  mean((cc_min[nu_ind] - beta1[nu_ind])^2)
  ########################## Cause 2 #####################################
  # Censor competing event
  y_train <- Surv(time = train$ftime, event = train$fstatus == 2)
  
  x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
  
  # Censor competing event
  y_test <- Surv(time = test$ftime, event = test$fstatus == 2)
  
  x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
  
  # Fit cause-specific cox model with glmnet on training set 
  cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7)
  
  # Fit on validation set 
  cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                        lambda = cox_mod$lambda.min)
  
  cc_min <- coef(cox_val_min)
  
  # Function to output variable selection performance metrics
  res_cox_min2 <- varsel_perc(cc_min, beta2)
  
  # let's calculate the bias for all the competitors as well (a task could be turning this one line into a function as well)
  cc_min_bias <- which(coef(cox_val_min) != 0)
  
  # MSE (bias)
  mean((cc_min[nu_ind] - beta2[nu_ind])^2)
  ########################## Fit PenCR model ##################################
  penCR = cv.glmnet.CR(data = train, family="cox", alpha= 0.7, standardize= TRUE,
                       nlambda = 20, t.BS = median(train$ftime), seed = 115, causeOfInt = 1,
                       nfold = 5)
  
  cc_min_penCR1 <- penCR$glmnet.fits$models$`Cause 1`$glmnet.res$lambda[penCR$min.index[1]]
  cc_min_penCR2 <- penCR$glmnet.fits$models$`Cause 2`$glmnet.res$lambda[penCR$min.index[2]]
  
  # Fit on validation set 
  penCR_val_min1 <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                           lambda = cc_min_penCR1)
  
  cc_min_penCR1 <- coef(penCR_val_min1)
  
  penCR_val_min2 <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                           lambda = cc_min_penCR2)
  
  
  cc_min_penCR2 <- coef(penCR_val_min2)
  
  res_pencr_min1 <- varsel_perc(cc_min_penCR1, beta1)
  res_pencr_min2 <- varsel_perc(cc_min_penCR2, beta2)
  
  # Free up memory for case-base
  rm(penCR)
  
  # Calculate MSE here as well (try and fill it out!)
  # MSE Cox model for cause of interest
  mean((cc_min_penCR1[nu_ind] - beta1[nu_ind])^2)
  # MSE Cox model for competing risk
  cc_min_bias_pen <- which(cc_min_penCR2 != 0)
  mean((cc_min_penCR2[nu_ind] - beta2[nu_ind])^2)
  
  ########################## Fit casebase model #################################
  # Test set 
  surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
  
  # Covariance matrix
  cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
  
  # Case-base dataset
  cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)
  
  # Train case-base model 
  cv.lambda <- mtool.multinom.cv(train, seed = 1, nfold = 5)
  
  # Case-base fits 
  # Lambda.min
  fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                             lambda = cv.lambda$lambda.min , alpha = 0.7, unpen_cov = 2)
  
  
  res_cb_min1 <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)
  
  res_cb_min2 <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 2], beta2)
  
  # Calculate MSE here as well!
  # MSE for casebase model for cause of interest
  fit_val_coef_1 <- fit_val_min$coefficients[1:eval(parse(text="p")), 1]
  mean((beta1[nu_ind] - fit_val_coef_1[nu_ind])^2)
  
  # MSE for casebase model for competing risk
  fit_val_coef_2 <- fit_val_min$coefficients[1:eval(parse(text="p")), 2]
  cb_min_bias <- which(fit_val_coef_2 != 0)
  mean((fit_val_coef_2[nu_ind] - beta2[nu_ind])^2)
  
  
  ########################## Fit casebase model with post LASSO#################################
  
  res <- multinom.post_enet(train, test)
  
  
  # Calculate MSE for this as well (fill in here)
  mean((res$coefficients[nu_ind]- beta1[nu_ind])^2)
  mean((res$coefficients[nu_ind] - beta2[nu_ind])^2)
  
  res_cb_post_lasso <- varsel_perc(res$coefficients, beta1)
  
  ###################################################################################
  Res <- rbind(res_cb_post_lasso, res_cb_min1, res_cb_min2, res_cox_min1, res_cox_min2, res_pencr_min1, res_pencr_min2, cen.prop)
  
  rownames(Res) <- c("casebase.post.lasso.lambda.min", "casebase.lambda.min_cause1", "casebase.lambda.min_cause2", "cox.lambda.min_cause1",
                     "cox.lambda.min_cause2", "pencr.lambda.mincause1", "pencr.lambda.mincause2", "cens.prop")
  
  Res
  
  }, simplify = FALSE)
  
  Results <- do.call(rbind, Results)
  
  Results
}

# res1 <- multinom.post_enet(train, test)
# res_cb_post_lasso <- varsel_perc(res$coefficients, beta1)

results_table = runCasebaseSim(400, 20, 5)




formatCaseBaseTable = function(sim_results) {
  # sim_results = results_table
  
  model_fit_labels = c("casebase.post.lasso.lambda.min", "casebase.lambda.min_cause1", "casebase.lambda.min_cause2", "cox.lambda.min_cause1",
                       "cox.lambda.min_cause2", "pencr.lambda.mincause1", "pencr.lambda.mincause2", "cens.prop")
  
  stat_names = colnames(as.data.frame(sim_results))
  
  num_models = length(model_fit_labels)
  rows = nrow(sim_results)
  # Add model name to model column
  model_list = c()
  for (i in c(1:rows)) {
    index = i %% num_models
    if (i %% num_models == 0)
      index = num_models
    curr_model = model_fit_labels[index]
    model_list = c(model_list, curr_model)
  }
  
  sim_results["Model"] = model_list
  
  formatted_table = data.frame(matrix(nrow = 0, ncol = length(sim_results)))
  col_names = colnames(as.data.frame(sim_results))
  colnames(formatted_table) = col_names

  for (i in model_fit_labels) {
    stat_table = sim_results[stat_names]
    average_stat_row = colMeans((sim_results %>% filter(Model == i))[1:ncol(sim_results) - 1])
    formatted_table = rbind(formatted_table, c(average_stat_row, i))
  }
  colnames(formatted_table) = col_names
  formatted_table = formatted_table %>% replace(is.na(.), NaN)  
  # formatted_table[, ncol(formatted_table) - 1] = lapply(formatted_table[, ncol(formatted_table) - 1], as.numeric)
  
  options(digits = 5)
  
  # formatted_table[, ncol(formatted_table) - 1] = round_df(formatted_table[, ncol(formatted_table) - 1], 3)
  
  # Make model column the first column of df
  formatted_table = formatted_table[,c(ncol(formatted_table), seq(1:ncol(formatted_table) - 1))]
  formatted_table = formatted_table[1:ncol(formatted_table) - 1]
  
  return(formatted_table)
}

summarizedSimCaseBaseTable = formatCaseBaseTable(results_table)






############ Sketch of function for post-LASSO (or post elastic net in this case) #########
# Look into ... argument to pass parameters from other functions because you want to pass cross-validation parameters
# multinom.post_enet <- function(train, test) {
#   # Train case-base model to get lambda.min
#   cv.lambda <- mtool.multinom.cv(train, seed = 1, nfold = 5)
#   # This fit (with lambda.min) needs to be de-biased
#   # Fit on test set 
#   # Covariance matrix
#   cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
#   
#   # Case-base dataset
#   surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
#   cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)
#   
#   # Case-base fits 
#   # Lambda.min
#   fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
#                              lambda = cv.lambda$lambda.min, alpha = 0.7, unpen_cov = 2)
#   
#   # Obtain all non-zero selected covariates across both classes
#   non_zero_coefs_cause1 <- which(fit_val_min$coefficients[1:eval(parse(text="p")), 1] != 0)
#   non_zero_coefs_cause2 <- which(fit_val_min$coefficients[1:eval(parse(text="p")), 2] != 0)
#   # Combine them 
#   non_zero_coefs <- union(non_zero_coefs_cause1, non_zero_coefs_cause2)
#   non_zero_coefs <- paste("X", non_zero_coefs , sep = "")
#   # Create new subsetted dataset 
#   testnew <- cbind(test[, c(colnames(test) %in% non_zero_coefs)], ftime = (test$ftime), fstatus = test$fstatus)
#   # Fit "OLS" (unparameterized multinomial model)
#   # For working of this function see: http://sahirbhatnagar.com/casebase/articles/competingRisk.html
#   model_cb <- fitSmoothHazard(fstatus ~. +log(ftime) -fstatus,
#                               data = testnew,
#                               ratio = 100,
#                               time = "ftime")
#   
#   # Only return estimated coefficients of covariates
#   exclude_coefs = c("(Intercept):1", "(Intercept):2", "ftime:1", "ftime:2", 
#                     "log(ftime):1", "log(ftime):2")
#   est_betas = coef(model_cb)[names(coef(model_cb))
#                              %in% exclude_coefs == FALSE]
#   all_coef_names = 
#   
#   for (l in c(1:length(coefs_all_lambdas))) {
#     if(is.na(fit_OLS_on_LASSO_subset$coefficients[j]))
#       break
#     if(colnames(coefs_all_lambdas)[l] == names(fit_OLS_on_LASSO_subset$coefficients[j])) {
#       coefs_all_lambdas[i, l] = fit_OLS_on_LASSO_subset$coefficients[j]
#       j = j + 1
#     }
#   }
#   
#   est_betas = est_betas[1:20]
#   
#   res <- list(coefficients = est_betas, lambda.min = cv.lambda$lambda.min,
#               lambdagrid = cv.lambda$lambdagrid)
#   
#   res
# }


################### Post-LASSO function #################################
############ Sketch of function for post-LASSO (or post elastic net in this case) #########
# Look into ... argument to pass parameters from other functions because you want to pass cross-validation parameters
multinom.post_enet <- function(fit_object, cause = 1) {
  # Obtain all non-zero selected covariates from cause 1 
  coef <- fit_val_min$coefficients[1:eval(parse(text="p")), 1]
  non_zero_coefs_cause1 <- which(fit_val_min$coefficients[1:eval(parse(text="p")), 1] != 0)
  non_zero_coefs <- paste("X", non_zero_coefs_cause1, sep = "")
  # Create new subsetted dataset 
  testnew <- cbind(test[, c(colnames(test) %in% non_zero_coefs)], ftime = (test$ftime), fstatus = test$fstatus)
  # Fit "OLS" (unparameterized multinomial model)
  model_cb <- fitSmoothHazard(fstatus ~. +log(ftime) -ftime -fstatus,
                              data = testnew,
                              ratio = 100, time = "ftime")
  coeffs <- matrix(coef(model_cb), nrow =  length(coef(model_cb))/2, byrow = TRUE)
  rownames(coeffs) <- c(colnames(testnew)[-length(colnames(testnew))], "Intercept")
  coef[non_zero_coefs_cause1] <- coeffs[1:length( non_zero_coefs_cause1)]
  return(list(coef_selected = coeffs, non_zero_coefs = non_zero_coefs_cause1, coefs_all = coef))
}




###########################################################
#' Cross-validation function for mtool 
multinom.relaxed_enet <- function(train, regularization = 'elastic-net', lambda_max = NULL, alpha = 1, nfold = 10, 
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
  lowest_deviance = .Machine$double.xmax
  best_fit = NULL
  min_lambda_index = 0
  for (j in lambdagrid) {
    curr_lambda_deviance = 0
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
      mult_deviance <- unlist(lapply(cvs_res, multi_deviance, cb_data = test_cv))
      all_deviances[, i] <- mult_deviance
      non_zero_coefs[, i] <-  unlist(lapply(cvs_res, function(x) {return(x$no_non_zero)}))
      cat("Completed Fold", i, "\n")
    }
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
              lambda.1se = lambda.1se, lambda.0.5se = lambda.0.5se, cv.se = cv_se, lambdagrid = lambdagrid, deviance_grid = all_deviances))
}
