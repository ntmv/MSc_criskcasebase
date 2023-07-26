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

regularization = 'elastic-net'; lambda_max = NULL; alpha = 1; nfold = 10; 
constant_covariates = 2; initial_max_grid = NULL; precision = 0.001; epsilon = .0001; grid_size = 100; plot = FALSE; 
ncores = parallelly::availableCores(); seed = NULL; train_ratio = 20


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
  
  current_deviance = .Machine$double.xmax
  best_fit = NULL
  min_lambda_index = 0

  #Perform 10 fold cross validation
  for(i in 1:nfold){
    #Segment your data by fold using the which() function 
    train_cv_cb <- cb_data_train[which(folds != i), ] #Set the training set
    # Create X and Y
    train_cv_cb <- list("time" = train_cv_cb$time,
                     "event_ind" = train_cv_cb$event_ind,
                     "covariates" = train_cv_cb[, grepl("covariates", names(train_cv_cb))],
                     "offset" = train_cv_cb$offset)
    # Standardize
    train_cv_cb$covariates <- as.data.frame(scale(train_cv_cb$covariates, center = T, scale = T))
    cvs_res <- mclapply(lambdagrid, 
                        function(lambda_val) {
                          fit_cbmodel(train_cv_cb, regularization = 'elastic-net',
                                      lambda = lambda_val, alpha = alpha, unpen_cov = constant_covariates)
                        }, mc.cores = ncores, mc.set.seed = seed)
    
    
    cvs_res_unlisted = unlist(lapply(cvs_res, function(x) {return(x$coefficients)}))
    
    # Obtain all non-zero selected covariates across both causes
    cvs_non_zero_coefs_cause1 <- mclapply(cvs_res_unlisted,
                        function(cb_fit_obj) {
                          which(cb_fit_obj[1:eval(parse(text="p")), 1] != 0)
                        }, mc.cores = ncores, mc.set.seed = seed)


    cvs_non_zero_coefs_cause2 <- mclapply(cvs_res_unlisted,
                                          function(cb_fit_obj) {
                                            which(cb_fit_obj[1:eval(parse(text="p")), 2] != 0)
                                          }, mc.cores = ncores, mc.set.seed = seed)
    
    train_unparameterized_cv <- cb_data_train[which(folds != i), ] #Set the training set
    # covariate_names = lapply(colnames(train_unparameterized_cv)[2:22], function(x){substring(text = x, first = 12)})
    # colnames(train_unparameterized_cv) = c(colnames(train_unparameterized_cv)[1], covariate_names, colnames(train_unparameterized_cv)[23])
    
    test_unparameterized_cv <- cb_data_train[which(folds == i), ] #Set the validation set
    
    # unparameterized_hazards = mclapply(cvs_non_zero_coefs_cause1, 
    #                                  function(list_cause1, list_cause2) {
    #                                    non_zero_coefs <- union(list_cause1, list_cause2)
    #                                    paste("X", non_zero_coefs , sep = "")
    #                                    trainnew <- cbind(train_cv[, c(colnames(train_cv) %in% non_zero_coefs)],
    #                                                      ftime = (train_cv$ftime), fstatus = train_cv$fstatus)
    #                                    fitSmoothHazard(fstatus ~. +log(ftime) -fstatus,
    #                                                                data = trainnew,
    #                                                                ratio = 100,
    #                                                                time = "ftime")
    #                                  },
    #                                  cvs_non_zero_coefs_cause2,
    #                                  mc.cores = ncores, mc.set.seed = seed)
    
    
    # TODO: Replace with mclapply (figure out how to use mclapply with nested lists)
    for (l in c(1:length(lambdagrid))) {
    # for (l in c(1:3)) {
      print(paste("CURRENT LAMBDA VALUE: ", as.character(lambdagrid[l])))
      non_zero_coefs <- union(cvs_non_zero_coefs_cause1[[l]], cvs_non_zero_coefs_cause2[[l]])
      non_zero_coefs <- paste("covariate.X", non_zero_coefs , sep = "")
      
      
      trainnew <- cbind(train_unparameterized_cv[, c(colnames(train_unparameterized_cv) %in% non_zero_coefs)], 
                        ftime = (train_unparameterized_cv$time), fstatus = train_unparameterized_cv$event_ind)


      model_cb <- fitSmoothHazard(fstatus ~. +log(ftime) -fstatus,
                                  data = trainnew,
                                  time = "ftime")
      
      all_deviances[l, i] = deviance(model_cb)
    
      
      exclude_coefs = c("(Intercept):1", "(Intercept):2", "ftime:1", "ftime:2",
                        "log(ftime):1", "log(ftime):2")
      est_betas = coef(model_cb)[names(coef(model_cb))
                                 %in% exclude_coefs == FALSE]
      non_zero_coefs[l, i] = est_betas
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
  rownames(all_deviances) <- lambdagrid
  rownames(non_zero_coefs) <- lambdagrid
  return(list(lambda.min = lambda.min, non_zero_coefs = non_zero_coefs, lambdagrid = lambdagrid, 
              cv_se = cv_se, deviance_grid = all_deviances))
}



############################### TEST CODE ######################
# res1 <- multinom.post_enet_old(train, test)
# res2 <- multinom.post_enet(train, test)
# res1_cb_post_lasso <- varsel_perc(res1$coefficients, beta1)
# res2_cb_post_lasso <- varsel_perc(res2$coefficients, beta1)

res_relaxed = multinom.relaxed_enet(train = train, seed = 2023)

results_table = runCasebaseSim(400, 20, 5)
summarizedSimCaseBaseTable = formatCaseBaseTable(results_table)







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
