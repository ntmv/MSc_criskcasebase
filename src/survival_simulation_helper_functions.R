#' ################################################# Survival simulation helper functions ##################################################


###########################################################
#' Runs simulation of hazard function fits

runCasebaseSim = function(n = 400, p = 20, N = 5, nfolds = 5, seed) {
  # n = 400
  # p = 20
  # N = 5
  # nfolds = 5
  # seed = 2023
  
  # tryCatch({
  # For each simulation you want to output the MSE
  Results <- replicate(N, {
    # Set seed
    # take the last five digits of the initial seed
    # seed = as.integer(Sys.time())
    # the_seed= seed %% 100000
    # set.seed(the_seed)
    
    set.seed(seed)
    
    num_true <- p/2
    beta1 <- c(rep(0, p))
    beta2 <- c(rep(0, p))
    nu_ind <- seq(num_true)
    # Here out of 20 predictors, 10 should be non-zero 
    beta1[nu_ind] <- c(rep(1, p/2), rep(0, p/2))
    beta2[nu_ind] <- c(rep(-1, p/2), rep(0, p/2))
    
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
    cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7, folds = nfolds)
    
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
    cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7, folds = nfolds)
    
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
                         nlambda = 20, t.BS = median(train$ftime), seed = seed, causeOfInt = 1,
                         nfold = nfolds)
    
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
    start_cv <- Sys.time()
    cv.lambda <- mtool.multinom.cv(train, seed = seed, nfold = nfolds, lambda_max =  0.9, alpha = 0.7)
    end_cv <- Sys.time()
    
    
    # Case-base fits 
    # Lambda.min
    fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                               lambda = cv.lambda$lambda.min , alpha = 0.7, unpen_cov = 2)
    
    
    res_cb_min1 <- varsel_perc(fit_val_min$coefficients[1:p, 1], beta1)
    
    res_cb_min2 <- varsel_perc(fit_val_min$coefficients[1:p, 2], beta2)
    
    # Calculate MSE here as well!
    # MSE for casebase model for cause of interest
    fit_val_coef_1 <- fit_val_min$coefficients[1:p, 1]
    mean((beta1[nu_ind] - fit_val_coef_1[nu_ind])^2)
    
    # MSE for casebase model for competing risk
    fit_val_coef_2 <- fit_val_min$coefficients[1:p, 2]
    cb_min_bias <- which(fit_val_coef_2 != 0)
    mean((fit_val_coef_2[nu_ind] - beta2[nu_ind])^2)
    
    
    ########################## Fit casebase model with post LASSO#################################
    
    start_post <- Sys.time()
    res <- multinom.post_enet_old(train, test, nfold = nfolds, seed = seed, lambda_max =  0.9, alpha = 0.7)
    end_post <- Sys.time()
    
    # Calculate MSE for this as well (fill in here)
    mean((res$coefficients[nu_ind, 1]- beta1[nu_ind])^2)
    mean((res$coefficients[nu_ind, 2] - beta2[nu_ind])^2)
    
    res_cb_post_lasso1 <- varsel_perc(res$coefficients[nu_ind, 1], beta1)
    res_cb_post_lasso2 <- varsel_perc(res$coefficients[nu_ind, 2], beta2)
    
    
    ########################## Fit casebase model with relaxed LASSO#################################
    
    start_relaxed <- Sys.time()
    res <- multinom.relaxed_enet(train, nfold = nfolds, seed = seed, lambda_max =  0.9, alpha = 0.7)
    end_relaxed <- Sys.time()
    print(res$lambda.min)
    
    model_final <- fitSmoothHazard(fstatus ~. +log(ftime) -fstatus,
                                   data = test,
                                   time = "ftime",
                                   lambda = res$lambda.min)
    
    all_coef_names_cause1 =  paste("X", seq(1:p), ":1", sep ="")
    all_coef_names_cause2 =  paste("X", seq(1:p), ":2", sep ="")
    
    exclude_coefs = c("(Intercept):1", "ftime:1",
                      "log(ftime):1", "(Intercept):2", "ftime:2",
                      "log(ftime):2")
    
    est_betas = coef(model_final)[!(names(coef(model_final)) %in% exclude_coefs)]
    est_betas_cause1 = est_betas[names(est_betas) %in% all_coef_names_cause1]
    est_betas_cause2 = est_betas[names(est_betas) %in% all_coef_names_cause2]
    
    
    # Calculate MSE for this as well (fill in here)
    mean((est_betas_cause1 - beta1[nu_ind])^2)
    mean((est_betas_cause2 - beta2[nu_ind])^2)
    
    res_cb_relaxed_lasso1 <- varsel_perc(est_betas_cause1, beta1)
    res_cb_relaxed_lasso2 <- varsel_perc(est_betas_cause2, beta2)
    
    ###################################################################################
    
    print(paste(round(end_cv - start_cv, 3), " taken to run multinomial cv"))
    print(paste(round(end_post - start_post, 3), " taken to run post LASSO"))
    print(paste(round(end_relaxed - start_relaxed, 3), " taken to run relaxed LASSO"))
    
    Res <- rbind(res_cb_relaxed_lasso1, res_cb_relaxed_lasso2, res_cb_post_lasso1, 
                 res_cb_post_lasso2, res_cb_min1, res_cb_min2, 
                 res_cox_min1, res_cox_min2, res_pencr_min1, res_pencr_min2, cen.prop)
    
    rownames(Res) <- c("casebase.relaxed.lasso.lambda.min_cause1", "casebase.relaxed.lasso.lambda.min_cause2", 
                       "casebase.post.lasso.lambda.min_cause1", "casebase.post.lasso.lambda.min_cause2", 
                       "casebase.lambda.min_cause1", 
                       "casebase.lambda.min_cause2", "cox.lambda.min_cause1",
                       "cox.lambda.min_cause2", "pencr.lambda.mincause1", "pencr.lambda.mincause2", "cens.prop")
    
    #     Res
    #     
    #   }, simplify = FALSE)
    #   
    # Results <- do.call(rbind, Results)
    # Results
    # # },
    # # 
    # # error = function(e) {
    # #   print(e)
  })
}


plot_post_relaxed.multinom <- function(cv_object) {
  nfold <- ncol(cv_object$deviance_grid)
  mean_dev <- rowMeans(cv_object$deviance_grid)
  row_stdev <- apply(cv_object$deviance_grid, 1, function(x) {sd(x)/sqrt(nfold)})
  plot.dat.p <- data.frame(lambdagrid = cv_object$lambdagrid, mean.dev = mean_dev, 
                           upper = mean_dev +row_stdev, lower = mean_dev - row_stdev)
  p <- ggplot(plot.dat.p, aes(log(lambdagrid), mean.dev)) + geom_point(colour = "red", size = 3) + theme_bw() + 
    geom_errorbar(aes(ymin= lower, ymax=upper), width=.2, position=position_dodge(0.05), colour = "grey") + 
    labs(x = "log(lambda)", y = "Multinomial Deviance")  + 
    geom_vline(xintercept = log(cv_object$lambda.min), linetype = "dotted", colour = "blue")
  
  return(p)
}





###########################################################
#' Formats resulting table from simulation of hazard function fits

formatCaseBaseTable = function(sim_results) {
  model_fit_labels = c("casebase.relaxed.lasso.lambda.min_cause1", "casebase.relaxed.lasso.lambda.min_cause2", 
                       "casebase.post.lasso.lambda.min_cause1", "casebase.post.lasso.lambda.min_cause2", 
                       "casebase.lambda.min_cause1", 
                       "casebase.lambda.min_cause2", "cox.lambda.min_cause1",
                       "cox.lambda.min_cause2", "pencr.lambda.mincause1", "pencr.lambda.mincause2", "cens.prop")
  
  # stat_names = colnames(as.data.frame(sim_results))
  
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
  
  stat_names = c("Model", "Sensitivity",	"Specificity",	"MCC",	"Coefficient_Bias")
  
  formatted_table = data.frame(matrix(nrow = 0, ncol = length(stat_names)))
  colnames(formatted_table) = stat_names
  stat_table = sim_results[stat_names]
  
  for (i in model_fit_labels) {
    average_stat_row = colMeans((stat_table %>% filter(Model == i))[2:ncol(stat_table)])
    formatted_table = rbind(formatted_table, c(i, average_stat_row))
  }
  colnames(formatted_table) = stat_names
  formatted_table = formatted_table %>% replace(is.na(.), NaN)
  
  options(digits = 5)
  
  return(formatted_table)
}






# # Penalized Multinomial Logistic Regression (LASSO) with decreased learning rate
# mtool.MNlogistic_new <- function(X, Y, offset, N_covariates,
#                              regularization = 'l1', transpose = F,
#                              lambda1, lambda2 = 0, lambda3 = 0,
#                              learning_rate = 1e-3, tolerance = 1e-4,
#                              niter_inner_mtplyr = 7, maxit = 100, ncores = -1,
#                              group_id, group_weights,
#                              groups, groups_var,
#                              own_variables, N_own_variables) {
#   ## Dimensions and checks
#   nx <- nrow(X)
#   
#   if (!is.vector(Y)) {Y <- as.vector(Y)}
#   ny <- length(Y)
#   
#   if (!is.vector(offset)) {offset <- as.vector(offset)}
#   noff <- length(offset)
#   
#   if (nx == ny & nx == noff) {
#     n <- nx
#   } else {
#     stop('X, Y and offset have different number of observations.')
#   }
#   
#   p <- ncol(X)
#   
#   K <- length(unique(Y)) - 1
#   
#   ## regularization
#   pen1 <- c("l0", "l1", "l2", "linf", "l2-not-squared",
#             "elastic-net", "fused-lasso",
#             "group-lasso-l2", "group-lasso-linf",
#             "sparse-group-lasso-l2", "sparse-group-lasso-linf",
#             "l1l2", "l1linf", "l1l2+l1", "l1linf+l1", "l1linf-row-column",
#             "trace-norm", "trace-norm-vec", "rank", "rank-vec", "none")
#   pen2 <- c("graph", "graph-ridge", "graph-l2", "multi-task-graph")
#   pen3 <- c("tree-l0", "tree-l2", "tree-linf", "multi-task-tree")
#   
#   if (regularization %in% pen1) { penalty <- 1 }
#   if (regularization %in% pen2) { penalty <- 2 }
#   if (regularization %in% pen3) { penalty <- 3 }
#   if (! regularization %in% c(pen1, pen2, pen3)) {
#     stop('The provided regularization is not supported.')
#   }
#   
#   ### check regularization-specific inputs
#   #### penalty = 1, call proximal(Flat), requires `group_id` in integer vector
#   if (penalty == 1) {
#     if (missing(group_id)) { group_id <- rep(0L, p) }
#     group_weights <- vector(mode = 'double')
#     groups <- matrix(NA)
#     groups_var <- matrix(NA)
#     own_variables <- vector(mode = 'integer')
#     N_own_variables <- vector(mode = 'integer')
#   }
#   
#   #### penalty = 2, call proximalGraph
#   #### requires `groups` and `groups_var` in integer matrices and `group_weights` in double vector
#   if (penalty == 2) {
#     if (missing(groups)) { stop('Required input `groups` is missing.') }
#     if (missing(groups_var)) { stop('Required input `groups_var` is missing.') }
#     if (missing(group_weights)) { stop('Required input `group_weights` is missing.') }
#     group_id <- rep(0L, p)
#     own_variables <- vector(mode = 'integer')
#     N_own_variables <- vector(mode = 'integer')
#   }
#   
#   #### penalty = 3, call proximalGraph
#   #### requires `own_variables` and `N_own_variables` in integer vectors, `group_weights` in double vector
#   #### and `groups` in integer matrix
#   if (penalty == 3) {
#     if (missing(groups)) { stop('Required input `groups` is missing.') }
#     if (missing(own_variables)) { stop('Required input `own_variables` is missing.') }
#     if (missing(N_own_variables)) { stop('Required input `N_own_variables` is missing.') }
#     if (missing(group_weights)) { stop('Required input `group_weights` is missing.') }
#     group_id <- rep(0L, p)
#     groups_var <- matrix(NA)
#   }
#   
#   ## call mtool main function
#   result <- MultinomLogistic(X = X, Y = Y, offset = offset, K = K, reg_p = p - N_covariates,
#                              penalty = penalty, regul = regularization, transpose = transpose,
#                              grp_id = group_id, etaG = group_weights,
#                              grp = groups, grpV = groups_var,
#                              own_var = own_variables, N_own_var = N_own_variables,
#                              lam1 = lambda1, lam2 = lambda2, lam3 = lambda3,
#                              learning_rate = learning_rate, tolerance = tolerance,
#                              niter_inner = niter_inner_mtplyr * nx, maxit = maxit,
#                              ncores = ncores)
#   nzc <- length(result$`Sparse Estimates`@i)
#   return(list(coefficients = result$`Sparse Estimates`,
#               no_non_zero = nzc))
# }



