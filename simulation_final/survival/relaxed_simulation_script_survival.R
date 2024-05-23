########################## Survival simulation for relaxed LASSO #######################
# Libraries 
library(casebase)
library(future.apply)
library(glmnet)
library(mtool)
library(parallel)
library(tictoc)
library(tidyverse)
library(foreach)
library(survival)
library(cmprsk)
library(glue)
library(pec)

source("src/final_relaxed_implementation.R")
source("src/fitting_functions.R")

# Setup
n <- 400
p <- 20
N <- 1
seed = 2023
beta1 <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, rep(0, p/2))
beta2 <-  c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, rep(0, p/2))

  # For each simulation you want to output the MSE
Results <- replicate(N, {
  tryCatch({
      
    # Simulate data
    sim.data <- cause_hazards_sim(n = n, p = p,  
                                  beta1 = beta1, beta2 = beta2, rate_cens = 0.15, 
                                  h1 = 0.55, h2 = 0.15, gamma1 = 1.5, gamma2 = 1.5)
    
    # Training-test split 
    train.index <- caret::createDataPartition(as.factor(sim.data$fstatus), p = 0.75, list = FALSE)
    train <- sim.data[train.index,]
    test <- sim.data[-train.index,]
    
    test_list = list("time" = test$ftime,
                     "event_ind" = test$fstatus,
                     "covariates" = as.matrix(test[, paste("X", as.character(seq(1:p)), sep = "")]),
                     "offset" = test$offset)  
      
      
    ############################ Coxnet model cause 1 #############################
    # Censor competing event
    y_train1 <- Surv(time = train$ftime, event = train$fstatus == 1)
    
    x_train1 <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
    
    # Censor competing event
    y_test1 <- Surv(time = test$ftime, event = test$fstatus == 1)
    
    x_test1 <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
    
    set.seed(seed)
    # Fit cause-specific cox model with glmnet on training set 
    cox_mod1 <- cv.glmnet(x = x_train1, y = y_train1, family = "cox", alpha = 1, nfolds = 5, relax = TRUE)
    
    # Fit on validation set 
    cox_val_min1 <- glmnet(x = x_test1, y = y_test1, family = "cox", alpha = 1, 
                           lambda = cox_mod1$lambda.min)
    
    cc_min1 <- coef(cox_val_min1)
    
    res_cox_min1 <- varsel_perc(cc_min1, beta1)
    
    ############################ Coxnet model cause 2 #############################
    # Censor competing event
    y_train2 <- Surv(time = train$ftime, event = train$fstatus == 2)
    
    x_train2 <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
    
    # Censor competing event
    y_test2 <- Surv(time = test$ftime, event = test$fstatus == 2)
    
    x_test2 <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
    
    # Fit cause-specific cox model with glmnet on training set 
    cox_mod2 <- cv.glmnet(x = x_train2, y = y_train2, family = "cox", alpha = 1, nfolds = 5, relax = TRUE)
    
    # Fit on validation set 
    cox_val_min2 <- glmnet(x = x_test2, y = y_test2, family = "cox", alpha = 1, 
                           lambda = cox_mod2$lambda.min)
    
    cc_min2 <- coef(cox_val_min2)
    
    res_cox_min2 <- varsel_perc(cc_min2, beta2)
    
    glmnet_cox_coefs = list("coefficients" = as.matrix(cbind(cc_min1, cc_min2)))
    
    
    # Test multi-deviance
    dev_cox <- multi_deviance_final(data = test_list, fit_object = glmnet_cox_coefs, cox = TRUE)
    
    ######################### casebase model with lassox ####################
    # Train case-base model
    lasso.lambda <- mtool.multinom.cv(train, seed = seed, alpha = 1, nfold = 5, train_ratio = 20)
    
    # Test set
    surv_obj_test <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
    
    # Covariance matrix
    cov_test <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
    
    # Case-base dataset
    cb_data_test <- create_cbDataset(surv_obj_test, as.matrix(cov_test), ratio = 10)
    
    # Case-base fit
    # Lambda.min
    ouptut <- capture.output(fit_lambda_min_lasso <- fit_model(cb_data_test, regularization = 'l1',
                               lambda = lasso.lambda$lambda.min , alpha = 1, unpen_cov = 2))
    
    coef_val1 <- fit_lambda_min_lasso$coefficients[1:p, 1]
    coef_val2 <- fit_lambda_min_lasso$coefficients[1:p, 2]
    
    res_lambda_min1_lasso <- varsel_perc(coef_val1, beta1)
    res_lambda_min2_lasso <- varsel_perc(coef_val2, beta2)
    
    
    # Test multi-deviance
    dev_cb_mtool_lasso <- multi_deviance_final(data = cb_data_test, fit_object = fit_lambda_min_lasso)
    
    
    
    ######################### casebase model (relaxed fit with no penalization) ####################
    # Train case-base model
    res_relaxed <- multinomial_casebase_relaxed_lasso(train, seed = seed, alpha = 1, gamma = 0,
                                           nfold = 5, train_ratio = 20)
    
    
    # Case-base relaxed fit
    # Lambda.min
    output <- capture.output(fit_lambda_min_relaxed_no_pen <- fit_model(cb_data_test, regularization = 'l1',
                                  lambda = res_relaxed$lambda.min , alpha = 1, unpen_cov = 2))
    
    coef_val1 <- fit_lambda_min_relaxed_no_pen$coefficients[1:p, 1]
    coef_val2 <- fit_lambda_min_relaxed_no_pen$coefficients[1:p, 2]
    
    res_lambda_min1_relaxed_no_pen <- varsel_perc(coef_val1, beta1)
    res_lambda_min2_relaxed_no_pen <- varsel_perc(coef_val2, beta2)
    
    
    # Test multi-deviance
    dev_cb_mtool_relaxed_no_pen <- multi_deviance_final(data = cb_data_test, fit_object = fit_lambda_min_relaxed_no_pen)
    
    
    ######################### casebase model (relaxed fit with small penalization) ####################
    # Train case-base model
    res_relaxed_small_pen <- multinomial_casebase_relaxed_lasso(train, seed = seed, alpha = 1, gamma = 0.00005,
                                                    nfold = 5, train_ratio = 20)
    
    
    # Case-base relaxed fit
    # Lambda.min
    output <- capture.output(fit_lambda_min_relaxed_small_pen <- fit_model(cb_data_test, regularization = 'l1',
                                              lambda = res_relaxed_small_pen$lambda.min , alpha = 1, unpen_cov = 2))
    
    coef_val1 <- fit_lambda_min_relaxed_small_pen$coefficients[1:p, 1]
    coef_val2 <- fit_lambda_min_relaxed_small_pen$coefficients[1:p, 2]
    
    res_lambda_min1_relaxed_small_pen <- varsel_perc(coef_val1, beta1)
    res_lambda_min2_relaxed_small_pen <- varsel_perc(coef_val2, beta2)
    
    
    # Test multi-deviance
    dev_cb_mtool_relaxed_small_pen <- multi_deviance_final(data = cb_data_test, fit_object = fit_lambda_min_relaxed_small_pen)
    
    
    ######################### casebase model (relaxed fit with sig penalization) ####################
    # Train case-base model
    res_relaxed_sig_pen <- multinomial_casebase_relaxed_lasso(train, seed = seed, alpha = 1, gamma = 0.001,
                                                          nfold = 5, train_ratio = 20)
    
    
    # Case-base relaxed fit
    # Lambda.min
    output <- capture.output(fit_lambda_min_relaxed_sig_pen <- fit_model(cb_data_test, regularization = 'l1',
                                                    lambda = res_relaxed_sig_pen$lambda.min , alpha = 1, unpen_cov = 2))
    
    coef_val1 <- fit_lambda_min_relaxed_sig_pen$coefficients[1:p, 1]
    coef_val2 <- fit_lambda_min_relaxed_sig_pen$coefficients[1:p, 2]
    
    res_lambda_min1_relaxed_sig_pen <- varsel_perc(coef_val1, beta1)
    res_lambda_min2_relaxed_sig_pen <- varsel_perc(coef_val2, beta2)
    
    
    # Test multi-deviance
    dev_cb_mtool_relaxed_sig_pen <- multi_deviance_final(data = cb_data_test, fit_object = fit_lambda_min_relaxed_sig_pen)
    
    
    
    #############################################################
    Res <- rbind(res_cox_min1, res_cox_min2, res_lambda_min1_lasso, res_lambda_min2_lasso,
                 res_lambda_min1_relaxed_no_pen, res_lambda_min2_relaxed_no_pen, res_lambda_min1_relaxed_small_pen, res_lambda_min2_relaxed_small_pen, 
                 res_lambda_min1_relaxed_sig_pen, res_lambda_min2_relaxed_sig_pen)
    
    rownames(Res) <- c("cox.lambda.min_cause1", "cox.lambda.min_cause2",
                       "casebase.lasso.lambda.min_cause1", "casebase.lasso.lambda.min_cause2",
                       "casebase.relaxed.lambda.min.no.pen_cause1", "casebase.relaxed.lambda.min.no.pen_cause2",
                       "casebase.relaxed.lambda.min.small.pen_cause1", "casebase.relaxed.lambda.min.small.pen_cause2",
                       "casebase.relaxed.lambda.min.sig.pen_cause1", "casebase.relaxed.lambda.min.sig.pen_cause2")
    
    
    devs <- rbind(dev_cb_mtool_lasso, dev_cb_mtool_relaxed_no_pen, dev_cb_mtool_relaxed_small_pen, dev_cb_mtool_relaxed_sig_pen)
    
    rownames(devs) <- c("casebase.lasso.lambda.min", "casebase.relaxed.lambda.min.no.pen",
                        "casebase.relaxed.lambda.min.small.pen", "casebase.relaxed.lambda.min.sig.pen")
    
    list(Res, devs)
  },
  
  error = function(e) {
    print(e)
    
    message("Error: ", conditionMessage(e))
    
    # Retrieve and print the call stack
    calls <- sys.calls()
    
    # Locate and print the offending call
    offending_call <- calls[[length(calls) - 1]]
    message("Offending line: ", deparse(offending_call))
    
    # Optional: Use traceback for more detailed information
    traceback()
    
    # Return NA or another indicator for the error
    return(list(Res, devs))
  })
  
}, simplify = FALSE)

print(Results[[1]][[1]])

sim_coefficient_results = data.frame(matrix(nrow = 0, ncol = ncol(Results[[1]][[1]])))
sim_deviance_results = data.frame(matrix(nrow = 0, ncol = ncol(Results[[1]][[2]])))
for(i in c(1:N)) {
  sim_coefficient_results = rbind(sim_coefficient_results, Results[[i]][[1]])
  sim_deviance_results = rbind(sim_deviance_results, Results[[i]][[2]])
}


write_csv(as.data.frame(sim_coefficient_results), file = glue("simulation_final/survival/single_run_results/coefficients/results_survival_{runif(1)}.csv"))
write_csv(as.data.frame(sim_deviance_results), file = glue("simulation_final/survival/single_run_results/deviances/deviances_survival_{runif(1)}.csv"))