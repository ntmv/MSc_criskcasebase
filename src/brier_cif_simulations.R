###################### Brier score and CIF simulations #########################################
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

# Test time-points for Brier score 
time_points <- sort(unique(test$ftime))

time_points <- seq(0, 6, 0.001)

####### Fit the following models: ##################
## 1)ICR
## 2) penCR
## 3) Oracle-Cox
## 4) Reference 
## 5) Oracle-SDH
## 6) CoxBoost (Fine-Gray)
## 7) Case-base


################ ############# Brier score for case base ############### ##################
# Combine scores
data_brier <- bind_rows(
  briercoxboost$Brier$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"),
  score_rdata$Brier$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"),
  # brierfinegray$AUC$score %>% 
  #   mutate(model = as.character(model)) %>% 
  #   filter(model != "Null model"),
  briercasebase$Brier$score %>%
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"),  
  # brierCoxKM$Brier$score %>% 
  #   mutate(model = as.character(model)) %>% 
  #   filter(model != "Null model")
) %>% 
  mutate(model = factor(model, levels = c("CoxBoost", "Pen. Cox", "Case-base")))

penCR = temp1$BootCvErr$cv.glmnet.CR)

iCS = two.i.CSlassos(data = train, nlambda = 100)

xx <- bind_rows(FineBoost = temp3$AppErr$iCoxBoost)

xx <- pivot_longer(xx, cols = c(1), values_to = "brier", names_to = "Model")

xx <- data.frame(cbind(time_points, xx))
  
ggplot(data = xx, aes(x = time_points, y = brier, col = Model)) +
  geom_line(size = 1) +
  xlab("Follow-up time (years)") +
  ylab("Brier Score") +
  labs(color = "Models") +
  theme_bw() + scale_colour_brewer(palette = "Dark2") 
################ ############# ############### ################## ############### ##################

# Penalized casebase
model_cb <- fitSmoothHazard(fstatus ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8
                            + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 +
                              X17 + X18 + X19 +X20,
                            ratio = 5, data = train, 
                            time = "ftime")

temp <-  pec(model_cb, formula=Hist(ftime,fstatus) ~1, data= test, cause = 1,
           times= time_points, start=NULL,exact=FALSE, cens.model = 'marginal')

ibs(temp,times= 1)

temp1 <- pec(penCR, formula=Hist(ftime,fstatus) ~ 1, data = test, cause = 1,
             times= time_points, start=NULL,exact=FALSE, cens.model = 'marginal')

ibs(temp1,times=1)

temp2 <- pec(iCS, formula=Hist(ftime,fstatus) ~ 1, data = test, cause = 1,
             times= time_points, start=NULL,exact=FALSE, cens.model = 'marginal')

ibs(temp2,times=1)

temp3 <- pec(cbfit, formula=Hist(ftime,fstatus) ~ 1, data = test, cause = 1,
             times= time_points, start=NULL,exact=FALSE, cens.model = 'marginal')

ibs(temp3,times=1)

temp$AppErr







###################################################################










# Fit case-base model 
# Convert to case-base dataset
surv_obj_train <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))

cov_train <- as.matrix(cbind(train[, c(grepl("X", colnames(train)))], time = log(train$ftime)))

# Create case-base dataset
cb_data_train <- create_cbDataset(surv_obj_train, cov_train, ratio = 10)

tic()
cv.alpha <- mtool.multinom.cv(cb_data_train, alpha = 1, nfold = 5)
toc()

cv.alpha

# Cross-validation plot 
p1 <- plot_cv.multinom(cv.alpha$deviance_grid, cv.alpha$lambdagrid, cv.alpha$lambda.min, cv.alpha$lambda.1se, nfold = 5)

# validation set
surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))

# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = cv.alpha$lambda.min , alpha = 1)

res_cb_min <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)

# Penalized casebase
model_cb <- fitSmoothHazard(fstatus ~. -ftime +log(ftime) -fstatus,
                            data = testnew,
                            ratio = 5,
                            time = "ftime")

coef_vector <- as.vector(t(fit_val_min$coefficients))

coef_vector_intercept <- coef_vector[319:320]

coef_vector <-  coef_vector[-c(319:320)]

coef_vector <- c(coef_vector_intercept, coef_vector)

names_coef <- names(coef(model_cb))

names_coef <- c("(Intercept):1", "(Intercept):2", "X1:1", "X1:2", )

model_cb@coefficients <-  coef_vector

names(model_cb@coefficients) <- names_coef 


# Brier score case-base
newx <- model.matrix(fstatus ~ . - ftime -fstatus,
                     data = test)[, -c(1)] 

briercasebase <- Score(list("Case-base" = model_cb),
                       data = test,
                       formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                       se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                       times = time_points)


brieriCS <- Score(list("iCS" = iCS),
                       data = test,
                       formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                       se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                       times = time_points)

brierpenCR <- Score(list("penCR" = penCR),
                  data = test,
                  formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                  se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                  times = time_points)

class(fit_val_min) <- "penalizedCompRisk"

brierpenCB <- Score(list("penCB" = fit_val_min),
                    data = test,
                    formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                    se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                    times = time_points)

#############################################################################
# Cox proportional hazards cause-specific model 
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

# Coefficients 
nonzero_coef_cox <- coef(cox_val_min, s = "lambda.min")

# Fitting a cox model using regular estimation, however we will not keep it.
# this is used more as an object place holder.
# We expect this not to converge
coxNet <- survival::coxph(Surv(time = ftime, event = fstatus == 1) ~ ., 
                          data = train, x = TRUE)

# The coefficients of this object will be replaced with the estimates from the
# original coxNet. Doing so makes it so that everything is invalid aside from
# the coefficients. In this case, all we need to estimate the absolute risk is
# the coefficients. Std. error would be incorrect here, if we were to draw error
# bars.
coxNet_coefnames <- names(coxNet$coefficients)
coxNet$coefficients <- as.vector(nonzero_coef_cox)
names(coxNet$coefficients) <- coxNet_coefnames

# Brier score for cox 
brierPenalized <- Score(list("Pen. Cox" = coxNet),
                        data = test,
                        formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                        se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                        times = time_points)



fit_csh <- CSC(Hist(ftime, fstatus) ~ 1,
               data = train,
               fitter = "cph")
fit_csc1 <- fit_csh$models$`Cause 1`
fit_csc2 <- fit_csh$models$`Cause 2`

# Overall performance measures ----------------
primary_event <- 1 # Set to 2 if cause 2 was of interest

# Development data
score_rdata <- Score(
  list("csh_development" = fit_csh),
  formula = Hist(ftime, fstatus) ~ 1,
  data = test,
  times = time_points,
  summary = c("ipa"),
  metrics = c("brier")
)

########################################################################
# Fit coxboost model
# Cross-validation to find optimal penalty parameter
optim.res <- iCoxBoost(Hist(ftime, fstatus)~., data = train, cause = 1, cv = TRUE)

cbfit <-  iCoxBoost(Hist(ftime, fstatus)~., data = test, cause = 1,
                    stepno = optim.res$cv.res$optimal.step)

# Brier score for coxboost
briercoxboost <- Score(list("CoxBoost" = cbfit),
                        data = test,
                        formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                        se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                        times = time_points)


cv.coxboostCR = function(data, failcode, maxstepno=200,  penalty, ...){
  data = data.frame(data)
  vars = colnames(data)[(!colnames(data) %in% c('time', 'status', 'ID','id'))]
  X = as.matrix(data[, vars])
  
  # subdistribution hazard with respect to failcode
  Event = rep(2, NROW(data$status))
  Event[data$status == failcode] <- 1
  Event[data$status == 0] <- 0
  y = Hist(data$time, Event)
  
  
  if (missing(penalty)) penalty <- sum(Event == 1) * (9)
  
  cv.res <- cv.CoxBoost(data$time, Event, X, maxstepno=maxstepno,
                        K=10, type="verweij", penalty = penalty, ...)
  
  cbfit <- CoxBoost(data$time, Event, X, stepno = cv.res$optimal.step, penalty = penalty, ...)
  
  out <- list('cv.res'=cv.res, 'cbfit'=cbfit, stepno = cv.res$optimal.step, 
              call = match.call(), 'response' = y)
  class(out) <- c("coxboostCR", 'CoxBoost')
  out
}


predictEventProb.coxboostCR = function(object, newdata, times, ...) {
  newdata = data.frame(newdata)
  vars = colnames(newdata)[(!colnames(newdata) %in% c('time', 'status','event', 'ID','id'))]
  newx = data.frame(newdata)
  newx = as.matrix(newx[, vars])
  
  p = predict(object$cbfit, newdata = newx, type = 'CIF', times = times)
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))  stop("Prediction failed")
  p
}

###############################################################################
fit.aj <- prodlim(Hist(ftime,fstatus)~1,data= train)

# Development data
score_aj <- Score(
  list("Aalen-Johnson" = fit.aj),
  formula = Hist(ftime, fstatus) ~ 1,
  data = test,
  times = time_points,
  summary = c("ipa"),
  metrics = c("brier")
)

library(timereg)

fit.fg <- comp.risk(Event(ftime, fstatus)~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 +
                      X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16
                    +X17 + X18 + X19 + X20, data = train, cause = 1)


# Development data
score_oracle <- Score(
  list("Aalen-Johnson" = fit.aj),
  formula = Hist(ftime, fstatus) ~ 1,
  data = test,
  times = time_points,
  summary = c("ipa"),
  metrics = c("brier")
)

########################## Plot Brier score ###############################
# Combine scores
data_brier <- bind_rows(
  briercoxboost$Brier$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"), 
  brieriCS$Brier$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"),
 # brierpenCB$Brier$score %>% 
 #   mutate(model = as.character(model)) %>% 
  #  filter(model != "Null model"),
 # score_oracle$Brier$score %>% 
  #  mutate(model = as.character(model)) %>% 
  #  filter(model != "Null model"),
 #brierpenCR$Brier$score %>% 
 #  mutate(model = as.character(model)) %>% 
  # filter(model != "Null model"),
 # brierPenalized$Brier$score %>% 
   # mutate(model = as.character(model)) %>% 
 #   filter(model != "Null model"),
 # score_rdata$Brier$score %>% 
 #   mutate(model = as.character(model)) %>% 
 #   filter(model != "Null model"),
  # brierfinegray$AUC$score %>% 
  #   mutate(model = as.character(model)) %>% 
  #   filter(model != "Null model"),
  briercasebase$Brier$score %>%
 mutate(model = as.character(model)) %>% 
 filter(model != "Null model"),  
 # brierCoxKM$Brier$score %>% 
 #   mutate(model = as.character(model)) %>% 
 #   filter(model != "Null model")
) %>% 
  mutate(model = factor(model, levels = c("CoxBoost", "iCR" , "Post-OLS Case-base")))

data_brier$model <- factor(data_brier$model)

levels(data_brier$model) <- c("Aalen-Johnson", "De-biased Case-base", "Boosted Fine-Gray", 
                              "iCR",  "Case-base")

levels(data_brier$model) <- c("De-biased Case-base",  "Boosted Fine-Gray",  "iCR")

cbPalette <- c("#56B4E9", "#D55E00", "#E69F00", "#009E73", "#CC79A7")

cbPalette <- c("#D55E00", "#56B4E9", "#E69F00")

data_brier$title <- "N = 400, p = 120, Mixture of Effects"

png(filename = "~/Desktop/brier_setting1.png", res = 300, height = 12, width = 20, units = "cm") 
ggplot(data = data_brier, aes(x = times, y = Brier, col = model)) +
  geom_line(size = 0.5) + 
  geom_line(data = filter(data_brier, model == "De-biased Case-base"), size = 1.5) + 
  geom_line(data = filter(data_brier, model == "Case-base"), size = 1.5) + 
  xlab("Follow-up time (years)") +
  ylab("Brier Score") +
  labs(color = "Models") +
  theme_bw() + scale_colour_brewer(palette = "Dark2") +
  scale_colour_manual(values=cbPalette) + labs(y = "Brier Score for Cause 1 Predictions") +
  facet_grid(. ~ title)
dev.off()

################## Plot CIF #########################################
# Estimate absolute risk curve
risk_cb <- absoluteRisk(
  object = model_cb, time = time_points,
  method = "numerical"
)

risk_cb <- risk_cb[-1, ]
risk_cb <- as.data.frame(rowMeans(risk_cb))

risk_fg <- as.data.frame(colMeans(risk_fg$P1))

# Estimate absolute risk curve
risk_cox <- survfit(coxNet, newdata = test)

coxNet <- coxph(Surv(time, status == 1) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
                  X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20,
                data = train)
risk_cox <- survival::survfit(coxNet, type = "breslow", 
                              newdata = test)

# Prepare data for coxph
bmtcrr_cox <- transform(train, 
                        id = seq_len(nrow(train)),
                        Status = factor(fstatus))

model_cox <- coxph(Surv(ftime, Status) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
                     X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 , data = bmtcrr_cox,
                   id = id, x = TRUE)

risk_cox <- survival::survfit(model_cox, type = "breslow", 
                              newdata = test)


risk_all <- dplyr::bind_rows(
  data.frame(
    Time = time_points,
    Method = "Case-base",
    Risk = risk_cb,
    stringsAsFactors = FALSE
  ),
 # data.frame(
 #   Time = risk_cox$time,
 #   Method = "Cox",
 #   Risk = risk_cox[,2]$pstate[,1,1],
 #   stringsAsFactors = FALSE
 # ),
 # data.frame(
 #   Time = time_point,
 #   Method = "Fine-Gray",
 #   Risk = risk_fg$`colMeans(risk_fg$P1)`,
  #  stringsAsFactors = FALSE
  )
) %>% 
  dplyr::filter(Time <= 60)

png(filename = "~/Desktop/cif_iidnonsparse.png", res = 300, height = 10, width = 15, units = "cm") 
ggplot(risk_all, aes(x = Time, y = Risk, colour = Method)) +
  # geom_line for smooth curve
  geom_line(data = dplyr::filter(risk_all, Method == "Case-base"), size = 1) +
  geom_step(data = dplyr::filter(risk_all, Method != "Case-base"), size = 1) +
  ylim(c(0, 1)) + theme(text = element_text(size = 30)) + 
  xlab("Time (in Years)") +
  ylab("Absolute Risk") + theme_bw()
dev.off()
