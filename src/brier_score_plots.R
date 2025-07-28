# libraries
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
library(CoxBoost)
library(riskRegression)


# Test time-points for Brier score 
time_points <- sort(unique(test$time))
################## Brier score casebase #########################
# Test set 
surv_obj_train <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))

# Covariance matrix
cov_val <- cbind(train[, c(grepl("X", colnames(train)))], time = log(train$ftime))

# Case-base dataset
cb_data_train <- create_cbDataset(surv_obj_train, as.matrix(cov_val), ratio = 10)

# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_train, regularization = 'elastic-net',
                           lambda = cv.lambda$lambda.min , alpha = 0.7, unpen_cov = 2)


res_cb_min1 <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)

non_zero_coefs_cause1 <- which(fit_val_min$coefficients[1:eval(parse(text="p")), 1] != 0)
non_zero_coefs <- paste("X", non_zero_coefs_cause1 , sep = "")
testnew <- cbind(test[, c(colnames(test) %in% non_zero_coefs)], ftime = (test$ftime), fstatus = test$fstatus)


# Fit unpenalized model
f <- fitSmoothHazard(fstatus ~ . + log(ftime) - fstatus, data = testnew, ratio = 100, time = "ftime")


briercasebase <- Score(list("Case-base" = model_cb),
                       data = test,
                       formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                       se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                       times = time_points, cause = 1)
############################# CoxBoost ########################################
optim.res <- iCoxBoost(Hist(ftime, fstatus)~., data = train, cause = 1, cv = TRUE)

cbfit <-  iCoxBoost(Hist(ftime, fstatus)~., data = train, cause = 1,
                    stepno = optim.res$cv.res$optimal.step)

# Brier score for coxboost
briercoxboost <- Score(list("CoxBoost" = cbfit),
                       data = test,
                       formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                       se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                       times = time_points, cause = 1)

########################## iCR ###############################################
# Prepare data for coxph
iCS = two.i.CSlassos(data = train, nlambda = 100)

brieriCS <- Score(list("iCS" = iCS),
                  data = test,
                  formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                  se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                  times = time_points, cause = 1)

###################### Aalen-Johnson #########################################
fit.aj <- prodlim(Hist(ftime,fstatus)~1,data= train)

# Aalen-Johnson 
score_aj <- Score(
  list("Aalen-Johnson" = fit.aj),
  formula = Hist(ftime, fstatus) ~ 1,
  data = test,
  times = time_points,
  summary = c("ipa"),
  metrics = c("brier"))
######################### Penalized casebase ####################################
class(fit_val_min) <- "penalizedCompRisk"

brierpenCB <- Score(list("penCB" = fit_val_min),
                    data = test,
                    formula = Hist(ftime, fstatus) ~ 1, summary = "ibs", 
                    se.fit = FALSE, metrics = "Brier", contrasts = FALSE, 
                    times = time_points, cause = 1)
########################## Plot brier score ######################################
# Combine scores
data_brier <- bind_rows(
  briercasebase$Brier$score %>% 
   mutate(model = as.character(model)) %>% 
    filter(model != "Null model"), 
  brierpenCB$Brier$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"), 
  score_aj$Brier$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"), 
  briercoxboost$Brier$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"), 
  brieriCS$Brier$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"), 
#  brierpenCR$Brier$score %>% 
#    mutate(model = as.character(model)) %>% 
#  filter(model != "Null model")
  ) 

# Convert into factor 
data_brier$model <- factor(data_brier$model)

levels(data_brier$model) <- c("Aalen-Johansen", "De-biased Case-base", "Boosted Fine-Gray", "iCR", "Case-base")

data_brier$title <- "N = 400, p = 1000, Non-proportional hazards"

cbPalette <- c("#808080", "#D55E00", "#CC79A7", "#56B4E9", "#E69F00")


png(filename = "~/Desktop/brier_setting1.png", res = 300, height = 12, width = 20, units = "cm") 
plot2 <- ggplot(data = data_brier, aes(x = times, y = Brier, col = model)) +
  geom_line(size = 0.5) + 
  geom_line(data = filter(data_brier, model == "De-biased Case-base"), size = 2) + 
  geom_line(data = filter(data_brier, model == "Case-base"), size = 1.5) + 
  xlab("Follow-up time (years)") +
  ylab("Brier Score") +
  labs(color = "Models") +
  theme_bw() + scale_colour_brewer(palette = "Dark2") +
  labs(y = "Brier Score for Cause 1 Predictions") +
  scale_colour_manual(values=cbPalette) + facet_grid(. ~ title)
dev.off()