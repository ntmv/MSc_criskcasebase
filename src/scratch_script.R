# Load the MASS package
library(MASS)

# Set the number of variables and covariate blocks
num_variables <- 16
num_blocks <- 4

# Set the number of variables per block
variables_per_block <- num_variables / num_blocks

# Set the correlation values for each covariate block
correlation_values <- c(0.5, 0.35, 0.05, 0.32)

# Initialize an empty covariance matrix
covariance_matrix <- matrix(0, nrow = num_variables, ncol = num_variables)

# Generate the covariance matrix with block correlations
for (i in 1:num_blocks) {
  start_index <- (i - 1) * variables_per_block + 1
  end_index <- i * variables_per_block
  covariance_matrix[start_index:end_index, start_index:end_index] <- correlation_values[i]
}

diag(covariance_matrix) <- rep(1, length(diag(covariance_matrix) ))


# Generate the correlated data
simulated_data <- mvrnorm(n = 100, mu = rep(0, num_variables), Sigma = covariance_matrix, empirical = TRUE)

# Calculate the correlation matrix of simulated data
simulated_correlation <- cor(simulated_data)

# Print the simulated data
print(simulated_data)

# Print the correlation matrix of simulated data

################################# Playing with CIF function ##########
data(bmtcrr)
time_points <- sort(unique(test$time))

model1 <- fitSmoothHazard(status ~ log(time) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
                            X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20, 
                          data = train, 
                          ratio = 100,
                          time = "time")
summary(model1)

# Estimate absolute risk curve
risk_cb <- absoluteRisk(
  object = model1, time = time_points,
  method = "numerical", newdata = test
)

risk_cb <- as.data.frame(rowMeans(risk_cb))

risk_cb <- as.data.frame(cbind(risk_cb[-1,], time_points))

ggplot(risk_all, aes(x = Time, y = Risk, colour = Method)) + geom_line() + theme_bw()


model_fg <- comp.risk(Event(time, status) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
                        X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20,
                      data = sim.data, cause = 1, model = "fg")

# Estimate CI curve
risk_fg <- predict(model_fg, newdata = test, times = time_points)


plot(risk_fg)
print(simulated_correlation)

risk_fg <- as.data.frame(rowMeans(risk_fg$RR))


risk_all <- dplyr::bind_rows(
  data.frame(
    Time = time_points,
    Method = "Case-base",
    Risk = risk_cb$V1,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Time = time_points,
    Method = "Nelson-Aalen Estimator",
    Risk = cif$`1 1`$est,
    stringsAsFactors = FALSE 
  )
  )

############################ Plot competitor models risk curves ###############
library(casebase)
library(future.apply)
library(glmnet)
library(mtool)
library(timereg)
library(parallel)
library(tictoc)
library(tidyverse)
library(riskRegression)
library(cmprsk)
library(survsim)

n <- 400
p <- 20
beta <- list(c(0.5, 0.5, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0.5), 
             c(0.2, 0.2, 0.2, 0, 0.2, 0, 0, 0.2, 0, 0.2, 0, 0, 0.2, 0, 0.2, 0, 0, 0.2, 0, 0.2))
dist.ev <- c("weibull", "weibull")
anc.ev <- c(0.8, 0.3)
beta0.ev <- c(0.1, 0.1)
x <- rep(list(c("normal", 0, 1)), p)

# Generating survival data 
sim.data <- crisk.sim(foltime = 4, dist.ev = dist.ev, 
                      anc.ev = anc.ev, beta0.ev = beta0.ev, beta0.cens = 0.05, anc.cens = 4, nsit = 2, 
                      beta = beta, x = x, n = n)

# fix status variable
sim.data$cause <- with(sim.data, ifelse(is.na(sim.data$cause), 0, sim.data$cause))
colnames(sim.data)[grepl("x", colnames(sim.data))]   <- paste0("X", seq_len(p))

# Format data
sim.data <- sim.data %>%
  dplyr::select(-nid, -status, -start, -stop, -z) %>%
  rename(status = cause)

# Average estimates of incidence and censoring rate
prop.table(table(sim.data$status))

# True cumulative incidence 
cif <- cuminc(ftime = sim.data$time, fstatus = sim.data$status)
cif


################### Cox-Regression Brier score #########################
# Split into training and validation sets (stratified)
train.index <- caret::createDataPartition(sim.data$status, p = 0.70, list = FALSE)
train <- sim.data[train.index,]
test <- sim.data[-train.index,]
######################### Cause-specific proportional hazards model ###############
# Absolute Risks
newx <- model.matrix(status ~ . - time,
                     data = test)[, -c(1)]
times <- sort(unique(test$time))

# Censor competing event
y_train <- Surv(time = train$time, event = train$status == 1)

x_train <- model.matrix(~ . -time -status, data = train)[, -1] 

# Censor competing event
y_test <- Surv(time = test$time, event = test$status == 1)

x_test <- model.matrix(~ . -time -status, data = test)[, -1] 

# Fit cause-specific cox model with glmnet on training set 
cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 1, standardize = TRUE)

nonzero_coef_cox <- coef(cox_val_min, s = "lambda.min")

# Creating a new dataset that only contains the covariates chosen through glmnet
# Taking the coefficient estimates for later use
nonzero_covariate_cox <- predict(cox_val_min, type = "nonzero", s = "lambda.min")
cleanCoxData <- as.data.frame(cbind(y_train, x_train[, nonzero_covariate_cox$X1]))

# Fitting a cox model using regular estimation, however we will not keep it.
# this is used more as an object place holder.
coxNet <- survival::coxph(Surv(time = time, event = status == 1) ~ ., 
                          data = train, x = TRUE)

# The coefficients of this object will be replaced with the estimates from the
# original coxNet. Doing so makes it so that everything is invalid aside from
# the coefficients. In this case, all we need to estimate the absolute risk is
# the coefficients. Std. error would be incorrect here, if we were to draw error
# bars.
coxNet_coefnames <- names(coxNet$coefficients)
coxNet$coefficients <- nonzero_coef_cox
names(coxNet$coefficients) <- coxNet_coefnames

# Absolute risks 
newx <- model.matrix(status ~ . - time -status,
                     data = test)[, -c(1)] 

times <- sort(unique(test$time))[c(1:97)]

# 2. Penalized models
brierPenalized <- Score(list("Pen. Cox" = coxNet),
                        data = test,
                        formula = Hist(time, status != 0) ~ 1, summary = NULL, 
                        se.fit = FALSE, metrics = "brier", contrasts = FALSE, 
                        times = times)

# Brier Cox 
briercox <- Score(list("Cox" = coxNet),
                        data = test,
                        formula = Hist(time, status == 1) ~ 1, summary = NULL, 
                        se.fit = FALSE, metrics = "AUC", contrasts = FALSE, 
                        times = times)

# 4. Kaplan-Meier
KM <- survival::coxph(Surv(time = time, event = status == 1) ~ 1, data = test,
                      x = TRUE)

# First fix NA coefficients, then compute Brier score
brierCoxKM <- Score(list("KM" = KM), data = test,
                    formula = Hist(time, status != 0) ~ 1, summary = NULL,
                    se.fit = FALSE, metrics = "AUC", contrasts = FALSE, 
                    times = times)


##################### Try plotting casebase Brier score #################
newx <- model.matrix(status ~ . - time -status,
                     data = test)[, -c(1)] 

newx1 <- model.matrix(status ~ .  -status,
                     data = test)[, -c(1)] 

fit_val_min
# Fit same model with fit smooth hazard
model_cb <- fitSmoothHazard(
  status ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + 
    X17 + X18 + X19 + X20 + log(time),
  data = train,
  ratio = 100,
  time = "time"
)

coef_model_cb <- coef(model_cb)

names_coef <- names(coef(model_cb))

coef_model_cb <- c(-1.10323798 , -0.858481455, 
                   0, -0.123129719, 0, -0.013260420, 0, 0, 0.04555622, 0, -0.17006286
                   , -0.237898485, 0, 0, 0, 0.070789893, 0, 0, 0, -0.180621619, 
                   -0.15629495, -0.478715058, 0, 0, 0.08266276, 0, 0, 0, 0, 0.005199849, 
                   0, -0.239810796, rep(0, 6), 0.04132416, 0, -0.16227504, -0.155466505, 
                   -1.04247451, -0.672779588)
names(coef_model_cb) <- names_coef

model_cb@coefficients <-  coef_model_cb

# Estimate absolute risk curve
risk_cb <- absoluteRisk(
  object = model_cb, time = time_points,
  method = "numerical", newdata = test
)

# Fine-Gray model 
library(timereg)
model_fg <- comp.risk(Event(time, status) ~ const(X1) + const(X2) + const(X3) + 
                        const(X4) + const(X5) + const(X6) + const(X7) + const(X8) + 
                        const(X9) + const(X10) + const(X11) + const(X12) + 
                        const(X13) + const(X14) + const(X15) + const(X16) + 
                        const(X17) + const(X18) + const(X19) + const(X20),
                      data = train, cause = 1, model = "fg")

time_points <- seq(0, 60, length.out = 50)
risk_fg <- predict(model_fg, test, times = time_points)

times <- sort(unique(test$time))
# Brier score 
brierfinegray <- Score(list("Fine-Gray" = model_fg),
                        data = cbind(subset(test, select = c(time, status)),
                                     as.data.frame(newx)),
                        formula = Hist(time, status != 0) ~ 1, summary = NULL, 
                        se.fit = FALSE, metrics = "AUC", contrasts = FALSE, 
                        times = times)

# Brier score case-base
briercasebase <- Score(list("Case-base" = model_cb),
                       data = cbind(subset(test, select = c(time, status)),
                                    as.data.frame(newx)),
                       formula = Hist(time, status != 0) ~ 1, summary = NULL, 
                       se.fit = FALSE, metrics = "AUC", contrasts = FALSE, 
                       times = times)
#######################################################################
# Plot Brier score 
# Combine scores
data_brier <- bind_rows(
# brierPenalized$Brier$score %>% 
 #  mutate(model = as.character(model)) %>% 
#   filter(model != "Null model"),
  brierfinegray$AUC$score %>% 
  mutate(model = as.character(model)) %>% 
  filter(model != "Null model"),
 briercasebase$AUC$score %>%
   mutate(model = as.character(model)) %>% 
  filter(model != "Null model"),  
 brierCoxKM$AUC$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model"),
  briercox$AUC$score %>% 
    mutate(model = as.character(model)) %>% 
    filter(model != "Null model")
) %>% 
  mutate(model = factor(model, levels = c("Fine-Gray", "Case-base", "KM","Cox")))

data_brier_mean <- data_brier %>%
  group_by(model) %>%
  summarize(mean_brier = mean(Brier))

write.csv(data_brier_mean, file = "~/Desktop/june15_images/brier_iidnonsparse.csv")

png(filename = "~/Desktop/AUC_nsparse.png", res = 300, height = 12, width = 20, units = "cm") 
 ggplot(data = data_brier, aes(x = times, y = AUC, col = model)) +
  geom_line() +
  xlab("Follow-up time (years)") +
  ylab("Brier score") +
  labs(color = "Models") +
  theme_bw() 
dev.off()
 
 
################################### Case-base Brier score for competing risks ################
 # We need this chunk until this method reaches the CRAN version of riskRegression
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
 
#############################################################
 time_point <- sort(unique(test$time))
 
 newx <- model.matrix(status ~ . - time -status,
                      data = test)[, -c(1)] 
 
 newx <- as.data.frame(newx)
 
 # get all x covariate values from test dataset
 times <- sort(unique(test$time)) # get all unique survival times (d.time) and sort in ascending order
 
 
 risk_fg <- predict(model_fg, test, times = time_point)
 
 # Estimate absolute risk curve
 risk_cb <- absoluteRisk(
   object = model_cb, time = time_point,
   method = "numerical", newdata = test
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
                        Status = factor(status))

model_cox <- coxph(Surv(time, Status) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
                     X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20,
                   data = bmtcrr_cox, id = id)

risk_cox <- survival::survfit(model_cox, type = "breslow", 
                              newdata = test)

 
 risk_all <- dplyr::bind_rows(
   data.frame(
     Time = time_point,
     Method = "Case-base",
     Risk = risk_cb,
     stringsAsFactors = FALSE
   ),
  data.frame(
    Time = risk_cox$time,
   Method = "Cox",
   Risk = risk_cox[,2]$pstate[,1,1],
   stringsAsFactors = FALSE
 ),
   data.frame(
     Time = time_point,
     Method = "Fine-Gray",
     Risk = risk_fg$`colMeans(risk_fg$P1)`,
     stringsAsFactors = FALSE
   )
 ) %>% 
   dplyr::filter(Time <= 60)
 
png(filename = "~/Desktop/cif_iidnonsparse.png", res = 300, height = 12, width = 20, units = "cm") 
 ggplot(risk_all, aes(x = Time, y = Risk, colour = Method)) +
   # geom_line for smooth curve
   geom_line(data = dplyr::filter(risk_all, Method == "Case-base")) +
   # geom_step for step function
   geom_step(data = dplyr::filter(risk_all, Method != "Case-base")) +
   ylim(c(0, 1)) +
 #  paper_gg_theme +  
   xlab("Time (in Months)") +
   ylab("Relapse risk") + theme_bw()
 dev.off()
 ############### Understanding survsim parameters #######################
 library("survsim")
 
 dat1 <- crisk.sim(n=1000, foltime=500, dist.ev=c("weibull","weibull"),
                   anc.ev=c(1, 1), beta0.ev=c(1, 1), dist.cens= "unif",
              anc.cens=600, beta0.cens=550, z=NULL, beta=list(c(1.5,1)),
                   x=list(c("bern", 0.5)), nsit=2)
 
 coxCause1 <- coxph(Surv(time, cause== 1) ~ x, data=dat1)
 
 WeibullC1 <-survreg(Surv(time, cause==2) ~ x, data=dat1,  dist="weibull")
 
                     