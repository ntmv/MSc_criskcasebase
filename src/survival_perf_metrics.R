##################### Prediction function for casebase penalized model ################
predict_CompRisk_penalized <- function(fit_object, newdata = NULL) {
  X <- as.matrix(cbind(newdata[, .SD, .SDcols = grep("X", colnames(newdata), value = TRUE)], 
                       ftime = newdata$ftime, Intercept = 1))
  preds <- as.matrix(X %*% fit_object$coefficients)
  colnames(preds) <- paste0("log(mu[,",
                            seq(2, length(3)),
                            "]/mu[,1])")
  return(preds)
}


########################## Survival analysis performance metrics #############
### Cause-specific cumulative incidence function for penalized case-base fit
absoluteRisk.penalized <- function(fit_object, time, newdata,
                                  nsamp = 100,
                                  addZero = TRUE) {
  typeEvents <- c(0, 1, 2)
  ###################################################
  # In competing risks, we can get a cumulative
  # incidence function using a nested double integral
  # subdensity f_j = lambda_j * Survival
  # F_j = P(T <= t, J = j : covariates) = int_0^t f_j
  ###################################################
  time_ordered <- unique(c(0, sort(time)))
  # Create array to store output
    output <- matrix(NA, ncol = nrow(newdata) + 1,
                     nrow = length(time_ordered))
    output[, 1] <- time_ordered
    colnames(output) <- c("time", rep("", nrow(newdata)))
    rownames(output) <- rep("", length(time_ordered))
    output[1, -1] <- 0
  # Compute subdensities
    # Compute points at which we evaluate integral
    knots <- seq(0, max(time_ordered),
                 length.out = (length(time_ordered) - 1) * nsamp)
    knots <- unique(sort(c(knots, time_ordered)))
    
    for (j in seq_len(nrow(newdata))) {
      # Extract current obs
      current_obs <- newdata[j, , drop = FALSE]
      # Create data.table for prediction
      newdata2 <- data.table::data.table(current_obs)
      newdata2 <- newdata2[rep(1, length(knots))]
      newdata2[, ftime := knots]
      newdata2[, "offset" := 0]
      # Compute all values for all hazards
      lambdas <- exp(predict_CompRisk_penalized(fit_object, newdata2))
      lambdas[which(lambdas %in% c(Inf, -Inf))] <- 0
      OverallLambda <- rowSums(lambdas)
      survFunction <- exp(-trap_int(knots, OverallLambda))
        # Only compute first subdensity
        subdensity <- lambdas[, 1] * survFunction
        pred <- trap_int(knots, subdensity)[knots %in% c(0, time)]
        output[, j + 1] <- pred
        } 
  return(output)
}

#########################################
# Streamlined version of pracma::cumtrapz
trap_int <- function(x, y) {
  x <- as.matrix(c(x))
  m <- length(x)
  y <- as.matrix(y)
  n <- ncol(y)
  dt <- kronecker(matrix(1, 1, n), 0.5 * diff(x))
  ct <- apply(dt * (y[1:(m - 1), ] + y[2:m, ]), 2, cumsum)
  return(rbind(0, ct))
}
############################################

