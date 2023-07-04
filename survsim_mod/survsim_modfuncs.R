######################### Survsim modified functions ###############################
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
# Survsim's function modified for block correlations
crisk.sim_block <-
  function (n, p, nblocks, foltime, dist.ev, anc.ev, beta0.ev, dist.cens = "weibull", 
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
    # Set the number of variables per block
    vpb <- p/nblocks
    # Set the correlation values for each covariate block
    correlation_values <- c(0.5, 0.35, 0.05, 0.32)
    # Initialize empty matrix
    correlation_matrix <- matrix(0, nrow = p, ncol = p)
    # Generate the covariance matrix with block correlations
    for (i in 1:nblocks) {
      start_index <- (i - 1) * vpb + 1
      end_index <- i * vpb
      correlation_matrix[start_index:end_index, start_index:end_index] <- correlation_values[i]
    }
    # Diagonal elements should be 1
    diag(correlation_matrix) <- rep(1, length(diag(correlation_matrix)))
    
    X <- mvtnorm::rmvnorm(n, sigma = correlation_matrix)
    
    for (i in 1:n) {
      eff <- mvtnorm::rmvnorm(1, sigma = correlation_matrix)
      
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