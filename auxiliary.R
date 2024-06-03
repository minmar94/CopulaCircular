# Packages ----------------------------------------------------------------
require(tidyverse)
require(magrittr)
require(circular)

# Auxiliary ---------------------------------------------------------------
# CDF of the Wrapped Cauchy 
pwrappedcauchy <- function(q, mu = circular(0), rho = exp(-1), from = NULL){
  
  q <- conversion.circular(q, units = "radians", zero = 0, 
                           rotation = "counter", modulo = "2pi")
  mu <- conversion.circular(mu, units = "radians", zero = 0, 
                            rotation = "counter", modulo = "2pi")
  if (is.null(from)) {
    from <- mu - pi
  }
  else {
    from <- conversion.circular(from, units = "radians", 
                                zero = 0, rotation = "counter", modulo = "2pi")
  }
  n <- length(q)
  mu <- (mu - from)%%(2 * pi)
  q <- (q - from)%%(2 * pi)
  intDwrappedcauchylRad <- function(q) {
    if (is.na(q)) {
      return(NA)
    }
    else {
      return(integrate(circular:::DwrappedcauchyRad, mu = mu, rho = rho, lower = 0, upper = q)$value)
    }
  }
  value <- sapply(X = q, FUN = intDwrappedcauchylRad)
  return(value)
}

# Quantile Function of the Wrapped Cauchy
qwrappedcauchy <- function(p, mu, rho, from = NULL, tol = .Machine$double.eps^(0.6)){
  epsilon <- 10 * .Machine$double.eps
  if (any(p > 1) | any(p < 0)) 
    stop("p must be in [0,1]")
  if (is.circular(mu)) {
    datacircularp <- circularp(mu)
  }
  else {
    datacircularp <- list(type = "angles", units = "radians", 
                          template = "none", modulo = "asis", zero = 0, rotation = "counter")
  }
  
  mu <- conversion.circular(mu, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  if (is.null(from)) {
    from <- mu - pi
  }
  else {
    from <- conversion.circular(from, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  }
  n <- length(p)
  mu <- (mu - from)%%(2 * pi)
  if (rho < 0 | rho > 1) 
    stop("rho must be between 0 and 1")
  
  zeroPwrappedcauchyRad <- function(x, p, mu, rho) {
    if (is.na(x)) {
      y <- NA
    }
    else {
      y <- integrate(circular:::DwrappedcauchyRad, mu = mu, rho = rho, lower = 0, upper = x)$value - p
    }
    return(y)
  }
  value <- rep(NA, length(p))
  sem <- options()$show.error.messages
  options(show.error.messages = FALSE)
  for (i in 1:length(p)) {
    res <- try(uniroot(zeroPwrappedcauchyRad, p = p[i], mu = mu, 
                       rho = rho, lower = 0, upper = 2 * pi - epsilon, tol = tol))
    if (is.list(res)) {
      value[i] <- res$root
    }
    else if (p[i] < 10 * epsilon) {
      value[i] <- 0
    }
    else if (p[i] > 1 - 10 * epsilon) {
      value[i] <- 2 * pi - epsilon
    }
  }
  options(show.error.messages = sem)
  value <- value + from
  return(value)
}


##### Inference for margins estimation #####
copula_estimate <- function(y, families){
  
  n <- nrow(y)
  J <- ncol(y)
  z <- matrix(0, ncol = J, nrow = n)
  ydens <- matrix(0, ncol = J, nrow = n)
  pars <- list()
  
  for(j in 1:J){
    pfun <- match.fun((paste0("p", families[j])))
    dfun <- match.fun((paste0("d", families[j])))
    if(families[j] %in% c("vonmises", "wrappednormal", "wrappedcauchy")){
      mle.fit <- match.fun((paste0("mle.", families[j])))
      fit <- mle.fit(atan2(sin(y[,j]), cos(y[,j])))
      pars[[j]] <- c(as.numeric(fit[2]), as.numeric(fit[3])) 
      mu <- fit$mu
      mu <- atan2(sin(mu), cos(mu))
      
      if(families[j] == "wrappedcauchy"){
        sigma2 <- fit$rho
        sigma2 <- (tanh(sigma2)+1)/2
      }else{
        sigma2 <- fit$kappa
      }
      z[,j] <- qnorm(pfun(y[,j], mu, sigma2))
      ydens[,j] <- dfun(y[,j], mu, sigma2)
    }
    if(families[j] %in% c("weibull", "gamma", "lnorm")){
      fit <- fitdistrplus::fitdist(y[,j], distr = families[j])
      mu <- fit$estimate[1]
      sigma2 <- fit$estimate[2]
      pars[[j]] <- fit$estimate
      z[,j] <- qnorm(pfun(y[,j], mu, sigma2))
      ydens[,j] <- dfun(y[,j], mu, sigma2)
    }
  }
  
  Omega <- cor(z)
  num <- mvtnorm::dmvnorm(z, sigma = Omega)
  # den <- -sum(apply(z, 2, dnorm))
  den <- eval(parse(text=paste0("dnorm(z[,", 1:J, "])", collapse = "*")))
  marg <- eval(parse(text=paste0("ydens[,", 1:J, "]", collapse = "*")))
  
  return(list(copula = num/den*marg,
              omega = Omega,
              pars = pars,
              z = z))
  
}

##### Density #####
copula_compute <- function(y, families, pars, Omega){
  
  n <- nrow(y)
  J <- ncol(y)
  z <- matrix(0, ncol = J, nrow = n)
  ydens <- matrix(0, ncol = J, nrow = n)
  
  for(j in 1:J){
    pfun <- match.fun((paste0("p", families[j])))
    dfun <- match.fun((paste0("d", families[j])))
    parcur <- pars[[j]]
    mu <- parcur[1]
    sigma2 <- parcur[2]
    if(families[j] %in% c("vonmises", "wrappednormal", "wrappedcauchy")){
      y[,j] <- atan2(sin(y[,j]), cos(y[,j]))
      z[,j] <- qnorm(pfun(y[,j], mu, sigma2, from = mu-pi))
      ydens[,j] <- dfun(y[,j], mu, sigma2)
    }
    if(families[j] %in% c("weibull", "gamma", "lnorm")){
      z[,j] <- qnorm(pfun(y[,j], mu, sigma2))
      ydens[,j] <- dfun(y[,j], mu, sigma2)
    }
  }
  
  num <- mvtnorm::dmvnorm(z, sigma = Omega)
  den <- eval(parse(text=paste0("dnorm(z[,", 1:J, "])", collapse = "*")))
  marg <- eval(parse(text=paste0("ydens[,", 1:J, "]", collapse = "*")))
  
  return(num/den*marg)
  
}

##### Likelihood #####

### KRIGING FOR CIRCULAR DATA WITH COPULAS ###
copula_lik_spatialkrig <- function(par, obs, Dist, covtype){
  
  n <- ncol(Dist)
  mu <- par[2]
  kappa <- exp(par[3])
  z <- qnorm(pvonmises(obs, mu, kappa)) # si assume identica distribuzione delle marginali
  I <- diag(1, n)
  ldens <- suppressWarnings(sum(dvonmises(obs, mu, kappa, log = T)))
  
  if(covtype == "exponential"){
    rho <- exp(par[1])
    Omega <- exp(-rho*Dist)
  }
  if(covtype == "gaussian"){
    rho <- par[1]
    Omega <- exp(-rho^2*Dist^2)
  }
  
  Omegamenouno <- solve(Omega)
  ldetomega <- log(det(Omegamenouno))
  ll <- .5 * ldetomega + .5*as.numeric(t(z)%*%(I-Omegamenouno)%*%z) + ldens
  # print(par)
  return(-ll)
  
}



### AR(1) ###
copula_lik_idar <- function(par, obs, families, modtype){
  
  J <- length(families)
  rhos <- c()
  
  if(modtype == "unstruct"){
    for(j in 1:(J*(J+1)/2 - J)){
      rhos[j] <- par[j]
    }
    dens_par <- par[-c(1:(J*(J+1)/2 - J))]
    Omega <- matrix(NA, nrow = J, ncol = J)
    Omega[lower.tri(Omega)] <- rhos
    Omega <- t(Omega)
    Omega[lower.tri(Omega)] <- rhos
    diag(Omega) <- 1
    # print(Omega)
  }
  
  # def Sigma
  if(modtype == "ar1"){
    # def Sigma
    dens_par <- par[-1]
    rho <- par[1]
    Omega.offd <- c(rho, rho^2, rho^3, rho^4, 
                    rho, rho^2, rho^3,
                    rho, rho^2, 
                    rho)
    Omega <- matrix(NA, J, J)
    Omega[lower.tri(Omega)] <- Omega.offd
    Omega <- t(Omega)
    Omega[lower.tri(Omega)] <- Omega.offd
    
    diag(Omega) <- 1
  }
  
  if(all(families == "vonmises")){
    pars_marg <- rep(paste0(paste0("dens_par[", 1, "], "), paste0("exp(dens_par[", 2, "])")), J)
  }else if(all(families == "wrappedcauchy")){
    pars_marg <- rep(paste0(paste0("dens_par[", 1, "], "), paste0("plogis(dens_par[", 2, "])")), J)
  }else{ #if(all(families == "weibull")){
    pars_marg <- rep(paste0(paste0("exp(dens_par[", 1, "]),"), paste0("exp(dens_par[", 2, "])")), J)
  }
  den_fam <- paste0("qnorm(p", families, "(obs[,", 1:J, "], ", pars_marg, "))")
  z <- suppressWarnings(sapply(den_fam, \(d) eval(parse(text = d))))
  
  cop <- "mvtnorm::dmvnorm(z, sigma = Omega, log = T)"
  
  den <- paste0("dnorm(z[,", 1:J, "], log = T)", collapse = " + ")
  den <- paste0("(", den, ")")
  names_marg <- paste0("d", families)
  
  marg <- paste0("log(", names_marg, "(obs[, ", 1:J, "], ", pars_marg, "))", collapse = "+")
  
  f <- eval(parse(text = paste0("sum(",paste0(cop, "-", den, " + ", marg), ")")))
  return(-f)
  
}

