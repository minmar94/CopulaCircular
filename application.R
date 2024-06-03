
# Packages ----------------------------------------------------------------
require(tidyverse)
require(magrittr)
require(circular)
require(spdep)
require(sf)

source("auxiliary.R")

# AR(1) -------------------------------------------------------------------
load("data/sandhoppers.RData")
y <- as.matrix(y)

### Fit
families <- rep("vonmises", 5)

# Inference for margins to initialize
fit_vm <- copula_estimate(y, families)

# Autoregressive - id
par <- c(fit_vm$omega[1,2], fit_vm$pars[[1]][1], log(fit_vm$pars[[1]][2]))
fit_hypertorus_optim_id <- nlminb(start = par, objective = copula_lik_idar, 
                                  obs = y, modtype = "ar1",
                                  families = families,  
                                  control = list(trace = 5))
# Hessian
hess <- numDeriv::hessian(copula_lik_idar, x = fit_hypertorus_optim_id$par, 
                          obs = y, modtype = "ar1",
                          families = families)
# Parameters and standard errors
fit_hypertorus_optim_id$par
diag(solve(hess)) %>% sqrt %>% round(3)


# Predictions
y1 <- seq(-pi, pi, length.out = 2000)
y1a <- seq(-pi, -0.1, length = 500)
y1b <- seq(-0.1, 0.1, length = 3000)
y1c <- seq(0.1, pi, length = 500)
y1 <- unique(sort(c(y1a, y1b, y1c)))

# z-scores with the estimated parameters
z1 <- qnorm(pvonmises(y1, fit_hypertorus_optim_id$par[2], exp(fit_hypertorus_optim_id$par[3])))
y2 <- qvonmises(pnorm(fit_hypertorus_optim_id$par[1]*z1), 
                fit_hypertorus_optim_id$par[2], 
                exp(fit_hypertorus_optim_id$par[3]))

q1 <- qnorm(0.025, mean = fit_hypertorus_optim_id$par[1]*z1, sd = sqrt(1 - fit_hypertorus_optim_id$par[1]^2))
q2 <- qnorm(0.975, mean = fit_hypertorus_optim_id$par[1]*z1, sd = sqrt(1 - fit_hypertorus_optim_id$par[1]^2))
# Confidence bands
yq1 <- qvonmises( pnorm(q1), fit_hypertorus_optim_id$par[2], exp(fit_hypertorus_optim_id$par[3]) )
yq2 <- qvonmises( pnorm(q2), fit_hypertorus_optim_id$par[2], exp(fit_hypertorus_optim_id$par[3]) )

# Kriging -----------------------------------------------------------------

load("data/seacurrent.RData")
iho <- st_read("data/iho.shp") %>% st_transform(crs = "+proj=utm +zone=32 +units=km")

Ds <- as.matrix(dist(scale(coord2_utm))) # Euclidean distances

#### Fit
# Inference for margins to initialize
marg_pars <- mle.vonmises(yobs)
mu <- marg_pars$mu
kappa <- marg_pars$kappa

# Compute the z-scores
z <- qnorm(pvonmises(yobs, mu, kappa))

# Joint maximization
fitkrig <- nlminb(start = c(runif(1, 0, 1), mu, log(kappa)), objective = copula_lik_spatialkrig, 
                  Dist = Ds, obs = yobs, covtype = "exponential",
                  control = list(trace = 1))
# Hessian
hess_krig <- numDeriv::hessian(copula_lik_spatialkrig, fitkrig$par, Dist = Ds, obs = yobs, covtype = "exponential")
# Standard errors
stderrs <- diag(solve(hess_krig)) %>% sqrt %>% round(4)

# Parameter estimates
fitkrig$par

# Delta method for the correlation which is parametrized on the log-scale
exp(fitkrig$par[1])
exp(fitkrig$par[1])*sqrt(solve(hess_krig)[1,1])

# Prediction at unobserved locations
# Build a grid of unobserved locations 
xutm <- seq(min(coord2_utm[,1])-1, max(coord2_utm[,1])+2, length.out = 25)
yutm <- seq(min(coord2_utm[,2])-.5, max(coord2_utm[,2])+2, length.out = 25)
xyutm <- expand.grid(xutm, yutm) %>% as.data.frame() %>% 
  st_as_sf(coords = c("Var1", "Var2"), crs = "+proj=utm +zone=32 +units=km")

ids <- st_intersects(xyutm, iho, sparse = F)
xyutm <- xyutm[ids[,1],] %>% st_coordinates()

xyutmall <- rbind(as.matrix(xyutm), as.matrix(coord2_utm))
# Distances
Ds <- as.matrix(dist(scale(xyutmall)))

# Compute z-scores with the estimated parameters
z <- qnorm(pvonmises(yobs, fitkrig$par[2], exp(fitkrig$par[3])))

# Function to compute predictions
mucondkrig <- function(rho, D, obs){
  D <- exp(-rho*D)
  nall <- ncol(D)
  nobs <- length(obs)
  npred <- nall-nobs
  sigma12 <- D[1:npred, (npred+1):nall, drop=F]
  sigma22 <- D[(npred+1):nall, (npred+1):nall]
  return(c(sigma12%*%solve(sigma22)%*%obs))
}

# Median
ypred <- as.numeric(qvonmises(pnorm(mucondkrig(exp(fitkrig$par[1]), Ds, z)), fitkrig$par[2], exp(fitkrig$par[3])))
# Quantiles
ypredq1 <- as.numeric(qvonmises(pnorm(mucondkrig(exp(fitkrig$par[1]-1.96*stderrs[1]), Ds, z)), fitkrig$par[2], exp(fitkrig$par[3])))
ypredq2 <- as.numeric(qvonmises(pnorm(mucondkrig(exp(fitkrig$par[1]+1.96*stderrs[1]), Ds, z)), fitkrig$par[2], exp(fitkrig$par[3])))

# Confidence arc length
r1 <- 1
errs <- r1*(acos(1 - sqrt((sin(ypredq1)-sin(ypredq2))^2 + (cos(ypredq1)-cos(ypredq2))^2)/(2*r1^2)))

