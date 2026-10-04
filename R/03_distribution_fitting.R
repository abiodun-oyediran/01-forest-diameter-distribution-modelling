# =============================================================================
# [REUSED] R/03_distribution_fitting.R
# Purpose : Define the custom PDFs/CDFs and fit 12 candidate distributions to
#           the DBH data by maximum likelihood.
# Inputs  : objects from R/01_data_cleaning.R (data, xmin, xmax, mean_data)
# Outputs : list `fits` (one element per successfully fitted model)
# Author  : Oyediran Abiodun
# Date    : 2026-10-04
# Notes   : [adapted from my script: custom PDF/CDF functions, fit_via_optim,
#           fit calls]. Two corrections, marked [FIX]: Log-Gamma definition (B1)
#           and the Burr 3P starting value (B6). See docs/methodology.md.
# =============================================================================

if(!exists("data") || !exists("xmin") || !exists("mean_data")) {
  stop("Run R/01_data_cleaning.R first (objects 'data', 'xmin', 'mean_data' not found).")
}

# ---------- CUSTOM PDF/CDF FUNCTIONS ----------
dbeta4 <- function(x, shape1, shape2, a, b, log = FALSE){
  t <- (x - a) / (b - a)
  valid <- (x >= a & x <= b)
  out <- rep(-Inf, length(x))
  if(any(valid)){
    ld <- dbeta(t[valid], shape1, shape2, log = TRUE) - log(b - a)
    out[valid] <- ld
  }
  if(!log) return(exp(out)) else return(out)
}

pbeta4 <- function(q, shape1, shape2, a, b){
  t <- (q - a) / (b - a)
  p <- ifelse(q < a, 0, ifelse(q > b, 1, pbeta(t, shape1, shape2)))
  return(p)
}

dburr4 <- function(x, c, k, scale, loc, log = FALSE){
  z <- (x - loc) / scale
  zpos <- x > loc
  out_log <- rep(-Inf, length(x))
  if(any(zpos)){
    out_log[zpos] <- log(c) + log(k) + (c - 1) * log(z[zpos]) - log(scale) - (k + 1) * log1p(z[zpos]^c)
  }
  if(!log) return(exp(out_log)) else return(out_log)
}

pburr4 <- function(q, c, k, scale, loc){
  z <- (q - loc) / scale
  p <- ifelse(q <= loc, 0, 1 - (1 + z^c)^(-k))
  return(p)
}

dburr3 <- function(x, c, k, scale, log = FALSE){
  dburr4(x, c=c, k=k, scale=scale, loc=0, log=log)
}

pburr3 <- function(q, c, k, scale){
  pburr4(q, c=c, k=k, scale=scale, loc=0)
}

dgamma3 <- function(x, shape, scale, loc, log = FALSE){
  pos <- x > loc
  out_log <- rep(-Inf, length(x))
  if(any(pos)){
    out_log[pos] <- dgamma(x[pos] - loc, shape=shape, scale=scale, log=TRUE)
  }
  if(!log) return(exp(out_log)) else return(out_log)
}

pgamma3 <- function(q, shape, scale, loc){
  ifelse(q <= loc, 0, pgamma(q - loc, shape=shape, scale=scale))
}

# [FIX B1] Log-Gamma: ln(X) ~ Gamma(alpha, scale), defined on X (DBH) itself.
# The source script defined the density of ln(Y) with Y ~ Gamma, fitted it to
# log(data) and then evaluated the CDF on raw DBH, so fit and tests disagreed.
dloggamma_custom <- function(x, alpha, scale, log = FALSE){
  pos <- x > 1
  out_log <- rep(-Inf, length(x))
  if(any(pos)){
    out_log[pos] <- dgamma(log(x[pos]), shape=alpha, scale=scale, log=TRUE) - log(x[pos])
  }
  if(!log) return(exp(out_log)) else return(out_log)
}

ploggamma_custom <- function(q, alpha, scale){
  ifelse(q <= 1, 0, pgamma(log(pmax(q, 1)), shape=alpha, scale=scale))
}

dloglogistic3 <- function(x, shape, scale, loc, log = FALSE){
  z <- (x - loc) / scale
  out_log <- rep(-Inf, length(x))
  pos <- x > loc
  if(any(pos)){
    out_log[pos] <- log(shape) - log(scale) + (shape - 1) * log(z[pos]) - 2 * log1p(z[pos]^shape)
  }
  if(!log) return(exp(out_log)) else return(out_log)
}

ploglogistic3 <- function(q, shape, scale, loc){
  z <- (q - loc) / scale
  p <- ifelse(q <= loc, 0, 1 / (1 + z^(-shape)))
  return(p)
}

dlnorm3 <- function(x, meanlog, sdlog, loc, log = FALSE){
  pos <- x > loc
  out_log <- rep(-Inf, length(x))
  if(any(pos)){
    out_log[pos] <- dlnorm(x[pos] - loc, meanlog=meanlog, sdlog=sdlog, log=TRUE)
  }
  if(!log) return(exp(out_log)) else return(out_log)
}

plnorm3 <- function(q, meanlog, sdlog, loc){
  ifelse(q <= loc, 0, plnorm(q - loc, meanlog=meanlog, sdlog=sdlog))
}

dweibull3 <- function(x, shape, scale, loc, log = FALSE){
  pos <- x > loc
  out_log <- rep(-Inf, length(x))
  if(any(pos)){
    out_log[pos] <- dweibull(x[pos] - loc, shape=shape, scale=scale, log=TRUE)
  }
  if(!log) return(exp(out_log)) else return(out_log)
}

pweibull3 <- function(q, shape, scale, loc){
  ifelse(q <= loc, 0, pweibull(q - loc, shape=shape, scale=scale))
}

# ---------- GENERIC MLE FITTING ----------
fit_via_optim <- function(dfunc, start, data, lower = NULL, upper = NULL, control = list()){
  par0 <- unlist(start)
  if(is.null(names(par0)) || any(names(par0) == "")) stop("start must be a named list")
  
  if(is.null(lower)) lower_vec <- rep(-Inf, length(par0)) else lower_vec <- rep(lower, length.out=length(par0))
  if(is.null(upper)) upper_vec <- rep(Inf, length(par0)) else upper_vec <- rep(upper, length.out=length(par0))
  
  neglog <- function(parvec){
    params <- as.list(parvec)
    names(params) <- names(par0)
    args <- c(list(x = data), params, list(log = TRUE))
    ld <- tryCatch(do.call(dfunc, args), error = function(e) rep(-Inf, length(data)))
    if(!is.numeric(ld) || length(ld) != length(data)) return(1e12)
    if(any(is.nan(ld)) || any(is.na(ld))) return(1e12)
    s <- sum(ld)
    if(!is.finite(s)) return(1e12)
    return(-s)
  }
  
  opt <- tryCatch(
    optim(par = par0, fn = neglog, method = "L-BFGS-B", lower = lower_vec, upper = upper_vec, 
          control = c(list(maxit = 1000), control)),
    error = function(e) e
  )
  
  if(inherits(opt, "error")){
    return(list(estimate = NULL, loglik = NA, convergence = 1, message = opt$message))
  } else {
    est <- setNames(opt$par, names(par0))
    return(list(estimate = est, loglik = -opt$value, convergence = opt$convergence, message = opt$message))
  }
}

# ---------- FIT ALL DISTRIBUTIONS ----------
fits <- list()

# Beta 4P
start_beta4 <- list(shape1 = 2, shape2 = 3, a = xmin - 1, b = xmax + 50)
fit_beta4 <- fit_via_optim(dfunc = dbeta4, start = start_beta4, data = data,
                           lower = c(0.1, 0.1, xmin - 100, xmax + 10),
                           upper = c(20, 20, xmin - 0.1, xmax + 500))
fits$Beta4 <- fit_beta4

# Burr 3P
# [FIX B6] start scale was the literal 46.347 (tuned to one dataset); now data-derived
start_burr3 <- list(c = 2, k = 1, scale = mean_data)
fit_burr3 <- fit_via_optim(dfunc = dburr3, start = start_burr3, data = data,
                           lower = c(0.1, 0.1, 1), upper = c(20, 20, 500))
fits$Burr3 <- fit_burr3

# Burr 4P
# [FIX B6] same data-derived start as Burr 3P (was the literal 46.347)
start_burr4 <- list(c = 2, k = 1, scale = mean_data, loc = xmin - 5)
fit_burr4 <- fit_via_optim(dfunc = dburr4, start = start_burr4, data = data,
                           lower = c(0.1, 0.1, 1, xmin - 50), 
                           upper = c(20, 20, 500, xmin - 0.1))
fits$Burr4 <- fit_burr4

# Gamma 2P
fit_gamma2 <- tryCatch(fitdist(data, "gamma", method="mle"), error=function(e) NULL)
fits$Gamma2 <- fit_gamma2

# Gamma 3P
start_gamma3 <- list(shape = 2, scale = mean_data/2, loc = xmin - 5)
fit_gamma3 <- fit_via_optim(dfunc = dgamma3, start = start_gamma3, data = data,
                            lower = c(0.1, 0.1, xmin - 50), 
                            upper = c(20, 100, xmin - 0.1))
fits$Gamma3 <- fit_gamma3

# Log-Gamma
# [FIX B1] fitted on DBH itself (source script fitted on log(data)). Bounds widened
# (alpha <= 1000, scale >= 0.001): on the log scale alpha is large and scale small,
# so the original bounds (0.1-20 for both) forced the estimate onto the boundary.
start_lgamma <- list(alpha = 2, scale = 1)
fit_lgamma <- fit_via_optim(dfunc = dloggamma_custom, start = start_lgamma, data = data,
                            lower = c(0.1, 0.001), upper = c(1000, 20))
fits$LogGamma <- fit_lgamma

# Log-Logistic 2P
start_llogis2 <- list(shape = 2, scale = mean_data)
fit_llogis2 <- fit_via_optim(dfunc = function(x, shape, scale, log=TRUE) 
  dloglogistic3(x, shape, scale, 0, log=log),
  start = start_llogis2, data = data,
  lower = c(0.1, 1), upper = c(20, 500))
fits$LogLogistic2 <- fit_llogis2

# Log-Logistic 3P
start_llogis3 <- list(shape = 2, scale = mean_data, loc = xmin - 5)
fit_llogis3 <- fit_via_optim(dfunc = dloglogistic3, start = start_llogis3, data = data,
                             lower = c(0.1, 1, xmin - 50), 
                             upper = c(20, 500, xmin - 0.1))
fits$LogLogistic3 <- fit_llogis3

# Lognormal 2P
fit_lnorm2 <- tryCatch(fitdist(data, "lnorm", method="mle"), error=function(e) NULL)
fits$Lognormal2 <- fit_lnorm2

# Lognormal 3P
start_lnorm3 <- list(meanlog = mean(log(data)), sdlog = sd(log(data)), loc = xmin - 5)
fit_lnorm3 <- fit_via_optim(dfunc = dlnorm3, start = start_lnorm3, data = data,
                            lower = c(-2, 0.1, xmin - 50), 
                            upper = c(5, 2, xmin - 0.1))
fits$Lognormal3 <- fit_lnorm3

# Weibull 2P
fit_weib2 <- tryCatch(fitdist(data, "weibull", method="mle"), error=function(e) NULL)
fits$Weibull2 <- fit_weib2

# Weibull 3P
start_weib3 <- list(shape = 2, scale = mean_data, loc = xmin - 5)
fit_weib3 <- fit_via_optim(dfunc = dweibull3, start = start_weib3, data = data,
                           lower = c(0.1, 1, xmin - 50), 
                           upper = c(20, 500, xmin - 0.1))
fits$Weibull3 <- fit_weib3

cat("\nDistribution fitting finished:", sum(!sapply(fits, function(f) is.null(f) || is.null(f$estimate))),
    "of 12 models returned parameter estimates.\n")
