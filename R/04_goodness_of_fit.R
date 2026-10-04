# =============================================================================
# [REUSED] R/04_goodness_of_fit.R
# Purpose : Compute KS / AD / CvM statistics, log-likelihood and AIC for each
#           fitted model, rank the models, and export parameter tables.
# Inputs  : objects from R/01 and R/03 (data, fits, out_tab_dir) and the
#           distribution functions defined in R/03
# Outputs : outputs/tables/goodness_of_fit_results.csv
#           outputs/tables/parameters_long_format.csv
#           outputs/tables/parameters_wide_format.csv
#           objects: gof, make_cdf_pdf, parameters_result
# Author  : Oyediran Abiodun
# Date    : 2026-10-04
# Notes   : [adapted from my script: GOF testing, ranking, parameter extraction]
#           Corrections marked [FIX]: removal of dead Johnson SB code (B2),
#           ranking of failed fits (B7). Output paths now point to out_tab_dir.
# =============================================================================

if(!exists("fits") || !exists("data") || !exists("out_tab_dir")) {
  stop("Run R/01_data_cleaning.R and R/03_distribution_fitting.R first.")
}

# ---------- GOF TESTING ----------
gof <- data.frame()

# Improved AD test function
ad_test_custom <- function(data, cdf_fun) {
  n <- length(data)
  if (n < 5) return(NA)  # AD test requires reasonable sample size
  
  # Sort data and compute empirical CDF
  x <- sort(data)
  F_emp <- (1:n) / n
  F_theo <- sapply(x, cdf_fun)
  
  # Avoid boundary issues
  F_theo <- pmin(pmax(F_theo, 1e-10), 1 - 1e-10)
  
  # Anderson-Darling statistic
  ad_stat <- -n - mean((2 * (1:n) - 1) * (log(F_theo) + log(1 - rev(F_theo))))
  
  return(ad_stat)
}

# [FIX B2] the Johnson SB branch was removed: no Johnson SB model was ever fitted
make_cdf_pdf <- function(model_name, fit_obj){
  if(is.null(fit_obj) || is.null(fit_obj$estimate)) return(NULL)
  
  switch(model_name,
         Beta4 = {
           est <- fit_obj$estimate
           pdf <- function(x) dbeta4(x, est["shape1"], est["shape2"], est["a"], est["b"])
           cdf <- function(q) pbeta4(q, est["shape1"], est["shape2"], est["a"], est["b"])
           list(pdf = pdf, cdf = cdf)
         },
         Burr3 = {
           est <- fit_obj$estimate
           pdf <- function(x) dburr3(x, est["c"], est["k"], est["scale"])
           cdf <- function(q) pburr3(q, est["c"], est["k"], est["scale"])
           list(pdf = pdf, cdf = cdf)
         },
         Burr4 = {
           est <- fit_obj$estimate
           pdf <- function(x) dburr4(x, est["c"], est["k"], est["scale"], est["loc"])
           cdf <- function(q) pburr4(q, est["c"], est["k"], est["scale"], est["loc"])
           list(pdf = pdf, cdf = cdf)
         },
         Gamma2 = {
           est <- fit_obj$estimate
           pdf <- function(x) dgamma(x, shape=est["shape"], rate=est["rate"])
           cdf <- function(q) pgamma(q, shape=est["shape"], rate=est["rate"])
           list(pdf = pdf, cdf = cdf)
         },
         Gamma3 = {
           est <- fit_obj$estimate
           pdf <- function(x) dgamma3(x, est["shape"], est["scale"], est["loc"])
           cdf <- function(q) pgamma3(q, est["shape"], est["scale"], est["loc"])
           list(pdf = pdf, cdf = cdf)
         },
         LogGamma = {
           est <- fit_obj$estimate
           pdf <- function(x) dloggamma_custom(x, est["alpha"], est["scale"])
           cdf <- function(q) ploggamma_custom(q, est["alpha"], est["scale"])
           list(pdf = pdf, cdf = cdf)
         },
         LogLogistic2 = {
           est <- fit_obj$estimate
           pdf <- function(x) dloglogistic3(x, est["shape"], est["scale"], 0)
           cdf <- function(q) ploglogistic3(q, est["shape"], est["scale"], 0)
           list(pdf = pdf, cdf = cdf)
         },
         LogLogistic3 = {
           est <- fit_obj$estimate
           pdf <- function(x) dloglogistic3(x, est["shape"], est["scale"], est["loc"])
           cdf <- function(q) ploglogistic3(q, est["shape"], est["scale"], est["loc"])
           list(pdf = pdf, cdf = cdf)
         },
         Lognormal2 = {
           est <- fit_obj$estimate
           pdf <- function(x) dlnorm(x, meanlog=est["meanlog"], sdlog=est["sdlog"])
           cdf <- function(q) plnorm(q, meanlog=est["meanlog"], sdlog=est["sdlog"])
           list(pdf = pdf, cdf = cdf)
         },
         Lognormal3 = {
           est <- fit_obj$estimate
           pdf <- function(x) dlnorm3(x, est["meanlog"], est["sdlog"], est["loc"])
           cdf <- function(q) plnorm3(q, est["meanlog"], est["sdlog"], est["loc"])
           list(pdf = pdf, cdf = cdf)
         },
         Weibull2 = {
           est <- fit_obj$estimate
           pdf <- function(x) dweibull(x, shape=est["shape"], scale=est["scale"])
           cdf <- function(q) pweibull(q, shape=est["shape"], scale=est["scale"])
           list(pdf = pdf, cdf = cdf)
         },
         Weibull3 = {
           est <- fit_obj$estimate
           pdf <- function(x) dweibull3(x, est["shape"], est["scale"], est["loc"])
           cdf <- function(q) pweibull3(q, est["shape"], est["scale"], est["loc"])
           list(pdf = pdf, cdf = cdf)
         },
         NULL)
}


models <- c("Beta4", "Burr3", "Burr4", "Gamma2", "Gamma3",
            "LogGamma", "LogLogistic2", "LogLogistic3", "Lognormal2", 
            "Lognormal3", "Weibull2", "Weibull3")

for(model in models) {
  
  fit_obj <- fits[[model]]
  
  if(is.null(fit_obj) || is.null(fit_obj$estimate)) {
    gof <- rbind(gof, data.frame(Model = model, KS = NA, AD = NA, CvM = NA, 
                                 LogLik = NA, AIC = NA, stringsAsFactors = FALSE))
    next
  }
  
  pfuns <- make_cdf_pdf(model, fit_obj)
  if(is.null(pfuns)) {
    gof <- rbind(gof, data.frame(Model = model, KS = NA, AD = NA, CvM = NA,
                                 LogLik = NA, AIC = NA, stringsAsFactors = FALSE))
    next
  }
  
  cdf_f <- pfuns$cdf
  
  # K-S test with better handling
  # (only the statistic is used; ties in rounded DBH make R warn about the p-value)
  ks_stat <- tryCatch({
    ks_res <- suppressWarnings(ks.test(data, cdf_f))
    as.numeric(ks_res$statistic)
  }, error = function(e) NA)
  
  # Anderson-Darling test
  ad_stat <- tryCatch({
    ad_test_custom(data, cdf_f)
  }, error = function(e) NA)
  
  # Cramer-von Mises
  cvm_stat <- tryCatch({
    as.numeric(goftest::cvm.test(data, null = cdf_f)$statistic)  # unnamed, keeps row names clean
  }, error = function(e) NA)
  
  # Log-Likelihood and AIC
  loglik <- if(!is.null(fit_obj$loglik)) fit_obj$loglik else NA
  k <- length(fit_obj$estimate)
  aic_val <- if(!is.na(loglik)) 2 * k - 2 * loglik else NA
  
  gof <- rbind(gof, data.frame(
    Model = model, 
    KS = ks_stat, 
    AD = ad_stat, 
    CvM = cvm_stat,
    LogLik = loglik,
    AIC = aic_val,
    stringsAsFactors = FALSE
  ))
}

# Sort by AD statistic (best first)
if(!all(is.na(gof$AD))) {
  gof <- gof[order(gof$AD, na.last = TRUE), ]
}

# Display results
cat("Goodness-of-Fit Results:\n")
print(gof)

# Identify best fitting models
if(!all(is.na(gof$AD))) {
  best_ad <- gof$Model[which.min(gof$AD)]
  cat("\nBest model by AD statistic:", best_ad, "\n")
}

if(!all(is.na(gof$AIC))) {
  valid_aic <- gof[!is.na(gof$AIC), ]
  if(nrow(valid_aic) > 0) {
    best_aic <- valid_aic$Model[which.min(valid_aic$AIC)]
    cat("Best model by AIC:", best_aic, "\n")
  }
}

if(!all(is.na(gof$KS))) {
  best_ks <- gof$Model[which.min(gof$KS)]
  cat("Best model by KS statistic:", best_ks, "\n")
}


# ---------- RANKING ----------
relative_rank <- function(vec){
  valid <- !is.na(vec)
  m <- sum(valid)
  ranks <- rep(NA, length(vec))
  if(m <= 1) return(rep(NA, length(vec)))
  Smin <- min(vec[valid]); Smax <- max(vec[valid])
  if(Smax == Smin) { ranks[valid] <- 1; return(ranks) }  # all statistics equal
  ranks[valid] <- 1 + (m - 1) * (vec[valid] - Smin) / (Smax - Smin)
  return(ranks)
}

gof$Rank_KS  <- relative_rank(gof$KS)
gof$Rank_AD  <- relative_rank(gof$AD)
gof$Rank_CvM <- relative_rank(gof$CvM)
# [FIX B7] rowSums(..., na.rm=TRUE) gave a rank sum of 0 (i.e. "best") to models
# with failed fits or missing statistics; such models now get NA and sort last.
rank_cols <- c("Rank_KS","Rank_AD","Rank_CvM")
gof$Rank_Sum <- ifelse(rowSums(is.na(gof[, rank_cols])) > 0, NA, rowSums(gof[, rank_cols]))
gof <- gof[order(gof$Rank_Sum, na.last = TRUE), ]

# Save results
write.csv(gof, file.path(out_tab_dir, "goodness_of_fit_results.csv"), row.names = FALSE)
cat("\nGoodness-of-fit results saved to '", file.path(out_tab_dir, "goodness_of_fit_results.csv"), "'\n", sep = "")



# ---------- EXTRACT PARAMETERS ----------
extract_parameters_table <- function(fits) {
  param_list <- list()
  
  for(model_name in names(fits)) {
    fit_obj <- fits[[model_name]]
    
    if(is.null(fit_obj) || is.null(fit_obj$estimate)) {
      cat("Skipping", model_name, "- NULL fit object\n")
      next
    }
    
    # Handle fitdist objects
    if(inherits(fit_obj, "fitdist")) {
      if(!is.null(fit_obj$estimate) && length(fit_obj$estimate) > 0) {
        params <- fit_obj$estimate
        param_df <- data.frame(
          Model = model_name, 
          Parameter = names(params),
          Value = as.numeric(params),
          stringsAsFactors = FALSE
        )
        param_list[[length(param_list) + 1]] <- param_df
        cat("✓ Added", model_name, "(fitdist object)\n")
      }
    } 
    # Handle custom fit objects
    else if(is.list(fit_obj) && !is.null(fit_obj$estimate) && length(fit_obj$estimate) > 0) {
      params <- fit_obj$estimate
      
      param_df <- data.frame(
        Model = model_name, 
        Parameter = names(params),
        Value = as.numeric(params),
        stringsAsFactors = FALSE
      )
      param_list[[length(param_list) + 1]] <- param_df
      cat("✓ Added", model_name, "(custom fit object)\n")
    }
  }
  
  if(length(param_list) > 0) {
    params_table <- do.call(rbind, param_list)
    
    # Convert to wide format
    params_wide <- tryCatch({
      wide_data <- reshape(params_table, 
                           idvar = "Model", 
                           timevar = "Parameter", 
                           direction = "wide")
      names(wide_data) <- gsub("Value.", "", names(wide_data))
      wide_data
    }, error = function(e) {
      cat("Error in wide format conversion:", e$message, "\n")
      return(NULL)
    })
    
    return(list(long = params_table, wide = params_wide))
  } else {
    return(NULL)
  }
}

# ---------- EXTRACT AND DISPLAY PARAMETERS ----------
cat("\n", strrep("=", 60), "\n")
cat("EXTRACTING PARAMETERS\n")
cat(strrep("=", 60), "\n")

parameters_result <- extract_parameters_table(fits)

if(!is.null(parameters_result)) {
  # Display results
  cat("\n", strrep("=", 60), "\n")
  cat("PARAMETERS TABLE (LONG FORMAT):\n")
  cat(strrep("=", 60), "\n")
  print(parameters_result$long)
  
  if(!is.null(parameters_result$wide)) {
    cat("\n", strrep("=", 60), "\n")
    cat("PARAMETERS TABLE (WIDE FORMAT):\n")
    cat(strrep("=", 60), "\n")
    print(parameters_result$wide)
    
    # Save to CSV files
    write.csv(parameters_result$long, file.path(out_tab_dir, "parameters_long_format.csv"), row.names = FALSE)
    write.csv(parameters_result$wide, file.path(out_tab_dir, "parameters_wide_format.csv"), row.names = FALSE)
    
    cat("\n✓ Parameters saved to:\n")
    cat("  - '", file.path(out_tab_dir, "parameters_long_format.csv"), "'\n", sep = "")
    cat("  - '", file.path(out_tab_dir, "parameters_wide_format.csv"), "'\n", sep = "")
  }
} else {
  cat("✗ No parameters extracted\n")
}
