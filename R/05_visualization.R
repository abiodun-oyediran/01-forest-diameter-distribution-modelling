# =============================================================================
# [REUSED] R/05_visualization.R
# Purpose : Draw the diagnostic PNG plots for the fitted distributions, the
#           2P Weibull / 2P Gamma histogram, and print the final summary.
# Inputs  : objects from R/01 - R/04 (fits, gof, data, make_cdf_pdf, out_fig_dir,
#           out_tab_dir)
# Outputs : outputs/figures/histogram_all_models.png
#           outputs/figures/histogram_top3_models.png
#           outputs/figures/qq_plot_best_model.png
#           outputs/figures/cdf_comparison_top3.png
#           outputs/figures/combined_diagnostic_plots.png
#           outputs/figures/histogram_weibull2_gamma2.png
# Author  : Oyediran Abiodun
# Date    : 2026-10-04
# Notes   : [adapted from my script: create_essential_png_plots and the final
#           summary; adapted from my graph snippet: Weibull/Gamma histogram].
#           Corrections marked [FIX]: dead Johnson SB branches removed (B2),
#           Q-Q panel state (B4), snippet objects (B5), summary text (B3),
#           failed fits excluded from rankings (B7), model legend (B9),
#           clipped Gamma curve (B10).
# =============================================================================

if(!exists("gof") || !exists("fits") || !exists("make_cdf_pdf") || !exists("out_fig_dir")) {
  stop("Run R/01 - R/04 first (objects 'gof', 'fits', 'make_cdf_pdf', 'out_fig_dir' not found).")
}

# ---------- ESSENTIAL PLOTS IN PNG FORMAT ----------
create_essential_png_plots <- function(fits, gof, data, top_n = 3, out_dir = out_fig_dir) {
  
  # [FIX B7] only models with a complete set of ranks can be ranked
  ranked_gof <- gof[!is.na(gof$Rank_Sum), ]
  if(nrow(ranked_gof) == 0) stop("No successfully fitted models are available to plot.")
  top_n <- min(top_n, nrow(ranked_gof))
  
  # Get top models
  top_models <- head(ranked_gof$Model, top_n)
  top_gof <- ranked_gof[seq_len(top_n), ]
  
  # Standard colors for better visibility
  colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")  # Standard matplotlib colors
  
  # 1. HISTOGRAM WITH FITTED DISTRIBUTIONS (ALL MODELS)
  png(file.path(out_dir, "histogram_all_models.png"), width = 12, height = 8, units = "in", res = 300)
  
  par(mar = c(5, 5, 4, 2) + 0.1, family = "sans")
  
  # Calculate histogram data
  hist_data <- hist(data, breaks = 25, plot = FALSE)
  x_range <- seq(0, max(data) * 1.1, length.out = 1000)
  
  # Create main plot
  plot(0, 0, type = "n", 
       xlim = c(0, max(data) * 1.1), 
       ylim = c(0, max(hist_data$density) * 1.3),
       xlab = "Diameter at Breast Height (cm)", 
       ylab = "Probability Density",
       main = "Diameter Distribution: Histogram with All Fitted Models",
       cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.4, 
       font.lab = 2, font.main = 2)
  
  # Add histogram
  rect(hist_data$breaks[-length(hist_data$breaks)], 0, 
       hist_data$breaks[-1], hist_data$density,
       col = "lightgray", border = "darkgray", lwd = 0.5)
  
  # Add ALL fitted distributions with different line types and colors
  all_models <- gof$Model[!is.na(gof$KS)]  # Only models with successful fits
  # [FIX B9] colours cycle through the palette, line type changes after each cycle,
  # so every model has a unique style that the legend can name
  model_cols <- colors[(seq_along(all_models) - 1) %% length(colors) + 1]
  model_lty  <- (seq_along(all_models) - 1) %/% length(colors) + 1
  
  for(i in seq_along(all_models)) {
    model <- all_models[i]
    
    # Get PDF function
    pdf_fun <- NULL
    pfuns <- make_cdf_pdf(model, fits[[model]])
    if(!is.null(pfuns)) pdf_fun <- pfuns$pdf
    
    if(!is.null(pdf_fun)) {
      tryCatch({
        y_vals <- sapply(x_range, pdf_fun)
        lines(x_range, y_vals, col = model_cols[i], 
              lwd = 2, lty = model_lty[i])
      }, error = function(e) NULL)
    }
  }
  
  # Add legend
  legend("topright", 
         legend = c("Empirical Data", all_models),
         fill = c("lightgray", rep(NA, length(all_models))),
         border = c("darkgray", rep(NA, length(all_models))),
         lty = c(NA, model_lty),
         col = c(NA, model_cols),
         lwd = c(NA, rep(2, length(all_models))),
         bg = "white", cex = 0.9, box.lwd = 1)
  
  grid()
  dev.off()
  cat("✓ Histogram with all models saved as:", file.path(out_dir, "histogram_all_models.png"), "\n")
  
  # 2. HISTOGRAM WITH TOP 3 MODELS ONLY
  png(file.path(out_dir, "histogram_top3_models.png"), width = 12, height = 8, units = "in", res = 300)
  
  par(mar = c(5, 5, 4, 2) + 0.1, family = "sans")
  
  # Create main plot for top 3 models
  plot(0, 0, type = "n", 
       xlim = c(0, max(data) * 1.1), 
       ylim = c(0, max(hist_data$density) * 1.3),
       xlab = "Diameter at Breast Height (cm)", 
       ylab = "Probability Density",
       main = "Top 3 Best-Fitting Distributions",
       cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.4, 
       font.lab = 2, font.main = 2)
  
  # Add histogram
  rect(hist_data$breaks[-length(hist_data$breaks)], 0, 
       hist_data$breaks[-1], hist_data$density,
       col = "lightgray", border = "darkgray", lwd = 0.5)
  
  # Add top 3 fitted distributions
  for(i in seq_along(top_models)) {
    model <- top_models[i]
    
    # Get PDF function
    pdf_fun <- NULL
    pfuns <- make_cdf_pdf(model, fits[[model]])
    if(!is.null(pfuns)) pdf_fun <- pfuns$pdf
    
    if(!is.null(pdf_fun)) {
      tryCatch({
        y_vals <- sapply(x_range, pdf_fun)
        lines(x_range, y_vals, col = colors[i], lwd = 3, lty = i)
      }, error = function(e) NULL)
    }
  }
  
  # Add detailed legend for top 3 models
  legend("topright", 
         legend = c("Empirical Data", 
                    paste0(top_models, " (Rank ", seq_along(top_models), ", KS = ", round(top_gof$KS, 4), ")")),
         fill = c("lightgray", rep(NA, length(top_models))),
         border = c("darkgray", rep(NA, length(top_models))),
         lty = c(NA, seq_along(top_models)),
         col = c(NA, colors[seq_along(top_models)]),
         lwd = c(NA, rep(3, length(top_models))),
         bg = "white", cex = 1.1, box.lwd = 1)
  
  grid()
  dev.off()
  cat("✓ Histogram with top 3 models saved as:", file.path(out_dir, "histogram_top3_models.png"), "\n")
  
  # 3. Q-Q PLOT (Best model only) - Using blue instead of red
  png(file.path(out_dir, "qq_plot_best_model.png"), width = 8, height = 8, units = "in", res = 300)
  
  par(mar = c(5, 5, 4, 2) + 0.1, family = "sans")
  
  best_model <- top_models[1]
  qq_ready <- FALSE   # [FIX B4] explicit flag used by panel B of the combined plot
  
  # Get CDF function for best model
  cdf_fun <- NULL
  pfuns <- make_cdf_pdf(best_model, fits[[best_model]])
  if(!is.null(pfuns)) cdf_fun <- pfuns$cdf
  
  if(!is.null(cdf_fun)) {
    # Calculate theoretical quantiles using quantile function approximation
    theoretical_quantiles <- sapply(ppoints(length(data)), function(p) {
      tryCatch({
        uniroot(function(x) cdf_fun(x) - p, 
                interval = c(0, max(data) * 3), 
                extendInt = "yes", tol = 1e-6)$root
      }, error = function(e) NA)
    })
    
    # Remove NA values
    valid_idx <- !is.na(theoretical_quantiles)
    theoretical_quantiles <- theoretical_quantiles[valid_idx]
    sample_quantiles <- sort(data)[valid_idx]
    
    if(length(theoretical_quantiles) > 0) {
      qq_ready <- TRUE
      # Create Q-Q plot with blue color scheme
      plot(theoretical_quantiles, sample_quantiles,
           main = paste("Q-Q Plot:", best_model, "(Best Fitting Model)"),
           xlab = "Theoretical Quantiles (cm)",
           ylab = "Sample Quantiles (cm)",
           pch = 21, bg = "#1f77b4", col = "darkblue", cex = 1.2,
           cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.4,
           font.lab = 2, font.main = 2)
      
      # Add reference line in blue
      abline(0, 1, col = "#1f77b4", lwd = 3, lty = 2)
      
      # Add confidence envelope
      n <- length(sample_quantiles)
      if(n > 10) {
        # Simple confidence bands based on order statistics
        k <- qnorm(ppoints(n))
        se <- (1 / dnorm(k)) * sqrt(ppoints(n) * (1 - ppoints(n)) / n)
        upper <- theoretical_quantiles + 1.96 * se * sd(theoretical_quantiles)
        lower <- theoretical_quantiles - 1.96 * se * sd(theoretical_quantiles)
        
        lines(sort(theoretical_quantiles), upper[order(theoretical_quantiles)], 
              col = "gray50", lwd = 1, lty = 3)
        lines(sort(theoretical_quantiles), lower[order(theoretical_quantiles)], 
              col = "gray50", lwd = 1, lty = 3)
      }
      
      legend("topleft", 
             legend = c("Q-Q Points", "Line of Perfect Fit", "95% Confidence Envelope"),
             pch = c(21, NA, NA),
             pt.bg = c("#1f77b4", NA, NA),
             col = c("darkblue", "#1f77b4", "gray50"),
             lty = c(NA, 2, 3),
             lwd = c(NA, 3, 1),
             bg = "white", cex = 1.1, box.lwd = 1)
      
      grid()
    }
  }
  
  dev.off()
  cat("✓ Q-Q plot saved as:", file.path(out_dir, "qq_plot_best_model.png"), "\n")
  
  # 4. CDF COMPARISON PLOT (Top 3 models)
  png(file.path(out_dir, "cdf_comparison_top3.png"), width = 10, height = 6, units = "in", res = 300)
  
  par(mar = c(5, 5, 4, 2) + 0.1, family = "sans")
  
  # Create empirical CDF
  ecdf_vals <- ecdf(data)
  x_vals <- seq(0, max(data) * 1.1, length.out = 500)
  emp_cdf <- ecdf_vals(x_vals)
  
  # Create main plot
  plot(x_vals, emp_cdf, type = "l", lwd = 3, col = "black",
       xlim = c(0, max(data) * 1.1), ylim = c(0, 1),
       xlab = "Diameter at Breast Height (cm)",
       ylab = "Cumulative Probability",
       main = "Cumulative Distribution Function: Top 3 Models",
       cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.4,
       font.lab = 2, font.main = 2)
  
  # Add theoretical CDFs for top 3 models
  for(i in seq_along(top_models)) {
    model <- top_models[i]
    
    # Get CDF function
    cdf_fun <- NULL
    pfuns <- make_cdf_pdf(model, fits[[model]])
    if(!is.null(pfuns)) cdf_fun <- pfuns$cdf
    
    if(!is.null(cdf_fun)) {
      tryCatch({
        theo_cdf <- sapply(x_vals, cdf_fun)
        lines(x_vals, theo_cdf, col = colors[i], lwd = 2.5, lty = i)
      }, error = function(e) NULL)
    }
  }
  
  # Add legend
  legend("bottomright", 
         legend = c("Empirical CDF", 
                    paste0(top_models, " (Rank ", seq_along(top_models), ")")),
         col = c("black", colors[seq_along(top_models)]),
         lty = c(1, seq_along(top_models)),
         lwd = c(3, rep(2.5, length(top_models))),
         bg = "white", cex = 1.1, box.lwd = 1)
  
  grid()
  dev.off()
  cat("✓ CDF comparison saved as:", file.path(out_dir, "cdf_comparison_top3.png"), "\n")
  
  # 5. COMBINED PLOT WITH ALL VISUALIZATIONS
  png(file.path(out_dir, "combined_diagnostic_plots.png"), width = 14, height = 10, units = "in", res = 300)
  
  # Set up 2x2 layout
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  par(mar = c(4.5, 4.5, 3, 1.5), family = "sans", cex = 0.9)
  
  # Subplot 1: Histogram with top 3 models
  hist_data <- hist(data, breaks = 25, plot = FALSE)
  plot(0, 0, type = "n", 
       xlim = c(0, max(data) * 1.1), 
       ylim = c(0, max(hist_data$density) * 1.3),
       xlab = "DBH (cm)", ylab = "Density",
       main = "A) Histogram with Top 3 Fitted Distributions",
       cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3, font.main = 2)
  
  rect(hist_data$breaks[-length(hist_data$breaks)], 0, 
       hist_data$breaks[-1], hist_data$density,
       col = "lightgray", border = "darkgray")
  
  for(i in seq_along(top_models)) {
    model <- top_models[i]
    
    pdf_fun <- NULL
    pfuns <- make_cdf_pdf(model, fits[[model]])
    if(!is.null(pfuns)) pdf_fun <- pfuns$pdf
    
    if(!is.null(pdf_fun)) {
      tryCatch({
        y_vals <- sapply(x_range, pdf_fun)
        lines(x_range, y_vals, col = colors[i], lwd = 2, lty = i)
      }, error = function(e) NULL)
    }
  }
  
  legend("topright", legend = c("Data", top_models),
         fill = c("lightgray", rep(NA, length(top_models))),
         border = c("darkgray", rep(NA, length(top_models))),
         lty = c(NA, seq_along(top_models)),
         col = c(NA, colors[seq_along(top_models)]), lwd = 2, cex = 0.9, bg = "white")
  
  # Subplot 2: Q-Q Plot (using blue colors)
  if(qq_ready) {
    plot(theoretical_quantiles, sample_quantiles,
         main = "B) Q-Q Plot (Best Model)",
         xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
         pch = 19, col = "#1f77b4", cex = 0.8,
         cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3, font.main = 2)
    abline(0, 1, col = "#1f77b4", lwd = 2, lty = 2)
    grid()
  } else {
    plot.new()   # keep the 2x2 layout when no Q-Q data are available
  }
  
  # Subplot 3: CDF Comparison
  plot(x_vals, emp_cdf, type = "l", lwd = 2, col = "black",
       xlab = "DBH (cm)", ylab = "Cumulative Probability",
       main = "C) CDF Comparison - Top 3 Models",
       cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3, font.main = 2)
  
  for(i in seq_along(top_models)) {
    model <- top_models[i]
    
    cdf_fun <- NULL
    pfuns <- make_cdf_pdf(model, fits[[model]])
    if(!is.null(pfuns)) cdf_fun <- pfuns$cdf
    
    if(!is.null(cdf_fun)) {
      tryCatch({
        theo_cdf <- sapply(x_vals, cdf_fun)
        lines(x_vals, theo_cdf, col = colors[i], lwd = 2, lty = i)
      }, error = function(e) NULL)
    }
  }
  
  legend("bottomright", legend = c("Empirical", top_models),
         col = c("black", colors[seq_along(top_models)]), lty = 1, lwd = 2, cex = 0.9, bg = "white")
  
  # Subplot 4: GOF Statistics for top models
  top_5 <- head(ranked_gof, 5)
  par(mar = c(7, 4.5, 3, 1.5))
  barpos <- barplot(top_5$Rank_Sum, names.arg = top_5$Model, las = 2,
                    main = "D) Goodness-of-Fit Ranking (Top 5)",
                    ylab = "Rank Sum (Lower = Better)",
                    col = colors[seq_len(nrow(top_5))], border = "black",
                    cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3, font.main = 2)
  text(barpos, top_5$Rank_Sum, round(top_5$Rank_Sum, 2), 
       pos = 3, cex = 0.8, font = 2)
  grid(NA, NULL)
  
  dev.off()
  cat("✓ Combined diagnostic plot saved as:", file.path(out_dir, "combined_diagnostic_plots.png"), "\n")
}

# ---------- CREATE THE ESSENTIAL PNG PLOTS ----------
cat("\n", strrep("=", 70), "\n")
cat("CREATING ESSENTIAL PNG PLOTS FOR RESEARCH\n")
cat(strrep("=", 70), "\n")

# close any half-written device if plotting fails
tryCatch(create_essential_png_plots(fits, gof, data, top_n = 3),
         error = function(e) {
           while(dev.cur() > 1) dev.off()
           stop("Plotting failed: ", conditionMessage(e))
         })

cat("\n", strrep("=", 70), "\n")
cat("PLOT FILES CREATED:\n")
cat(strrep("=", 70), "\n")
cat("1.", file.path(out_fig_dir, "histogram_all_models.png"), "- All fitted distributions\n")
cat("2.", file.path(out_fig_dir, "histogram_top3_models.png"), "- Top 3 distributions only\n")
cat("3.", file.path(out_fig_dir, "qq_plot_best_model.png"), "- Q-Q plot validation (blue color scheme)\n") 
cat("4.", file.path(out_fig_dir, "cdf_comparison_top3.png"), "- Cumulative distribution comparison\n")
cat("5.", file.path(out_fig_dir, "combined_diagnostic_plots.png"), "- All diagnostic plots combined\n")
cat("\nAll plots are high-resolution (300 DPI) PNG files suitable for publications.\n")

# ---------- 2P WEIBULL vs 2P GAMMA HISTOGRAM ----------
# [adapted from my graph snippet] [FIX B5] the snippet used objects that do not
# exist in this pipeline (d, dW2, dGamma, fitw2p, fitgammap). They are replaced
# by `data` and the Weibull2 / Gamma2 fits already estimated in R/03.
if(!is.null(fits$Weibull2) && !is.null(fits$Gamma2)) {
  png(file.path(out_fig_dir, "histogram_weibull2_gamma2.png"), width = 8, height = 6, units = "in", res = 300)
  x <- seq(min(data), max(data), 1)
  y_weibull <- dweibull(x, shape = fits$Weibull2$estimate["shape"], scale = fits$Weibull2$estimate["scale"])
  y_gamma   <- dgamma(x, shape = fits$Gamma2$estimate["shape"], rate = fits$Gamma2$estimate["rate"])
  # [FIX B10] y-axis sized from the histogram and both curves; before, the Gamma
  # peak could be clipped at the top of the plot
  hist(data, prob = TRUE, xlab = "Dbh class (cm)", ylab = "Relative frequency of trees", main = "",
       ylim = c(0, 1.05 * max(hist(data, plot = FALSE)$density, y_weibull, y_gamma)))
  lines(x, y_weibull, lty = 1, lwd = 2)
  lines(x, y_gamma, lty = 2, lwd = 2)
  legend("topright", c("2P Weibull", "Gamma 2p"), lty = c(1, 2), lwd = c(2, 2))
  dev.off()
  cat("✓ Weibull / Gamma histogram saved as:", file.path(out_fig_dir, "histogram_weibull2_gamma2.png"), "\n")
} else {
  cat("Weibull 2P / Gamma 2P histogram skipped: one of the two fits is unavailable.\n")
}

# ---------- FINAL SUMMARY ----------
# [FIX B3] the source summary named files that the script never created
ranked_gof_final <- gof[!is.na(gof$Rank_Sum), ]
cat("\n=== FINAL SUMMARY ===\n")
if(nrow(ranked_gof_final) > 0) {
  cat("Best fitting distribution:", ranked_gof_final$Model[1], "\n")
  cat("KS statistic:", round(ranked_gof_final$KS[1], 4), "\n")
  cat("Rank Sum:", round(ranked_gof_final$Rank_Sum[1], 2), "\n\n")
} else {
  cat("No model produced a complete set of goodness-of-fit statistics.\n\n")
}

cat("Output files created:\n")
cat("1.", file.path(out_tab_dir, "diameter_class_table.csv"), "- diameter class frequencies\n")
cat("2.", file.path(out_tab_dir, "goodness_of_fit_results.csv"), "- GOF statistics table\n")
cat("3.", file.path(out_tab_dir, "parameters_long_format.csv"), "and parameters_wide_format.csv - parameter estimates\n")
cat("4. PNG figures in", out_fig_dir, "\n")

cat("\nAnalysis complete!\n")
