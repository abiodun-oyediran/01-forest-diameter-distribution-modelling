# Methodology

This document describes what the code in `R/` does. It does not describe any
analysis beyond the scripts in this repository.

## 1. Pipeline

```
data/raw/dbh_sample.csv
        │
R/01_data_cleaning.R        read CSV, drop invalid DBH values, summary values
        │
R/02_diameter_classes.R     diameter-class frequency table (10 cm classes)
        │
R/03_distribution_fitting.R maximum-likelihood fit of 12 candidate distributions
        │
R/04_goodness_of_fit.R      KS / AD / CvM / AIC, relative ranking, parameter tables
        │
R/05_visualization.R        PNG diagnostics and final summary
```

`run_all.R` sources the five scripts in order in one R session. Scripts share
objects through the global environment (`data`, `fits`, `gof`, ...), so each
script checks that its prerequisites exist and stops with a message otherwise.

## 2. Data and cleaning

The input is a single column of tree diameters at breast height (DBH, cm).
`R/01_data_cleaning.R` coerces the column to numeric and keeps finite values
greater than zero; the number of records removed is printed. It stops if fewer
than five valid values remain. `xmin`, `xmax`, the mean and the standard
deviation are computed after cleaning and are used for starting values and
parameter bounds.

## 3. Diameter classes

`R/02_diameter_classes.R` assigns each tree to a class of width
`DBH_CLASS_WIDTH` cm (default 10, classes are `[lower, upper)`) and writes the
count and relative frequency per class to `outputs/tables/diameter_class_table.csv`.
This step is descriptive. The distribution fitting in step 4 uses the individual
DBH values, not the class counts.

## 4. Candidate distributions

| Model | Parameters (k) | Definition |
|---|---|---|
| Beta4 | shape1, shape2, a, b (4) | Beta on [a, b] |
| Burr3 | c, k, scale (3) | Burr XII, location 0 |
| Burr4 | c, k, scale, loc (4) | Burr XII with location |
| Gamma2 | shape, rate (2) | Gamma |
| Gamma3 | shape, scale, loc (3) | Gamma shifted by loc |
| LogGamma | alpha, scale (2) | ln(DBH) ~ Gamma(alpha, scale), defined for DBH > 1 |
| LogLogistic2 | shape, scale (2) | Log-logistic, location 0 |
| LogLogistic3 | shape, scale, loc (3) | Log-logistic with location |
| Lognormal2 | meanlog, sdlog (2) | Lognormal |
| Lognormal3 | meanlog, sdlog, loc (3) | Lognormal shifted by loc |
| Weibull2 | shape, scale (2) | Weibull |
| Weibull3 | shape, scale, loc (3) | Weibull shifted by loc |

## 5. Parameter estimation

* **2-parameter Gamma, Lognormal and Weibull** are fitted with
  `fitdistrplus::fitdist(..., method = "mle")`.
* **All other models** are fitted with `fit_via_optim()`, which minimises the
  negative log-likelihood with `optim(method = "L-BFGS-B", maxit = 1000)` under
  box constraints. Invalid or non-finite likelihood values return a large
  penalty (1e12). A fit that errors returns no estimate and is reported as `NA`
  in the goodness-of-fit table.

Bounds used by `fit_via_optim()` (`xmin` = smallest DBH, `xmax` = largest DBH):

| Model | Bounds |
|---|---|
| Beta4 | shape1, shape2 in [0.1, 20]; a in [xmin − 100, xmin − 0.1]; b in [xmax + 10, xmax + 500] |
| Burr3 | c, k in [0.1, 20]; scale in [1, 500] |
| Burr4 | c, k in [0.1, 20]; scale in [1, 500]; loc in [xmin − 50, xmin − 0.1] |
| Gamma3 | shape in [0.1, 20]; scale in [0.1, 100]; loc in [xmin − 50, xmin − 0.1] |
| LogGamma | alpha in [0.1, 1000]; scale in [0.001, 20] |
| LogLogistic2 | shape in [0.1, 20]; scale in [1, 500] |
| LogLogistic3 | shape in [0.1, 20]; scale in [1, 500]; loc in [xmin − 50, xmin − 0.1] |
| Lognormal3 | meanlog in [−2, 5]; sdlog in [0.1, 2]; loc in [xmin − 50, xmin − 0.1] |
| Weibull3 | shape in [0.1, 20]; scale in [1, 500]; loc in [xmin − 50, xmin − 0.1] |

Starting values are in `R/03_distribution_fitting.R`. Scale starting values for
Burr, Log-logistic and Weibull models use the sample mean (half the mean for
Gamma3); location starting values are `xmin − 5` (Beta4 starts at a = `xmin − 1`,
b = `xmax + 50`).

## 6. Goodness of fit

For each model with parameter estimates, `R/04_goodness_of_fit.R` computes:

* **KS**: the Kolmogorov–Smirnov statistic `D` from `ks.test(data, cdf)`. Only the
  statistic is used.
* **AD**: the Anderson–Darling statistic
  `A² = −n − (1/n) Σ (2i − 1) [ln F(x₍ᵢ₎) + ln(1 − F(x₍ₙ₊₁₋ᵢ₎))]`, with `F`
  clipped to [1e-10, 1 − 1e-10].
* **CvM**: the Cramér–von Mises statistic from `goftest::cvm.test`.
* **LogLik** and **AIC** = `2k − 2·logLik`, with `k` the number of estimated parameters.

## 7. Ranking

Each of KS, AD and CvM is rescaled across the valid models to a relative rank
between 1 (smallest statistic) and `m` (largest), where `m` is the number of
valid models:

`rank = 1 + (m − 1) · (S − S_min) / (S_max − S_min)`

The **Rank_Sum** is the sum of the three relative ranks; lower is better.
Models missing any of the three statistics receive `NA` and are sorted last and
excluded from the plots. Rankings are by Rank_Sum; AD, AIC and KS winners are
also printed separately.

## 8. Diagnostics produced

| Figure | Content |
|---|---|
| `histogram_all_models.png` | Histogram with all fitted densities |
| `histogram_top3_models.png` | Histogram with the three best-ranked densities |
| `qq_plot_best_model.png` | Q-Q plot for the best-ranked model (quantiles by numerical inversion of its CDF) |
| `cdf_comparison_top3.png` | Empirical CDF against the CDFs of the top three models |
| `combined_diagnostic_plots.png` | 2 × 2 panel: histogram, Q-Q, CDF, rank-sum bar chart |
| `histogram_weibull2_gamma2.png` | Histogram with the 2-parameter Weibull and Gamma fits |

## 9. Implementation notes and corrections

This repository was built from an exploratory single-script analysis. The
following changes were made so the code is correct and reproducible; everything
else keeps the original logic, function names and variable names.

| ID | Issue in the exploratory script | Change |
|---|---|---|
| B1 | Log-Gamma was fitted to `log(data)` but its CDF was evaluated on raw DBH, so fit, GOF statistics and plots were inconsistent; its log-likelihood was also not comparable with other models | Defined on DBH itself (`ln X ~ Gamma`), fitted to `data`; bounds widened to alpha ≤ 1000, scale ≥ 0.001 |
| B2 | Johnson SB branches referenced `djsb` (not defined) and a model that was never fitted | Dead Johnson SB code and the unused `pjsb` removed |
| B3 | Final summary listed files that were never written | Summary lists the files actually produced |
| B4 | Combined plot panel B depended on `cdf_fun`, overwritten by an earlier loop | Explicit `qq_ready` flag; blank panel keeps the layout if no Q-Q data |
| B5 | Weibull/Gamma histogram snippet used undefined objects (`d`, `dW2`, `dGamma`, `fitw2p`, `fitgammap`) | Rewritten to use `data` and the existing Weibull2 / Gamma2 fits |
| B6 | Burr starting scale hard-coded as 46.347 | Starting scale = sample mean |
| B7 | `rowSums(..., na.rm = TRUE)` gave failed fits a rank sum of 0 (best) | Rank sum is `NA` unless all three ranks exist; such models sort last |
| B8 | `cvm.test` statistic carried a name that leaked into data-frame row names | Converted with `as.numeric` |
| B9 | All-models plot had a legend that did not name the models | Legend names every model; colour and line type are unique per model |
| — | `file.choose()` and working-directory outputs | Input path and output folders via `DBH_CSV`, `OUT_FIG_DIR`, `OUT_TAB_DIR` |
| — | Eight packages loaded, six unused | Only `fitdistrplus` and `goftest` are loaded |
| — | KS ties warning | `suppressWarnings` around `ks.test` (p-value is not used) |

## 10. Assumptions and caveats

* The goodness-of-fit statistics are computed on the same data used for
  estimation and no p-values are reported; the statistics are for ranking
  candidates, not for formal hypothesis tests.
* Rank sums depend on the set of candidate models: adding or removing a model
  changes every relative rank.
* Several location parameters (and Beta4 `a`, `shape2`) can sit on their
  bounds, especially when the data are truncated at an inventory threshold.
