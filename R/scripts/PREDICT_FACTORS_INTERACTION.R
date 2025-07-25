##
## ---- Introduction ----
##

# This script demonstrates how to calculate individual coefficients,
# standard errors, and p-values for each factor level (combination) from a
# partitioned GLS t-test.
#
# It also compares this method against the "predict.lm" function for
# linear regression models, to demonstrate that it is correct.

##
## ---- Setup example data----
##

# Load the libraries
library(tidyverse)
library(remotePARTS)

# set a random seed for replicability
set.seed(1)

# Load the data
data(ndvi_AK10000)

# Define the two factors, from the alaska dataset
ndvi_AK10000 <- ndvi_AK10000 %>%
  mutate(
    # land class, with unused levels dropped
    F1 = factor(land),
    # latitude group
    F2 = findInterval(lat, quantile(lat, na.rm = TRUE), all.inside = TRUE),
    F2 = factor(F2)
  )

# get the levels of these factors
F1_lvls <- with(ndvi_AK10000, factor(levels(F1), levels = levels(F1)))
F2_lvls <- with(ndvi_AK10000, factor(levels(F2), levels = levels(F2)))
nF1 <- length(F1_lvls)
nF2 <- length(F2_lvls)

# randomly sample the pixels
pm = sample_partitions(nrow(ndvi_AK10000), partsize = 2000, npart = NA)

##
## ---- Conduct partitioned GLS ----
##

# fit the GLS, with an interaction between two factors
partGLS <- fitGLS_partition(
  formula = CLS_coef ~ 0 + F1 * F2, partmat = pm, data = ndvi_AK10000,
  save.GLS = TRUE, ncores = 6
)

##
## ---- Correlated t-test ----
##

# Calculate the overall test statistics for each (original) coefficient:
ttest_results <- t.test(partGLS) # t-test results
coefs <- ttest_results$p.t[, "Est"] # just the coefficient estimates
vcov <- ttest_results$covar_coef # just the variance-covariance matrix
vars <- diag(vcov) # just the variances

# extract the degrees of freedom for the t-test
deg_free <- partGLS$overall$dfs[2]

##
## ---- Adjust t-test results (i.e., "predict" at each level) ----
##

# create an empty table to fill
updated_stats <- tibble()

factor_one_name = "F1" # change to the name of the 1st factor (i.e., "pft")
factor_two_name = "F2" # change to the name of the 2nd factor (i.e., "env_zone)

# loop through each factor level to alter them
for (i in seq_len(nF1)) for (j in seq_len(nF2)) {

  ## label the current levels of each factors that we're comparing
  f1_lvl <- F1_lvls[i]
  f2_lvl <- F2_lvls[j]

  ## get the coefficient names
  f1_nm = paste0(factor_one_name, f1_lvl)
  f2_nm = paste0(factor_two_name, f2_lvl)
  f1f2_nm = paste(f1_nm, f2_nm, sep = ":")

  ## base-line coefficient/variance value
  est = coefs[f1_nm]
  v = vars[f1_nm]

  ## add factors together, for higher levels of the second factor
  if (j > 1) {
    est = est + coefs[f2_nm]
    v = v + vars[f2_nm] + 2*vcov[f1_nm, f2_nm]
  }

  ## add interactions, for relevant levels
  if (i > 1 & j > 1) {
    est = est + coefs[f1f2_nm]
    v = v + vars[f1f2_nm] + 2*(vcov[f1_nm, f1f2_nm] + vcov[f2_nm, f1f2_nm])
  }

  ## calculate final stats
  adj_se = sqrt(v)
  adj_tval = est / adj_se
  adj_pval = 2 * pt(abs(adj_tval), df = deg_free, lower.tail = F)

  ## collect the results in a table
  tab = tibble(
    F1 = f1_lvl, F2 = f2_lvl, adj_coef = est, adj_se = adj_se,
    t_stat = adj_tval,
    adj_pval = adj_pval
  )

  ## add the table to the output
  updated_stats <- bind_rows(updated_stats, tab)
}

# calculate 95% critical t-value (for confidence intervals)
t_crit <- qt((1 - 0.95) / 2, df = deg_free, lower.tail = FALSE)

# add confidence interval to table
updated_stats <- updated_stats %>%
  mutate(
    CI_lwr = adj_coef - (t_crit * adj_se),
    CI_upr = adj_coef + (t_crit * adj_se)
  )

# plot the results
updated_stats %>%
  ggplot(aes(x = F1, y = adj_coef, col = F2, fill = adj_pval <= 0.05)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_pointrange(
    aes(ymin = CI_lwr, ymax = CI_upr),
    position = position_dodge(0.5), shape = 21
  ) +
  theme_bw() +
  scale_fill_manual(breaks = c(TRUE, FALSE), values = c("black", "white")) +
  labs(
    x = "F1 (land cover)", col = "F2 (latitude quartile)",
    y = "Adjusted coefficient (± 2 SE)"
  )

##
## ---- Prove that it works with simple linear regression ----
##

# fit the linear model
fm <- lm(CLS_coef ~ 0 + F1*F2, data = ndvi_AK10000)

# extract model stats
lm_coefs <- coefficients(fm)
lm_vcov <- summary(fm)$cov.unscaled
lm_vars <- diag(lm_vcov)

# get the table of levels, to predict with the model
new_dat <- expand_grid(F1 = F1_lvls, F2 = F2_lvls) %>% arrange(F2, F1)

# extract the degrees of freedom for the t-test
lm_deg_free <- fm$df.residual

# get the critical t-value
lm_t_crit <- qt((1 - 0.95) / 2, df = lm_deg_free, lower.tail = FALSE)

# predict with the model to the factor levels:
preds <- predict(
  fm, newdata = new_dat, interval = "confidence", level = 0.95, se.fit = TRUE
)

# build the table of expectations, which we will test against
expected <- cbind(
  new_dat, pred_coef = preds$fit[, "fit"], pred_se = preds$se.fit/preds$residual.scale
)

# perform the method on the linear regression output
lm_updated_stats <- tibble()

for (i in seq_len(nF1)) for (j in seq_len(nF2)) {
  f1_lvl <- F1_lvls[i]
  f2_lvl <- F2_lvls[j]

  f1_nm = paste0(factor_one_name, f1_lvl)
  f2_nm = paste0(factor_two_name, f2_lvl)
  f1f2_nm = paste(f1_nm, f2_nm, sep = ":")

  est = lm_coefs[f1_nm]
  v = lm_vars[f1_nm]

  if (j > 1) {
    est = est + lm_coefs[f2_nm]
    v = v + lm_vars[f2_nm] + 2*lm_vcov[f1_nm, f2_nm]
  }

  if (i > 1 & j > 1) {
    est = est + lm_coefs[f1f2_nm]
    v = v + lm_vars[f1f2_nm] + 2*(lm_vcov[f1_nm, f1f2_nm] + lm_vcov[f2_nm, f1f2_nm])
  }

  adj_se = sqrt(v)
  adj_tval = est / adj_se
  adj_pval = 2 * pt(abs(adj_tval), df = lm_deg_free, lower.tail = F)

  tab = tibble(
    F1 = f1_lvl, F2 = f2_lvl, adj_coef = est, adj_se = adj_se,
    t_stat = adj_tval,
    adj_pval = adj_pval
  )

  lm_updated_stats <- bind_rows(lm_updated_stats, tab)
}

# combine these adjusted stats with the output from predict()
check_tab <- lm_updated_stats %>%
  select(F1:adj_se) %>%
  left_join(expected)

# Are all coefficients the same?
all.equal(check_tab$adj_coef, check_tab$pred_coef, check.attributes = FALSE)

# Are all standard errors the same?
all.equal(check_tab$adj_se, check_tab$pred_se, check.attributes = FALSE)
