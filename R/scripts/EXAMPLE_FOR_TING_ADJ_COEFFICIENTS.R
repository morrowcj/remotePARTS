## Load the libraries
library(tidyverse)
library(remotePARTS)

# set a random seed for repeatability
set.seed(1)

## Load the data
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

# Setup altered coefficient table to fill
coef_table <- expand_grid(F1 = F1_lvls, F2 = F2_lvls) %>%
  arrange(F2, F1)

# randomly sample the pixels
pm = sample_partitions(nrow(ndvi_AK10000), partsize = 2000, npart = NA)

# partition the data set
df_list <- lapply(
  seq_len(ncol(pm)),
  function(i){
    indx = pm[, i]
    df = ndvi_AK10000[indx, ]
    return(df)
  }
)

# fit the spatial covariance matrix for each partition (arbitrary range)
V_list <- lapply(
  df_list,
  function(df){
    covar_exp(distm_scaled(cbind(df$lng, df$lat)), range = 0.01)
  }
)

# fit the GLS for each partition, saving key components
GLS_list <- lapply(
  seq_len(ncol(pm)),
  function(i){
    df = df_list[[i]] # data partition
    V = V_list[[i]] # covariance matrix
    GLS <- fitGLS(
      CLS_coef ~ 0 + F1 + F2, data = df, V = V,
      save.invchol = TRUE, save.xx = TRUE
    )
    return(GLS)
  }
)

# extract coefficient estimates from each GLS
coef_list <- lapply(GLS_list, function(x) x$coefficients)

# extract coefficient covariances from each GLS
vcov_list <- lapply(GLS_list, function(x) x$covar_coef)

# get the degrees of freedom, for each partition (should be identical)
deg_frees <- lapply(GLS_list, function(x) x$df_t)

# calculate new coefficient tables for each partition
new_list <- lapply(
  seq_len(ncol(pm)),
  function(x){
    # get the important components for this partition
    coefs <- coef_list[[x]]
    vcov <- vcov_list[[x]]
    vars <- diag(vcov)
    deg_free <- deg_frees[[x]]


    coef_table <- tibble() # empty table to fill

    for (i in seq_len(nF1)) for (j in seq_len(nF2)) {
      # label the levels that we're comparing
      f1_lvl <- F1_lvls[i]
      f2_lvl <- F2_lvls[j]

      # get the coef values for each factor
      f1_coef <- coefs[i]
      f2_coef <- if_else(j > 1, coefs[nF1 + j - 1], 0)

      # get the variance for each factor
      var_f1 <- vars[i]
      var_f2 <- if_else(j > 1, vars[nF1 + j - 1], 0)

      # get the covariance between the two factors
      cov_f1f2 <- ifelse(j > 1, vcov[i, nF1 + j - 1], 0)

      # get the join coefficient and joint variance (addition)
      coef_f1f2 <- f1_coef + f2_coef
      var_f1f2 <- var_f1 + var_f2 + 2*cov_f1f2
      se_f1f2 <- sqrt(var_f1f2)

      tab = tibble(
        partition = x,
        F1 = f1_lvl, F2 = f2_lvl, adj_coef = coef_f1f2, adj_se = se_f1f2,
        t_stat = adj_coef / adj_se,
        pval = 2 * pt(abs(t_stat), df = deg_free, lower.tail = F)
      )

      coef_table <- bind_rows(coef_table, tab)
    }
    return(coef_table)
  }
)

# ## fit spatial covariance matrix
# V = covar_exp(distm_scaled(cbind(df$lng, df$lat)), range = .01)
#
# ## Fit the GLS
# (GLS <- fitGLS(df$CLS_coef ~ 0 + F1 + F2, data = df, V = V, save.invchol = TRUE, save.xx = TRUE))

# ## variance-covariance matrix
# (vcov <- GLS$covar_coef)
# vars <- diag(vcov)
#
# # get degrees of freedom
# deg_free <- GLS$df_t
#
# ## coefficients
# coefs <- GLS$coefficients

# calculate individual coefficients
# coef_table <- tibble()
#
# for (i in seq_len(nF1)) for (j in seq_len(nF2)) {
#   # label the levels that we're comparing
#   f1_lvl <- F1_lvls[i]
#   f2_lvl <- F2_lvls[j]
#
#   # get the coef values for each factor
#   f1_coef <- coefs[i]
#   f2_coef <- if_else(j > 1, coefs[nF1 + j - 1], 0)
#
#   # get the variance for each factor
#   var_f1 <- vars[i]
#   var_f2 <- if_else(j > 1, vars[nF1 + j - 1], 0)
#
#   # get the covariance between the two factors
#   cov_f1f2 <- ifelse(j > 1, vcov[i, nF1 + j - 1], 0)
#
#   # get the join coefficient and joint variance (addition)
#   coef_f1f2 <- f1_coef + f2_coef
#   var_f1f2 <- var_f1 + var_f2 - 2*cov_f1f2
#   se_f1f2 <- sqrt(var_f1f2)
#
#   tab = tibble(
#     F1 = f1_lvl, F2 = f2_lvl, adj_coef = coef_f1f2, adj_se = se_f1f2,
#     t_stat = adj_coef / adj_se,
#     pval = 2 * pt(abs(t_stat), df = deg_free, lower.tail = F)
#   )
#
#   coef_table <- bind_rows(coef_table, tab)
# }

# Plot the adjusted coefficients
pd = position_dodge(0.5)
new_list %>%
  bind_rows() %>%
  ggplot(aes(x = partition, y = adj_coef, col = F2)) +
  facet_wrap(~F1, labeller = "label_both") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_errorbar(
    aes(ymin = adj_coef - 2*adj_se, ymax = adj_coef + 2*adj_se),
    width = 0, position = pd, size = 1
  ) +
  geom_point(position = pd, shape = 21, fill = "white", size = 2) +
  theme_bw() +
  labs(y = "Adjusted coefficient (± 2 SE)", x = "Partition", col = "F2:")


## ---- calculate partitioned GLS ----

partGLS <- fitGLS_partition(
  formula = CLS_coef ~ 0 + F1 + F2, partmat = pm, data = ndvi_AK10000, save.GLS = TRUE
)

ttest_results <- t.test(partGLS)
coef_tab <- ttest_results$p.t
coefs <- coef_tab[, "Est"]
vcov <- ttest_results$covar_coef
vars <- diag(vcov)

coef_table <- tibble() # empty table to fill

for (i in seq_len(nF1)) for (j in seq_len(nF2)) {
  # label the levels that we're comparing
  f1_lvl <- F1_lvls[i]
  f2_lvl <- F2_lvls[j]

  # get the coef values for each factor
  f1_coef <- coefs[i]
  f2_coef <- if_else(j > 1, coefs[nF1 + j - 1], 0)

  # get the variance for each factor
  var_f1 <- vars[i]
  var_f2 <- if_else(j > 1, vars[nF1 + j - 1], 0)

  # get the covariance between the two factors
  cov_f1f2 <- ifelse(j > 1, vcov[i, nF1 + j - 1], 0)

  # get the join coefficient and joint variance (addition)
  coef_f1f2 <- f1_coef + f2_coef
  var_f1f2 <- var_f1 + var_f2 + 2*cov_f1f2
  se_f1f2 <- sqrt(var_f1f2)

  tab = tibble(
    F1 = f1_lvl, F2 = f2_lvl, adj_coef = coef_f1f2, adj_se = se_f1f2,
    t_stat = adj_coef / adj_se,
    pval = 2 * pt(abs(t_stat), df = deg_free, lower.tail = F)
  )

  coef_table <- bind_rows(coef_table, tab)
}

coef_table

coef_table %>%
  ggplot(aes(x = F1, y = adj_coef, col = F2)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_pointrange(
    aes(ymin = adj_coef - 2*adj_se, ymax = adj_coef + 2*adj_se),
    position = position_dodge(0.5)
  ) +
  theme_bw() +
  labs(
    y = "Adjusted coefficient (± 2 SE)", x = "F1 (landclass)",
    col = "F2 (lat group)"
  )


## ---- With interaction ----

inx_partGLS <- fitGLS_partition(
  formula = CLS_coef ~ 0 + F1 * F2, partmat = pm, data = ndvi_AK10000, save.GLS = TRUE
)

inx_ttest_results <- t.test(inx_partGLS)
inx_coefs <- inx_ttest_results$p.t
inx_vcov <- inx_ttest_results$covar_coef

# CLS_coef ~ 0 + F1 + F2, data = df, V = V,
# save.invchol = TRUE, save.xx = TRUE

## ---- Correlated T-test ----

# # get a table of the new coefficient estimates
# new_coef_tab <- new_list %>% bind_rows() %>%
#   select(partition:adj_coef) %>%
#   pivot_wider(
#     values_from = adj_coef, names_from = partition, names_prefix = "part"
#   ) %>%
#   # calculate the coefficient means across all partitions
#   mutate(avg = rowMeans(across(starts_with("part"))))
#
# # calculate the (rcoefficients)
# # within GLS_partition, this is normally calculated by crosspart_GLS().
#
# # Vcoefij = Wi %*% (QZ::H(xxi) %*% Rij %*% xxj) %*% GZ::H(Wj)
# # rcoefficient = SEi * Vcoefij *SEi
#
# est_rcoefs <- new_coef_tab %>% select(starts_with("part")) %>% t() %>% cor()
#
# # properly format the arguments for part_ttest
# mean_coefs <- new_coef_tab$avg # average coefficient estimates
# vcov_array <- simplify2array(vcov_list) # array of covariances for each partition
#
#
# t_dfs <- remotePARTS:::calc_dfpart(
#   partsize = nrow(pm), p = nrow(new_coef_tab), p0 = 1
# )
#
# df2 <- t_df2
#
# remotePARTS:::part_ttest(
#   coefs = mean_coefs, part.covar_coef = vcov_array, rcoefficients = est_rcoefs,
#   df2 = df2, npart = ncol(pm)
# )
#
# # t.test.partGLS <- function(x, ...){
# #   part_ttest(
# #     coefs = x$overall$coefficients,
# #     part.covar_coef = x$part$covar_coef,
# #     rcoefficients = x$overall$rcoefficients,
# #     df2 = x$overall$dfs[2],
# #     npart = x$overall$partdims["npart"]
# #   )
# # }



