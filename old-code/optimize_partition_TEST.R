library(remotePARTS)
set.seed(456)
run.sims = FALSE # change this to run the loop

df <- read.csv("data-raw/AK_ndvi_common-land.csv")
n.pix = nrow(df)
n.part = 15 # 15
part.size = 2000 # 2000
parts = sample_partitions(npix = n.pix, partsize = part.size, npart = n.part)

## parameters
spcor = NA
# spcor = c(.5)
nugget = NA
# nugget = .2
pars.start = c(r = 0.1, a = 1,  nug = 1e-9)
V.meth = "exponential-power"

# iterate through all partitions:
if (run.sims){
  pars.iter <-  data.frame(NULL)
  # progress bar
  pb = txtProgressBar(min = 0, max = ncol(parts), style = 3)
  for(p in 1:ncol(parts)){
    setTxtProgressBar(pb = pb, value = p)

    # subset data
    df.tmp = df[parts[, p], ]
    loc.tmp = df.tmp[, c("lng", "lat")]

    ## optimize_GLS method ----
    # estimate params
    opt = optimize_GLS(formula = cls.coef ~ 0 + land, data = df.tmp,
                       pars.start = pars.start,
                       coords = loc.tmp,
                       V.meth = V.meth,
                       verbose = F, ret.GLS = TRUE)
    # collect some data
    tmp.opt = data.frame(r = opt$spatial.pars["r"],
                         a = opt$spatial.pars["a"],
                         nug = opt$spatial.pars["nug"],
                         B1 = opt$GLS$betahat[1],
                         B2 = opt$GLS$betahat[2],
                         B3 = opt$GLS$betahat[3],
                         p1 = opt$GLS$pval.t[1],
                         p2 = opt$GLS$pval.t[2],
                         p3 = opt$GLS$pval.t[3],
                         partition = p,
                         type = "optim")

    ## old method - equivalent to fit_spatialcor() ----
    # fit time-series
    X.tmp = as.matrix(df.tmp[, -c(1:7)])
    ar = fitAR.map(X = X.tmp, t = 1:ncol(X.tmp), ret_resid = TRUE)
    # get residuals and reformat
    resids = ar$resids
    cor.resid = cor(t(resids))
    vec_cor = cor.resid[upper.tri(cor.resid, diag = TRUE)]
    # distance matrix and reformat
    D.tmp <- geosphere::distm(loc.tmp)
    vec_dist = scales::rescale((D.tmp[upper.tri(D.tmp, diag = TRUE)]),
                               to = c(0,1))
    # combine into dataframe
    W = data.frame(dist = vec_dist, cor = vec_cor)
    # fot the residuals and estimate the parameters
    r.fit = try(nls(cor ~ exp(-(dist/r)^a), data = W,
                    start = list(r = .1, a = 1)))
    sp.est = try(coef(r.fit))
    # rescale to the data
    sp.est["r.scaled"] = try(sp.est["r"] * max(D.tmp))
    # fit the variance matrix
    V.tmp = try(fitV(D.tmp, spatialcor = c(r = sp.est["r.scaled"], a = sp.est["a"]), "exponential-power"))
    # fit the GLS
    GLS = GLS_worker(y = df.tmp$cls.coef,
                     X = model.matrix(~ 0 + land, data = df.tmp),
                     V = V.tmp, X0 = model.matrix(~1, data = df.tmp))
    # collect data
    tmp.res = data.frame(r = sp.est["r"],
                         a = sp.est["a"],
                         nug = GLS$nugget,
                         B1 = GLS$betahat[1],
                         B2 = GLS$betahat[2],
                         B3 = GLS$betahat[3],
                         p1 = GLS$pval.t[1],
                         p2 = GLS$pval.t[2],
                         p3 = GLS$pval.t[3],
                         partition = p,
                         type = "resid")

    ## Stack ----
    pars.iter = rbind(pars.iter, tmp.opt, tmp.res)
  }

  close(pb)

  # # fix stupid r screw up
  # pars.iter = pars.iter %>% mutate(r.scaled = r)
  # pars.iter.bckp = pars.iter
  #
  # for(i in 1:ncol(parts)){
  #   print(i)
  #   df.tmp = df[parts[, p], ]
  #   loc.tmp = df.tmp[, c("lng", "lat")]
  #   D.tmp = geosphere::distm(loc.tmp)
  #   max.D = max(D.tmp)
  #
  #   pars.iter[pars.iter$partition == i & pars.iter$type == "resid", "r"] = pars.iter[pars.iter$partition == i & pars.iter$type == "resid", "r.scaled"]/max.D
  # }

  save(pars.iter, file = "old-code/spatial-parameter-estimation-comparison_2Kpx.RData")
  usethis::use_data(pars.iter, internal = TRUE)
} else {
  load("old-code/spatial-parameter-estimation-comparison_2Kpx.RData")
}

# Plots ----
library(ggplot2);library(dplyr)

# histogram of spatial parameters
pars.iter %>%
  reshape2::melt(id.vars = c("partition", "type")) %>%
  filter(variable %in% c("r", "a", "nug")) %>%
  ggplot(aes(x = value, fill = type)) +
  geom_histogram(position = "dodge") +
  facet_wrap(~ variable, scales = "free", ncol = 1)

# line plot of spatial parameters
pars.iter %>% reshape2::melt(id.vars = c("partition", "type")) %>%
  filter(variable %in% c("r", "a", "nug")) %>%
  ggplot(aes(x = partition, y = value, group = type, col = type, shape = type)) +
  geom_point() +
  geom_line(aes(lty = type)) +
  theme(legend.position = "bottom") +
  labs(lty = NULL, col = NULL, shape = NULL) +
  facet_wrap(~ variable, scales = "free", ncol = 1)

# line plot of coefficeint estimates
pars.iter %>% reshape2::melt(id.vars = c("partition", "type")) %>%
  filter(variable %in% c("B1", "B2", "B3")) %>%
  ggplot(aes(x = partition, y = value, group = type, shape = type, col = type)) +
  geom_point() +
  geom_line(aes(lty = type)) +
  theme(legend.position = "bottom") +
  labs(lty = NULL, shape = NULL, col = NULL) +
  facet_wrap(~ variable, scales = "free", ncol = 1)


# line plot of p-values
pars.iter %>% reshape2::melt(id.vars = c("partition", "type")) %>%
  filter(variable %in% c("p1", "p2", "p3")) %>%
  ggplot(aes(x = partition, y = -log10(value), group = type, color = type, shape = type)) +
  geom_hline(yintercept = -log10(.05), col = "blue", alpha = .2) +
  geom_point() +
  geom_line(aes(lty = type)) +
  theme(legend.position = "bottom") +
  labs(lty = NULL, shape = NULL, y = "-log10(p)", color = NULL) +
  facet_wrap(~ variable, scales = "free", ncol = 1)

# get the semivariograms from these data ----

spcor.func = function(d, r, a, nug){
  return((1 - nug)*exp(-(d/r)^a))
}

semivar.table = NULL
for(i in 1:nrow(pars.iter)){
  pars = as.data.frame(pars.iter[i, ])
  d = seq(0, 1, length.out = 1000)
  r = pars$r
  a = pars$a
  nug = pars$nug
  var = spcor.func(d = d, r = r, a = a, nug = nug)
  semivar = 1-var
  tmp = data.frame(d = d, var = var, semivar = semivar,
                   r = r, a = a, nug = nug,
                   partition = pars$partition, type = pars$type
  )
  semivar.table <- rbind(semivar.table, tmp)
}

ggplot(data = semivar.table, aes(x = d, y = semivar, group = partition)) +
  geom_path(aes(col = r)) +
  geom_point(aes(x = 0, y = nug, col = r), size = 2, shape = "_") +
  facet_rep_wrap(~ type) +
  labs(x = "distance", y = "semivariance") +
  scale_color_viridis_c(option = "plasma", begin = 0, end = 0.75)

