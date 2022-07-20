# brms::mcmc_plot(bfit, variable = "^b_X", regex = TRUE, type = "areas")

simresrng <- readRDS(file = "simresrng.rds")
SIM_ITER <- 4L
.Random.seed <- simresrng[[SIM_ITER]]

sim_dat_etc <- dataconstructor()
refm_fit <- fit_ref(dat = sim_dat_etc$dat, fml = sim_dat_etc$fml)
seed_vs <- sample.int(.Machine$integer.max, 1)

# debug(projpred:::.get_p_clust)
# debug(brms:::link_cumulative)
# debug(brms:::link)
projpred_aug <- run_projpred(refm_fit, dat_indep = sim_dat_etc$dat_indep,
                             seed = seed_vs, nclusters_pred = 50,
                             nterms_max = 5)
dat_indep_lat <- sim_dat_etc$dat_indep
dat_indep_lat$projpredY <- colMeans(
  rstantools::posterior_linpred(refm_fit, newdata = dat_indep_lat)
)
projpred_lat <- run_projpred(refm_fit, dat_indep = dat_indep_lat,
                             seed = seed_vs, latent_proj = TRUE,
                             nclusters_pred = 50, nterms_max = 5)
