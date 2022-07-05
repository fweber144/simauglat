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
aug = run_projpred(refm_fit, seed = seed_vs, nclusters_pred = 50, nterms_max = 5)
lat = run_projpred(refm_fit, seed = seed_vs, latent_proj = TRUE, nclusters_pred = 50, nterms_max = 5)
