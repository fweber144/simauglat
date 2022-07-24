#_______________________________________________________________________________
# Simulation study for comparing augmented-data and latent projection in
# projpred
#_______________________________________________________________________________

# Setup -------------------------------------------------------------------

## Installation (only required once) --------------------------------------

### With more checks:
# devtools::install_github("fweber144/brms", "projpred_latent")
# devtools::install_github("fweber144/projpred", "augdat_latent")
###
### With less checks:
# remotes::install_github("fweber144/brms", "projpred_latent")
# remotes::install_github("fweber144/projpred", "augdat_latent")
###

## Timestamp --------------------------------------------------------------

cat("\n-----\n")
cat("Timestamp at the beginning of the script:\n")
print(Sys.time())
cat("-----\n")

## Global options ---------------------------------------------------------

# For source()-ing this script:
warn_orig_glob <- options(warn = 1)

# options(projpred.warn_augdat_experimental = FALSE)

## User options -----------------------------------------------------------

only_init_fit <- F

if (!only_init_fit) {
  nsim <- 50 # The number of simulation iterations
  # par_type <- "doSeq"
  # par_type <- "doParallel"
  par_type <- "doMPI"
  if (par_type %in% c("doParallel")) {
    # The number of CPU cores (more generally, "workers") for the simulation:
    ncores <- 8 # parallel::detectCores(logical = FALSE)
    if (nsim < ncores) {
      warning("Increasing `nsim` to `ncores = ", ncores, "`.")
      nsim <- ncores
    }
  }
}

nobsv <- 100L
nobsv_indep <- nobsv
ncat <- 5L
yunq <- paste0("ycat", seq_len(ncat))
nthres <- ncat - 1L
npreds_tot <- 50L
p0 <- as.integer(ceiling(npreds_tot * 0.10))
source("rh_sigma_tilde.R")
sigti <- calc_sigma_tilde(ncats = ncat)
seed_glob <- 856824715

cat("-----\np0:\n")
print(p0)
cat("-----\n")

cat("-----\nsigma_tilde:\n")
print(sigti)
cat("-----\n")

## Prepare simulation -----------------------------------------------------

if (!only_init_fit) {
  library(foreach)
  # library(iterators)
  library(doRNG)
  if (par_type %in% c("doParallel")) {
    library(doParallel)
    registerDoParallel(ncores)
  } else if (par_type %in% c("doMPI")) {
    library(doMPI)
    cl_obj <- startMPIcluster()
    registerDoMPI(cl_obj)
  } else if (par_type %in% c("doSeq")) {
    registerDoSEQ()
  } else {
    stop("Unknown `par_type`.")
  }
  cat("-----\nWorker information:\n")
  print(getDoParName())
  print(getDoParWorkers())
  cat("-----\n")
}

set.seed(seed_glob)

# Functions (part 1) ------------------------------------------------------

# For generating regression coefficients according to the regularized horseshoe
# (RH) prior (note: scale_global = tau_0):
rhorseshoe <- function(
    npreds_RH,
    df_lambdas = 1,
    df_global = 1, scale_global = par_ratio * sigma_tilde / sqrt(N),
    par_ratio, sigma_tilde, N,
    df_slab = 4, scale_slab = 2
) {
  csq <- 1 / rgamma(1, shape = df_slab / 2, rate = scale_slab^2 * df_slab / 2)
  tau <- abs(scale_global * rt(1, df = df_global) + 0)
  lambdas <- abs(1 * rt(npreds_RH, df = df_lambdas) + 0)
  lambdassq <- lambdas^2
  lambdas_tilde <- sqrt(csq * lambdassq / (csq + tau^2 * lambdassq))
  betas <- rnorm(npreds_RH, mean = 0, sd = tau * lambdas_tilde)
  return(betas)
}

# For simulating data in a single simulation iteration:
dataconstructor <- function() {

  ## Definitions ------------------------------------------------------------

  nobsv_sim <- nobsv + nobsv_indep

  # The intercepts at centered predictors ("Intercept"s in brms parlance, not
  # "b_Intercept"s); note that when switching the option, the prior in the
  # data-fitting model usually should be adjusted:
  ### Option 1:
  # thres <- seq(-1.5, 1.5, length.out = nthres)
  ###
  ### Option 2:
  # thres <- sort(rnorm(nthres))
  ###
  ### Option 3:
  # thres <- sort(2.5 * rt(nthres, df = 3) + 0)
  ###
  ### Option 4:
  thres <- qnorm(seq_len(nthres) / ncat)
  ###
  # print(diff(c(0, pnorm(thres), 1)))

  npreds_cont <- npreds_tot
  coefs_cont <- rhorseshoe(npreds_cont,
                           par_ratio = p0 / (npreds_cont - p0),
                           sigma_tilde = sigti,
                           N = nobsv)

  npreds_grCP <- 0L # 1L
  ngrCP <- integer() # c(3L)
  stopifnot(identical(length(ngrCP), npreds_grCP))
  if (npreds_grCP > 0) {
    coefs_grCP <- list(
      c(0, seq(-2, 2, length.out = ngrCP[1] - 1L))
    )
    stopifnot(identical(length(coefs_grCP), npreds_grCP))
  }

  npreds_grPP <- 0L # 1L
  if (!exists("nobsv_per_grPP")) {
    # Number of observations per group (only roughly; some groups might get a
    # few observations more):
    nobsv_per_grPP <- integer() # c(2L)
  }
  ngrPP <- c(nobsv %/% nobsv_per_grPP)
  stopifnot(identical(length(ngrPP), npreds_grPP))
  if (npreds_grPP > 0) {
    coefs_grPP <- list(
      list("icpt" = rnorm(ngrPP[1], sd = 1.5)) # , "Xcont1" = <...>
    )
    stopifnot(identical(length(coefs_grPP), npreds_grPP))
  }

  ## Continuous predictors --------------------------------------------------

  dat_sim <- setNames(replicate(npreds_cont, {
    rnorm(nobsv_sim)
  }, simplify = FALSE), paste0("Xcont", seq_len(npreds_cont)))
  dat_sim <- as.data.frame(dat_sim)
  # Start constructing the linear predictor:
  eta <- drop(as.matrix(dat_sim) %*% coefs_cont)

  ## Completely pooled (CP) categorical predictors --------------------------

  if (npreds_grCP == 1) {
    dat_sim$XgrCP1 <- sample(
      gl(n = ngrCP[1], k = floor(nobsv_sim / ngrCP[1]), length = nobsv_sim,
         labels = paste0("gr", seq_len(ngrCP[1])))
    )
    # Continue constructing the linear predictor:
    eta <- eta + coefs_grCP[[1]][dat_sim$XgrCP1]
  } else if (npreds_grCP > 1) {
    stop("This value of `npreds_grCP` is currently not supported.")
  }

  ## Partially pooled (PP) categorical predictors ---------------------------

  if (npreds_grPP == 1) {
    dat_sim$XgrPP1 <- sample(
      gl(n = ngrPP[1], k = floor(nobsv_sim / ngrPP[1]), length = nobsv_sim,
         labels = paste0("gr", seq_len(ngrPP[1])))
    )
    # Continue constructing the linear predictor:
    eta <- eta + coefs_grPP[[1]]$icpt[dat_sim$XgrPP1]
    # + coefs_grPP[[1]]$Xcont1[dat_sim$XgrPP1] * dat_sim$Xcont1
  } else if (npreds_grPP > 1) {
    stop("This value of `npreds_grPP` is currently not supported.")
  }

  ## Construct "epred" ------------------------------------------------------

  thres_eta <- sapply(thres, function(thres_k) {
    thres_k - eta
  })
  # Shouldn't be necessary (just to be safe): Emulate a single posterior draw:
  dim(thres_eta) <- c(1L, nobsv_sim, nthres)
  yprobs <- brms:::inv_link_cumulative(thres_eta, link = "logit")
  # Because of emulating a single posterior draw above:
  dim(yprobs) <- c(nobsv_sim, ncat)

  ## Response ---------------------------------------------------------------

  dat_sim$Y <- factor(sapply(seq_len(nobsv_sim), function(i) {
    sample(yunq, size = 1, prob = yprobs[i, ])
  }), ordered = TRUE)

  ## Formula ----------------------------------------------------------------

  voutc <- "Y"
  vpreds <- grep("^Xcont|^XgrCP", names(dat_sim), value = TRUE)
  vpreds_PP <- grep("^XgrPP", names(dat_sim), value = TRUE)
  if (length(vpreds_PP) > 0L) {
    if (!all(lengths(coefs_grPP) == 1)) {
      stop("Group-level slopes are currently not supported.")
    }
    vpreds_PP <- paste0("(1 | ", vpreds_PP, ")")
  }

  fml_sim <- as.formula(paste(
    voutc, "~", paste(c(vpreds, vpreds_PP), collapse = " + ")
  ))

  ## Output -----------------------------------------------------------------

  return(list(true_coefs_cont = coefs_cont,
              true_PPEs = if (npreds_grPP > 0) coefs_grPP[[1]]$icpt else NULL,
              dat = dat_sim[1:nobsv, , drop = FALSE],
              fml = fml_sim,
              dat_indep = dat_sim[(nobsv + 1):nobsv_sim, , drop = FALSE]))
}

# Initial reference model fit ---------------------------------------------

if (only_init_fit) {
  # Fit the brms reference model once, so it only needs to be updated (and not
  # recompiled) during the simulation. MCMC convergence diagnostics are only
  # checked here, not during the simulation.
  sim_dat_etc <- dataconstructor()
  # Check that all response categories are present in the initial model fit:
  stopifnot(identical(levels(sim_dat_etc$dat$Y), yunq))
  options(mc.cores = parallel::detectCores(logical = FALSE))
  # options(cmdstanr_write_stan_file_dir = getwd())
  ### Needed because the computing cluster complains about `p0` not found:
  par_ratio_sigti <- p0 / (npreds_tot - p0) * sigti
  cat("-----\npar_ratio_sigti:\n")
  print(par_ratio_sigti, digits = 9)
  cat("-----\n")
  stopifnot(all.equal(par_ratio_sigti, 0.117657042, tolerance = 1e-9))
  ###
  bfit <- brms::brm(
    formula = sim_dat_etc$fml,
    data = sim_dat_etc$dat,
    family = brms::cumulative(link = "probit"),
    prior = brms::prior(horseshoe(par_ratio = 0.117657042)) +
      brms::prior(normal(0, 2.5), class = "Intercept"),
    ### For backend = "rstan":
    control = list(adapt_delta = 0.99), # , max_treedepth = 15L
    init_r = 1,
    ###
    ### For backend = "cmdstanr":
    # backend = "cmdstanr",
    # adapt_delta = 0.99,
    # # max_treedepth = 15L,
    # init = 1,
    ###
    silent = 0,
    file = "bfit",
    file_refit = "on_change",
    seed = sample.int(.Machine$integer.max, 1)
  )
  # Check MCMC diagnostics:
  source("check_MCMC_diagn.R")
  # debugonce(check_MCMC_diagn)
  MCMC_diagn <- check_MCMC_diagn(
    C_stanfit = bfit$fit,
    # exclude_NAs = TRUE,
    pars = "disc",
    include = FALSE
  )
  cat("\n-----\n")
  cat("Are all MCMC diagnostics OK?:\n")
  print(MCMC_diagn$all_OK)
  cat("-----\n")

  ## Teardown

  # Reset global options:
  options(warn = warn_orig_glob$warn)

  # Timestamp:
  cat("\n-----\n")
  cat("Timestamp at the end of the script:\n")
  print(Sys.time())
  cat("-----\n")

  # Exit code:
  cat("\nExit code: 0\n")
  if (interactive() && !isatty(stdout())) {
    beepr::beep(2)
  }
  stop("Exiting cleanly (only throwing an error to stop source()-ing the ",
       "script).")
} else {
  bfit <- readRDS("bfit.rds")
}

# Functions (part 2) ------------------------------------------------------

# For fitting (updating) the reference model:
fit_ref <- function(dat, fml) {
  return(update(
    bfit,
    formula. = fml,
    newdata = dat,
    cores = 1,
    silent = 2,
    refresh = 0,
    ### For backend = "rstan":
    init_r = 1,
    ###
    ### For backend = "cmdstanr":
    # adapt_delta = 0.99,
    # # max_treedepth = 15L,
    # init = 1,
    ###
    seed = sample.int(.Machine$integer.max, 1)
  ))
}

# For running projpred's variable selection:
run_projpred <- function(refm_fit, dat_indep, ...) {
  d_indep <- list(
    data = dat_indep,
    offset = rep(0, nrow(dat_indep)),
    weights = rep(1, nrow(dat_indep)),
    y = if (!is.null(dat_indep$projpredY)) dat_indep$projpredY else dat_indep$Y
  )
  time_bef <- Sys.time()
  vs <- projpred::varsel(refm_fit, d_test = d_indep, method = "forward", ...)
  time_aft <- Sys.time()
  # Currently, we need to use the internal projpred function .tabulate_stats()
  # to obtain the reference model's performance:
  stats_man <- projpred:::.tabulate_stats(vs, stats = "mlpd")
  return(list(
    time_vs = as.numeric(time_aft - time_bef, units = "mins"),
    refstat = stats_man$value[stats_man$size == Inf],
    # plot_obj = plot(vs, deltas = TRUE, stats = "mlpd"),
    sgg_size = projpred::suggest_size(vs, stat = "mlpd"),
    smmry = summary(vs, deltas = TRUE, stats = "mlpd",
                    type = c("mean", "se", "lower", "upper"))$selection,
    soltrms = projpred::solution_terms(vs)
  ))
}

sim_runner <- function(...) {
  foreach(
    sit = seq_len(nsim), # Needs package "iterators": icount(nsim),
    .inorder = FALSE,
    # .packages = c("brms", "projpred"), # , "rstanarm"
    .export = c("rhorseshoe", "dataconstructor", "fit_ref", "run_projpred",
                "nobsv", "nobsv_indep", "ncat", "yunq", "nthres", "npreds_tot",
                "p0", "sigti", "bfit"),
    # .noexport = c("<object_name>"),
    .options.snow = list(attachExportEnv = TRUE)
  ) %dorng% {
    cat("\nSimulation iteration: ", sit, "\n", sep = "")
    Rseed <- .Random.seed
    sim_dat_etc <- dataconstructor()
    refm_fit <- fit_ref(dat = sim_dat_etc$dat, fml = sim_dat_etc$fml)
    seed_vs <- sample.int(.Machine$integer.max, 1)
    projpred_aug <- try(run_projpred(refm_fit,
                                     dat_indep = sim_dat_etc$dat_indep,
                                     seed = seed_vs,
                                     ...),
                        silent = TRUE)
    dat_indep_lat <- sim_dat_etc$dat_indep
    dat_indep_lat$projpredY <- colMeans(
      rstantools::posterior_linpred(refm_fit, newdata = dat_indep_lat)
    )
    projpred_lat <- try(run_projpred(refm_fit,
                                     dat_indep = dat_indep_lat,
                                     seed = seed_vs,
                                     latent_proj = TRUE,
                                     ...),
                        silent = TRUE)
    if (inherits(projpred_aug, "try-error")) {
      dot_args <- list(...)
      save(sit, refm_fit, seed_vs, dot_args, Rseed, file = "failed.rda")
      stop("The augmented-data projpred run failed in simulation iteration ",
           sit, ". Error message: \"", attr(projpred_aug, "condition")$message,
           "\". Objects for replicating this failure were saved to ",
           "\"failed.rda\". Use `loaded_objs <- load(\"failed.rda\")` to ",
           "restore it (including `.Random.seed`).")
    }
    if (inherits(projpred_lat, "try-error")) {
      dot_args <- list(...)
      save(sit, refm_fit, seed_vs, dot_args, Rseed, file = "failed.rda")
      stop("The latent projpred run failed in simulation iteration ",
           sit, ". Error message: \"", attr(projpred_lat, "condition")$message,
           "\". Objects for replicating this failure were saved to ",
           "\"failed.rda\". Use `loaded_objs <- load(\"failed.rda\")` to ",
           "restore it (including `.Random.seed`).")
    }
    return(list(
      aug = projpred_aug,
      lat = projpred_lat,
      true_coefs_cont = sim_dat_etc$true_coefs_cont,
      true_PPEs = sim_dat_etc$true_PPEs
    ))
  }
}

# Run simulation ----------------------------------------------------------

print(system.time({
  simres <- sim_runner(
    ### For a faster (experimental) run:
    nclusters_pred = 50,
    nterms_max = 5
    ###
  )
}))
saveRDS(simres, file = "simres.rds") # simres <- readRDS(file = "simres.rds")

# Post-processing ---------------------------------------------------------

if (!dir.exists("figs")) dir.create("figs")

## True completely pooled coefficients ------------------------------------
## (i.e., the draws from the regularized horseshoe distribution)

true_coefs_cont <- unlist(lapply(simres, "[[", "true_coefs_cont"))
true_coefs_cont <- data.frame(coef = true_coefs_cont)
cat("\n-----\n")
cat("Proportion of (completely pooled) coefficient draws with absolute value >",
    "0.5 (across all simulation iterations and all `npreds_tot` coefficient",
    "draws within each simulation iteration):\n")
print(proportions(table(abs(true_coefs_cont$coef) > 0.5, useNA = "ifany")))
cat("-----\n")
print(ggplot2::ggplot(data = true_coefs_cont,
                      mapping = ggplot2::aes(x = coef)) +
        ggplot2::geom_histogram(bins = 40)) # + ggplot2::geom_density()
ggplot2::ggsave(file.path("figs", "true_coefs_cont.pdf"),
                width = 7, height = 7 * 0.618)

## Timing -----------------------------------------------------------------

mins_vs <- do.call(rbind, lapply(seq_along(simres), function(sim_idx) {
  return(data.frame(sim_idx = sim_idx,
                    t_aug = simres[[sim_idx]]$aug$time_vs,
                    t_lat = simres[[sim_idx]]$lat$time_vs))
}))
mins_vs <- reshape(
  mins_vs,
  direction = "long",
  v.names = "minutes",
  varying = list("minutes" = grep("^t_", names(mins_vs), value = TRUE)),
  timevar = "prj_meth",
  times = c("aug", "lat"),
  idvar = "sim_idx_ch",
  sep = "_"
)
stopifnot(identical(mins_vs$sim_idx_ch, mins_vs$sim_idx))
mins_vs$sim_idx_ch <- NULL
print(ggplot2::ggplot(data = mins_vs,
                      mapping = ggplot2::aes(x = prj_meth, y = minutes)) +
        ggplot2::geom_boxplot()) # + ggplot2::geom_violin()
ggplot2::ggsave(file.path("figs", "time.pdf"),
                width = 7, height = 7 * 0.618)

## True partially pooled effects ------------------------------------------

if (!all(sapply(lapply(simres, "[[", "true_PPEs"), is.null))) {
  true_PPEs <- do.call(cbind, lapply(simres, "[[", "true_PPEs"))
  dimnames(true_PPEs) <- list(
    grp = paste0("grp", seq_len(nrow(true_PPEs))),
    simiter = paste0("simiter", seq_len(ncol(true_PPEs)))
  )
  stopifnot(!all(lengths(apply(true_PPEs, 1, unique, simplify = FALSE)) == 1))
  cat("\n-----\n")
  cat("Quartiles of the true partially pooled effects (across",
      "all simulation iterations):\n")
  print(quantile(true_PPEs))
  cat("-----\n")
  if (nsim <= 10) {
    cat("\n-----\n")
    cat("Empirical SDs of the true partially pooled effects (for",
        "all simulation iterations):\n")
    print(apply(true_PPEs, 2, sd))
    cat("-----\n")
  } else {
    cat("\n-----\n")
    cat("Quartiles of the empirical SDs of the true partially pooled effects",
        "(across all simulation iterations):\n")
    print(quantile(apply(true_PPEs, 2, sd)))
    cat("-----\n")
  }
}

## Solution paths ---------------------------------------------------------

same_solpths <- sapply(simres, function(simres_i) {
  identical(simres_i$aug$soltrms,
            simres_i$lat$soltrms)
})
cat("\n-----\n")
cat("Differing solution paths:\n")
for (sim_idx in seq_along(simres)[!same_solpths]) {
  cat("\n---\n")
  cat("Simulation iteration: ", sim_idx, "\n", sep = "")
  cat("Solution path of the augmented-data variable selection:\n",
      paste(simres[[sim_idx]]$aug$soltrms, collapse = ", "),
      "\n", sep = "")
  cat("Solution path of the latent variable selection:\n",
      paste(simres[[sim_idx]]$lat$soltrms, collapse = ", "),
      "\n", sep = "")
  cat("---\n")
}
cat("-----\n")

## Model size selection plots ---------------------------------------------

### Currently can't be used since returning the output of plot.vsel() leads to
### memory issues when running the simulation in parallel:
# plotter_single <- function(sim_idx, prj_meth) {
#   if (prj_meth == "aug") {
#     title_raw <- "Augmented-data"
#   } else if (prj_meth == "lat") {
#     title_raw <- "Latent"
#   }
#   title_pretty <- paste0(title_raw, "; simulation iteration ", sim_idx)
#   stopifnot(inherits(simres[[sim_idx]][[prj_meth]]$plot_obj, "ggplot"))
#   print(simres[[sim_idx]][[prj_meth]]$plot_obj +
#           ggplot2::labs(title = title_pretty))
#   return(invisible(TRUE))
# }
# plotter_sep <- function(sim_idx) {
#   aug_succ <- plotter_single(sim_idx = sim_idx, prj_meth = "aug")
#   lat_succ <- plotter_single(sim_idx = sim_idx, prj_meth = "lat")
#   return(invisible(aug_succ && lat_succ))
# }
# sep_succs <- lapply(seq_along(simres), plotter_sep)
# stopifnot(all(unlist(sep_succs)))
###
plotter_ovrlay <- function(prj_meth) {
  if (prj_meth == "aug") {
    title_raw <- "Augmented-data"
  } else if (prj_meth == "lat") {
    title_raw <- "Latent"
  }
  y_chr <- setdiff(names(simres[[1L]][[prj_meth]]$smmry),
                   c("solution_terms", "se", "lower", "upper", "size"))
  plotdat_comm <- do.call(rbind, lapply(seq_along(simres), function(sim_idx) {
    cbind(sim_idx = sim_idx,
          simres[[sim_idx]][[prj_meth]]$smmry[, c("size", y_chr)])
  }))
  print(ggplot2::ggplot(data = plotdat_comm,
                        mapping = ggplot2::aes_string(x = "size",
                                                      y = y_chr,
                                                      group = "sim_idx",
                                                      alpha = I(0.4))) +
          ggplot2::geom_hline(yintercept = 0,
                              color = "firebrick",
                              linetype = "dashed") +
          ### Only required when using the extended suggest_size() heuristics
          ### from branch `elpd4` of repo `fweber144/projpred`:
          # ggplot2::geom_hline(yintercept = -4 / nobsv_indep,
          #                     color = "dodgerblue",
          #                     linetype = "dotdash") +
          ###
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::labs(title = title_raw))
  ggplot2::ggsave(file.path("figs", paste0(y_chr, "_", prj_meth, ".pdf")),
                  width = 7, height = 7 * 0.618)
  return(invisible(TRUE))
}
plotter_com <- function() {
  aug_succ <- plotter_ovrlay(prj_meth = "aug")
  lat_succ <- plotter_ovrlay(prj_meth = "lat")
  return(invisible(aug_succ && lat_succ))
}
com_succs <- plotter_com()
stopifnot(com_succs)

## Suggested sizes --------------------------------------------------------

sgger_size <- function(sim_idx) {
  return(c(
    sgg_size_aug = simres[[sim_idx]]$aug$sgg_size,
    sgg_size_lat = simres[[sim_idx]]$lat$sgg_size
  ))
}
sgg_sizes <- sapply(seq_along(simres), sgger_size)
if (nsim <= 10) {
  print(sgg_sizes)
}
if (anyNA(sgg_sizes)) {
  warning("Found suggested sizes which are `NA`.")
}
sgg_sizes_lat_minus_aug <- apply(sgg_sizes, 2, function(x) {
  x["sgg_size_lat"] - x["sgg_size_aug"]
})
cat("\n-----\n")
cat("Differences of the suggested sizes (latent minus augmented-data):\n")
sgg_sizes_tab <- table(sgg_sizes_lat_minus_aug, useNA = "always")
print(sgg_sizes_tab)
print(proportions(sgg_sizes_tab))
cat("-----\n")
xlab_long <- "Difference of the suggested sizes (latent minus augmented-data)"
print(ggplot2::qplot(factor(sgg_sizes_lat_minus_aug),
                     geom = "bar",
                     xlab = xlab_long) +
        ggplot2::scale_y_continuous(breaks = scales::breaks_pretty()))
ggplot2::ggsave(file.path("figs", "sgg_sizes_diff.pdf"),
                width = 7, height = 7 * 0.618)

# doRNG -------------------------------------------------------------------

# Information from package "doRNG" about the seeds used:
simresrng <- attr(simres, "rng")
saveRDS(simresrng, file = "simresrng.rds")
stopifnot(identical(class(simresrng), "list"))
stopifnot(length(simresrng) == nsim)
# These are random seeds of length 7 each (see `?RNGkind`, sections "Details"
# (subpart "L'Ecuyer-CMRG") and "Value"):
stopifnot(identical(lengths(simresrng), rep(7L, nsim)))
simresrngv <- attr(simres, "doRNG_version")
saveRDS(simresrngv, file = "simresrngv.rds")

# Teardown ----------------------------------------------------------------

## Session info -----------------------------------------------------------

# Needed to list these packages and their dependencies in the session info:
suppressPackageStartupMessages({
  library(brms)
  library(rstan)
  library(StanHeaders)
  library(cmdstanr)
  library(rstanarm)
  library(projpred)
})
if (!isNamespaceLoaded("rmarkdown")) loadNamespace("rmarkdown")
if (!isNamespaceLoaded("yaml")) loadNamespace("yaml")

# sessioninfo::session_info(to_file = TRUE)
sink(file = paste0("session_info_", par_type, ".txt"))
print(sessioninfo::session_info())
cat("\n-----\n")
cat("CmdStan version:\n")
print(cmdstanr::cmdstan_version(error_on_NA = FALSE))
cat("-----\n")
sink()

## Free resources ---------------------------------------------------------

if (par_type %in% c("doParallel")) {
  stopImplicitCluster()
} else if (par_type %in% c("doMPI")) {
  ### Ideally, the following code should be used (at least according to the
  ### "doMPI" vignette). However, `closeCluster(cl_obj)` as well as `mpi.quit()`
  ### (as well as simply `q("no")`) seem to cause the execution to get stuck (so
  ### that the R session is not quit):
  # closeCluster(cl_obj)
  # mpi.quit()
  ###
  mpi.finalize()
  stop("Quitting cleanly (this is just a dummy error since the R session ",
       "somehow can't be quitted using `q(\"no\")` here).")
}

## Reset global options ---------------------------------------------------

options(warn = warn_orig_glob$warn)

## Timestamp --------------------------------------------------------------

cat("\n-----\n")
cat("Timestamp at the end of the script:\n")
print(Sys.time())
cat("-----\n")

## Exit -------------------------------------------------------------------

cat("\nExit code: 0\n")
if (interactive() && !isatty(stdout())) {
  beepr::beep(2)
}
