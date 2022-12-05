#_______________________________________________________________________________
# simauglat: Simulation study comparing augmented-data and latent projection in
#            projpred
# Copyright (C) 2022  Frank Weber
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
# options(projpred.warn_latent_experimental = FALSE)

## User options -----------------------------------------------------------

only_init_fit <- F

if (!only_init_fit) {
  nsim <- 50 # The number of simulation iterations
  par_type <- "auto"
  # par_type <- "doSeq"
  # par_type <- "doParallel"
  # par_type <- "doMPI"
  if (par_type == "auto") {
    if (nsim <= 1) {
      par_type <- "doSeq"
    } else if (interactive() && !isatty(stdout())) {
      par_type <- "doParallel"
    } else {
      par_type <- "doMPI"
    }
  }
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
link_str <- "probit"
nthres <- ncat - 1L
npreds_tot <- 50L
p0 <- as.integer(ceiling(npreds_tot * 0.20))
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
  # data-fitting model should usually be adjusted:
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
                           N = nobsv,
                           df_slab = 100, scale_slab = 1)

  npreds_grPL <- 0L # 1L
  ngrPL <- integer() # c(3L)
  stopifnot(identical(length(ngrPL), npreds_grPL))
  if (npreds_grPL > 0) {
    coefs_grPL <- list(
      c(0, seq(-2, 2, length.out = ngrPL[1] - 1L))
    )
    stopifnot(identical(length(coefs_grPL), npreds_grPL))
  }

  npreds_grGL <- 0L # 1L
  if (!exists("nobsv_per_grGL")) {
    # Number of observations per group (only roughly; some groups might get a
    # few observations more):
    nobsv_per_grGL <- integer() # c(2L)
  }
  ngrGL <- c(nobsv %/% nobsv_per_grGL)
  stopifnot(identical(length(ngrGL), npreds_grGL))
  if (npreds_grGL > 0) {
    coefs_grGL <- list(
      list("icpt" = rnorm(ngrGL[1], sd = 0.5)) # , "Xcont1" = <...>
    )
    stopifnot(identical(length(coefs_grGL), npreds_grGL))
  }

  ## Continuous predictors --------------------------------------------------

  dat_sim <- setNames(replicate(npreds_cont, {
    rnorm(nobsv_sim)
  }, simplify = FALSE), paste0("Xcont", seq_len(npreds_cont)))
  dat_sim <- as.data.frame(dat_sim)
  # Start constructing the linear predictor:
  eta <- drop(as.matrix(dat_sim) %*% coefs_cont)

  ## Population-level (PL) categorical predictors ---------------------------

  if (npreds_grPL == 1) {
    dat_sim$XgrPL1 <- sample(
      gl(n = ngrPL[1], k = floor(nobsv_sim / ngrPL[1]), length = nobsv_sim,
         labels = paste0("gr", seq_len(ngrPL[1])))
    )
    # Continue constructing the linear predictor:
    eta <- eta + coefs_grPL[[1]][dat_sim$XgrPL1]
  } else if (npreds_grPL > 1) {
    stop("This value of `npreds_grPL` is currently not supported.")
  }

  ## Group-level (GL) categorical predictors --------------------------------

  if (npreds_grGL == 1) {
    dat_sim$XgrGL1 <- sample(
      gl(n = ngrGL[1], k = floor(nobsv_sim / ngrGL[1]), length = nobsv_sim,
         labels = paste0("gr", seq_len(ngrGL[1])))
    )
    # Continue constructing the linear predictor:
    eta <- eta + coefs_grGL[[1]]$icpt[dat_sim$XgrGL1]
    # + coefs_grGL[[1]]$Xcont1[dat_sim$XgrGL1] * dat_sim$Xcont1
  } else if (npreds_grGL > 1) {
    stop("This value of `npreds_grGL` is currently not supported.")
  }

  ## Construct "epred" ------------------------------------------------------

  thres_eta <- sapply(thres, function(thres_k) {
    thres_k - eta
  })
  # Shouldn't be necessary (just to be safe): Emulate a single posterior draw:
  dim(thres_eta) <- c(1L, nobsv_sim, nthres)
  yprobs <- brms:::inv_link_cumulative(thres_eta, link = link_str)
  # Because of emulating a single posterior draw above:
  dim(yprobs) <- c(nobsv_sim, ncat)

  ## Response ---------------------------------------------------------------

  dat_sim$Y <- factor(sapply(seq_len(nobsv_sim), function(i) {
    sample(yunq, size = 1, prob = yprobs[i, ])
  }), ordered = TRUE)

  ## Formula ----------------------------------------------------------------

  voutc <- "Y"
  vpreds <- grep("^Xcont|^XgrPL", names(dat_sim), value = TRUE)
  vpreds_GL <- grep("^XgrGL", names(dat_sim), value = TRUE)
  if (length(vpreds_GL) > 0L) {
    if (!all(lengths(coefs_grGL) == 1)) {
      stop("Group-level slopes are currently not supported.")
    }
    vpreds_GL <- paste0("(1 | ", vpreds_GL, ")")
  }

  fml_sim <- as.formula(paste(
    voutc, "~", paste(c(vpreds, vpreds_GL), collapse = " + ")
  ))

  ## Check response categories ----------------------------------------------

  dat_train <- dat_sim[1:nobsv, , drop = FALSE]
  dat_test <- dat_sim[(nobsv + 1):nobsv_sim, , drop = FALSE]
  # Check that all response categories are present in training and test dataset:
  stopifnot(identical(levels(droplevels(dat_train$Y)), yunq) &&
              identical(levels(droplevels(dat_test$Y)), yunq))

  ## Output -----------------------------------------------------------------

  return(list(true_coefs_cont = coefs_cont,
              true_GLEs = if (npreds_grGL > 0) coefs_grGL[[1]]$icpt else NULL,
              dat = dat_train,
              fml = fml_sim,
              dat_indep = dat_test))
}

# Initial reference model fit ---------------------------------------------

if (only_init_fit) {
  # Fit the brms reference model once, so it only needs to be updated (and not
  # recompiled) during the simulation. MCMC convergence diagnostics are only
  # checked here, not during the simulation.
  sim_dat_etc <- dataconstructor()
  options(mc.cores = parallel::detectCores(logical = FALSE))
  if (packageVersion("cmdstanr") >= "0.5.3") {
    options(cmdstanr_write_stan_file_dir = ".")
  } else {
    options(cmdstanr_write_stan_file_dir = getwd())
  }
  ### Needed because the computing cluster complains about `p0` not found:
  par_ratio_sigti <- p0 / (npreds_tot - p0) * sigti
  cat("-----\npar_ratio_sigti:\n")
  print(par_ratio_sigti, digits = 9)
  cat("-----\n")
  stopifnot(all.equal(par_ratio_sigti, 0.264728344, tolerance = 1e-8))
  ###
  bfit <- brms::brm(
    formula = sim_dat_etc$fml,
    data = sim_dat_etc$dat,
    family = brms::cumulative(link = link_str),
    prior = brms::prior(horseshoe(df_slab = 100, scale_slab = 1,
                                  par_ratio = 0.264728344)) +
      brms::prior(normal(0, 2.5), class = "Intercept"),
    ### For backend = "rstan":
    # control = list(adapt_delta = 0.99), # , max_treedepth = 15L
    # init_r = 1,
    ###
    ### For backend = "cmdstanr":
    backend = "cmdstanr",
    adapt_delta = 0.99,
    # max_treedepth = 15L,
    init = 1,
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
    # init_r = 1,
    ###
    ### For backend = "cmdstanr":
    adapt_delta = 0.99,
    # max_treedepth = 15L,
    init = 1,
    ###
    seed = sample.int(.Machine$integer.max, 1)
  ))
}

# For running projpred's variable selection:
run_projpred <- function(refm_fit, dat_indep, latent = FALSE, ...) {
  d_indep <- list(
    data = dat_indep,
    offset = rep(0, nrow(dat_indep)),
    weights = rep(1, nrow(dat_indep)),
    y = if (latent) rep(NA, nrow(dat_indep)) else dat_indep$Y
  )
  if (latent) {
    d_indep$yOrig <- dat_indep$Y
  }
  time_bef <- Sys.time()
  vs <- projpred::varsel(refm_fit, d_test = d_indep, method = "forward",
                         latent = latent, ...)
  time_aft <- Sys.time()
  out_projpred <- list(
    time_vs = as.numeric(time_aft - time_bef, units = "mins"),
    soltrms = projpred::solution_terms(vs)
  )
  if (latent) {
    respOrig_vals <- c(TRUE) # , FALSE
  } else {
    respOrig_vals <- TRUE
  }
  respOrig_vals <- setNames(nm = respOrig_vals)
  names(respOrig_vals) <- paste0("respOrig_", names(respOrig_vals))
  out_projpred <- c(out_projpred, lapply(respOrig_vals, function(respOrig_val) {
    # Currently, we need to use the internal projpred function .tabulate_stats()
    # to obtain the reference model's performance:
    stats_man <- projpred:::.tabulate_stats(vs, stats = "mlpd",
                                            respOrig = respOrig_val)
    refsmms <- vs$summaries$ref
    if (latent && respOrig_val) refsmms <- refsmms$Orig
    return(list(
      refsmms = refsmms,
      refstat = stats_man$value[stats_man$size == Inf],
      # plot_obj = plot(vs, deltas = TRUE, stats = "mlpd",
      #                 respOrig = respOrig_val),
      sgg_size = projpred::suggest_size(vs, stat = "mlpd",
                                        respOrig = respOrig_val),
      smmry = summary(vs, deltas = TRUE, stats = "mlpd",
                      type = c("mean", "se", "lower", "upper"),
                      respOrig = respOrig_val)$selection
    ))
  }))
  return(out_projpred)
}

sim_runner <- function(...) {
  foreach(
    sit = seq_len(nsim), # Needs package "iterators": icount(nsim),
    # .packages = c("brms", "projpred"), # , "rstanarm"
    .export = c("rhorseshoe", "dataconstructor", "fit_ref", "run_projpred",
                "nobsv", "nobsv_indep", "ncat", "yunq", "link_str", "nthres",
                "npreds_tot", "p0", "sigti", "bfit"),
    # .noexport = c("<object_name>"),
    .options.snow = list(attachExportEnv = TRUE)
  ) %dorng% {
    cat("\nSimulation iteration: ", sit, "\n", sep = "")
    Rseed <- .Random.seed
    sim_dat_etc <- try(dataconstructor(), silent = TRUE)
    dat_failed <- inherits(sim_dat_etc, "try-error")
    if (dat_failed) {
      save(sit, Rseed, file = "failed.rda")
      stop("The data generation failed in simulation iteration ", sit, ". ",
           "Error message: \"", attr(sim_dat_etc, "condition")$message, "\". ",
           "Objects for replicating this failure were saved to ",
           "\"failed.rda\". Use `loaded_objs <- load(\"failed.rda\")` to ",
           "restore it.")
    }
    refm_fit <- fit_ref(dat = sim_dat_etc$dat, fml = sim_dat_etc$fml)
    seed_vs <- sample.int(.Machine$integer.max, 1)
    projpred_aug <- try(run_projpred(refm_fit,
                                     dat_indep = sim_dat_etc$dat_indep,
                                     seed = seed_vs,
                                     ...),
                        silent = TRUE)
    projpred_lat <- try(run_projpred(refm_fit,
                                     dat_indep = sim_dat_etc$dat_indep,
                                     seed = seed_vs,
                                     latent = TRUE,
                                     ...),
                        silent = TRUE)
    aug_failed <- inherits(projpred_aug, "try-error")
    lat_failed <- inherits(projpred_lat, "try-error")
    if (aug_failed || lat_failed) {
      dot_args <- list(...)
      save(sit, refm_fit, seed_vs, dot_args, Rseed, file = "failed.rda")
      mssgs_failed <- NULL
      if (aug_failed) {
        which_failed <- "augmented-data"
        mssgs_failed <- c(mssgs_failed, attr(projpred_aug, "condition")$message)
      } else if (lat_failed) {
        which_failed <- "latent"
        mssgs_failed <- c(mssgs_failed, attr(projpred_lat, "condition")$message)
      }
      if (aug_failed && lat_failed) {
        which_failed <- "both"
      }
      stop("The following of the two projpred runs (augmented-data, latent) ",
           "failed in simulation iteration ", sit, ": ", which_failed, ". ",
           "Error message(s for augmented-data and latent run, respectively): ",
           "\"", paste(mssgs_failed, collapse = "\", \""),
           "\". Objects for replicating this failure were saved to ",
           "\"failed.rda\". Use `loaded_objs <- load(\"failed.rda\")` to ",
           "restore it.")
    }
    stopifnot(all.equal(projpred_aug$respOrig_TRUE$refsmms,
                        projpred_lat$respOrig_TRUE$refsmms))
    projpred_aug$respOrig_TRUE$refsmms <- NULL
    projpred_lat$respOrig_TRUE$refsmms <- NULL
    # projpred_lat$respOrig_FALSE$refsmms <- NULL
    return(list(
      aug = projpred_aug,
      lat = projpred_lat,
      true_coefs_cont = sim_dat_etc$true_coefs_cont,
      true_GLEs = sim_dat_etc$true_GLEs
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

source("gg_to_tikz/tikzpicture-template.R")
ggsave_cust <- function(fname_no_ext, plot = ggplot2::last_plot(),
                        width = 6, height = width * 0.618,
                        timestamp = FALSE, verbose = FALSE, ...) {
  ggplot2::ggsave(filename = paste0(fname_no_ext, ".pdf"), plot = plot,
                  width = width, height = height)
  save_tikz_plot(plot = plot, filename = paste0(fname_no_ext, ".tex"),
                 width = width, height = height,
                 timestamp = timestamp, verbose = verbose, ...)
  return(invisible(TRUE))
}

if (!dir.exists("figs")) dir.create("figs")

## True population-level coefficients -------------------------------------
## (i.e., the draws from the regularized horseshoe distribution)

true_coefs_cont <- unlist(lapply(simres, "[[", "true_coefs_cont"))
true_coefs_cont <- data.frame(coef = true_coefs_cont)
cat("\n-----\n")
cat("Proportion of (population-level) coefficient draws with absolute value >",
    "0.5 (across all simulation iterations and all `npreds_tot` coefficient",
    "draws within each simulation iteration):\n")
print(proportions(table(abs(true_coefs_cont$coef) > 0.5, useNA = "ifany")))
cat("-----\n")
gg_true_coefs_cont <- ggplot2::ggplot(data = true_coefs_cont,
                                      mapping = ggplot2::aes(x = coef)) +
  ggplot2::geom_histogram(bins = 40) # + ggplot2::geom_density()
ggsave_cust(file.path("figs", "true_coefs_cont"))

## True group-level effects -----------------------------------------------

if (!all(sapply(lapply(simres, "[[", "true_GLEs"), is.null))) {
  true_GLEs <- do.call(cbind, lapply(simres, "[[", "true_GLEs"))
  dimnames(true_GLEs) <- list(
    grp = paste0("grp", seq_len(nrow(true_GLEs))),
    simiter = paste0("simiter", seq_len(ncol(true_GLEs)))
  )
  stopifnot(!all(lengths(apply(true_GLEs, 1, unique, simplify = FALSE)) == 1))
  cat("\n-----\n")
  cat("Quartiles of the true group-level effects (across",
      "all simulation iterations):\n")
  print(quantile(true_GLEs))
  cat("-----\n")
  if (nsim <= 10) {
    cat("\n-----\n")
    cat("Empirical SDs of the true group-level effects (for",
        "all simulation iterations):\n")
    print(apply(true_GLEs, 2, sd))
    cat("-----\n")
  } else {
    cat("\n-----\n")
    cat("Quartiles of the empirical SDs of the true group-level effects",
        "(across all simulation iterations):\n")
    print(quantile(apply(true_GLEs, 2, sd)))
    cat("-----\n")
  }
}

## Runtime ----------------------------------------------------------------

mins_vs <- do.call(rbind, lapply(seq_along(simres), function(sim_idx) {
  return(data.frame(sim_idx = sim_idx,
                    t_aug = simres[[sim_idx]]$aug$time_vs,
                    t_lat = simres[[sim_idx]]$lat$time_vs))
}))
mins_vs <- reshape(
  mins_vs,
  direction = "long",
  v.names = "Runtime [min]",
  varying = list("Runtime [min]" = grep("^t_", names(mins_vs), value = TRUE)),
  timevar = "prj_meth",
  times = c("Augmented-data", "Latent"),
  idvar = "sim_idx_ch",
  sep = "_"
)
stopifnot(identical(mins_vs$sim_idx_ch, mins_vs$sim_idx))
mins_vs$sim_idx_ch <- NULL
gg_time <- ggplot2::ggplot(
  data = mins_vs,
  mapping = ggplot2::aes(x = prj_meth, y = `Runtime [min]`)
) +
  ggplot2::geom_boxplot() + # ggplot2::geom_violin() +
  ggplot2::labs(x = "Projection method") +
  ggplot2::coord_cartesian(ylim = c(0, NA)) # + ggplot2::coord_flip()
ggsave_cust(file.path("figs", "time"),
            width = 0.5 * 6, height = 0.75 * 6 * 0.618)

## Solution paths ---------------------------------------------------------

same_solpths <- sapply(simres, function(simres_i) {
  identical(simres_i$aug$soltrms,
            simres_i$lat$soltrms)
})
cat("\n-----\n")
cat("Differing solution paths:\n")
cat("There are", sum(!same_solpths), "simulation iterations with the solution",
    "path differing between augmented-data and latent projection. The first 10",
    "(if there are less, then only those) in detail:\n")
for (sim_idx in head(seq_along(simres)[!same_solpths], 10)) {
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

plotter_ovrlay <- function(prj_meth, eval_scale = "response",
                           ylim_full = NULL, ylim_eb = NULL) {
  if (prj_meth == "aug") {
    title_gg <- "Augmented-data"
    stopifnot(eval_scale == "response")
    respOrig_nm <- paste0("respOrig_", TRUE)
  } else if (prj_meth == "lat") {
    title_gg <- "Latent"
    respOrig_nm <- paste0("respOrig_", eval_scale == "response")
  }
  title_gg <- paste0(title_gg, " (evaluation scale: ", eval_scale, ")")
  y_chr <- setdiff(names(simres[[1L]][[prj_meth]][[respOrig_nm]]$smmry),
                   c("solution_terms", "se", "lower", "upper", "size"))
  stopifnot(length(y_chr) == 1)
  plotdat <- do.call(rbind, lapply(seq_along(simres), function(sim_idx) {
    cbind(sim_idx = sim_idx,
          simres[[sim_idx]][[prj_meth]][[respOrig_nm]]$smmry[
            c("size", "se", y_chr)
          ])
  }))
  ylab <- paste0("$\\Delta\\mathrm{", toupper(y_chr), "}_{\\mathrm{", prj_meth,
                 "}}$")
  ggobj <- ggplot2::ggplot(data = plotdat,
                           mapping = ggplot2::aes(x = size,
                                                  y = .data[[y_chr]],
                                                  group = sim_idx,
                                                  alpha = I(0.4))) +
    ggplot2::geom_hline(yintercept = 0,
                        color = "firebrick",
                        linetype = "dashed") +
    ### Only applicable when using the extended suggest_size() heuristics by
    ### setting argument `thres_elpd` to `4`:
    # ggplot2::geom_hline(yintercept = -4 / nobsv_indep,
    #                     color = "dodgerblue",
    #                     linetype = "dotdash") +
    ###
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(
      # title = title_gg,
      x = "Submodel size",
      # y = bquote(Delta*.(toupper(y_chr)))
      y = ylab
    )
  fnm_base <- paste(y_chr, prj_meth, eval_scale, sep = "_")
  ggsave_cust(file.path("figs", fnm_base))
  ggobj_full <- ggobj +
    ggplot2::coord_cartesian(ylim = ylim_full)
  ggsave_cust(file.path("figs", paste0(fnm_base, "_full")))
  ggobj_zoom <- ggobj +
    ggplot2::coord_cartesian(ylim = c(-0.75, 0.05))
  ggsave_cust(file.path("figs", paste0(fnm_base, "_zoom")))
  nsub_eb <- 5L
  sub_idxs <- seq_len(nsub_eb) # sample.int(length(simres), size = nsub_eb)
  plotdat_sub <- plotdat[plotdat$sim_idx %in% sub_idxs, ]
  ggobj_eb <- ggplot2::ggplot(data = plotdat_sub,
                              mapping = ggplot2::aes(x = size,
                                                     y = .data[[y_chr]],
                                                     ymin = .data[[y_chr]] - se,
                                                     ymax = .data[[y_chr]] + se,
                                                     group = sim_idx,
                                                     alpha = I(0.4))) +
    ggplot2::geom_hline(yintercept = 0,
                        color = "firebrick",
                        linetype = "dashed") +
    ### Only applicable when using the extended suggest_size() heuristics by
    ### setting argument `thres_elpd` to `4`:
    # ggplot2::geom_hline(yintercept = -4 / nobsv_indep,
    #                     color = "dodgerblue",
    #                     linetype = "dotdash") +
    ###
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar() +
    ggplot2::labs(
      # title = title_gg,
      x = "Submodel size",
      # y = bquote(Delta*.(toupper(y_chr)))
      y = ylab
    ) +
    ggplot2::coord_cartesian(ylim = ylim_eb)
  ggsave_cust(file.path("figs",
                        paste(y_chr, prj_meth, eval_scale, "eb", sep = "_")))
  ggobj_se <- ggplot2::ggplot(data = plotdat,
                              mapping = ggplot2::aes(x = factor(size),
                                                     y = se,
                                                     alpha = I(0.4))) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(width = 0.25, height = 0) +
    ggplot2::labs(
      x = "Submodel size",
      # y = bquote(SE(Delta*.(toupper(y_chr))))
      y = paste0("SE($\\Delta$", toupper(y_chr), ")")
    )
  ggsave_cust(file.path("figs",
                        paste(y_chr, prj_meth, eval_scale, "se", sep = "_")))
  return(list(succ_ind = TRUE, ggobj = ggobj, ggobj_full = ggobj_full,
              ggobj_zoom = ggobj_zoom, ggobj_eb = ggobj_eb,
              ggobj_se = ggobj_se))
}
comm_lat <- plotter_ovrlay(prj_meth = "lat")
# comm_lat_nonOrig <- plotter_ovrlay(prj_meth = "lat", eval_scale = "latent")
ylim_lat <- ggplot2::ggplot_build(
  comm_lat$ggobj
)$layout$panel_scales_y[[1]]$range$range
ylim_lat_eb <- ggplot2::ggplot_build(
  comm_lat$ggobj_eb
)$layout$panel_scales_y[[1]]$range$range
comm_aug <- plotter_ovrlay(prj_meth = "aug", ylim_full = ylim_lat,
                           ylim_eb = ylim_lat_eb)
library(patchwork)
gg_aug_lat <- comm_aug$ggobj_full / comm_lat$ggobj_full
ggsave_cust(file.path("figs", "aug_lat"), height = 2 * 6 * 0.618)

plotter_ovrlay_diff <- function(eval_scale = "response") {
  stopifnot(eval_scale == "response")
  respOrig_nm_aug <- paste0("respOrig_", TRUE)
  respOrig_nm_lat <- paste0("respOrig_", eval_scale == "response")
  title_gg <- "Performance diff."
  title_gg <- paste0(title_gg, " (evaluation scale: ", eval_scale, ")")

  # Check that the reference model (performance) is the same, so that the
  # difference of the Delta MLPDs can be interpreted as the difference of the
  # submodel MLPDs (i.e., the reference model MLPD cancels out):
  refstats_aug <- do.call(c, lapply(seq_along(simres), function(sim_idx) {
    simres[[sim_idx]]$aug[[respOrig_nm_aug]]$refstat
  }))
  refstats_lat <- do.call(c, lapply(seq_along(simres), function(sim_idx) {
    simres[[sim_idx]]$lat[[respOrig_nm_lat]]$refstat
  }))
  stopifnot(all.equal(refstats_aug, refstats_lat,
                      tolerance = .Machine$double.eps))
  refstats <- refstats_aug

  smmry_nms <- names(simres[[1L]]$aug[[respOrig_nm_aug]]$smmry)
  stopifnot(identical(smmry_nms,
                      names(simres[[1L]]$lat[[respOrig_nm_lat]]$smmry)))
  y_chr <- setdiff(smmry_nms,
                   c("solution_terms", "se", "lower", "upper", "size"))
  stopifnot(length(y_chr) == 1)
  y_chr_aug <- paste(y_chr, "aug", sep = "_")
  y_chr_lat <- paste(y_chr, "lat", sep = "_")
  plotdat <- do.call(rbind, lapply(seq_along(simres), function(sim_idx) {
    smmry_aug <- simres[[sim_idx]]$aug[[respOrig_nm_aug]]$smmry
    smmry_lat <- simres[[sim_idx]]$lat[[respOrig_nm_lat]]$smmry
    stopifnot(identical(smmry_aug["size"], smmry_lat["size"]))
    cbind(sim_idx = sim_idx, smmry_aug["size"],
          setNames(smmry_aug[c(y_chr, "se")], c(y_chr_aug, "se_aug")),
          setNames(smmry_lat[c(y_chr, "se")], c(y_chr_lat, "se_lat")))
  }))
  y_chr_diff <- paste("diff", y_chr, sep = "_")
  plotdat[[y_chr_diff]] <- plotdat[[y_chr_lat]] - plotdat[[y_chr_aug]]
  plotdat$diff_se <- plotdat[["se_lat"]] - plotdat[["se_aug"]]
  ggobj <- ggplot2::ggplot(data = plotdat,
                           mapping = ggplot2::aes(x = size,
                                                  y = .data[[y_chr_diff]],
                                                  group = sim_idx,
                                                  alpha = I(0.4))) +
    ggplot2::geom_hline(yintercept = 0,
                        color = "gray30",
                        linetype = "dotted") +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(
      # title = title_gg,
      x = "Submodel size",
      y = bquote(.(toupper(y_chr))[lat] - .(toupper(y_chr))[aug])
    )
  fnm_base <- paste(y_chr_diff, eval_scale, sep = "_")
  ggsave_cust(file.path("figs", fnm_base))
  ggobj_se <- ggplot2::ggplot(data = plotdat,
                              mapping = ggplot2::aes(x = factor(size),
                                                     y = diff_se,
                                                     alpha = I(0.4))) +
    ggplot2::geom_hline(yintercept = 0,
                        color = "gray30",
                        linetype = "dotted") +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(width = 0.25, height = 0) +
    ggplot2::labs(
      x = "Submodel size",
      y = paste0(
        "SE($\\Delta\\mathrm{", toupper(y_chr), "}_{\\mathrm{lat}}$) - ",
        "SE($\\Delta\\mathrm{", toupper(y_chr), "}_{\\mathrm{aug}}$)"
      )
    )
  ggsave_cust(file.path("figs", paste("diff_se", eval_scale, sep = "_")))
  return(list(succ_ind = TRUE, ggobj = ggobj, ggobj_se = ggobj_se,
              refstats = refstats, q_diff_se = quantile(plotdat$diff_se)))
}
diff_out <- plotter_ovrlay_diff()
cat("\n-----\n")
cat("Quartiles of the reference model's performance statistic (across all",
    "simulation iterations):\n")
print(quantile(diff_out$refstats))
cat("-----\n")
cat("\n-----\n")
cat("Quartiles of the SE difference (across all simulation iterations and all",
    "submodel sizes):\n")
print(diff_out$q_diff_se)
cat("-----\n")

## Suggested sizes --------------------------------------------------------

sgger_size <- function(sim_idx, eval_scale_lat = "response") {
  respOrig_nm_aug <- paste0("respOrig_", TRUE)
  respOrig_nm_lat <- paste0("respOrig_", eval_scale_lat == "response")
  return(c(
    sgg_size_aug = simres[[sim_idx]]$aug[[respOrig_nm_aug]]$sgg_size,
    sgg_size_lat = simres[[sim_idx]]$lat[[respOrig_nm_lat]]$sgg_size
  ))
}
for (eval_scale_lat_val in c("response")) { # , "latent"
  cat("\n----------\n")
  cat("Evaluation scale for the latent projection (CAUTION: always using ",
      "response scale for the augmented-data projection): ", eval_scale_lat_val,
      "\n", sep = "")
  sgg_sizes <- sapply(seq_along(simres), sgger_size,
                      eval_scale_lat = eval_scale_lat_val)
  cat("\n-----\n")
  cat("Suggested sizes (printing only the first 10 simulation iterations):\n")
  print(sgg_sizes[, head(seq_len(ncol(sgg_sizes)), 10), drop = FALSE])
  cat("-----\n")
  if (anyNA(sgg_sizes)) {
    warning("Found suggested sizes which are `NA`.")
  }
  sgg_sizes_lat_minus_aug <- apply(sgg_sizes, 2, function(x) {
    x["sgg_size_lat"] - x["sgg_size_aug"]
  })
  sgg_sizes_lat_minus_aug <- factor(sgg_sizes_lat_minus_aug)
  sgg_sizes_NA <- apply(sgg_sizes, 2, function(x) {
    if (!is.na(x["sgg_size_lat"]) && is.na(x["sgg_size_aug"])) {
      return("NA (aug.)")
    } else if (is.na(x["sgg_size_lat"]) && !is.na(x["sgg_size_aug"])) {
      return("NA (lat.)")
    } else if (is.na(x["sgg_size_lat"]) && is.na(x["sgg_size_aug"])) {
      return("NA (both)")
    } else {
      return(NA)
    }
  })
  sgg_sizes_NA <- factor(sgg_sizes_NA,
                         levels = c("NA (aug.)", "NA (lat.)", "NA (both)"))
  stopifnot(identical(is.na(sgg_sizes_lat_minus_aug), !is.na(sgg_sizes_NA)))
  sgg_sizes_lat_minus_aug <- factor(
    sgg_sizes_lat_minus_aug,
    levels = union(levels(sgg_sizes_lat_minus_aug), levels(sgg_sizes_NA))
  )
  sgg_sizes_lat_minus_aug[is.na(sgg_sizes_lat_minus_aug)] <- sgg_sizes_NA[
    !is.na(sgg_sizes_NA)
  ]
  cat("\n-----\n")
  cat("Differences of the suggested sizes (latent minus augmented-data):\n")
  sgg_sizes_tab <- table(sgg_sizes_lat_minus_aug, useNA = "ifany")
  print(sgg_sizes_tab)
  print(proportions(sgg_sizes_tab))
  cat("-----\n")
  xlab_long <- bquote(G[lat] - G[aug])
  title_gg <- "Suggested size diff."
  title_gg <- paste0(title_gg, " (evaluation scale for the latent projection: ",
                     eval_scale_lat_val, ")")
  gg_sgg_sizes_diff <- ggplot2::ggplot(
    data = data.frame(sgg_sizes_lat_minus_aug),
    mapping = ggplot2::aes(x = sgg_sizes_lat_minus_aug)
  ) +
    ggplot2::geom_bar() +
    ggplot2::labs(
      # title = title_gg,
      x = xlab_long,
      y = paste0("Number of simulation iterations (total: ", nsim, ")")
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty())
  ggsave_cust(
    file.path("figs", paste0("sgg_sizes_diff_", eval_scale_lat_val, "Lat"))
  )
  cat("----------\n")
}

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
