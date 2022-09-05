# MCMC diagnostics (convergence (stationarity, mixing)) -------------------

check_MCMC_diagn <- function(
    C_stanfit,
    n_chains_spec = 4L,
    C_pars = character(),
    C_selPars = character(),
    overview = FALSE,
    HMC_auto = FALSE,
    sampler_pars = FALSE,
    exclude_NAs = FALSE,
    ESS_classical = FALSE, # Deactivated by default, but often makes sense.
    AC = FALSE,
    rstan_plots = FALSE,
    bayesplot_MCMC = FALSE,
    trace_chain_highlight = 1,
    bayesplot_HMC = FALSE,
    bayesplot_parallelCoord = FALSE,
    new_plots = FALSE, # Deactivated by default, but often makes sense.
    ...
) {
  ## Preparations -----------------------------------------------------------

  if (!("data.table" %in% loadedNamespaces() &&
        "data.table" %in% .packages())) {
    suppressPackageStartupMessages({
      library(data.table)
    })
  }

  out_list <- list()

  C_draws_arr <- as.array(C_stanfit, ...)
  n_drawsPerChain <- dim(C_draws_arr)[1]
  n_chains <- dim(C_draws_arr)[2] # n_chains <- C_stanfit@sim$chains
  C_draws_mat <- as.matrix(C_stanfit, ...)
  n_draws <- nrow(C_draws_mat)
  stopifnot(identical(n_draws, n_chains * n_drawsPerChain))

  if (identical(length(C_pars), 0L)) {
    C_pars <- colnames(C_draws_mat)
  }
  if (identical(length(C_selPars), 0L)) {
    C_selPars <- C_pars # C_selPars <- C_pars[seq_len(min(2, length(C_pars)))]
  }

  # Check that the mode of the resulting "stanfit" object is the "normal" mode
  # (0L), i.e. neither test gradient mode (1L) nor error mode (2L):
  stopifnot(identical(C_stanfit@mode, 0L))

  # Check the number of chains in the output:
  if (n_chains < n_chains_spec) {
    warning("At least one chain exited with an error.",
            "The posterior results should not be used.")
  }

  ## HMC-specific diagnostics -----------------------------------------------

  if (isTRUE(HMC_auto)) {
    rstan::check_hmc_diagnostics(C_stanfit)
  }

  C_div <- rstan::get_num_divergent(C_stanfit)
  C_div_OK <- identical(C_div, 0L)

  C_tree <- rstan::get_num_max_treedepth(C_stanfit)
  C_tree_OK <- identical(C_tree, 0L)

  C_EBFMI <- rstan::get_bfmi(C_stanfit)
  C_EBFMI_OK <- all(C_EBFMI >= 0.3)

  out_list <- c(out_list,
                list("div" = C_div,
                     "div_OK" = C_div_OK,
                     "tree" = C_tree,
                     "tree_OK" = C_tree_OK,
                     "EBFMI" = C_EBFMI,
                     "EBFMI_OK" = C_EBFMI_OK))

  if (isTRUE(sampler_pars)) {
    C_sampler_params <- rstan::get_sampler_params(C_stanfit, inc_warmup = FALSE)
    C_sampler_params <- lapply(
      seq_along(C_sampler_params),
      function(idx_chain) {
        data.table(C_sampler_params[[idx_chain]],
                   "chain" = paste0("chain", idx_chain))
      }
    )
    C_sampler_params <- rbindlist(C_sampler_params, fill = TRUE)

    C_sampler_params_smmry <- C_sampler_params[
      ,
      .("mean_accept_stat" = mean(accept_stat__),
        "max_treedepth" = max(treedepth__),
        "sum_divergent" = sum(divergent__),
        "median_stepsize" = median(stepsize__),
        "median_leapfrog" = median(n_leapfrog__),
        "max_leapfrog" = max(n_leapfrog__),
        "median_energy" = mean(energy__)),
      by = chain
    ]

    out_list <- c(out_list,
                  list("sampler_params" = C_sampler_params,
                       "sampler_params_smmry" = C_sampler_params_smmry))
  }

  ## General MCMC diagnostics -----------------------------------------------

  ### Bulk-ESS --------------------------------------------------------------
  ### (cf. Vehtari et al., 2021, DOI: 10.1214/20-BA1221)

  C_essBulk <- apply(C_draws_arr, MARGIN = 3, FUN = rstan::ess_bulk)
  if (isTRUE(exclude_NAs)) {
    C_essBulk_OK <- all(C_essBulk > 100 * n_chains)
    if (any(is.na(C_essBulk))) {
      warning("\"C_essBulk\" contains at least one NA.")
      C_essBulk_OK <- all(C_essBulk > 100 * n_chains, na.rm = TRUE)
    }
  } else {
    if (any(is.na(C_essBulk))) {
      C_essBulk_OK <- FALSE
    } else {
      C_essBulk_OK <- all(C_essBulk > 100 * n_chains)
    }
  }

  # Bulk-ESS ratio to total number of (post-warmup) draws:
  C_essBulkRatio <- C_essBulk / n_draws

  out_list <- c(out_list,
                list("essBulk" = C_essBulk,
                     "essBulk_OK" = C_essBulk_OK,
                     "essBulkRatio" = C_essBulkRatio))

  ### "New" R-hat -----------------------------------------------------------
  ### (cf. Vehtari et al., 2021, DOI: 10.1214/20-BA1221)

  C_rhat <- apply(C_draws_arr, MARGIN = 3, FUN = rstan::Rhat)
  if (isTRUE(exclude_NAs)) {
    C_rhat_OK <- all(C_rhat < 1.01)
    if (any(is.na(C_rhat))) {
      warning("\"C_rhat\" contains at least one NA.")
      C_rhat_OK <- all(C_rhat < 1.01, na.rm = TRUE)
    }
  } else {
    if (any(is.na(C_rhat))) {
      C_rhat_OK <- FALSE
    } else {
      C_rhat_OK <- all(C_rhat < 1.01)
    }
  }

  out_list <- c(out_list,
                list("rhat" = C_rhat,
                     "rhat_OK" = C_rhat_OK))

  ### Tail-ESS --------------------------------------------------------------
  ### (cf. Vehtari et al., 2021, DOI: 10.1214/20-BA1221)

  C_essTail <- apply(C_draws_arr, MARGIN = 3, FUN = rstan::ess_tail)
  if (isTRUE(exclude_NAs)) {
    C_essTail_OK <- all(C_essTail > 100 * n_chains)
    if (any(is.na(C_essTail))) {
      warning("\"C_essTail\" contains at least one NA.")
      C_essTail_OK <- all(C_essTail > 100 * n_chains, na.rm = TRUE)
    }
  } else {
    if (any(is.na(C_essTail))) {
      C_essTail_OK <- FALSE
    } else {
      C_essTail_OK <- all(C_essTail > 100 * n_chains)
    }
  }

  # Tail-ESS ratio to total number of (post-warmup) draws:
  C_essTailRatio <- C_essTail / n_draws

  out_list <- c(out_list,
                list("essTail" = C_essTail,
                     "essTail_OK" = C_essTail_OK,
                     "essTailRatio" = C_essTailRatio))

  ### ESS -------------------------------------------------------------------
  ### (only necessary if posterior mean is of interest which is not advised for
  ### heavy-tailed priors with infinite mean, especially if the data is weakly
  ### identifying)

  if (isTRUE(ESS_classical)) {
    if (!("rstan" %in% loadedNamespaces() && "rstan" %in% .packages())) {
      suppressPackageStartupMessages({
        library(rstan)
      })
      on.exit({
        detach("package:rstan")
      })
    }
    C_smmry <- summary(C_stanfit, pars = C_pars)$summary
    C_ess <- C_smmry[, "n_eff"]
    if (isTRUE(exclude_NAs)) {
      C_ess_OK <- all(C_ess > 100 * n_chains)
      if (any(is.na(C_ess))) {
        warning("\"C_ess\" contains at least one NA.")
        C_ess_OK <- all(C_ess > 100 * n_chains, na.rm = TRUE)
      }
    } else {
      if (any(is.na(C_ess))) {
        C_ess_OK <- FALSE
      } else {
        C_ess_OK <- all(C_ess > 100 * n_chains)
      }
    }

    # ESS ratio to total number of (post-warmup) draws:
    C_essRatio <- C_ess / n_draws

    out_list <- c(out_list,
                  list("ess" = C_ess,
                       "ess_OK" = C_ess_OK,
                       "essRatio" = C_essRatio))
  }

  ### Autocorrelation -------------------------------------------------------

  if (isTRUE(AC)) {
    C_ac_obj <- rstan::stan_ac(C_stanfit, pars = C_pars, ncol = 1, lags = 10) # , separate_chains = TRUE
    C_ac <- as.data.table(C_ac_obj$data)
    C_ac <- C_ac[lag != 0L, ]
    C_ac_max <- C_ac[, .("ac_max" = max(abs(ac))), parameters]

    out_list <- c(out_list,
                  list("ac" = C_ac,
                       "ac_max" = C_ac_max))
  }

  ## Overall check for all MCMC diagnostics ---------------------------------

  C_OKs <- c("C_div_OK" = C_div_OK, "C_tree_OK" = C_tree_OK,
             "C_EBFMI_OK" = C_EBFMI_OK, "C_essBulk_OK" = C_essBulk_OK,
             "C_rhat_OK" = C_rhat_OK, "C_essTail_OK" = C_essTail_OK)
  if (exists("C_ess_OK")) {
    C_OKs <- c(C_OKs, "C_ess_OK" = C_ess_OK)
  }
  # Improve readability of the names:
  names(C_OKs) <- sub("^C_", "", names(C_OKs))
  stopifnot(identical(
    names(C_OKs),
    names(unlist(out_list[grep("_OK$", names(out_list))]))
  ))
  names(C_OKs) <- sub("_OK$", "", names(C_OKs))
  C_all_OK <- all(C_OKs)

  if (!C_all_OK) {
    warning("At least one MCMC diagnostic is worrying. This should be ",
            "inspected (see the output from this function for details). ",
            "In general, this indicates that the posterior results should not ",
            "be used. The concerned diagnostic(s) is/are: ",
            paste(names(C_OKs)[!C_OKs], collapse = ", "))
  }

  ## Overview ---------------------------------------------------------------

  if (isTRUE(overview)) {
    if (!("rstan" %in% loadedNamespaces() && "rstan" %in% .packages())) {
      suppressPackageStartupMessages({
        library(rstan)
      })
      on.exit({
        detach("package:rstan")
      })
    }

    ### rstan -----------------------------------------------------------------

    C_smmry <- summary(C_stanfit, pars = C_pars)$summary
    C_mon_old <- rstan::monitor(C_stanfit)

    ### "New" diagnostics -----------------------------------------------------
    ### (cf. Vehtari et al., 2021, DOI: 10.1214/20-BA1221)

    devtools::source_url("https://raw.githubusercontent.com/avehtari/rhat_ess/master/code/monitornew.R")
    devtools::source_url("https://raw.githubusercontent.com/avehtari/rhat_ess/master/code/monitorplot.R")

    C_mon <- monitor(C_stanfit)
    C_monex <- monitor_extra(C_stanfit)

    out_list <- c(out_list,
                  list("smmry" = C_smmry,
                       "mon_old" = C_mon_old,
                       "mon" = C_mon,
                       "monex" = C_monex))
  }

  ## Diagnostic plots -------------------------------------------------------

  ### rstan -----------------------------------------------------------------

  if (isTRUE(rstan_plots)) {
    # NOTE: As may be seen from rstan:::pairs.stanfit(), for 4 chains, using
    # "condition = NULL" means to plot chains 1 and 2 in the lower panel and
    # chains 3 and 4 in the upper panel.
    print(pairs(C_stanfit, pars = C_selPars, condition = NULL))
    print(pairs(C_stanfit, pars = C_selPars, condition = "accept_stat__"))
    print(pairs(C_stanfit, pars = C_selPars, condition = "divergent__"))

    print(rstan::traceplot(C_stanfit, pars = C_selPars, inc_warmup = TRUE,
                           nrow = 2))

    for (idx in seq_along(C_selPars)) {
      print(rstan::stan_par(C_stanfit, par = C_selPars[idx], chain = 0))
    }
  }

  ### bayesplot -------------------------------------------------------------

  if (isTRUE(bayesplot_MCMC) || isTRUE(bayesplot_HMC) ||
      isTRUE(bayesplot_parallelCoord)) {
    bayesplot::color_scheme_set("blue")
    # Extract log-posterior and NUTS statistics:
    C_np <- bayesplot::nuts_params(C_stanfit) # , inc_warmup = TRUE
  }

  #### General MCMC plots ---------------------------------------------------

  if (isTRUE(bayesplot_MCMC)) {
    print(bayesplot::mcmc_trace(C_stanfit, pars = C_selPars, np = C_np))
    print(bayesplot::mcmc_trace_highlight(C_stanfit, pars = C_selPars,
                                          highlight = trace_chain_highlight))
    print(bayesplot::mcmc_rank_hist(C_stanfit, pars = C_selPars))
    print(bayesplot::mcmc_rank_overlay(C_stanfit, pars = C_selPars))
  }

  #### HMC-specific plots ---------------------------------------------------

  if (isTRUE(bayesplot_HMC)) {
    bayesplot::color_scheme_set("darkgray")
    C_lp <- log_posterior(C_stanfit) # , inc_warmup = TRUE
    # available_mcmc(pattern = "_nuts_")
    print(bayesplot::mcmc_nuts_acceptance(x = C_np, lp = C_lp))
    print(bayesplot::mcmc_nuts_divergence(x = C_np, lp = C_lp))
    print(bayesplot::mcmc_nuts_energy(x = C_np)) # , merge_chains = TRUE
    print(bayesplot::mcmc_nuts_stepsize(x = C_np, lp = C_lp))
    print(bayesplot::mcmc_nuts_treedepth(x = C_np, lp = C_lp))
    bayesplot::color_scheme_set("blue")
  }

  #### Parallel coordinates plot --------------------------------------------
  #### (useful for "debugging" divergent transitions)

  ### WARNING: Might take long!:
  if (isTRUE(bayesplot_parallelCoord)) {
    print(bayesplot::mcmc_parcoord(C_draws_arr, np = C_np))
  }
  ###

  ### "New" diagnostic plots ------------------------------------------------
  ### (cf. Vehtari et al., 2021, DOI: 10.1214/20-BA1221)

  if (isTRUE(new_plots)) {
    if (!exists("plot_local_ess") ||
        !exists("plot_quantile_ess") ||
        !exists("plot_change_ess")) { #  || !exists("mcmc_hist_r_scale")
      devtools::source_url("https://raw.githubusercontent.com/avehtari/rhat_ess/master/code/monitornew.R")
      devtools::source_url("https://raw.githubusercontent.com/avehtari/rhat_ess/master/code/monitorplot.R")
    }

    if (!("rstan" %in% loadedNamespaces() && "rstan" %in% .packages())) {
      suppressPackageStartupMessages({
        library(rstan)
      })
      on.exit({
        detach("package:rstan")
      })
    }
    for (selPar_i in C_selPars) {
      print(
        plot_local_ess(fit = C_stanfit, par = selPar_i, nalpha = 20) # , rank = FALSE
      )
      print(
        plot_quantile_ess(fit = C_stanfit, par = selPar_i, nalpha = 40) # , rank = FALSE
      )
      print(
        plot_change_ess(fit = C_stanfit, par = selPar_i)
      )
      ### Not needed because bayesplot::mcmc_rank_hist() and
      ### bayesplot::mcmc_rank_overlay() now exist:
      # print(
      #   mcmc_hist_r_scale(C_draws_arr[, , selPar_i])
      # )
      ###
    }

    ### Not used, but defined in the file sourced above:
    # plot_ranknorm()
    # plot_rhat()
    # plot_ess()
    ###
  }

  ## Output -----------------------------------------------------------------

  out_list <- c(out_list,
                list("all_OK" = C_all_OK))

  return(out_list)
}
