# Derivative of the probability density function of the standard normal
# distribution (i.e., the second derivative of the cumulative distribution
# function of the standard normal distribution):
dnorm_deriv <- function(x) {
  x_out <- rep(NaN, length(x))
  x_is_fin <- is.finite(x)
  x_out[!x_is_fin] <- 0
  x_out[x_is_fin] <- -x[x_is_fin] * dnorm(x[x_is_fin])
  return(x_out)
}

# Derivative of the probability density function of the standard logistic
# distribution (i.e., the second derivative of the cumulative distribution
# function of the standard logistic distribution):
dlogis_deriv <- function(x) {
  x_out <- rep(NaN, length(x))
  x_is_fin <- is.finite(x)
  x_out[!x_is_fin] <- 0
  exp_minus_x <- exp(-x[x_is_fin])
  # Note: For better numerical accuracy, we could use `expm1(-x[x_is_fin])`
  # instead of `(exp_minus_x - 1)` in the following:
  x_out[x_is_fin] <- exp_minus_x * (exp_minus_x - 1) / (1 + exp_minus_x)^3
  return(x_out)
}

# Function for calculating the (approximate) variance of the Gaussian
# pseudo-observations (`\tilde{\sigma}^2` from @piironen_sparsity_2017) in a
# "cumulative" family for a given response value `y \in {1, ..., C}` with `C`
# denoting the number of response categories:
calc_sigma2_tilde_y <- function(y,
                                ilink_fun = pnorm,
                                ilink_deriv_fun = dnorm,
                                ilink_deriv2_fun = dnorm_deriv,
                                linpred = 0,
                                thres_vec = NULL,
                                ncats = 5) {
  nthres <- ncats - 1L
  if (!is.null(thres_vec)) {
    stopifnot(length(thres_vec) == nthres)
  } else {
    thres_vec <- qnorm(seq_len(nthres) / ncats)
  }
  thres_vec <- c(-Inf, thres_vec, Inf)
  y <- y + 1L

  diff_orig <- ilink_fun(thres_vec[y] - linpred) -
    ilink_fun(thres_vec[y - 1L] - linpred)
  diff_deriv <- ilink_deriv_fun(thres_vec[y] - linpred) -
    ilink_deriv_fun(thres_vec[y - 1L] - linpred)
  diff_deriv2 <- ilink_deriv2_fun(thres_vec[y] - linpred) -
    ilink_deriv2_fun(thres_vec[y - 1L] - linpred)

  L_deriv2 <- (diff_deriv2 * diff_orig - diff_deriv^2) / diff_orig^2

  return(-1 / L_deriv2)
}

# Function for aggregating the (approximate) variances of the Gaussian
# pseudo-observations across all response values `y \in {1, ..., C}` with `C`
# denoting the number of response categories:
calc_sigma_tilde <- function(agg_type = "geometric_mean", ncats = 5, ...) {
  sigma2s_tilde <- do.call(c, lapply(
    seq_len(ncats),
    calc_sigma2_tilde_y,
    ncats = ncats,
    ...
  ))
  if (agg_type == "simple_mean") {
    # Root of the mean of the variances:
    return(sqrt(mean(sigma2s_tilde)))
  } else if (agg_type == "harmonic_mean") {
    # Root of the harmonic mean of the variances:
    return(sqrt(1 / mean(1 / sigma2s_tilde)))
  } else if (agg_type == "geometric_mean") {
    # Root of the geometric mean of the variances
    # (= exp of the mean of the log standard deviations):
    return(sqrt(exp(mean(log(sigma2s_tilde)))))
  } else {
    stop("Unexpected `agg_type`.")
  }
}
