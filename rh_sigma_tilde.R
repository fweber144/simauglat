# Derivative of the probability density function of the standard normal
# distribution (i.e., the second derivative of the cumulative distribution
# function of the standard normal distribution):
dnorm_deriv <- function(x) {
  x_out <- rep(NaN, length(x))
  x_is_fin <- is.finite(x)
  x_out[!x_is_fin] <- 0
  x_out[x_is_fin] <- -x[x_is_fin] / sqrt(2 * pi) * exp(-0.5 * x[x_is_fin]^2)
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

print(sapply(seq_len(2), calc_sigma2_tilde_y, ncats = 2))
print(sapply(seq_len(3), calc_sigma2_tilde_y, ncats = 3))
print(sapply(seq_len(4), calc_sigma2_tilde_y, ncats = 4))
print(sapply(seq_len(5), calc_sigma2_tilde_y, ncats = 5))

# Function for aggregating the (approximate) variances of the Gaussian
# pseudo-observations across all response values `y \in {1, ..., C}` with `C`
# denoting the number of response categories:
calc_sigma_tilde <- function(agg_type = "type3", ncats = 5, ...) {
  sigma2s_tilde <- do.call(c, lapply(
    seq_len(ncats),
    calc_sigma2_tilde_y,
    ncats = ncats,
    ...
  ))
  if (agg_type == "type1") {
    # Root of the mean of the variances:
    return(sqrt(mean(sigma2s_tilde)))
  } else if (agg_type == "type2") {
    # Root of the harmonic mean of the variances:
    return(sqrt(1 / mean(1 / sigma2s_tilde)))
  } else if (agg_type == "type3") {
    # Root of the geometric mean of the variances
    # (= exp of the mean of the log standard deviations):
    return(sqrt(exp(mean(log(sigma2s_tilde)))))
  } else {
    stop("Unexpected `agg_type`.")
  }
}

print(calc_sigma_tilde(agg_type = "type1"))
print(calc_sigma_tilde(agg_type = "type2"))
print(calc_sigma_tilde(agg_type = "type3"))

print(round(calc_sigma_tilde(ncats = 2), 3))
## --> Gives: 1.253
print(round(calc_sigma_tilde(ncats = 3), 3))
## --> Gives: 1.127
print(round(calc_sigma_tilde(ncats = 4), 3))
## --> Gives: 1.082
print(round(calc_sigma_tilde(ncats = 5), 3))
## --> Gives: 1.059

# For checking, use the Bernoulli family with the logit link which has a value
# of `\tilde{\sigma}^2 = 4` in @piironen_sparsity_2017:
brnll_logit_ch <- sapply(seq_len(2),
                         calc_sigma2_tilde_y,
                         ilink_fun = plogis,
                         ilink_deriv_fun = dlogis,
                         ilink_deriv2_fun = dlogis_deriv,
                         ncats = 2)
print(brnll_logit_ch)
stopifnot(all.equal(brnll_logit_ch, rep(4, 2), tolerance = .Machine$double.eps))
