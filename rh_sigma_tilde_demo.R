source("rh_sigma_tilde.R")

print(sapply(seq_len(2), calc_sigma2_tilde_y, ncats = 2))
print(sapply(seq_len(3), calc_sigma2_tilde_y, ncats = 3))
print(sapply(seq_len(4), calc_sigma2_tilde_y, ncats = 4))
print(sapply(seq_len(5), calc_sigma2_tilde_y, ncats = 5))

print(calc_sigma_tilde(agg_type = "simple_mean"))
print(calc_sigma_tilde(agg_type = "harmonic_mean"))
print(calc_sigma_tilde(agg_type = "geometric_mean"))

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
