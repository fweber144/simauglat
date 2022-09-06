extreme_coef_prespec <- 3
exp(extreme_coef_prespec)
## --> Ca. 20.
df_slab <- 15
scale_slab <- 1
( extreme_quantile <- pt(extreme_coef_prespec / scale_slab, df = df_slab) )
## --> Ca. 0.995, so 99% between this value and its negative counterpart.
( extreme_coef <- scale_slab * qt(0.995, df = df_slab) )
exp(extreme_coef)
exp(-extreme_coef)
