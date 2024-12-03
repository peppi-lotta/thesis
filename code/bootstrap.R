source("./contingency_table.R")
source("./par_package.R")

bootstrap <- function(
  data,
  exposure_col,
  outcome_col,
  interval = 0.95
) {
  ct <- create_2x2_contingency_table(data, exposure_col,outcome_col)
  sample_count <- 10000
  p_11 <- ct$x_1e1d / ct$n
  p_10 <- ct$x_1e0d / ct$n
  p_01 <- ct$x_0e1d / ct$n
  p_00 <- ct$x_0e0d / ct$n
  n <- ct$n

  # Calculate the population attributable risk
  bootstrap_par <- calculate_par(
    c(ct$x_1e1d, ct$x_1e0d, ct$x_0e1d, ct$x_0e0d, ct$n)
  )

  # Initialize samples dataframe
  samples <- data.frame()

  for (i in 1:sample_count) {
    p <- t(rmultinom(1, n, c(p_11, p_10, p_01, p_00)))
    samples[i, "x_e1d1"] <- p[1, 1]
    samples[i, "x_e1d0"] <- p[1, 2]
    samples[i, "x_e0d1"] <- p[1, 3]
    samples[i, "x_e0d0"] <- p[1, 4]
  }

  par_samples <- apply(samples, 1, calculate_par)

  # Calculate the confidence interval
  confidence_interval <- quantile(
    par_samples,
    c(
      (1 - interval) / 2,
      1 - (1 - interval) / 2
    )
  )

  return(
    list(
      contingency_table = ct$contingency_table,
      bootstrap_par = bootstrap_par,
      confidence_interval = confidence_interval,
      par_samples = par_samples
    )
  )
}