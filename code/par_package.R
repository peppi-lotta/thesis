library(MCMCpack)

######################################################################
#' This function creates a 2x2 contingency table from
#' the specified exposure and outcome columns in the provided dataset.
#'
#' @param data Dataset.
#' @param exposure_col String, name of the exposure column in the dataset.
#' @param outcome_col String , name of the outcome column in the dataset.
#'
#' @return List containing:
#'   x_0e0d: The count of rows where both exposure and outcome are 0.
#'   x_0e1d: The count of rows where exposure is 0 and outcome is 1.
#'   x_1e0d: The count of rows where exposure is 1 and outcome is 0.
#'   x_1e1d: The count of rows where both exposure and outcome are 1.
#'   n: The total number of observations.
#'   contingency_table: The 2x2 contingency table with row and column totals.
#' }
#'
#' @examples
#' # Example usage:
#' data <- data.frame(exposure = c(0, 0, 1, 1, 0, 1, 0, 1),
#'                    outcome = c(0, 1, 0, 1, 0, 1, 1, 0))
#' result <- create_2x2_contingency_table(data, "exposure", "outcome")
#'
#' @export
create_2x2_contingency_table <- function(data, exposure_col, outcome_col) {
# Create a 2x2 table of chosen columns with dimension names
  table <- table(
    data[[exposure_col]],
    data[[outcome_col]],
    dnn = c(
      exposure_col,
      outcome_col
    )
)

# Count the number of each possible combination of exposure and outcome
  x_0e0d <- sum(data[[exposure_col]] == 0 & data[[outcome_col]] == 0)
  x_0e1d <- sum(data[[exposure_col]] == 0 & data[[outcome_col]] == 1)
  x_1e0d <- sum(data[[exposure_col]] == 1 & data[[outcome_col]] == 0)
  x_1e1d <- sum(data[[exposure_col]] == 1 & data[[outcome_col]] == 1)
  n    <- sum(table)

# Add row and column totals to the contingency table
  contingency_table <- addmargins(table)

  return(
    list(
        x_0e0d = x_0e0d,
        x_0e1d = x_0e1d,
        x_1e0d = x_1e0d,
        x_1e1d = x_1e1d,
        n = n, 
        contingency_table = contingency_table
        )
    )
}
######################################################################
######################################################################
#'
#' This function calculates the population attributable risk (PAR)
#' for a binary exposure and outcome.
#' PAR = P(D+) - P(D+| E-) = (a + c)/(a + b + c + d) - (c + d)/c
#'
#' @param a The count of rows where both exposure and outcome are 1.
#' @param b The count of rows where exposure is 1 and outcome is 0.
#' @param c The count of rows where exposure is 0 and outcome is 1.
#' @param d The count of rows where both exposure and outcome are 0.
#'
#' @return The population attributable risk (PAR).
#'
#' @examples
#' # Example usage:
#' PAR <- calculate_PAR(1, 2, 3, 4)
#'
#' @export
calculate_par <- function(x) {
  a <- as.numeric(x[1])
  b <- as.numeric(x[2])
  c <- as.numeric(x[3])
  d <- as.numeric(x[4])

  if (c + d == 0 || a + b + c + d == 0) {
    return(0)
  }
  par <- (a + c) / (a + b + c + d) - c / (c + d)
  return(par)
}
######################################################################
######################################################################
#' This function calculates the population attributable risk's (PAR) credibility
#' interval using a Bayesian approach.
#' The function samples from the Dirichlet distribution to estimate the
#' posterior distribution.
#'
#' @param ct A 2x2 contingency table.
#' @param interval Interval for the confidence interval (default is 0.95).
#' Possible values are between 0 and 1.
#' @param prior A list of prior values for the Dirichlet distribution.
#' Default is c(1, 1, 1, 1).
#' List position coresponds to the following:
#'  1. Prior for 0e0d
#'  2. Prior for 0e1d
#'  3. Prior for 1e0d
#'  4. Prior for 1e1d
#'
#' @return List containing:
#' par_samples: The population attributable risk samples from the Dirichlet
#' distribution.
#' confidence_interval: The confidence interval for the population attributable
#' risk.
#'
#' @examples
#' # Example usage:
#' ct <- create_2x2_contingency_table(data, "exposure", "outcome")
#' calculate_bayesian_ci(ct$x_1e1d, ct$x_1e0d, ct$x_0e1d, ct$x_0e0d, ct$n, 0.95)
#'
#' @export
calculate_bayesian_ci <- function(
  x,
  interval = 0.95,
  prior = c(1, 1, 1, 1)
) {
  a <- as.numeric(x[1])
  b <- as.numeric(x[2])
  c <- as.numeric(x[3])
  d <- as.numeric(x[4])
  n <- as.numeric(x[5])

  # Sample from the Dirichlet distribution 10,000 times
  # + prior is to add the prior distribution's effect
  samples <- rdirichlet(
    10000,
    c(a + as.numeric(prior[1]),
      b + as.numeric(prior[2]),
      c + as.numeric(prior[3]),
      d + as.numeric(prior[4]),
      n
    )
  )
  samples <- apply(samples, 2, function(x) x * n)

  par_samples <- apply(samples, 1, calculate_par)

  # Calculate the confidence interval
  confidence_interval <- quantile(
    par_samples,
    c(
      (1 - interval) / 2,
      1 - (1 - interval) / 2
    )
  )
  return(list(
    par_samples = par_samples,
    confidence_interval = confidence_interval
  ))
}
######################################################################
######################################################################
calculate_bootstrap_ci <- function(
  x,
  interval = 0.95
) {
  a <- as.numeric(x[1])
  b <- as.numeric(x[2])
  c <- as.numeric(x[3])
  d <- as.numeric(x[4])
  n <- as.numeric(x[5])

  sample_count <- 10000
  p_11 <- a / n
  p_10 <- b / n
  p_01 <- c / n
  p_00 <- d / n

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
      confidence_interval = confidence_interval,
      par_samples = par_samples
    )
  )
}
