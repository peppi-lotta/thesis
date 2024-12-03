library(MCMCpack)
source("./contingency_table.R")
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
#' @param n The total number of observations.
#'
#' @return The population attributable risk (PAR).
#'
#' @examples
#' # Example usage:
#' PAR <- calculate_PAR(1, 2, 3, 4, 10)
#'
#' @export
calculate_par <- function(x) {
  a <- x[1]
  b <- x[2]
  c <- x[3]
  d <- x[4]
  par <- (a + c) / (a + b + c + d) - c / (c + d)
  return(par)
}
######################################################################
#' This function print the Population Attributable Risk (PAR)
#' and its specified confidence interval
#' based on a 2x2 contingency table created from the provided data.
#'
#' @param data dataset.
#' @param exposure_col String, name of the column representing
#' the exposure variable.
#' @param outcome_col String, name of the column representing
#' the outcome variable.
#' @param interval Interval for the confidence interval (default is 0.95).
#' Possible values are between 0 and 1.
#' @param prior A list of prior values for the Dirichlet distribution.
#' Default is c(1, 1, 1, 1).
#' List position coresponds to the following:
#'   1. Prior for 0e0d
#'   2. Prior for 0e1d
#'   3. Prior for 1e0d
#'   4. Prior for 1e1d
#'
#' @return List containing:
#'  par: The population attributable risk.
#'  confidence_interval: The confidence interval for the population attributable
#' risk.
#'
#' par_samples: The population attributable risk samples from the Dirichlet
#' distribution.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(exposure = sample(c(0, 1), 100, replace = TRUE),
#'                    outcome = sample(c(0, 1), 100, replace = TRUE))
#' calculate_par_and_ci(data, "exposure", "outcome", 0.5)
#' }
calculate_par_and_ci_single <- function(
  data,
  exposure_col,
  outcome_col,
  interval = 0.95,
  prior = list(1, 1, 1, 1)
) {
  # Create a 2x2 contingency table
  ct <- create_2x2_contingency_table(data, exposure_col, outcome_col)

  # Calculate the population attributable risk
  bayesian_par <- calculate_par(
    c(ct$x_1e1d, ct$x_1e0d, ct$x_0e1d, ct$x_0e0d, ct$n)
  )

  # Sample from the Dirichlet distribution 10,000 times
  # + prior is to add the prior distribution's effect
  samples <- rdirichlet(
    10000,
    c(
      ct$x_1e1d + unlist(prior)[1],
      ct$x_1e0d + unlist(prior)[2],
      ct$x_0e1d + unlist(prior)[3],
      ct$x_0e0d + unlist(prior)[4],
      ct$n
    )
  )
  samples <- apply(samples, 2, function(x) x * ct$n)

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
      bayesian_par = bayesian_par,
      confidence_interval = confidence_interval,
      par_samples = par_samples
    )
  )
}
