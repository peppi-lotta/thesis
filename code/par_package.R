######################################################################
#' This function extracts values from the specified exposure and outcome
#' columns in the provided dataset. This is done by counting the number of
#' each possible combination of exposure and outcome
#'
#' @param data Dataset.
#' @param exposure_col String, name of the exposure column in the dataset.
#' @param outcome_col String , name of the outcome column in the dataset.
#'
#' @return List containing:
#'   x_1e1d: The count of rows where both exposure and outcome are 1.
#'   x_1e0d: The count of rows where exposure is 1 and outcome is 0.
#'   x_0e1d: The count of rows where exposure is 0 and outcome is 1.
#'   x_0e0d: The count of rows where both exposure and outcome are 0.
#'
#' @examples
#' # Example usage:
#' data <- data.frame(exposure = c(0, 0, 1, 1, 0, 1, 0, 1),
#'                    outcome = c(0, 1, 0, 1, 0, 1, 1, 0))
#' result <- extract_abcd(data, "exposure", "outcome")
#'
#' @export
extract_abcd <- function(data, exposure_col, outcome_col) {
  x_0e0d <- sum(data[[exposure_col]] == 0 & data[[outcome_col]] == 0)
  x_0e1d <- sum(data[[exposure_col]] == 0 & data[[outcome_col]] == 1)
  x_1e0d <- sum(data[[exposure_col]] == 1 & data[[outcome_col]] == 0)
  x_1e1d <- sum(data[[exposure_col]] == 1 & data[[outcome_col]] == 1)

  return(c(
    x_1e1d = x_1e1d,
    x_1e0d = x_1e0d,
    x_0e1d = x_0e1d,
    x_0e0d = x_0e0d
  ))
}

######################################################################
######################################################################
#'
#' This function calculates the population attributable risk (PAR)
#' for a binary exposure and outcome.
#' PAR = P(D+) - P(D+| E-) = (a + c)/(a + b + c + d) - (c + d)/c
#' @param x A vector containing the values of a, b, c, and d
#' in this order. Where
#' a: The count of rows where both exposure and outcome are 1.
#' b: The count of rows where exposure is 1 and outcome is 0.
#' c: The count of rows where exposure is 0 and outcome is 1.
#' d: The count of rows where both exposure and outcome are 0.
#'
#' @return The population attributable risk (PAR).
#'
#' @examples
#' # Example usage:
#' x <- c(1, 2, 3, 4)
#' PAR <- calculate_par(x)
#' @export
calculate_par <- function(x) {
  x <- as.numeric(x)
  a <- x[1]
  b <- x[2]
  c <- x[3]
  d <- x[4]

  if (c + d == 0 || a + b + c + d == 0) {
    return(0)
  }
  par <- (a + c) / (a + b + c + d) - c / (c + d)
  return(par)
}

######################################################################
######################################################################
calculate_paf <- function(x) {
  a <- as.numeric(x[1])
  b <- as.numeric(x[2])
  c <- as.numeric(x[3])
  d <- as.numeric(x[4])

  pos_d_count <- (a + c) / (a + b + c + d)

  if (c + d == 0 || a + b + c + d == 0) {
    return(0)
  }
  paf <- ( pos_d_count - c / (c + d) ) / pos_d_count
  return(paf)
}

######################################################################
######################################################################
#' This function calculates the population attributable risk's (PAR) credibility
#' interval using a Bayesian approach.
#' The function samples from the Dirichlet distribution to estimate the
#' posterior distribution.
#'
#' @importFrom MCMCpack rdirichlet
#'
#' @param type String, type of the measure to calculate. Possible values are
#' "par" and "paf".
#' @param x A vector containing the values of a, b, c, and d
#' in this order. Where
#' a: The count of rows where both exposure and outcome are 1.
#' b: The count of rows where exposure is 1 and outcome is 0.
#' c: The count of rows where exposure is 0 and outcome is 1.
#' d: The count of rows where both exposure and outcome are 0.
#' @param interval Interval for the confidence interval (default is 0.95).
#' Possible values are between 0 and 1.
#' @param prior A list of prior values for the Dirichlet distribution.
#' Default is c(1, 1, 1, 1).
#' List position coresponds to the following:
#'  1. Prior for a
#'  2. Prior for b
#'  3. Prior for c
#'  4. Prior for d
#' @param sample_count Number of samples to draw from the Dirichlet
#' distribution.
#'
#' @return Matrix containing the lower and upper bounds of the credibility
#' interval. The first column contains the lower bound and the second column
#' contains the upper bound.
#'
#' @examples
#' # Example usage:
#' x <- extract_abcd(data, "exposure_col_name", "outcome_col_name")
#' calculate_bayesian_ci("par", x, 0.99, c(1, 1, 0.001, 0.001), 5000)
#'
#' @export
calculate_bayesian_ci <- function(
  type,
  x,
  interval = 0.95,
  prior = c(1, 1, 1, 1),
  sample_count = 10000
) {
  x <- as.numeric(x)
  a <- x[1]
  b <- x[2]
  c <- x[3]
  d <- x[4]
  n <- a + b + c + d
  prior <- as.numeric(prior)

  # Sample from the Dirichlet distribution 10,000 times
  # + prior is to add the prior distribution's effect
  samples <- MCMCpack::rdirichlet(
    sample_count,
    c(a + prior[1],
      b + prior[2],
      c + prior[3],
      d + prior[4],
      n
    )
  )
  samples <- apply(samples, 2, function(x) x * n)

  if (type == "par") {
    par_samples <- apply(samples, 1, calculate_par)
  } else if (type == "paf") {
    par_samples <- apply(samples, 1, calculate_paf)
  } else {
    stop("Invalid type. Please use 'par' or 'paf'")
  }

  # Calculate the confidence interval
  confidence_interval <- quantile(
    par_samples,
    c(
      (1 - interval) / 2,
      1 - (1 - interval) / 2
    )
  )
  return(matrix(c(
    confidence_interval[1],
    confidence_interval[2]
  )))
}

######################################################################
######################################################################
#' This function calculates the population attributable risk's (PAR) credibility
#' interval using a Bootstrap
#' The function samples from the Multinomial distribution to estimate the
#' posterior distribution.
#'
#' @param type String, type of the measure to calculate. Possible values are
#' "par" and "paf".
#' @param x A vector containing the values of a, b, c, and d
#' in this order. Where
#' a: The count of rows where both exposure and outcome are 1.
#' b: The count of rows where exposure is 1 and outcome is 0.
#' c: The count of rows where exposure is 0 and outcome is 1.
#' d: The count of rows where both exposure and outcome are 0.
#' @param interval Interval for the confidence interval (default is 0.95).
#' Possible values are between 0 and 1.
#' @param sample_count Number of samples to draw from the Dirichlet
#' distribution.
#'
#' @return Matrix containing the lower and upper bounds of the credibility
#' interval. The first column contains the lower bound and the second column
#' contains the upper bound.
#'
#' @examples
#' # Example usage:
#' x <- extract_abcd(data, "exposure", "outcome")
#' calculate_bootstrap_ci("par", x, 0.99, 5000)
#'
#' @export
calculate_bootstrap_ci <- function(
  type,
  x,
  interval = 0.95,
  sample_count = 10000
) {
  x <- as.numeric(x)
  a <- x[1]
  b <- x[2]
  c <- x[3]
  d <- x[4]
  n <- a + b + c + d

  p_11 <- a / n
  p_10 <- b / n
  p_01 <- c / n
  p_00 <- d / n

  # Initialize samples dataframe
  samples <- t(rmultinom(sample_count, n, c(p_11, p_10, p_01, p_00)))

  if (type == "par") {
    par_samples <- apply(samples, 1, calculate_par)
  } else if (type == "paf") {
    par_samples <- apply(samples, 1, calculate_paf)
  } else {
    stop("Invalid type. Please use 'par' or 'paf'")
  }

  # Calculate the confidence interval
  confidence_interval <- quantile(
    par_samples,
    c(
      (1 - interval) / 2,
      1 - (1 - interval) / 2
    )
  )

  return(matrix(c(confidence_interval[1], confidence_interval[2])))
}

######################################################################
######################################################################
#' Compile all functions
#' This function compiles all functions in the package. This is useful for
#' optimizing the performance of the functions.
#' @importFrom compiler cmpfun
#' @export
compile_all <- function() {
  calculate_bayesian_ci <- compiler::cmpfun(calculate_bayesian_ci)
  calculate_bootstrap_ci <- compiler::cmpfun(calculate_bootstrap_ci)
  calculate_par <- compiler::cmpfun(calculate_par)
  calculate_paf <- compiler::cmpfun(calculate_paf)
  extract_abcd <- compiler::cmpfun(extract_abcd)
}
