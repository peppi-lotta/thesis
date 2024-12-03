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
