#' Permutation test for local and global association measures
#'
#' Permutation test: statistical significance of local and
#' global association measures
#'
#' @param x \code{\link[zebu]{lassie}} S3 object.
#'
#' @param nb number of resampling iterations.
#'
#' @param group list of column names specifying which columns
#' should be permuted together. This is useful for the multivariate case,
#' for example, when there is many dependent variables and one
#' independent variable. By default, permutes all columns separately.
#'
#' @param p_adjust multiple testing correction method.
#' (see \code{\link[stats]{p.adjust.methods}} for a list of methods).
#'
#' @return \code{permtest} returns an S3 object of \link[base]{class}
#' \code{\link[zebu]{lassie}} and \code{\link[zebu]{permtest}}.
#' Adds the following to the lassie object \code{x}:
#' \itemize{
#' \item global_p: global association p-value.
#' \item local_p: array of local association p-values.
#' \item global_perm: numeric global association values obtained with permutations.
#' \item local_perm: matrix local association values obtained with permutations. Column number correspond to positions in local association array after converting to numeric (e.g. local_perm[, 1] corresponds to local[1]).
#' \item perm_params: parameters used when calling permtest (nb and p_adjust).
#' }
#'
#' @seealso \code{\link[zebu]{lassie}}
#'
#' @examples
#'
#' # Calling lassie on cars dataset
#' las <- lassie(cars, continuous = colnames(cars))
#'
#' # Permutation test using default settings
#' permtest(las, nb = 30) # keep resampling low for example
#'
#' @export
#' @import data.table
#'
permtest <- function(x,
                     nb = 1000L,
                     group = as.list(colnames(x$data$pp)),
                     p_adjust = "BH") {

  # Estimate p-value
  p_value = function(observed_value, permuted_values) {
    sum(abs(observed_value) < abs(permuted_values)) / length(permuted_values)
  }

  # Check arguments ----
  if (! "lassie" %in% class(x)) {
    stop("Invalid 'x' argument: must be a 'lassie' object")
  }

  if (! is.numeric(nb) || nb <= 1) {
    stop("Invalid 'nb' argument: must a numeric")
  }

  if (! p_adjust %in% stats::p.adjust.methods) {
    stop(paste("Invalid 'p_adjust' argument: methods supported are", stats::p.adjust.methods))
  }

  if (! is.list(group)) {
    stop("Invalid 'group' argument: must be a list of characters that correspond to colnames.")
  }

  # Compute general variables ----
  dimensions <- dim(x$prob$observed)  # Dimensions
  dim_names <- dimnames(x$prob$observed)  # Dimension names
  measure <- x$lassie_params[["measure"]]  # Get association measure used
  perm_params <- list(nb = nb, p_adjust = p_adjust) # Save parameters in list

  # Convert group argument integers to characters refering to colnames
  if (any(! unlist(group) %in% colnames(x$data$pp))) {
    stop("Invalid 'group' argument: must be a list of characters that correspond colnames.")
  }

  # Resampling ----
  permutations <- permtest_rcpp(x, nb, group)

  # Compute p-values ----
  # Compute global p-value
  global_perm <- sapply(permutations, `[[`, "global")
  global_p <- p_value(x$global, global_perm)

  # Format local measure results to matrix: rows are local measure cells and
  # columns are resampling iterations
  local_perm <- sapply(permutations, `[[`, "local")

  # browser()

  # Compute local p-value for each cell
  local_p <- sapply(seq_len(length(x$local)), function(i) {

    # Retrieve observed value and permutation vector
    local_i <- x$local[i]
    local_perm_i <- local_perm[i, ]

    # If observed value is non-numerical (NA, NaN of Inf), return NA for p-value
    if (is.infinite(local_i) | is.na(local_i) | is.nan(local_i)) {
      return(NA_real_)
    }

    # Remove NA and NaNs values from permutations
    local_perm_i <- local_perm_i[! (is.na(local_perm_i) | is.nan(local_perm_i))]

    # Replace Inf values by the largest integer which can be represented
    inf <- is.infinite(local_perm_i)
    local_perm_i[inf] <- sign(local_perm_i[inf]) * .Machine$integer.max

    # Estimate p-value by intersection
    p_value(local_i, local_perm_i)
  })

  # Multiple test correction
  local_p <- stats::p.adjust(local_p, method = p_adjust)

  # Format to array
  local_p <- array(local_p, dim = dimensions, dimnames = dim_names)

  x[["global_p"]] <- global_p
  x[["local_p"]] <- local_p
  x[["global_perm"]] <- global_perm
  x[["local_perm"]] <- local_perm
  x[["perm_params"]] <- perm_params
  structure(x, class = c(class(x), "permtest"))
}
