#' Chi-squared test
#'
#' Chi-squared test: statistical significance of (global) chi-squared statistic
#' and (local) chi-squared residuals
#'
#' @param x \code{\link[zebu]{lassie}} S3 object.
#'
#' @param p_adjust multiple testing correction method.
#' (see \code{\link[stats]{p.adjust.methods}} for a list of methods).
#'
#' @return \code{chisqtest} returns an S3 object of \link[base]{class}
#' \code{\link[zebu]{lassie}} and \code{\link[zebu]{chisqtest}}.
#' Adds the following to the lassie object \code{x}:
#' \itemize{
#' \item global_p: global association p-value.
#' \item local_p: array of local association p-values.
#' }
#'
#' @seealso \code{\link[zebu]{lassie}}
#'
#' @examples
#'
#' # Calling lassie on cars dataset
#' las <- lassie(cars, continuous = colnames(cars), measure = "chisq")
#'
#' # Permutation test using default settings
#' chisqtest(las)
#'
#' @export
#' @import data.table
#'
chisqtest <- function(x, p_adjust = "BH") {

  # Check arguments ----
  if (! "lassie" %in% class(x)) {
    stop("Invalid 'x' argument: must be a 'lassie' object")
  }

  if (! p_adjust %in% stats::p.adjust.methods) {
    stop(paste("Invalid 'p_adjust' argument: methods supported are", stats::p.adjust.methods))
  }

  # Compute general variables ----
  dimensions <- dim(x$prob$observed)  # Dimensions
  dim_names <- dimnames(x$prob$observed)  # Dimension names
  measure <- x$lassie_params[["measure"]]  # Get association measure used

  if (measure != "chisq" | length(dimensions) != 2) {
    stop(paste("This function is only valid for the chi-squared method and two dimensions."))
  }
  # Compute p-values ----
  global_p <- stats::pchisq(x$global, df = prod(dimensions - 1), lower.tail = FALSE)
  local_p <- 2 * (1 - stats::pnorm(abs(x$local)))

  # Multiple test correction
  local_p <- stats::p.adjust(local_p, method = p_adjust)

  # Format to array
  local_p <- array(local_p, dim = dimensions, dimnames = dim_names)

  x[["global_p"]] <- global_p
  x[["local_p"]] <- local_p
  structure(x, class = c(class(x), "chisqtest"))
}
