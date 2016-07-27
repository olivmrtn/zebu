#' Estimate marginal and multivariate probabilities
#'
#' Maximum-likelihood estimation of marginal and multivariate observed and expected independence probabilities. Marginal probability refers to probability of each factor per individual column. Multivariate probability refer to cross-classifying factors for all columns.
#'
#' @param x data.frame or matrix.
#'
#' @return List containing the following values:
#'   \itemize{
#'   \item margins: a list of marginal probabilities. Names correspond to colnames(x).
#'   \item observed: observed multivariate probability array.
#'   \item expected: expected multivariate probability array
#'   }
#'
#' @example /inst/examples/lassie.R
#'
#' @export
estimate_prob <- function(x) {

  # Check if matrix or data frame
  if (!is.data.frame(x)) {
    x <- tryCatch(data.frame(x), error = function(c) {
      stop("Invalid 'x' argument: must be a 'data.frame' or 'matrix'")
    })
  }
  # Check if there are at least two variables
  if (ncol(x) < 2) {
    stop("Invalid 'x' argument: needs to have at least two columns.")
  }

  cnames <- colnames(x)
  prob <- list()

  # Compute univariate observed probabilities
  prob[["margins"]] <- lapply(cnames, function(i) {
    v <- x[, i]
    dim_names <- list()
    lev <- levels(v)
    if (is.null(lev))
      lev <- unique(v)
    dim_names[[i]] <- lev
    array(prop.table(table(v, useNA = "ifany")), dimnames = dim_names)
  })
  names(prob[["margins"]]) <- cnames

  # Compute multivariate observed probabilities
  observed <- as.array(prop.table(table(x, useNA = "ifany")))
  prob[["observed"]] <- observed

  # Compute expected probabilities if events were independent from each other
  expected <- prob[["margins"]][[1]]
  for (i in cnames[2:length(cnames)]) {
    expected <- outer(expected, prob[["margins"]][[i]])
  }
  expected <- array(expected, dim = dim(observed), dimnames = dimnames(observed))
  prob[["expected"]] <- expected

  prob
}
