#' @title Local Association Measures
#'
#' @description Subroutines called by \code{\link[zebu]{lassie}} to compute
#' local and global association measures from a list of probabilities.
#'
#' @param x list of probabilities as outputted by \code{\link[zebu]{estimate_prob}}.
#' @param measure name of measure to be used:
#' \itemize{
#' \item'z': Ducher's 'z'.
#' \item'pmi': Pointwise mutual information (in bits).
#' \item'npmi': Normalized pointwise mutual information.
#' }
#'
#' @return List containing the following values:
#' \itemize{
#' \item local: local association array (may contain NA, NaN and Inf values).
#' \item global: global association numeric value.
#' }
#'
#' @details
#' \itemize{
#' \item \code{local_association(x, measure = 'z')} is equivalent to
#' \code{duchers_z(x)}.
#' \item \code{local_association(x, measure = 'pmi')} is equivalent to
#' \code{pmi(x)}.
#' \item \code{local_association(x, measure = 'npmi')} is equivalent to
#' \code{npmi(x)} and \code{pmi(x, normalize = TRUE)}.
#' }
#'
#' @seealso \code{\link[zebu]{lassie}}
#'
#' @example /inst/examples/lassie.R
#'
#' @name local_association
#'
#' @export
#'
local_association <- function(x,
                              measure) {
  if (measure == "z") {
    duchers_z(x)
  } else if (measure == "pmi") {
      pmi(x, FALSE)
  } else if (measure == "npmi") {
      pmi(x, TRUE)
  } else {
    measure_error()
  }
}

#' @rdname local_association
#' @export
#'
duchers_z <- function(x) {

  expected <- x$expected
  observed <- x$observed
  margins <- expand.grid(x$margins)

  # Compute minimal margins probabilities
  minimum <- array(apply(margins, 1, min), dim = dim(expected), dimnames = dimnames(expected))

  # Compute global and local Ducher's Z
  local <- array((observed - expected), dim = dim(expected), dimnames = dimnames(expected))
  pos <- which(local > 0, arr.ind = TRUE)
  neg <- which(local < 0, arr.ind = TRUE)
  local[pos] <- local[pos]/(minimum - expected)[pos]
  local[neg] <- local[neg]/expected[neg]
  global <- sum(observed * local)

  lam <- list()
  lam[["local"]] <- local
  lam[["global"]] <- global
  lam
}

#' @rdname local_association
#' @param normalize Normalizes pointwise mutual information when calling \code{pmi}
#' @export
pmi <- function(x, normalize = FALSE) {

  expected <- x$expected
  observed <- x$observed

  # Compute pointwise multi-information
  local <- log2(observed/expected)

  # Normalize
  if (normalize) {
    pos <- which(local > 0, arr.ind = TRUE)
    neg <- which(local < 0, arr.ind = TRUE)
    margins <- expand.grid(x$margins)
    prod_max <- apply(margins, 1, function(i) prod(sort(i)[2:length(i)]))
    prod_max <- array(prod_max, dim = dim(expected), dimnames = dimnames(expected))
    local[pos] <- local[pos]/-log2(prod_max[pos])
    local[neg] <- local[neg]/-log2(observed[neg])
    local[which(is.infinite(local) | is.nan(local), arr.ind = TRUE)] <- -1
  }

  # Global
  not_na_inf <- which(!is.infinite(local) & !is.na(local), arr.ind = TRUE)
  global <- sum(observed[not_na_inf] * local[not_na_inf])

  lam <- list()
  lam[["local"]] <- local
  lam[["global"]] <- global
  lam
}

#' @rdname local_association
#' @export
npmi <- function(x) {
  pmi(x, normalize = TRUE)
}
