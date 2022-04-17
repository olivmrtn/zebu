#' @title Preprocess data
#'
#' @description Subroutine called by \code{\link[zebu]{lassie}}. Discretizes, subsets and remove missing data from a data.frame.
#'
#' @param x data.frame or matrix.
#'
#' @param select optional vector of column numbers or column names
#' specifying a subset of data to be used. By default, uses all columns.
#'
#' @param continuous optional vector of column numbers or column names specifying
#' continuous variables that should be discretized.
#' By default, assumes that every variable is categorical.
#'
#' @param breaks numeric vector or list passed on to \code{\link[base]{cut}} to discretize
#' continuous variables. When a numeric vector is specified, break points are
#' applied to all continuous variables. In order to specify variable-specific
#' breaks, lists are used. List names identify variables and list values identify
#' breaks. List names are column names (not numbers). If a continuous
#' variable has no specified breaks, then \code{default_breaks} will be applied.
#'
#' @param default_breaks default break points for discretizations.
#' Same syntax as in \code{\link[base]{cut}}.
#'
#' @return List containing the following values:
#'   \itemize{
#'   \item raw: raw subsetted data.frame
#'   \item pp: discretized, subsetted and complete data.frame
#'   \item select
#'   \item continuous
#'   \item breaks
#'   \item default_breaks
#'   }
#'
#' @example /inst/examples/lassie.R
#'
#' @export
#'
preprocess <- function(x,
                       select,
                       continuous,
                       breaks,
                       default_breaks = 4) {

  # Try to convert 'x' to data.frame
  x <- tryCatch(as.data.table(x), error = function(c) {
    stop("Invalid 'x' argument: must be convertible to a data.frame")
  })

  # Be sure to have colnames
  if (is.null(colnames(x))) {
    colnames(x) <- paste("V", seq_len(ncol(x)))
  }

  # Handle missing arguments
  if (missing(select) || is.null(select)) {
    select <- colnames(x)
  }

  if (missing(continuous) || is.null(continuous)) {
    continuous <- character(0L)
  }

  # Handle missing 'breaks' arguments
  if (missing(breaks) || is.null(breaks)) {
    breaks <- default_breaks
  }

  # Convert numeric columns to characters
  if (all(is.numeric(select))) {
    select <- colnames(x)[select]
  }

  if (all(is.numeric(continuous))) {
    continuous <- colnames(x)[continuous]
  }

  # Check if 'continuous' variables not in 'subset'
  if (any(! continuous %in% select)) {
    stop("Invalid 'continous' argument: contains variables names not in 'select' argument")
  }

  # Categorical variables are variables that are not continuous
  categorical <- setdiff(select, continuous)

  # Breaks define as list
  if (is.list(breaks)) {
    if (is.null(names(breaks))) {
      stop("Invalid 'breaks' argument: breaks does not have any names. Must correspond to colnames in 'x'")

      # Check if breaks names are in continuous
    } else if (any(!names(breaks) %in% select)) {
      stop("Invalid 'breaks' argument: contains variable names not in 'continuous'")
    }

    # Breaks specified as numeric vector
  } else if (all((is.numeric(breaks) && all(! is.na(breaks))))) {
    breaks <- rep(list(breaks), length(continuous))
    names(breaks) <- continuous

  } else {
    stop("Invalid 'breaks' argument: must be numeric vector or list (see help for format)")
  }

  # Subset data
  remove_cols <- colnames(x)[! colnames(x) %in% select]
  if (length(remove_cols) > 0) {
    set(x, j = remove_cols, value = NULL)
  }

  # Save 'x' data.frame to 'raw'
  raw <- as.data.frame(x)

  # Discretize continuous variables according to breaks
  # browser()
  for (j in select) {

    v <- x[[j]]

    # Discretize if continuous
    if (j %in% continuous) {
      b <- breaks[[j]]
      if (any(is.null(b) | is.na(b))) {
        b <- default_breaks
      }
      if (! is.numeric(v)) {
        v <- as.numeric(v)
      }
      v <- cut(v, breaks = b, include.lowest = TRUE)
    }

    v <- as.factor(v)
    # Check if discretizations worked
    if (is.null(v)) {
      stop(paste("Could not discretize", j))
    }

    data.table::set(x, j = j, value = v)
  }

  # Remove NAs
  x <- stats::na.omit(x)
  x <- as.data.frame(x)

  # Check if rows are left
  if (nrow(x) == 0) {
    stop("Too much missing data. Preprocessed data.frame contains 0 rows.")
  }

  list(raw = raw,
       pp = x,
       select = select,
       continuous = continuous,
       breaks = breaks,
       default_breaks = default_breaks)
}
