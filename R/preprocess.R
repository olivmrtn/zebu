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
preprocess <- function(x, select, continuous, breaks, default_breaks = 4) {

  # Try to convert 'x' to data.frame
  x <- tryCatch(as.data.frame(x), error = function(c) {
    stop("Invalid 'x' argument: must be convertible to a data.frame")
  })

  # Save 'x' data.frame to 'raw'
  raw <- x

  # Handle missing 'select' argument
  if (missing(select) || is.null(select)) {
    select <- colnames(x)
  }

  # Handle missing 'breaks' arguments
  if (missing(breaks) || is.null(breaks)) {
    breaks <- default_breaks
  }

  # Subset 'x' data.frame
  x <- tryCatch(subset(raw, select = select), error = function(c) {
    stop("Invalid 'select' argument: needs to be a vector of column numbers or column names")
  })

  # If 'select' argument is specified as column numbers, convert to variable names
  select <- colnames(x)

  # Handle missing 'continuous' argument
  if (missing(continuous) || is.null(continuous)) {
    continuous <- NULL
  }

  # If 'continuous' argument is specified as column numbers, convert to variable names
  continuous <- tryCatch(colnames(subset(raw, select = continuous)), error = function(c) {
    stop("Invalid 'continuous' argument: needs to be a vector of column numbers or column names")
  })

  # Check if 'continuous' variables not in 'subset'
  if (any(! continuous %in% select)) {
    stop("Invalid 'continous' argument: contains variables names not in 'select' argument")
  }

  # Save 'x' data.frame to 'raw' so that raw is subsetted
  raw <- x

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

  # Discretize continuous variables according to breaks
  x <- as.data.frame(lapply(select, function(i) {
    v <- x[, i, drop = TRUE]

    # Discretize if continuous
    if (i %in% continuous) {
      b <- breaks[[i]]
      if (is.null(b) || is.na(b)) {
        b <- default_breaks
      }
      v <- cut(as.numeric(v), breaks = b, include.lowest = TRUE) }

    # Convert to character if categorical
    else {
      v <- as.character(v)
    }
    # Check if discretizations worked
    if (is.null(v)) {
      stop(paste("Could not discretize", i))
    }
    # Return factor
    factor(v)
  }))
  colnames(x) <- select

  # Remove missing values
  x <- stats::na.omit(x)

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
