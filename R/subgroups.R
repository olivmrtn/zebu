#' @title Local Association Subgroup Analysis
#'
#' @description Identifies if the local association between variables (named associated variables)
#' is dependent on the value of an another variable (named interacting variable).
#' Associated variables are specified by \code{las}. Interacting variable(s)
#' values are specified by \code{x}.
#’
#' @details Associated variables events are recoded into a subgroup variable according to local
#' association values (and eventually significance) into 'positive', 'negative' and 'independent'.
#' This is specified by the \code{thresholds}, \code{significance} and \code{alpha} arguments.
#' The local (and global) association between the new subgroup variable
#' and the interacting variable is then estimated using \code{\link[zebu]{lassie}}.
#'
#' @inheritParams lassie
#'
#' @param las \code{\link[zebu]{lassie}} S3 object. Corresponds to associated variables.
#' @param x data.frame or matrix. Corresponds to interacting variable(s) specified by \code{select}.
#' @param select optional vector of column numbers or column names specifying a subset of data to be used.
#' By default, uses all colnames in \code{x} except those in \code{las} object.
#' @param thresholds vector specifying respectively the negative and the positive
#' association threshold. Local association values between these thresholds are considered independent.
#' Values not contained in this range are classified as independent.
#' @param significance optional logical value specifying if only non-significant local association
#' values should be considered as independent.
#' Only available if \code{las} is also a \code{\link[zebu]{permtest}} object.
#' @param alpha alpha error level. Local association with p-values above this value are considered
#' as independent. Only available if \code{las} is also a \code{\link[zebu]{permtest}} object.
#'
#' @return An instance of S3 \link[base]{class} \code{\link[zebu]{lassie}}.
#'
#' @seealso Significance can be accessed using a permutation test: \code{\link[zebu]{permtest}}.
#'
#' @examples
#' # In this example, we will use the zebu 'trial' dataset.
#' # See vignette example for more detailed explanation
#'
#' # 'trial' corresponds to a simulated clinical trial where patient recovery
#' # is dependent on drug intake ('drug') and resistance status ('resistance').
#' # Patient recovery is monitored by a biomarker (continuous variable from 0 to 1)
#' # Patients with post-treatment biomarker ('postbiom') above 0.7 is have recovered.
#'
#' # Load 'trial' dataset
#' data(trial)
#'
#' # Compute the association between drug intake and patient recovery
#' las <- lassie(trial,
#'               select = c("drug", "postbiom"),
#'               continuous = c("postbiom"),
#'               breaks = c(0, 0.7, 1))
#'
#' # Permuation test
#' # Access significance of global and local association
#' las <- permtest(las)
#'
#' # Global association between drug intake and recovery but not for all patients
#' # Being in the drug group is locally independent of having not recovered
#' print(las)
#'
#' # Local association subgroup analysis
#' sub <- subgroups(las, trial, select = "resistance", alpha = 0.01)
#'
#' # Variable 'resistance' explains differences between sensitive and resistance patients
#' print(sub)
#'
#' @export
#'
subgroups <- function(las,
                      x,
                      select,
                      continuous,
                      breaks,
                      default_breaks = 4,
                      thresholds = c(-0.05, 0.05),
                      significance,
                      alpha = 0.01) {

  # Check if 'las' is lassie object
  if (!"lassie" %in% class(las)) {
    stop("Invalid 'las' argument: must be a 'lassie' object")
  }

  # Check if 'x' is data.frame
  x <- tryCatch(data.frame(x), error = function(c) {
    stop("Invalid 'x' argument: must be convertible to a data.frame")
  })

  # Check if raw 'las' data.frame has the same number of rows than 'x' data.frame
  if (nrow(x) != nrow(las$data$raw)) {
    stop("Invalid 'las' or 'x' argument: raw 'las' data.frame must have the same number of rows as 'x'. Check if you are not using two different data.frames")
  }

  # Save 'x' data.frame to 'raw'
  raw <- x

  # Handle missing 'select' argument. Select all variables except those in 'las'
  if (missing(select) || is.null(select)) {
    select <- colnames(raw)[! colnames(raw) %in% las$lassie_params$var]
  }

  # Subset 'x' data.frame
  x <- tryCatch(subset(raw, select = select), error = function(c) {
    stop("Invalid 'select' argument: needs to be a vector of column numbers or column names")
  })

  # If 'select' argument is specified as column numbers, convert to variable names
  select <- colnames(x)

  # Check if 'select' variables are not in 'las'
  if (any(select %in% las$lassie_params$var)) {
    stop("Invalid 'select' argument: variable names in 'select' cannot correspond to variable names in 'las'")
  }

  # Handle missing 'continuous' argument
  if (missing(continuous) || is.null(continuous)) {
    continuous <- NULL
  }

  # If 'continuous' argument is specified as column numbers, convert to variable names
  continuous <- tryCatch(colnames(subset(raw, select = continuous)), error = function(c) {
    stop("Invalid 'continuous' argument: needs to be a vector of column numbers or column names")
  })

  # Handle missing 'breaks' argument
  if (missing(breaks)) {
    breaks <- NULL
  }

  # Check thresholds values are numeric and in decreasing order
  if (!is.numeric(thresholds) || length(thresholds) != 2) {
    stop("Invalid 'tresholds' argument: must be two numeric values")
  }
  if (thresholds[1] > thresholds[2]) {
    thresholds <- sort(thresholds, decreasing = FALSE)
  }

  # Handle missing 'significance argument
  if ((missing(significance) || is.null(significance))) {
    if ("permtest" %in% class(x)) significance <- TRUE
    else significance <- FALSE
  }

  # Check if permtest was run before if significance is TRUE
  if (! "permtest" %in% class(las) && significance) {
    warning("Invalid 'significance' argument: run permtest to take into account significance of values. Resuming subgroup analysis with 'significance' set to FALSE")
    significance <- FALSE
  }

  # Check alpha level of significance
  if (significance && (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)) {
    stop("Invalid 'alpha' argument: must be a numeric between 0 and 1 non included")
  }

  # Set measure to the one used in 'las'
  measure <- las$lassie_params$measure

  # Extract local association values and melt array to data.frame 'df'
  df <- reshape2::melt(las$local, value.name = "local")

  # Define positive and negative local association subgroups vectors
  pos <- df$local > thresholds[2]
  neg <- df$local < thresholds[1]

  # Take into account p-values for subgroup definition
  if (significance) {
    local_p <- reshape2::melt(las$local_p, value.name = "local_p")
    df <- plyr::join(df, local_p, by = las$lassie_params$var)
    pos <- pos & df$local_p <= alpha
    neg <- neg & df$local_p <= alpha
  }

  # Define independent local association subgroup vector
  ind <- !pos & !neg

  # Make local association subgroups factor and add to 'df' data.frame
  df$subgroups <- NA_character_
  assoc_levels <- c()

  if (sum(pos) > 0) {
    df$subgroups[pos] <- "Positive"
    assoc_levels <- c(assoc_levels, "Positive")
  }

  if (sum(ind) > 0) {
    df$subgroups[ind] <- "Independent"
    assoc_levels <- c(assoc_levels, "Independent")
  }

  if (sum(neg) > 0) {
    df$subgroups[neg] <- "Negative"
    assoc_levels <- c(assoc_levels, "Negative")
  }

  # Check if at least two subgroups were defined
  if (length(assoc_levels) < 2) {
    stop("Analysis requires at least two local association subgroups. Try setting less stringent 'thresholds' or 'alpha' level, or setting 'significance' to FALSE.")
  }

  # Transform subgroup character into factor
  df$subgroups <- factor(df$subgroups, levels = assoc_levels)

  # Merge 'df' subgroups data.frame with raw data from 'las'
  df <- plyr::join(las$data$pp, df, by = las$lassie_params$var)

  # Merge local association subgroups with 'x' data.frame for lassie function
  na <- stats::na.action(las$data$pp) # find values that were removed by na.omit
  x <- subset(x, subset = ! 1:nrow(x) %in% na) # remove those rows for 'x' data.frame
  df <- cbind(df$subgroups, x)

  # Change name from 'subgroups' to variable names separated by '_'
  new_varname <- paste(las$lassie_params$var, collapse="_")
  colnames(df)[1] <- new_varname
  select <- c(select, new_varname) # add to select

  # Run lassie
  lassie(df,
         select = select,
         continuous = continuous,
         measure = measure,
         breaks = breaks,
         default_breaks = default_breaks)
}
