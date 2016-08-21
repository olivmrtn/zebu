#' @title Local Association Measures
#'
#' @description Estimates local (and global) association measures: Ducher's Z and pointwise mutual
#' information and normalized pointwise mutual information.
#'
#' @inheritParams preprocess
#' @inheritParams local_association
#'
#' @return An instance of S3 \link[base]{class} \code{\link[zebu]{lassie}} with
#'   the following objects:
#'   \itemize{
#'   \item data: raw and preprocessed data.frames (see \link[zebu]{preprocess}).
#'   \item prob probability arrays (see \link[zebu]{estimate_prob}).
#'   \item global global association (see \link[zebu]{local_association}).
#'   \item local local association arrays (see \link[zebu]{local_association}).
#'   \item lassie_params parameters used in lassie.
#'   }
#'
#' @seealso Results can be visualized using \code{\link[zebu]{plot.lassie}} and
#' \code{\link[zebu]{print.lassie}}
#' methods. \code{\link[zebu]{plot.lassie}} is only available in the bivariate case and returns
#' a tile plot representing the probability or local association measure matrix.
#' \code{\link[zebu]{print.lassie}} shows an array or a data.frame.
#'
#' Results can be saved using \code{\link[zebu]{write.lassie}}.
#'
#' The \code{\link[zebu]{permtest}} function accesses the significance of local and global
#'association values using p-values estimated by permutations.
#'
#' The \code{\link[zebu]{subgroups}} function identifies if the
#' association between variables is dependent on the value of another variable.
#'
#' @examples
#' # In this example, we will use the 'mtcars' dataset
#'
#' # Selecting a subset of mtcars.
#' # Takes column names or numbers.
#' # If nothing was specified, all variables would have been used.
#' select <- c('mpg', 'cyl') # or select <- c(1, 2)
#'
#' # Specifying 'mpg' as a continuous variables using column numbers
#' # Takes column names or numbers.
#' # If nothing was specified, all variables would have been used.
#' continuous <- 'mpg' # or continuous <- 1
#'
#' # How should breaks be specified?
#' # Specifying equal-width discretization with 5 bins for all continuous variables ('mpg')
#' # breaks <- 5
#'
#' # Specifying user-defined breakpoints for all continuous variables.
#' # breaks <- c(10, 15, 25, 30)
#'
#' # Same thing but only for 'mpg'.
#' # Here both notations are equivalent because 'mpg' is the only continuous variable.
#' # This notation is useful if you wish to specify different break points for different variables
#' # breaks <- list('mpg' = 5)
#' # breaks <- list('mpg' = c(10, 15, 25, 30))
#'
#' # Calling lassie
#' # Not specifying breaks means that the value in default_breaks (4) will be used.
#' las <- lassie(mtcars, select = c(1, 2), continuous = 1)
#'
#' # Print local association to console as an array
#' print(las)
#'
#' # Print local association and probabilities
#' # Here only rows having a positive local association are printed
#' # The data.frame is also sorted by observed probability
#' print(las, type = 'df', range = c(0, 1), what_sort = 'obs')
#'
#' # Plot results as heatmap
#' plot(las)
#'
#' # Plot observed probabilities using different colours
#' plot(las, what_x = 'obs', low = 'white', mid = 'grey', high = 'black', text_colour = 'red')
#'
#' # Write results to text file
#' write.lassie(las, file = 'test.csv')
#'
#' # Retrieve results
#' lassie_df <- read.table('test.csv', sep = ',', header = TRUE)
#'
#' @export
#'
lassie <- function(x,
                   select,
                   continuous,
                   breaks,
                   measure = "z",
                   default_breaks = 4) {

  # Preprocess data (handle missing data and discretize continuous variables)
  pre <- preprocess(x,
                    select = select,
                    continuous = continuous,
                    breaks = breaks,
                    default_breaks = default_breaks)

  # Compute probabilities
  prob <- estimate_prob(pre$pp)

  # Compute local and global association
  lam <- local_association(prob, measure)
  global <- lam$global
  local <- lam$local

  # Retrieve variable and measure names
  lassie_params <- list(var = colnames(pre$pp),
                        measure = measure,
                        select = pre$select,
                        continuous = pre$continuous,
                        breaks = pre$breaks,
                        default_breaks = pre$default_breaks)

  # Return lassie and 'measure' S3 object
  structure(list(data = pre[c("raw", "pp")],
                 prob = prob,
                 global = global,
                 local = local,
                 lassie_params = lassie_params),
            class = c("lassie"))
}

#' @title Print a lassie object
#'
#' @description Print a \code{\link[zebu]{lassie}} object as an array or a data.frame.
#'
#' @inheritParams lassie_get
#' @inheritParams format.lassie
#'
#' @param type print style: 'array' for array or 'df' for data.frame.
#' @param ... other arguments passed on to methods. Not currently used.
#'
#' @seealso \code{\link[zebu]{lassie}}, \code{\link[zebu]{permtest}}
#'
#' @export
#'
print.lassie <- function(x,
                         type,
                         what_x,
                         range,
                         what_range,
                         what_sort,
                         decreasing,
                         na.rm,
                         ...) {

  # Check arguments
  if (missing(type) || is.null(type)) {
    if (ncol(x$data$pp) == 2) {
      type <- "array"
    } else {
      type <- "df"
    }
  } else if (!type %in% c("array", "df")) {
    stop("Invalid 'type' argument: choose from\n 'array': Array\n 'df' Melt array into a data frame")
  }

  if (missing(range)) {
    range <- NULL
  }

  # Format for array
  if (type == "array") {
    if (missing(what_x) || is.null(what_x)) {
      what_x <- "local"
    }
    print_x <- lassie_get(x, what_x)

    if (!is.null(range)) {
      range_x <- lassie_get(x, what_range)
      in_range <- !is.na(range_x) & range_x >= range[1] & range_x <= range[2]
      print_x[which(!in_range, arr.ind = TRUE)] <- NA_real_
    }

    # Format for data.frame
  } else if (type == "df") {
    if (missing(what_x)) {
      what_x <- NULL
    }
    if (missing(what_range)) {
      what_range <- NULL
    }
    if (missing(what_sort)) {
      what_sort <- NULL
    }
    if (missing(decreasing)) {
      decreasing <- NULL
    }
    if (missing(na.rm)) {
      na.rm <- NULL
    }

    print_x <- format.lassie(x,
                             what_x = what_x,
                             range = range,
                             what_range = what_range,
                             what_sort = what_sort,
                             decreasing = decreasing,
                             na.rm = na.rm)
  }

  # Make header (measure name and global association value)
  header <- paste0("Measure: ", measure_name(x), "\nGlobal: ", x$global)
  if (! is.null(x$global_p)) {
    header <- paste0(header, " (p-value=", x$global_p, ")")
  }
  header <- paste0(header, "\n")

  # Print to console
  cat(header)
  print(print_x)
  invisible(x)
}

#' @title Format a lassie object
#'
#' @description Formats a \code{\link[zebu]{lassie}} object for printing to console
#' (see \code{\link[zebu]{print.lassie}}) and for writing to a file
#' (see \code{\link[zebu]{write.lassie}}). Melts probability or local association
#' measure arrays into a data.frame.
#'
#' @inheritParams lassie_get
#' @param range range of values to be retained (vector of two numeric values).
#' @param what_range character specifying what value \code{range} refers to
#' (same options as \code{what_x}).
#' By default, takes the first value in \code{what_x}.
#' @param what_sort character specifying according to which values should \code{x} be sorted
#' (same options as \code{what_x}).
#' By default, takes the first value in \code{what_x}.
#' @param decreasing logical value specifying sort order.
#' @param na.rm logical value indicating whether NA values should be stripped.
#' @param ... other arguments passed on to methods. Not currently used.
#'
#' @seealso \code{\link[zebu]{lassie}}
#'
#' @export
#'
format.lassie <- function(x,
                          what_x,
                          range,
                          what_range,
                          what_sort,
                          decreasing,
                          na.rm, ...) {

  # Check arguments
  if (missing(what_x) || is.null(what_x)) {
    what_x <- c("local", "obs", "exp")
    if ("permtest" %in% class(x))
      what_x <- c(what_x, "local_p")
  }

  if (missing(na.rm) || is.null(na.rm)) {
    na.rm <- FALSE
  }

  # Melt array to data frame
  local <- reshape2::melt(x$local)
  format_x <- local[, 1:(ncol(local) - 1)]
  for (i in what_x) {
    format_x[i] <- reshape2::melt(lassie_get(x, i))$value
  }

  # Retain only rows that in the specified range
  if (!missing(range) && !is.null(range)) {
    if (missing(what_range) || is.null(what_range)) {
      what_range <- what_x[1]
    } else if (any(!what_range %in% what_x)) {
      stop("Invalid 'what_range' argument: values must also be present in 'what_x' argument")
    }
    if (length(what_range) > 1) {
      stop("Invalid 'what_range' argument: can only handle one variable name")
    }
    if (missing(na.rm)) {
      na.rm <- FALSE
    }

    range_x <- lassie_get(x, what_range)
    in_range <- !is.na(range_x) & range_x >= range[1] & range_x <= range[2]
    format_x <- format_x[in_range, ]
  }

  # Sort data frame
  if (missing(what_sort) || is.null(what_sort)) {
    what_sort <- what_x[1]
  } else if (any(!what_sort %in% colnames(format_x))) {
    stop(paste("Invalid 'what_sort' argument: Choose from:",
               paste(colnames(format_x), collapse = ", ")))
  }
  index <- do.call("order", subset(format_x, select = what_sort))

  # Reverse sort order
  if (missing(decreasing) || is.null(decreasing))
    decreasing <- TRUE
  if (decreasing) {
    index <- rev(index)
  }

  format_x <- format_x[index, ]
  if (na.rm) {
    format_x <- stats::na.omit(format_x)
  }
  format_x
}

#' @title Plot a lassie object
#'
#' @description Plots a \code{\link[zebu]{lassie}} object as a tile plot using
#' the ggplot2 package. Only available for bivariate association.
#' @inheritParams lassie_get
#'
#' @param digits integer indicating the number of decimal places.
#' @param low colour for low end of the gradient.
#' @param mid colour for midpoint of the gradient.
#' @param high colour for high end of the gradient.
#' @param na colour for NA values.
#' @param text_colour colour of text inside cells.
#' @param text_size integer indicating text size inside cells.
#' @param limits limits of gradient.
#' @param midpoint midpoint of gradient.
#' @param ... other arguments passed on to methods. Not currently used.
#'
#' @seealso \code{\link[zebu]{lassie}}
#'
#' @export
#'
plot.lassie <- function(x,
                        what_x = "local",
                        digits = 3,
                        low = "blue",
                        mid = "white",
                        high = "red",
                        na = "purple",
                        text_colour = "black",
                        text_size,
                        limits,
                        midpoint,
                        ...) {

  # Check if lassie object is bivariate
  if (ncol(x$data$pp) != 2) {
    stop("Cannot make tile plot for more than 2 variables.\n Use 'print()' function to visualize results.")
  }

  # Define text size in cells
  if (missing(text_size)) {
    if (what_x == "local" & "permtest" %in% class(x)) {
      text_size <- 4
    } else {
      text_size <- 5
    }
  }

  # Get values from lassie object
  measure <- x$lassie_params$measure
  tile_value <- lassie_get(x, what_x)
  tile_value[is.infinite(tile_value) | is.nan(tile_value)] <- NA
  tile_text <- as.character(round(tile_value, digits))

  # Define limits of gradients
  if (missing(limits) || is.null(limits)) {
    if (what_x %in% c("local")) {
      midpoint <- 0
      if (measure %in% c("z", "npmi")) {
        limits <- c(-1, 1)
      } else {
        # Pointwise Mutual Information
        limits <- range(tile_value, na.rm = TRUE)
      }
    } else {
      # Probabilities
      limits <- c(0, 1)
      midpoint <- 0.5
    }
  }

  # Define midpoint of gradient
  if (missing(midpoint) || is.null(midpoint)) {
    if (what_x %in% c("local")) {
      midpoint <- 0
    } else {
      midpoint <- 0.5
    }
  }

  # Format values to be displayed in cells
  if (what_x == "obs" | what_x == "exp") {
    tile_text <- paste0(tile_text, "\n", round(nrow(x$data$pp) * tile_value, digits))

  } else if (what_x == "local" & "permtest" %in% class(x)) {
    tile_text <- paste0(tile_text, "\n", format(x$local_p, digits = digits, scientific = TRUE))
  }

  # Format values in data.frame for ggplot2
  ggdf <- data.frame(reshape2::melt(tile_value), tile_text, stringsAsFactors = FALSE)
  colnames(ggdf) <- c("Rows", "Columns", "Value", "Text")
  ggdf$Text[is.na(tile_value)] <- ""
  axis_labels <- names(dimnames(tile_value))

  # Make title
  title <- paste0(measure_name(x), "\n", round(x$global, digits))
  if (! is.null(x$global_p)) {
    title <- paste0(title, " (p-value=", x$global_p, ")")
  }

  # Plot as tiles
  ggp <- ggplot2::ggplot(ggdf, ggplot2::aes_string(x = "factor(Columns)",
                                                   y = "factor(Rows)",
                                                   fill = "Value")) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes_string(label = "Text"), size = text_size, colour = text_colour) +
    ggplot2::scale_fill_gradient2(low = low,
                                  mid = mid,
                                  high = high,
                                  na.value = na,
                                  midpoint = midpoint,
                                  limits = limits) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(axis_labels[2]) +
    ggplot2::ylab(axis_labels[1]) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  ggp
}

#' @title Write a lassie object
#'
#' @description Writes \code{\link[zebu]{lassie}} object to a file in a table structured format.
#'
#' @inheritParams utils::write.table
#'
#' @param x \code{\link[zebu]{lassie}} S3 object.
#' @param file character string naming a file.
#' @param ... other arguments passed on to write.table.
#'
#' @seealso \code{\link[zebu]{lassie}}, \code{\link[zebu]{permtest}}
#'
#' @export
#'
write.lassie <- function(x,
                         file,
                         sep = ",",
                         dec = ".",
                         col.names = TRUE,
                         row.names = FALSE,
                         quote = TRUE,
                         ...) {

  if (missing(file) || is.null(file)) {
    stop("Invalid 'file' argument: please use a character string naming a file ")
  }
  if (sep == dec) {
    stop("Column separator ('sep') cannot be the same as decimal separator ('dec')")
  }

  com <- generate_comments(x)
  writeLines(text = com, con = file)

  # Suppress warnings because append is set to TRUE
  suppressWarnings(utils::write.table(format.lassie(x), file = file, append = TRUE, ...))
}
