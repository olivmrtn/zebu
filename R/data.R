#' Resistance to drug treatment
#'
#' Simulated clinical trial where patient recovery is dependent on drug intake and resistance status.
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{drug}{binary variable (placebo, drug), did patient receive drug}
#'   \item{resistance}{binary variable (sensitive, resistant), is patient resistance to drug}
#'   \item{prebiom}{continuous variable between 0 and 1, biomarker that represents health status of patient before treatment; healthy patients have values around 0.6}
#'   \item{postbiom}{continuous variable between 0 and 1, biomarker that represents health status of patient after treatment; healthy patients have values above 0.6}}
"trial"
