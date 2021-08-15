#' sim_hidetify_data for the hidetify package
#'
#' Simulated toy data that comes with the hidetify package.
#' The data is generated according to model I of the single detection paper
#' where only the response variable is contaminated.
#' The data set has 1002 variables. The first variable y is the original response variable.
#' The second variable youtlier is the contaminated response variable, where the first 10 observations 
#' are contaminated. The rest of the variables (1000) are the predictors. 
#' 
#'
#' @docType data
#'
#' @usage data(sim_hidetify_data)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{y}{The original response variable}
#'  \item{youtlier}{The contaminated response variable, where the first 10 observations are contaminated}
#'  \item{x1}{The first predictor}
#'  \item{x1000}{The last predictor. There is 1000 predictors in total.}
#' }
#' @references This data set was artificially created for the hidetify package.
#' @keywords datasets
#' @examples
#'
#' data(sim_hidetify_data)
#' head(sim_hidetify_data)
#'
"sim_hidetify_data"