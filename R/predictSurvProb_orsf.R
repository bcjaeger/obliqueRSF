
#' Compute predictions using an oblique random survival forest.
#' @param object An object fitted using the ORSF function.
#' @param newdata A data frame containing observations to predict.
#' @param times A vector of times in the range of the response variable, e.g. times when the response is a survival object, at which to return the survival probabilities.
#' @param ... Additional parameters for predict.ost
#' @return A matrix of survival probabilities containing 1 row for each observation and 1 column for each value in times.


predictSurvProb.orsf <- function(object,newdata,times,...){
  predict.orsf(object=object,newdata=newdata,times=times,internal=FALSE,...)
}
