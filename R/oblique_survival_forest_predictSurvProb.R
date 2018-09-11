
#' Compute predictions using an oblique random survival forest.
#' @importFrom pec predictSurvProb
#' @param object A fitted model from which to extract predicted survival probabilities
#' @param newdata A data frame containing predictor variable combinations for which to compute predicted survival probabilities.
#' @param times A vector of times in the range of the response variable, e.g. times when the response is a survival object, at which to return the survival probabilities.
#' @param ... Additional arguments that are passed on to the current method.
#' @return A matrix of survival probabilities containing 1 row for each observation and 1 column for each value in times.
#' @export
#' @examples
#' 
#' data("pbc",package='survival')
#' pbc$status[pbc$status>=1]=pbc$status[pbc$status>=1]-1
#' pbc$id=NULL
#' fctrs<-c('trt','ascites','spiders','edema','hepato','stage')
#' for(f in fctrs)pbc[[f]]=as.factor(pbc[[f]])
#' pbc=na.omit(pbc)
#'
#' orsf=ORSF(data=pbc)
#' times=seq(365, 365*4,length.out = 10)
#'
#' predict(orsf,newdata=pbc[c(1:5),],times=times)

predictSurvProb.orsf <- function(object, newdata, times, ...){

  predict(object,newdata,times,...)

}
