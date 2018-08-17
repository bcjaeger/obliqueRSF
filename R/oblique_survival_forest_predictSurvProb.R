
#' Compute predictions using an oblique random survival forest.
#' @importFrom pec predictSurvProb
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
#' orsf=ORSF(data=pbc,ntree=10,minsplit=75)
#' times=seq(365, 365*2,length.out = 250)
#'
#' prds=pec::predictSurvProb(orsf,newdata=pbc,times=times)

predictSurvProb.orsf <- function(object, newdata, times, ...){

  predict(object,newdata,times,...)

}
