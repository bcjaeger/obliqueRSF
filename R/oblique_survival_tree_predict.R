
#' Compute predictions using an oblique survival tree
#' @importFrom stats predict
#' @param object The oblique survival tree used to compute predictions.
#' @param newdata A data frame containing observations to predict.
#' @param times A vector of times in the range of the response variable, e.g. times when the response is a survival object, at which to return the survival probabilities.
#' @param internal This argument is only used when computing out of bag error for oblique random survival forests and does not require manual supervision.
#' @param ... Other arguments passed to or from other functions.
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
#' ost=OST(data=pbc, minsplit = 125)
#'
#' #One and 2 year survival probability
#' times=seq(365, 365*2,length.out = 250)
#'
#' prds=predict(ost,newdata=pbc,times=times)

predict.ost=function(object,newdata,times,internal=FALSE,...){

  newdata=data.table::data.table(newdata)
  missing_data <- apply(newdata,2,function(x) any(is.na(x)))
  
  if(any(missing_data)){
    cat("performing imputation with missForest:")
    imp_data=suppressWarnings(missForest::missForest(xmis=newdata))
    newdata=imp_data$ximp
  }
  
  for(i in names(newdata)){
    
    ordered_fac = all(c("ordered", "factor")%in%class(newdata[[i]]))
    if(ordered_fac) newdata[[i]]=as.numeric(newdata[[i]])
    
  }
  
  newx = stats::model.matrix(~.,data=newdata)[,-1L]
  
  predict_ost(object,newx=newx,times)
  

}
