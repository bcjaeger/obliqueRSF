
#' Compute predictions using an oblique random survival forest.
#' @importFrom stats predict
#' @param object An object fitted using the ORSF function.
#' @param newdata A data frame containing observations to predict.
#' @param times A vector of times in the range of the response variable, e.g. times when the response is a survival object, at which to return the survival probabilities.
#' @param ... Other arguments passed to or from other functions.
#' @return A matrix of survival probabilities containing 1 row for each observation and 1 column for each value in times.
#' @export
#' @examples
#' data("pbc",package='survival')
#' pbc$status[pbc$status>=1]=pbc$status[pbc$status>=1]-1
#' pbc$id=NULL
#' fctrs<-c('trt','ascites','spiders','edema','hepato','stage')
#' for(f in fctrs)pbc[[f]]=as.factor(pbc[[f]])
#' pbc=na.omit(pbc)
#'
#' orsf=ORSF(data=pbc,ntree=5)
#' times=seq(365, 365*4,length.out = 10)
#'
#' predict(orsf,newdata=pbc[c(1:5),],times=times)

predict.orsf <- function(object, newdata, times, ...){


  missing_data <- apply(newdata,2,function(x) any(is.na(x)))
  use_imputation=any(missing_data)

  if(use_imputation){
    cat("performing imputation with missForest:")
    imp_data=suppressWarnings(missForest::missForest(xmis=newdata))
    newdata=imp_data$ximp
  }

  for(i in names(newdata)){

    ordered_fac = all(c("ordered", "factor")%in%class(newdata[[i]]))
    if(ordered_fac) newdata[[i]]=as.numeric(newdata[[i]])

  }

  if(class(newdata)[1]!='matrix'){
    newdata = stats::model.matrix(~.,data=newdata)[,-1L, drop = FALSE]
  }

  predict_orsf(object$forest,newx=newdata,times=times)

}
