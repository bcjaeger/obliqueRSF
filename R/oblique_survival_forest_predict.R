
#' Compute predictions using an oblique random survival forest.
#' @importFrom stats predict
#' @param object An object fitted using the ORSF function.
#' @param newdata A data frame containing observations to predict.
#' @param times A vector of times in the range of the response variable, e.g. times when the response is a survival object, at which to return the survival probabilities.
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
#' orsf=ORSF(data=pbc,ntree=20,minsplit=75)
#' times=seq(365, 365*2,length.out = 250)
#'
#' prds=predict(orsf,newdata=pbc,times=times)

predict.orsf <- function(object, newdata, times, ...){

  #newdata=randomForestSRC::impute(data=newdata,nodesize=10,nsplit=10,nimpute=5)

  # apply the prediction function above to each survival tree in object
  lst = purrr::map(object$forest,
    predict,
    newdata=newdata,
    times=times)

  # create an array of prediction matrices for each tree
  arr <- do.call(cbind, lst)
  arr <- array(arr, dim=c(dim(lst[[1]]), length(lst)))

  # average predictions over each tree,
  # remove missing values as needed
  apply(arr, c(1, 2), mean, na.rm = TRUE)

}


