
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

  # retrieve information about events and times from current tree
  time = object$time
  status = object$status
  ftrs = object$features

  # turn data into a data table if needed
  if(!('data.table'%in%class(newdata))){
    newdata=data.table::data.table(newdata)
  }

  # create data matrix (dmat) using designated features
  dmat = stats::model.matrix(~.,data=newdata[,ftrs,with=F])

  # loop_over indicates which observations to generate predictions for
  loop_over <- 1:nrow(dmat)
  pred_srv = matrix(0, nrow=nrow(dmat), ncol=length(times))

  # this logic is used when computing out-of-bag estimates for ORSF
  if(internal){
    # initialize indicator for in-bag observations
    newdata$inbag=0
    # set inbag=1 for observations with IDs in the current tree's data
    newdata$inbag[newdata$orsf_id%in%object$bootstrap_ids]=1
    # remove in-bag observations from loop_over
    loop_over <- loop_over[newdata$inbag==0]
    # set pred_srv to missing at in-bag observations
    pred_srv[-loop_over,]=NA
  }

  for(i in loop_over){

    # Start at the root of the tree
    current_node = object$nodes$R

    while(!current_node$leaf){

      # compute linear predictor at the current node
      wt = as.numeric(current_node$bwts %*% dmat[i,current_node$bvrs])

      # go to child 1 if the linear predictor is <= cut point, child 2 o.w.
      go_to_node =  if(wt <= current_node$cut_pnt){
        current_node$children[1]
      } else {
        current_node$children[2]
      }

      # reset current node and repeat, unless current node is a leaf
      current_node=object$nodes[[go_to_node]]

    }

    # use step function to generate predicted survival curve from
    # the observed data in the leaf node
    pred_srv[i,] <- current_node$srv(times)

  }

  pred_srv

}
