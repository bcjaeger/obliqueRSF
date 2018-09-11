

#' Grow an oblique random survival forest (ORSF)
#' @param x an ORSF object (i.e. the object returned from the ORSF function)
#' @param ... additional arguments passed to print 
#' @return A printed summary of the oblique random survival forest.
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
#' print(orsf)

print.orsf <- function(x,...){
  cat('\nOblique Random Survival Forest: ')
  print(x$call)
  print(x$oob_error)
  
}