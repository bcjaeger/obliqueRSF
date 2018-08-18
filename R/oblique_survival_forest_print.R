

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
#' orsf=ORSF(data=pbc,ntree=15,minsplit=75)
#' print(orsf)
#' 

print.orsf <- function(x,...){
  cat('\nOblique Random Survival Forest: ')
  print(x$call)
  cat('\nOut of bag error estimates: (lower is better )\n')
  cat('  integrated Brier score:', 
      format(round(x$oob_perr[length(x$oob_perr)],5),nsmall=5),'\n')
  cat('  integrated concordance:', 
      format(round(x$oob_cerr[length(x$oob_cerr)],5),nsmall=5),'\n')
  
}