
#' Grow an oblique random survival forest (ORSF)
#' @param data The data used to grow the forest.
#' @param alpha The elastic net mixing parameter. A value of 1 gives the lasso penalty, and a value of 0 gives the ridge penalty. If multiple values of alpha are given, then a penalized model is fit using each alpha value prior to splitting a node.
#' @param ntree The number of trees to grow.
#' @param time A character value indicating the name of the column in the data that measures time.
#' @param status A character value indicating the name of the column in the data that measures participant status. A value of zero indicates censoring and a value of 1 indicates that the event occurred.
#' @param eval_times A numeric vector holding the time values where ORSF out-of-bag predictions should be computed and evaluated.
#' @param features A character vector giving the names of columns in the data set that will be used as features. If NULL, then all of the variables in the data apart from the time and status variable are treated as features. None of these names should contain special characters or spaces.
#' @param min_events_to_split_node The minimum number of events required to split a node.
#' @param min_obs_to_split_node The minimum number of observations required to split a node.
#' @param min_obs_in_leaf_node The minimum number of observations in child nodes.
#' @param min_events_in_leaf_node The minimum number of events in child nodes.
#' @param nsplit The number of random cut-points assessed for each variable.
#' @param max_pval_to_split_node The maximum p-value corresponding to the log-rank test for splitting a node. If the p-value exceeds this cut-point, the node will not be split.
#' @param mtry Number of variables randomly selected as candidates for splitting a node. The default is the square root of the number of features.
#' @param dfmax Maximum number of variables used in a linear combination for node splitting.
#' @param use.cv if TRUE, cross-validation is used to identify optimal values of lambda, a hyper-parameter in penalized regression. if FALSE, a set of candidate lambda values are used. The set of candidate lambda values is built by picking the maximum value of lambda such that the penalized regression model has k degrees of freedom, where k is between 1 and mtry. 
#' @param verbose If verbose=TRUE, then the ORSF function will print output to console while it grows the tree.
#' @param random_seed If a number is given, then that number is used as a random seed prior to growing the forest. Use this seed to replicate a forest if needed.  
#' @return An oblique random survival forest.
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

ORSF <- function(data,
                 alpha=0.50,
                 ntree=100,
                 time='time',
                 status='status',
                 eval_times=NULL,
                 features=NULL,
                 min_events_to_split_node=10,
                 min_obs_to_split_node=10,
                 min_obs_in_leaf_node=3,
                 min_events_in_leaf_node=10,
                 nsplit=15,
                 max_pval_to_split_node=0.05,
                 mtry=ceiling(sqrt(ncol(data)-2)),
                 dfmax=mtry,
                 use.cv=FALSE,
                 verbose=TRUE,
                 random_seed=NULL){
 
  missing_data <- apply(data,2,function(x) any(is.na(x)))
  use_imputation=any(missing_data)

  if(use_imputation){
    cat("\nperforming imputation with missForest:\n")
    imp_data=suppressWarnings(missForest::missForest(xmis=data))
    data=imp_data$ximp
  }
  
  if(is.factor(data[[status]])){
    data[[status]] = as.numeric(data[[status]])
  }
  
  if(any(sort(unique(data[[status]]))!=c(0,1))){
    data[[status]][data[[status]]==min(data[[status]])]=0
    data[[status]][data[[status]]==max(data[[status]])]=1
  }
  
  for(i in names(data)){
    ordered_fac = all(c("ordered", "factor")%in%class(data[[i]]))
    if(ordered_fac) data[[i]]=as.numeric(data[[i]])
  }
  
  if(is.null(features)){
    features = setdiff(names(data),c(time,status))
  } 
  
  #dmat=data.matrix(data[,features])
  data=dplyr::arrange(data,!!rlang::sym(time))
  dmat=model.matrix(~.,data=data[,features])[,-1L]
  time=data[[time]]
  status=data[[status]]
  orsf_ids=1:nrow(data)
  
  srvR <- function(time_indx, status_indx){
    s=survival::survfit(
      survival::Surv(time_indx, status_indx) ~ 1,
      type = "kaplan-meier",se.fit=FALSE)
    list(times=s$time,
         probs=s$surv)
  }
  
  bootR <- function(mat,time,status,inb){
    inb=inb+1
    list(mat=mat[inb,],
         time=time[inb],
         status=status[inb])
  }
  
  netR <- function(dmat,time,status,indx,cols,dfmax,alpha){
    
    indx=indx+1
    cols=cols+1
    
    out = purrr::map(alpha, .f=function(a){
      fit <- suppressWarnings(glmnet::glmnet(
        dmat[indx,cols],
        survival::Surv(time[indx],status[indx]),
        family="cox",
        alpha=a,
        dfmax=dfmax))
      dfs=unique(fit$df)
      dfs=dfs[dfs>=1]
      if(length(dfs)>=1){
        out_indx=sapply(dfs, function(s) min(which(fit$df==s)))
        as.matrix(fit$beta[,out_indx])
      } else {
        matrix(rep(0,dfmax),ncol=1)
      }
    }) 
    
    purrr::reduce(out, cbind)
    
  }
  
  cv.netR <- function(dmat,time,status,indx,cols,dfmax,alpha){
    indx=indx+1
    cols=cols+1
    out = purrr::map(alpha, .f=function(a){
      cv.fit <- suppressWarnings(glmnet::cv.glmnet(
        dmat[indx,cols],
        survival::Surv(time[indx],status[indx]),
        family="cox", keep = FALSE, grouped = TRUE,
        alpha=a, nfolds=min(5,length(indx)), dfmax=dfmax))
      as.matrix(cbind(
        coef(cv.fit, cv.fit$lambda[min(which(cv.fit$glmnet.fit$df>0))]),
        coef(cv.fit,'lambda.1se'),
        coef(cv.fit,'lambda.min')
      ))
    })
    purrr::reduce(out, cbind)
  }
  
  fevalR <- function(prd,time,status,eval.times){
    
    ntimes=length(eval.times)
    conc=pec::cindex(matrix(prd[,ntimes],ncol=1),
                     formula=Surv(time,status)~1,
                     eval.times=eval.times[ntimes],
                     data=data.frame(time=time,status=status))
    
    conc=conc$AppCindex$matrix
    
    intbs=suppressMessages(pec::crps(pec::pec(
      list(ORSF=prd),times=eval.times,exact=F,
      start=eval.times[1],maxtime=eval.times[length(eval.times)],
      formula=Surv(time,status)~1,
      data=data.frame(time=time,status=status))))
    
    list(concordance=conc,integrated_briscr=intbs)
    
  }
  
  if(!is.null(random_seed)){
    set.seed(random_seed)
  }
  
  if(is.null(eval_times)){
    eval_times=seq(min(time[status==1]),max(time[status==1]),length.out=50)
  } 
  
  orsf=ORSFcpp(dmat=dmat,
               features=colnames(dmat),
               alpha=alpha,
               time=time,
               status=status,
               eval_times=eval_times,
               min_events_to_split_node=min_events_to_split_node,
               min_obs_to_split_node=min_obs_to_split_node,
               min_obs_in_leaf_node=min_obs_in_leaf_node,
               min_events_in_leaf_node=min_events_in_leaf_node,
               mtry=mtry,
               dfmax=dfmax,
               nsplit=nsplit,
               ntree=ntree,
               mincriterion=qchisq(1-max_pval_to_split_node,df=1),
               verbose=verbose,
               surv_KM_Rfun=srvR,
               bootstrap_Rfun=bootR,
               glmnet_Rfun=if(use.cv){cv.netR} else {netR},
               forest_eval_Rfun=fevalR)
  
  output=structure(
    list(forest = orsf$forest,
         oob_times = eval_times,
         oob_preds = orsf$oob_preds,
         oob_error = orsf$oob_error,
         data=data,
         call = match.call()),
    class = "orsf")
  
  output$imputation_used= if(use_imputation) TRUE else FALSE 

  return(output)
  
}

