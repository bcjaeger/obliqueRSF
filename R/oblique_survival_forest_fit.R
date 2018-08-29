
#' Grow an oblique random survival forest (ORSF)
#' @param data The data used to grow the forest.
#' @param time A character value indicating the name of the column in the data that measures time.
#' @param status A character value indicating the name of the column in the data that measures participant status. A value of zero indicates censoring and a value of 1 indicates that the event occurred.
#' @param features A character vector giving the names of columns in the data set that will be used as features. If NULL, then all of the variables in the data apart from the time and status variable are treated as features. None of these names should contain special characters or spaces.
#' @param verbose If verbose=TRUE, then the ORSF function will print output to console while it grows the tree.
#' @param ntree The number of trees to grow.
#' @param nmin The minimum number of events required to split a node.
#' @param nsplit The number of random cut-points assessed for each variable.
#' @param minsplit The minimum number of observations required to split a node.
#' @param mtry Number of variables randomly selected as candidates for splitting a node. The default is ceiling(sqrt(p)), where p is the number of features.
#' @param nmin_leaf The minimum number of events required for a leaf node.
#' @param alpha The elastic net mixing parameter. A value of 1 gives the lasso penalty, and a value of 0 gives the ridge penalty.
#' @param use.cv if TRUE, then cross-validation is used prior to splitting each node in order to identify optimal values of lambda, a hyper-parameter in penalized regression.
#' @param maxdf_lincombo This argument is ignored if cross-validation is applied, i.e. if use.cv=TRUE. The maximum number of degrees of freedom used by penalized regression models.
#' @param tree.err (Not yet implemented) if TRUE, then out-of-bag error is computed after each new tree is grown.
#' @param importance (Not yet implemented) if TRUE, variable importance is computed using axis based random survival forests.
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
#' orsf=ORSF(data=pbc,ntree=15,minsplit=75)


ORSF <- function(data,
                 time='time',
                 status='status',
                 features=NULL,
                 nmin=5,
                 nsplit=10,
                 minsplit=30,
                 mtry=ceiling(sqrt(ncol(data)-2)),
                 nmin_leaf=1,
                 alpha=0.01,
                 use.cv=TRUE,
                 maxdf_lincombo=if(use.cv) NULL else 5,
                 ntree=100,
                 verbose=TRUE,
                 importance=FALSE,
                 tree.err = FALSE){

  .integral <- function(xvec,yvec){
    integral=stats::integrate(stats::approxfun(xvec,yvec),
                       min(xvec),max(xvec),
                       subdivisions=2000)
    integral$value / (max(xvec)-min(xvec))
  }

  Hist=function(...){
    prodlim::Hist(...)
  }

  Surv=function(...){
    survival::Surv(...)
  }
  
  if(ntree<20 & tree.err) stop('ntree should be >20 if tree.err is TRUE')
  
  missing_data <- apply(data,2,function(x) any(is.na(x)))

  if(any(missing_data)){
    cat("\nperforming imputation with missForest:\n")
    imp_data=suppressWarnings(missForest::missForest(xmis=data))
    data=imp_data$ximp
  }
  
  data$orsf_id=1:nrow(data)
  data=data.table::data.table(data,key='orsf_id')
  boot_ids=unique(data$orsf_id)
  
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

  features = setdiff(names(data),c(time,status,'orsf_id'))
  forest = vector(mode='list',length=ntree); i = 1; error_count=0

  while(i <= ntree){

    if(error_count>5) stop('Too many lumberjacks')

    if(verbose){
      cat('\n Growing tree no.', i, '.')
    }

    inb = sample(data$orsf_id, nrow(data), replace=TRUE)
    oob = setdiff(data$orsf_id,inb)
    inbd = data[inb,]
    oobd = data[oob,]

    new_tree <- try(OST(data=inbd,
                        time=time,
                        status=status,
                        features=features,
                        nmin=nmin,
                        nsplit=nsplit,
                        minsplit=minsplit,
                        mtry=mtry,
                        nmin_leaf=nmin_leaf,
                        alpha=alpha,
                        boot_ids=inb,
                        verbose=verbose,
                        tree_lab=i,
                        data.preprocessed=TRUE),
      silent = FALSE)

    if(class(new_tree)!='try-error'){

      if(length(new_tree$nodes)>=3){

        forest[[i]]=new_tree; i=i+1

      }

      if(verbose) cat('| done |', length(new_tree$nodes), 'nodes grown')

    } else {
      # lumberjacks!
      error_count<-error_count+1
    }

  }

  names(forest)=paste0('Tree',1:ntree)

  times=unique(sort(data[[time]][data[[status]]==1]))
  
  dmat = stats::model.matrix(~.,data=data)[,-1L]
  
  oob_prd = predict_orsf(forest=forest,
                         newx=dmat,
                         times=times)
  
  forest_eval<-function(tree_num){
    
    if(verbose) cat('\nEvaluating predictions using', tree_num, 'trees')
    
    oob_prd = predict_orsf(forest=forest[1:tree_num],
                           newx=dmat,
                           times=times)
    
    sfrm=stats::as.formula(paste0("Surv(",time,',',status,')~1'))
    
    oob_perr = suppressMessages(
      pec::pec(oob_prd,
               data=data,
               times=times[-length(times)],
               exact=FALSE,formula=sfrm,
               cens.model="cox"))
    
    oob_perr=pec::crps(oob_perr)
    oob_perr=oob_perr[2,1]
    
    oob_cstats = suppressMessages(
      pec::cindex(oob_prd,
                  data = data,
                  eval.times=times,
                  formula = sfrm,
                  cens.model = 'cox'))
    
    oob_cstats=oob_cstats$AppCindex
    oob_cstats=data.frame(cbind(oob_cstats[[1]],times))
    names(oob_cstats)=c('cstat','time')
    
    oob_cerr = 1-with(oob_cstats,.integral(xvec=time,yvec=cstat))
    
    c(oob_perr=oob_perr,oob_cerr=oob_cerr)
    
  }

  if(tree.err){
    
    eval_indx=seq(15,ntree,by=5)
    tree_err=purrr::map(eval_indx,~forest_eval(tree_num=.))
    tree_err=purrr::reduce(tree_err,rbind)
    tree_err=data.frame(tree_err)
    rownames(tree_err)=NULL
    tree_err$trees=eval_indx
    
    oob_perr = tree_err$oob_perr[length(tree_err$oob_perr)]
    oob_cerr = tree_err$oob_cerr[length(tree_err$oob_cerr)]
    
  } else {
    
    tree_err=forest_eval(ntree)
    oob_perr = tree_err['oob_perr']
    oob_cerr = tree_err['oob_cerr']
    
  }
  
  if(verbose) cat('\n')

  structure(
    list(forest = forest,
         tree_err= if(tree.err) tree_err else NULL,
         oob_perr = oob_perr,
         oob_cerr = oob_cerr,
         call = match.call()
    ),
    class = "orsf")

}

# if(importance){
#
#   vimp = suppressWarnings(purrr::map(features,.f=function(v){
#
#     tmpF = forest
#
#     for(i in 1:length(tmpF)){
#       for(j in 1:length(tmpF[[i]]$nodes))
#       if('bvrs' %in% names(tmpF[[i]]$nodes[[j]])){
#         if(v %in% tmpF[[i]]$nodes[[j]]$bvrs){
#           tmpF[[i]]$nodes[[j]]$bwts[v]=0
#         }
#       }
#     }
#
#     oob_lst=purrr::map(tmpF, predict.internal_tree,newdata=data,times=times)
#     arr=oob_lst %>%
#       reduce(cbind)%>%
#       array(dim=c(dim(oob_lst[[1]]),length(oob_lst)))
#
#     oob_prd=apply(arr, c(1, 2), mean, na.rm = TRUE)
#     all_nan=apply(oob_prd,1,function(x) all(is.nan(x)))
#
#     sfrm=paste0("Surv(",time,',',status,')~1')%>%stats::as.formula()
#
#     shuffled_perr=suppressMessages(
#       pec::pec(oob_prd[!all_nan,],
#                data=data[!all_nan,],
#                times=times[-length(times)],
#                exact=FALSE,formula=sfrm,
#                cens.model="cox") %>%
#         pec::crps() %>% magrittr::extract(2,1))
#
#     shuffled_cerr = 1 - suppressMessages(
#       pec::cindex(oob_prd[!all_nan,],
#                   data = data[!all_nan,],
#                   eval.times=times,
#                   formula = sfrm,
#                   cens.model = 'cox')) %>%
#       use_series('AppCindex') %>%
#       magrittr::extract2(1) %>%
#       cbind(times) %>% data.frame() %>%
#       set_names(c('cstat','time')) %>%
#       mutate(diff=c(0,diff(time))) %>%
#       dplyr::summarise(sum(cstat*diff)/max(time)) %>%
#       as.numeric()
#
#     data.frame(ibris_vimp=shuffled_perr-oob_perr[length(oob_perr)],
#                cstat_vimp=shuffled_cerr-oob_cerr[length(oob_cerr)])
#
#   }) %>%
#     reduce(rbind) %>%
#     data.frame() %>%
#     set_names(c('ibris_vimp','cstat_vimp')) %>%
#     mutate(variable=pvars) %>%
#     left_join(elastic_importance,by='variable')) %>%
#     mutate(comb_vimp=scales::rescale(ibris_vimp,to=c(0,1))-1+elast_vimp) %>%
#     dplyr::arrange(desc(comb_vimp))
#
# } else {
#
#   vimp=NULL
#
# }