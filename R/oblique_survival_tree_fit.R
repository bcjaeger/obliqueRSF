

#' Grow an oblique survival tree
#' @param data The data used to grow the survival tree
#' @param time A character value indicating the name of the column in the data that measures time.
#' @param  status A character value indicating the name of the column in the data that measures participant status. A value of zero indicates censoring and a value of 1 indicates that the event occurred.
#' @param features A character vector giving the names of columns in the data set that will be used as features. If NULL, then all of the variables in the data apart from the time and status variable are treated as features. None of these names should contain special characters or spaces
#' @param verbose If verbose=TRUE, then the OST function will print output to console while it grows the tree.
#' @param nmin The minimum number of events required to split a node.
#' @param nsplit The number of random cut-points assessed for each variable.
#' @param minsplit The minimum number of observations required to split a node.
#' @param mtry Number of variables randomly selected as candidates for splitting a node. The default is ceiling(sqrt(p)), where p is the number of features.
#' @param nmin_leaf The minimum number of events required for a leaf node.
#' @param alpha The elastic net mixing parameter. A value of 1 gives the lasso penalty, and a value of 0 gives the ridge penalty.
#' @param use.cv if TRUE, then cross-validation is used prior to splitting each node in order to identify optimal values of lambda, a hyper-parameter in penalized regression.
#' @param maxdf_lincombo This argument is ignored if cross-validation is applied, i.e. if use.cv=TRUE. The maximum number of degrees of freedom used by penalized regression models.
#' @param boot_ids This argument is only used when building an oblique random survival forest and does not require manual supervision.
#' @param data.preprocessed This argument is only used when building an oblique random survival forest and does not require manual supervision.
#' @param tree_lab This argument is only used when building an oblique random survival forest and does not require manual supervision.
#' @return An oblique survival tree
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

# data=pbc
# time='time'
# status='status'
# features=setdiff(names(data),c('time','status'))
# verbose=TRUE
# nmin=5
# nsplit=10
# minsplit=30
# mtry=sqrt(length(features))
# nmin_leaf=1
# alpha=0.01
# use.cv=TRUE
# maxdf_lincombo=NULL
# boot_ids=NULL
# data.preprocessed=FALSE
# tree_lab=NULL
# tree.err=TRUE

OST <- function(data, time='time', status='status',features=NULL,verbose=TRUE,
  nmin=5,nsplit=10,minsplit=30,mtry=ceiling(sqrt(ncol(data)-2)),nmin_leaf=1,
  alpha=0.01,use.cv=TRUE,maxdf_lincombo=if(use.cv) NULL else 5,
  boot_ids=NULL,data.preprocessed=FALSE,tree_lab=NULL){

  is.error<-function(x) class(x)[1] == 'try-error'

  cat_mode <- function(cat_var){
    tt <- table(cat_var)
    names(tt)[which.max(tt)]
  }

  soft_sample <- function(x,size,...){
    sample(x,size=min(length(x),size),replace=FALSE,...)
  }

  if(nmin_leaf==0) stop("Minimum number of events in leaf nodes must be > 0")
  if(length(unique(data[[status]]))>2)
    stop('Status should only have 2 values')

  if(!('orsf_id'%in%names(data))){
    data$orsf_id=1:nrow(data)
  }

  if(!('data.table'%in%class(data))){
    data=data.table::data.table(data,key='orsf_id')
  }

  if(!data.preprocessed){
    
    data=data.table::data.table(data,key='orsf_id')
    boot_ids=unique(data$orsf_id)
    
    if(is.factor(data[[status]])){
      data[[status]] = as.numeric(data[[status]])
    }
    
    if(any(sort(unique(data[[status]]))!=c(0,1))){

      data[[status]][data[[status]]==min(data[[status]])]=0
      data[[status]][data[[status]]==max(data[[status]])]=1

    }
  }  
  
  if(is.null(features)){
    .ftrs = setdiff(names(data),c(time,status,'orsf_id'))
  } else {
    .ftrs = features
  }

  nodes = list(
    R = structure(
      list(name='R',
           leaf=TRUE,
           dont_split=FALSE,
           parent='root',
           children=NULL,
           depth = 0,
           nsubs=length(boot_ids),
           nevents=sum(data[[status]])),
      class = 'ost_node'))

  nodes_to_grow <- nodes$R$name
  node_ids  <- rep(nodes$R$name,nrow(data))
  dmat=stats::model.matrix(~.+0,data=data)

  #parent=nodes[['R']]
  current_node=nodes_to_grow[1]

  while(length(nodes_to_grow)>0){

    for(current_node in nodes_to_grow){

      # indx indicates which ids are in the current node
        indx=node_ids==nodes[[current_node]]$name

        ..ftrs=.ftrs

        # remove predictors that are constant for current ids
        for(i in ..ftrs){
          if(length(unique(data[[i]][indx]))==1){
            ..ftrs=..ftrs[-which(..ftrs==i)]
          }
        }

        nevnt = sum(data[indx,status,with=F])
        nrows = length(unique(data$orsf_id[indx]))
        ncols = length(..ftrs)

        # Split using lasso regression
        if(nrows >= minsplit & nevnt >= nmin & ncols >=1){

          split_candidates <- soft_sample(..ftrs, mtry)

          if(use.cv){

            fit = suppressWarnings(try(glmnet::cv.glmnet(
              stats::model.matrix(
                ~.+1,data=data[indx,split_candidates,with=FALSE]),
              survival::Surv(data[[time]][indx],data[[status]][indx]),
              family="cox", keep=FALSE, grouped=TRUE,
              alpha=alpha, nfolds=min(10,nrows)),silent = TRUE))

            if(is.error(fit)){

              if(verbose) cat("x")
              nodes[[current_node]]$dont_split=TRUE

            } else {

              if(any(fit$glmnet.fit$df>0)){
                lambda.one<-suppressWarnings(
                  fit$glmnet.fit$lambda[min(which(fit$glmnet.fit$df>0))]
                )
              } else {
                lambda.one=NULL
              }
              
              if(length(fit$glmnet.fit$df)>=2){
                lambda.val<- suppressWarnings(
                  purrr::map_dbl((1:max(fit$glmnet.fit$df)),.f=function(x){
                    fit$glmnet.fit$lambda[min(which(fit$glmnet.fit$df>=x))]})
                )
                lambda.val = lambda.val[lambda.val>=fit$lambda.min]
                lambda.val = lambda.val[lambda.val<=fit$lambda.1se]
                
                if(!is.null(lambda.one)) lambda.val = c(lambda.one,lambda.val)
              } else {
                if(!is.null(lambda.one)){
                  lambda.val = c(lambda.one)
                } else {
                  class(fit)='try-error'
                }
              }
            }
          } else {

            if (is.null(maxdf_lincombo)) dfmax=5 else dfmax = maxdf_lincombo

            fit <- try(suppressWarnings(glmnet::glmnet(
              stats::model.matrix(
                ~.+1,data=data[indx,split_candidates,with=FALSE]),
              survival::Surv(data[[time]][indx],data[[status]][indx]),
              family="cox",
              alpha=alpha,dfmax=dfmax)))

            if(is.error(fit)){

              if(verbose) cat("x")
              nodes[[current_node]]$dont_split=TRUE

            } else {

              lambda.val<-suppressWarnings(
                purrr::map_dbl(1:dfmax,~fit$lambda[min(which(fit$df>=.))])
              )
              
            }

          }

          if(!is.error(fit)){

            if(verbose) cat('-')

            betas = as.matrix(
              purrr::reduce(purrr::map(lambda.val,~coef(fit,s=.)), cbind)
              )

            if(all(betas==0)){
              if(verbose) cat("0")
              nodes[[current_node]]$dont_split=TRUE
            } else {
              for(col in 1:ncol(betas)){
                if(all(betas[,col]==0)) betas=betas[,-col]
              }
            }

            if(is.null(ncol(betas))){

              beta=betas
              bvrs <- names(beta)[beta!=0]
              bwts <- beta[beta!=0]
              split_vars=list(dmat[indx,bvrs] %*% matrix(bwts,ncol=1))

            } else {

              split_vars <- purrr::map(1:ncol(betas),.f=function(col){
                beta=betas[,col]
                bvrs <- names(beta)[beta!=0]
                bwts <- beta[beta!=0]
                dmat[indx,bvrs] %*% matrix(bwts,ncol=1)
                
              })
            }

            # log-rank splits
            lrs_out = purrr::map(split_vars, function(cut_var){

              uvals=unique(cut_var)
              cut_pnts=soft_sample(uvals,size=nsplit)
              lrstats=purrr::map_dbl(cut_pnts,.f=function(cut_pnt){

                tmp = rep(1,nrow(cut_var))
                tmp[cut_var <= cut_pnt] = 0

                if(sum(tmp==0)<nmin | sum(tmp==1)<nmin) return(0)

                sdiff = try(survival::survdiff(
                  survival::Surv(data[[time]][indx],data[[status]][indx])~tmp),
                  silent=TRUE)

                if(is.error(sdiff)){

                  return(0)

                } else {

                  pv = 1-stats::pchisq(sdiff$chisq,df=1)

                  if(pv > 0.01){
                    return(0)
                  } else {
                    return(sdiff$chisq)
                  }

                }

              })

              c(cut_pnt=cut_pnts[which.max(lrstats)],lrank_val=max(lrstats))

            })

            lrs_out = purrr::reduce(lrs_out, rbind)

            if(is.null(dim(lrs_out))){

              all_zero=lrs_out['lrank_val']==0

              if(all_zero){
                if(verbose)cat('0')
                nodes[[current_node]]$dont_split=TRUE
              }

              if(is.null(dim(betas))){
                beta=betas
              } else {
                beta=betas[,1]
              }

              cut_pnt=lrs_out['cut_pnt']

            } else {

              all_zero=all(lrs_out[,'lrank_val']==0)

              if(all_zero){
                if(verbose)cat('0')
                nodes[[current_node]]$dont_split=TRUE
              }

              win_indx=which.max(lrs_out[,'lrank_val'])
              beta=betas[,win_indx]
              cut_pnt=lrs_out[win_indx,'cut_pnt']

            }

            nodes[[current_node]]$bwts=beta[beta!=0]
            nodes[[current_node]]$bvrs=names(nodes[[current_node]]$bwts)
            nodes[[current_node]]$cut_pnt=cut_pnt

            prds = dmat[,nodes[[current_node]]$bvrs] %*%
              matrix(nodes[[current_node]]$bwts,ncol=1)
            lwr  = prds <= cut_pnt
            upr  = prds >  cut_pnt

          } else {

            lwr=upr=0

          }

          if(sum(lwr)<nmin_leaf | sum(upr)<nmin_leaf |
             nodes[[current_node]]$dont_split){

            if(verbose) cat('M')
            nodes[[current_node]]$dont_split=TRUE

          } else {

            if(verbose)  cat('.')

            new_names=paste0(nodes[[current_node]]$name,0:1)
            node_ids[indx & lwr] = new_names[1]
            node_ids[indx & upr] = new_names[2]

            nodes[[new_names[1]]]=structure(
              list(name = new_names[1],
                   leaf = TRUE,
                   dont_split=FALSE,
                   parent = current_node,
                   children = NULL,
                   depth = nodes[[current_node]]$depth+1,
                   nsubs = sum(node_ids==new_names[1]),
                   nevents = sum(node_ids==new_names[1]&data[[status]]==1)),
              class = 'ost_node')

            nodes[[new_names[2]]]=structure(
              list(name = new_names[2],
                   leaf = TRUE,
                   dont_split=FALSE,
                   parent = current_node,
                   children = NULL,
                   depth = nodes[[current_node]]$depth+1,
                   nsubs = sum(node_ids==new_names[2]),
                   nevents = sum(node_ids==new_names[2]&data[[status]]==1)),
              class = 'ost_node')

            nodes[[current_node]]$leaf=FALSE
            nodes[[current_node]]$dont_split=TRUE
            nodes[[current_node]]$children=new_names

          }

          nodes_to_grow = purrr::map_lgl(nodes, function(x){
            x$leaf & !x$dont_split & x$nsubs>=minsplit & x$nevents>=nmin
          })

          nodes_to_grow = names(which(nodes_to_grow))

        } else {
          if(verbose) cat("n")
          nodes[[current_node]]$dont_split=TRUE
          nodes_to_grow=nodes_to_grow[-which(nodes_to_grow==current_node)]
        }

    }

  }

  leaves <- purrr::map_lgl(nodes, .f=~.$leaf)
  leaves <- names(which(leaves))

  for(current_leaf in leaves){

    srv <- try(
      survival::survfit(survival::Surv(
        data[[time]][node_ids==current_leaf],
        data[[status]][node_ids==current_leaf]) ~ 1,
        type = "kaplan-meier",se.fit=FALSE))

    if(!is.error(srv)){

      nodes[[current_leaf]]$times=srv$time
      nodes[[current_leaf]]$probs=srv$surv
      
      nodes[[current_leaf]]$nevents=sum(srv$n.event)
      nodes[[current_leaf]]$ncensrd=sum(srv$n.censor)

    } else {

      print(data[node_ids==current_leaf,])
      stop("Error in survfit")

    }

    nodes[[current_leaf]]$vals <-
      as.data.frame(purrr::map(.ftrs,.f=function(ftr){
        if(is.factor(data[[ftr]])){
          cat_mode(data[[ftr]][node_ids==current_leaf])
        } else {
          mean(data[[ftr]][node_ids==current_leaf])
        }
      }))

    names(nodes[[current_leaf]]$vals)<-.ftrs

  }

  #purrr::map_dbl(times,~tmp$srv.surv[which.max(tmp$srv.time[tmp$srv.time<=.])])

  structure(
    list(time          = time,
         status        = status,
         features      = .ftrs,
         nodes         = nodes,
         bootstrap_ids = boot_ids),
    class = "ost")

}

