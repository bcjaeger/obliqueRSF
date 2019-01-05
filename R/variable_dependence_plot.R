
#' Plot variable dependence using an oblique random survival forest
#' @param object an ORSF object (i.e. object returned from the ORSF function)
#' @param xvar a string giving the name of the x-axis variable
#' @param include.hist if true, a histogram showing the distribution of values for the x-axis variable will be included at the bottom of the plot.
#' @param include.points if true, the predictions for each observation are plotted along with a smoothed population estimate. Note that points are always included if xvar is categorical.
#' @param ptsize only relevant if include.points = TRUE. The size of the points in the plot are determined by this numeric value.
#' @param ytype String. Use 'event' if you would like to plot the probability of the event, and 'nonevent' if you prefer to plot the probability of a non-event.
#' @param event_lab string that describes the event
#' @param nonevent_lab string that describes a non-event.
#' @param fvar (optional) a string indicating a variable to facet the plot with
#' @param flab the labels to be printed describing the facet variable. For a facet variable with k categories, flab should be a vector with k labels, given in the same order as the levels of the facet variable.
#' @param time_units the unit of time, e.g. days, since baseline.
#' @param xlab the label to be printed describing the x-axis variable
#' @param xvar_units the unit of measurement for the x-axis variable. For example, age is usually measured in years.
#' @param xlvls a character vector giving the labels that correspond to categorical xvar. This does not need to be specified if xvar is continuous. 
#' @param sub_times the times you would like to plot predicted values for. If left unspecified, the ORSF function will use all of the times in oob_times.
#' @param se.show if true, standard errors of the population estimate will be included in the plot.
#' @return A ggplot2 object
#' @export
#' @examples
#'\dontrun{
#' data("pbc",package='survival')
#' pbc$status[pbc$status>=1]=pbc$status[pbc$status>=1]-1
#' pbc$time=pbc$time/365.25
#' pbc$id=NULL
#' fctrs<-c('trt','ascites','spiders','edema','hepato','stage')
#' for(f in fctrs)pbc[[f]]=as.factor(pbc[[f]])
#' pbc=na.omit(pbc)
#'
#' orsf=ORSF(data=pbc, eval_time=5, ntree=30)
#' 
#' vdplot(object=orsf, xvar='bili', xlab='Bilirubin', xvar_units='mg/dl')
#'}

# object=orsf
# xvar='sbp'
# include.hist=FALSE
# include.points=TRUE
# ytype='event'
# event_lab='CHD or stroke'
# nonevent_lab='survival'
# fvar=NULL
# flab=NULL
# time_units='year'
# xlab='Systolic blood pressure'
# xvar_units='mm Hg'
# xlvls=NULL
# sub_times=10
# se.show=FALSE
# ptsize=0.75

vdplot <- function(object,
                   xvar,
                   include.hist=TRUE,
                   include.points=FALSE,
                   ptsize=0.75,
                   ytype="nonevent",
                   event_lab="death",
                   nonevent_lab="survival",
                   fvar=NULL,
                   flab=NULL,
                   time_units="years",
                   xlab=xvar,
                   xvar_units=NULL,
                   xlvls=NULL,
                   sub_times=NULL,
                   se.show=FALSE){
  
  times=object$oob_times
  prds=object$oob_preds
  
  if(!is.null(sub_times)){
    prds=prds[,which(times%in%sub_times)]
    times=times[times%in%sub_times]
  }
  
  ntimes=length(times)
  nxvals=length(unique(object$data[[xvar]]))
  
  
  if(!is.null(fvar)){
    if(!is.factor(object$data[[fvar]])){
      stop("facet variable must be a factor") 
    }
  }
  
  if(is.null(flab) & !is.null(fvar)){
    flab=paste(fvar,'=',levels(object$data[[fvar]]))
  }
  
  if(nxvals<4 & !is.factor(object$data[[xvar]])){
    warning("xvar has < 4 unique values but is not a factor")
  }
  
  
  if(is.null(fvar)){
    ggdat=cbind(prds,object$data[[xvar]])
    colnames(ggdat)=c(paste0("t",1:ntimes),'xvar')
    ggdat=data.frame(ggdat)
    ggdat=tidyr::gather(ggdat,key='time',value='pred',-xvar)
  } else {
    ggdat=data.frame(cbind(prds,object$data[[xvar]],object$data[[fvar]]))
    names(ggdat)=c(paste0("t",1:ntimes),'xvar','fvar')
    ggdat=tidyr::gather(ggdat,key='time',value='pred',-xvar,-fvar)
    ggdat=dplyr::mutate(ggdat,fvar=factor(fvar,labels=flab))
  }
  
  ggdat=na.omit(ggdat)
  
  if(is.factor(object$data[[xvar]])){
    if(is.null(xlvls)) xlvls=levels(object$data[[xvar]])
    if(is.null(xlab)) xlab=xvar
    ggdat$xvar=factor(ggdat$xvar,labels=xlvls)
    if(is.null(xvar_units)) xvar_units=levels(object$data[[xvar]])
  } 
  
  if(ytype=='event'){
    ggdat$pred=1-ggdat$pred
  }
  
  if(ntimes>1){
    
    ggdat=mutate(ggdat,time=factor(time,
                                   levels=paste0("t",1:ntimes),
                                   labels=paste(times)))
    
    if(is.factor(object$data[[xvar]])){

      p=ggplot(ggdat,aes_string(x='time',y='pred',col='xvar'))+
        geom_point(position=position_jitterdodge(dodge.width = 1/2),
                   size=1/3, alpha=1/4, shape=15)+
        geom_pointrange(stat = "summary",
                        fun.ymin = function(z) {quantile(z,0.25)},
                        fun.ymax = function(z) {quantile(z,0.75)},
                        fun.y = median,
                        size=1,shape=16,
                        position=position_dodge(width=1/2))+
        scale_color_viridis_d(end=0.75)+
        labs(y=paste0("Probability of ",
                      ifelse(ytype=='event',event_lab,nonevent_lab)),
             x=paste0("Time since baseline, ",time_units),
             col=xlab)+
        theme_Publication()+
        scale_x_discrete(labels = 1:ntimes)+
        theme(legend.position = 'right',
              legend.direction = 'vertical')
      
    } else {
      
      time_units=paste(",",time_units)
      color_label = paste0("Time since \nbaseline",time_units)
      
      
      p=ggplot(ggdat,aes_string(x='xvar',y='pred',col='time'))+
        scale_color_viridis_d()
      
      if(include.hist&!is.factor(object$data[[xvar]])){
        ..density..=NULL
        p=p+geom_histogram(
          aes(x=xvar,y=scales::rescale(
            ..density..,to=c(0,max(min(prds),1/3)))),
          inherit.aes=F,fill='grey50',alpha=1/3,col='white')
      }
      
      if(include.points){
        p=p+geom_point(size=ptsize, color='grey50',alpha=1/6,
                       position=position_jitter())
      }
      
       p=p+geom_smooth(se=se.show,method='gam',formula=y~s(x),size=1.2)+
        labs(y=paste0("Probability of ",
                      ifelse(ytype=='event',event_lab,nonevent_lab)),
             x=paste0(xlab,', ',xvar_units),
             col=paste0("Time since \nbaseline",time_units))+
        theme_Publication()+
        theme(legend.position = 'right',
              legend.direction = 'vertical')
      
      if(ntimes>5) p = p+theme(legend.position='')
  }
    
  } else if (ntimes==1) {
    
    if(is.factor(object$data[[xvar]])){
      
      p=ggplot(ggdat,aes_string(x='xvar',y='pred'))+
        geom_boxplot()+
        labs(y=paste0(times,'-',time_units,
                      ifelse(ytype=='event',' risk',' probability'),
                      ' of ', ifelse(ytype=='event',event_lab,nonevent_lab)),
             x=paste(xlab))+
        scale_x_discrete(labels = xvar_units)+
        theme_Publication()
      
    } else {
     
      p=ggplot(ggdat,aes_string(x='xvar',y='pred'))
      
      if(include.hist&!is.factor(object$data[[xvar]])){
        ..density..=NULL
        p=p+geom_histogram(
          aes(x=xvar,y=scales::rescale(
            ..density..,to=c(0,max(min(prds),1/3)))),
          inherit.aes=F,fill='grey50',alpha=1/3,col='white')
      }
      
      if(include.points){
        p=p+geom_point(size=ptsize, color='grey50',alpha=1/6,
                       position=position_jitter())
      }
      
      p=p+stat_smooth(se=se.show,method='gam',formula=y~s(x),
                    size=1.4,col='red', fill='grey80')+
        labs(y=paste0(times,'-',time_units,
                      ifelse(ytype=='event',' risk',' probability'),
                      ' of ', ifelse(ytype=='event',event_lab,nonevent_lab)),
             x=paste0(xlab,', ', xvar_units))+
        theme_Publication()
       
    }
    
  }
  
  
  if(is.null(fvar)){
    p+scale_y_continuous(labels=scales::percent)
  } else {
    p+scale_y_continuous(labels=scales::percent)+facet_wrap(~fvar)
  }
  
}



