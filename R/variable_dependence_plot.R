
#' Plot variable dependence using an oblique random survival forest
#' @param object an ORSF object (i.e. object returned from the ORSF function)
#' @param xvar a string giving the name of the x-axis variable
#' @param times a vector of times to compute predicted survival probabilities
#' @param include.hist if true, a histogram showing the distribution of values for the x-axis variable will be included at the bottom of the plot.
#' @param nonevent_lab the label that describes a non-event.
#' @param fvar (optional) a string indicating a variable to facet the plot with
#' @param flab the labels to be printed describing the facet variable. For a facet variable with k categories, flab should be a vector with k labels, given in the same order as the levels of the facet variable.
#' @param data a dataframe to be used for predictions and x-values in the plot
#' @param time_units the unit of time, e.g. days, since baseline.
#' @param xlab the label to be printed describing the x-axis variable
#' @param xvar_units the unit of measurement for the x-axis variable. For example, age is usually measured in years.
#' @return A ggplot2 object
#' @export
#' @examples
#'
#' data("pbc",package='survival')
#' pbc$status[pbc$status>=1]=pbc$status[pbc$status>=1]-1
#' pbc$time=pbc$time/365.25
#' pbc$id=NULL
#' fctrs<-c('trt','ascites','spiders','edema','hepato','stage')
#' for(f in fctrs)pbc[[f]]=as.factor(pbc[[f]])
#' pbc=na.omit(pbc)
#'
#' orsf=ORSF(data=pbc)
#' prds=predict(orsf,newdata=pbc,times=3)
#' 
#' vdplot(object=prds,
#'        xvar='copper',
#'        times=3,
#'        data=pbc,
#'        fvar = 'ascites',
#'        flab = c('No','Yes'),
#'        time_units='years',
#'        xlab='serum bilirubin',
#'        xvar_units = '(mg/dl)',
#'        nonevent_lab="survival",
#'        include.hist=TRUE)
#'


# object=orsf
# xvar='sex'
# times=1:10
# data=pbc
# fvar = 'ascites'
# flab = NULL
# time_units='years'
# xlab=xvar
# xvar_units = NULL
# nonevent_lab="survival"
# include.hist=TRUE

vdplot <- function(object,xvar,times,data,include.hist=TRUE,
                   nonevent_lab="survival",
                   fvar=NULL,flab=NULL,time_units="years",
                   xlab=xvar,xvar_units=NULL,xlvls=NULL){
  
  ntimes=length(times)
  nxvals=length(unique(data[[xvar]]))
  
  if(!is.null(fvar)){
    if(!is.factor(data[[fvar]])){
      stop("facet variable must be a factor") 
    }
  }
  
  if(is.null(flab) & !is.null(fvar)){
    flab=paste(fvar,'=',levels(data[[fvar]]))
  }
  
  if(nxvals<4 & !is.factor(data[[xvar]])){
    warning("xvar has < 4 unique values but is not a factor")
  }
  
  if(class(object)[1]=='orsf'){
    prds=predictSurvProb(object,newdata=data,times=times)
  } else {
    prds=object
  }
  
  if(is.null(fvar)){
    ggdat=cbind(prds,data[[xvar]])
    colnames(ggdat)=c(letters[1:ntimes],'xvar')
    ggdat=data.frame(ggdat)
    ggdat=tidyr::gather(ggdat,key='time',value='pred',-xvar)
  } else {
    ggdat=data.frame(cbind(prds,data[[xvar]],data[[fvar]]))
    names(ggdat)=c(letters[1:ntimes],'xvar','fvar')
    ggdat=tidyr::gather(ggdat,key='time',value='pred',-xvar,-fvar)
    ggdat=dplyr::mutate(ggdat,fvar=factor(fvar,labels=flab))
  }
  
  ggdat=na.omit(ggdat)
  
  if(is.factor(data[[xvar]])){
    if(is.null(xlab)) xlab=levels(data[[xvar]])
    ggdat$xvar=factor(ggdat$xvar,labels=xlab)
    if(is.null(xvar_units)) xvar_units=levels(data[[xvar]])
  } 
  
  if(ntimes>1){
    
    ggdat=mutate(ggdat,time=factor(time,
                                   levels=letters[1:ntimes],
                                   labels=paste(times)))
    
    if(is.factor(data[[xvar]])){

      p=ggplot(ggdat,aes_string(x='time',y='pred',col='xvar'))+
        geom_point(position=position_jitterdodge(dodge.width = 1/2),
                   size=3/4, alpha=1/3)+
        stat_summary(fun.data = mean_cl_boot,size=1, 
                     position=position_dodge(width=1/2))+
        scale_color_viridis_d(end=0.75)+
        labs(y=paste0("Probability of ",nonevent_lab),
             x=paste0("Time since baseline, ",time_units),
             col=xvar)+
        theme_Publication()+
        scale_x_discrete(labels = 1:ntimes)+
        theme(legend.position = 'right',
              legend.direction = 'vertical')
      
    } else {
      
      time_units=paste(",",time_units)
      color_label = paste0("Time since \nbaseline",time_units)
      
      
      p=ggplot(ggdat,aes_string(x='xvar',y='pred',col='time'))+
        scale_color_viridis_d()+
        geom_smooth(se=F,method='gam',formula=y~s(x),size=1.2)+
        labs(y=paste0("Probability of ",nonevent_lab),
             x=paste(xlab,xvar_units,collapse=', '),
             col=paste0("Time since \nbaseline",time_units))+
        theme_Publication()+
        theme(legend.position = 'right',
              legend.direction = 'vertical')
  }
    
  } else if (ntimes==1) {
    
    if(is.factor(data[[xvar]])){
      
      p=ggplot(ggdat,aes_string(x='xvar',y='pred'))+
        geom_boxplot()+
        labs(y=paste0("Probability of ",nonevent_lab,
                      " at ",times," ",time_units),
             x=paste(xlab))+
        scale_x_discrete(labels = xvar_units)+
        theme_Publication()
      
    } else {
     
      p=ggplot(ggdat,aes_string(x='xvar',y='pred'))+
        geom_smooth(se=T,method='gam',formula=y~s(x),size=1.2,col='red')+
        labs(y=paste0("Probability of ",nonevent_lab,
                      " at ",times," ",time_units),
             x=paste0(xlab,', ', xvar_units))+
        theme_Publication()
       
    }
    
  }
  
  if(include.hist&!is.factor(data[[xvar]])){
    ..density..=NULL
    p=p+geom_histogram(
      aes(x=xvar,y=scales::rescale(
        ..density..,to=c(0,max(min(prds),1/3)))),
      inherit.aes=F,fill='grey50',alpha=1/3,col='white')
  }
  
  return(if(is.null(fvar)) p else p+facet_wrap(~fvar))
  
}


theme_Publication <- function(base_size=16) {
  
  (theme_foundation(base_size=base_size)+ 
     theme(plot.title = element_text(face = "bold",
                                     size = rel(1), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.border = element_rect(colour = 'black'),
           legend.key.size = unit(3,"line"),
           legend.key.height = unit(3,"line"),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size=rel(1)),
           axis.line = element_blank(), 
           axis.ticks = element_line(),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.title = element_text(face="italic"),
           legend.text = element_text(size=rel(1)),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="black",fill="#f0f0f0"),
           strip.text = element_text(size=rel(1), face="bold")
     ))
  
}

