
#' Plot partial variable dependence using an oblique random survival forest
#' @param object an ORSF object (i.e. object returned from the ORSF function)
#' @param xvar a string giving the name of the x-axis variable
#' @param xlab the label to be printed describing the x-axis variable
#' @param xvar_units the unit of measurement for the x-axis variable. For example, age is usually measured in years.
#' @param xvals a vector containing the values that partial dependence will be computed with.
#' @param nxpts instead of specifying xvals, you can specify how many points on the x-axis you would like to plot predicted responses for, and a set of nxpts equally spaced percentile values from the distribution of xvar will be used.
#' @param ytype String. Use 'event' if you would like to plot the probability of the event, and  'nonevent' if you prefer to plot the probability of a non-event.
#' @param event_lab string that describes the event
#' @param nonevent_lab string that describes a non-event.
#' @param fvar a string indicating a variable to facet the plot with
#' @param flab a label describing the facet variable.
#' @param flvls the labels to be printed describing the facet variable. For a facet variable with k categories, flab should be a vector with k labels, given in the same order as the levels of the facet variable.
#' @param time_units the unit of time, e.g. days, since baseline.
#' @param xlvls A character vector with descriptions of each category in the x-variable. This is only relevant if x is categorical.
#' @param sub_times a vector of times to compute predicted survival probabilities. Note that the eval_times from the ORSF object are used to compute predictions, and sub_times must be a subset of those times.
#' @param separate_panels true or false. If true, the plot will display predictions in two separate panels, determined by the facet variable.
#' @param color_palette Palette to use for colors in the figure. Options are Diverging (BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral), Qualitative (Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3), Sequential (Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd), and viridis.
#' @return A ggplot2 object showing partial dependence according to the oblique random survival forest object.
#' @export
#' @examples
#' data("pbc",package='survival')
#' pbc$status[pbc$status>=1]=pbc$status[pbc$status>=1]-1
#' pbc$time=pbc$time/365.25
#' pbc$id=NULL
#' fctrs<-c('trt','ascites','spiders','edema','hepato','stage')
#' for(f in fctrs)pbc[[f]]=as.factor(pbc[[f]])
#' pbc=na.omit(pbc)
#'
#' orsf=ORSF(data=pbc, eval_time=1:10)
#' 
#' pdplot(object=orsf, xvar='bili', xlab='Bilirubin', 
#'        xvar_units='mg/dl', sub_times=10)
#'

pdplot <- function(object,
                   xvar,
                   xlab=NULL,
                   xvar_units=NULL,
                   xvals=NULL,
                   nxpts=10,
                   ytype=c("nonevent","event"),
                   event_lab='death',
                   nonevent_lab='survival',
                   fvar=NULL,
                   flab=NULL,
                   flvls=NULL,
                   time_units='years',
                   xlvls=NULL,
                   sub_times=NULL,
                   separate_panels=TRUE,
                   color_palette='Dark2'){
  
  x=pred_sd=pred_mean=pred_se=ymin=ymax=NULL
  
  if(!is.null(sub_times)){
    if(!all(sub_times%in%object$oob_times)){
      stop("sub_times must be a subset of eval_times (specified in ORSF call)")
    }
  }
  
  xvar_is_factor <- is.factor(object$data[[xvar]])
  if(is.null(xlab)) xlab = xvar
  
  if(xvar_is_factor){
    if(is.null(xvals)){
      xvals=levels(object$data[[xvar]])
    }
  } else {
    if(is.null(xvals)){
      #xvals=sort(unique(object$data[[xvar]]))
      xvals=unique(quantile(
        object$data[[xvar]],probs=seq(0.10,0.90,length.out=nxpts)))
    }
  }
  
  times=object$oob_times
  prds=object$oob_preds
  
  if(!is.null(sub_times)){
    prds=prds[,which(times%in%sub_times)]
    times=times[times%in%sub_times]
  }
  
  nx = length(xvals)
  nt = length(times)
  indx=1:nt
  
  if(!is.null(fvar)){
    gg_fvar = levels(object$data[[fvar]])
  } else {
    gg_fvar = 0
  }
  
  ggdat=expand.grid(x=xvals,fvar=gg_fvar,time=times,pred_mean=0,pred_sd=0)
  
  if(is.null(fvar)){
    ggdat = dplyr::arrange(ggdat, x)
    ggdat$fvar=NULL
  } else {
    ggdat = dplyr::arrange(ggdat,x,fvar)
    if(is.null(flvls)) flvls=paste(fvar,'=',levels(object$data[[fvar]]))
  }
  
  mat=model.matrix(~., data=object$data)
  
  if(is.null(fvar)){
    for(i in 1:nx){
      if(xvar_is_factor){
        if(i==1){
          mat[,paste0(xvar,xvals[-1])]=0
        } else {
          mat[,paste0(xvar,xvals[i])]=1
        }
      } else {
        mat[,xvar]=xvals[i]
      }
      prds=predict(object,newdata=mat,times=times)
      ggdat[indx,'pred_mean']=apply(prds,2,mean)
      ggdat[indx,'pred_sd']=apply(prds,2,sd)
      indx=indx+nt
    }
  } else {
    for(i in 1:nx){
      if(xvar_is_factor){
        if(i==1){
          mat[,paste0(xvar,xvals[-1])]=0
        } else {
          mat[,paste0(xvar,xvals[i])]=1
        }
      } else {
        mat[,xvar]=xvals[i]
      }
      for(f in gg_fvar){
        if(f==gg_fvar[1]){
          mat[,paste0(fvar,gg_fvar[-1])]=0
        } else {
          mat[,paste0(fvar,f)]=1
        }
        prds=predict(object,newdata=mat,times=times)
        ggdat[indx,'pred_mean']=apply(prds,2,mean)
        ggdat[indx,'pred_sd']=apply(prds,2,sd)
        indx=indx+nt
      }
    }
  }
  
  if(ytype[1]=='event'){
    ggdat$pred_mean=1-ggdat$pred_mean
  }
  
  ggdat <- mutate(ggdat,
                  pred_se=pred_sd/sqrt(nrow(mat)),
                  ymin=pred_mean-1.96*pred_se,
                  ymax=pred_mean+1.96*pred_se)
  
  if(!is.null(fvar)){
    ggdat$fvar=factor(ggdat$fvar,levels=gg_fvar,labels=flvls)
  }
  
  if(nt>1){
    # Plot will display predictions from >1 time point
    if(xvar_is_factor){
      # Categorical variable
      dwidth=max(times)/nt
      
      if(!is.null(xlvls)) ggdat$x=factor(ggdat$x,labels=xlvls)

      p=ggplot(ggdat,aes(x=time,y=pred_mean,fill=x,
                         label=paste0(100*round(pred_mean,2),'%')))+
        geom_bar(stat='identity',position='dodge')+
        geom_text(aes(col=x), vjust=-1/2,position=position_dodge(width=dwidth))+
        labs(y=paste0("Probability of ",
                      ifelse(ytype[1]=='event',event_lab,nonevent_lab)),
             x=ifelse(is.null(time_units),'Time since baseline',
                      paste0("Time since baseline, ",time_units)),
             fill=xlab)+
        theme_Publication()+guides(col=FALSE)+
        scale_x_continuous(breaks=times)+
        theme(legend.position = 'right',
              legend.direction = 'vertical')
      
      if(color_palette=='viridis'){
        p=p+scale_color_viridis_d()+
          scale_fill_viridis_d()
      } else {
        p=p+scale_color_brewer(palette=color_palette)+
          scale_fill_brewer(palette=color_palette)
      }
      
      
    } else {
      # Continuous variable
      p=ggplot(ggdat,aes(x=x,y=pred_mean,col=factor(time)))+
        geom_line(size=1.2)+
        labs(y=paste0("Probability of ",
                      ifelse(ytype[1]=='event',event_lab,nonevent_lab),
                      ifelse(is.null(time_units),'',
                             paste('',times, time_units, 'after baseline'))),
             x=ifelse(is.null(xvar_units),xlab,paste0(xlab,', ',xvar_units)),
             col=ifelse(is.null(time_units),
                        'Time since \nbaseline',
                        paste0("Time since \nbaseline, ",time_units)))+
        theme_Publication()+
        theme(legend.position = 'right',
              legend.direction = 'vertical')
      
      if(nt>8 | color_palette=='viridis'){
        p=p+scale_color_viridis_d()
      } else {
        p=p+scale_color_brewer(palette=color_palette)
      }
      
    }
    
  } else if(nt==1) {
    # Plot will display predictions from 1 time point
    if(xvar_is_factor){
      # Categorical variable
      dwidth=0.9
      if(!is.null(xlvls)){
        ggdat$x=factor(ggdat$x,levels=xvals,labels=xlvls)
      }
      
      if(separate_panels){
        warning("Must use > 1 time if you prefer to have separate panels")
        separate_panels=FALSE
      }
      
      p=ggplot(ggdat,aes(x=x,y=pred_mean,fill=fvar,
                         label=paste0(100*round(pred_mean,2),'%')))+
        geom_bar(stat='identity',position='dodge')+
        geom_text(aes(col=fvar), vjust=-1/2, size=8,
                  position=position_dodge(width=dwidth))+
        labs(y=paste0("Probability of ",
                      ifelse(ytype[1]=='event',event_lab,nonevent_lab),
                      ifelse(is.null(time_units),'',
                             paste('',times, time_units, 'after baseline'))),
             x=xlab,
             fill=ifelse(is.null(fvar),fvar,ifelse(is.null(flab),fvar,flab)))+
        theme_Publication()+guides(col=FALSE)+
        coord_cartesian(ylim=c(0,max(ggdat$pred_mean*1.15)))+
        theme(legend.position = 'right',
              legend.direction = 'vertical')
      
      if(color_palette=='viridis'){
        p=p+scale_color_viridis_d()+
          scale_fill_viridis_d()
      } else {
        p=p+scale_color_brewer(palette=color_palette)+
          scale_fill_brewer(palette=color_palette)
      }
      
    } else {
      # Continuous variable
      if(separate_panels){
        p=ggplot(ggdat,aes(x=x, y=pred_mean, ymin=ymin, ymax=ymax, group=time))+
          geom_errorbar(col='grey60',width=0.1)+
          geom_line(col='grey60')+
          geom_point(col='red',size=2)+
          labs(y=paste0("Probability of ",
                        ifelse(ytype[1]=='event',event_lab,nonevent_lab),
                        ifelse(is.null(time_units),'',
                               paste('',times, time_units, 'after baseline'))),
               x=ifelse(is.null(xvar_units),xlab,paste0(xlab,', ',xvar_units)))+
          theme_Publication()
      } else {
        p=ggplot(ggdat,aes(x=x, y=pred_mean, ymin=ymin, ymax=ymax, 
                           group=fvar, col=fvar))+
          geom_errorbar(width=0.1,
                        position=position_dodge2(width=0.1))+
          geom_line(position=position_dodge2(width=0.1))+
          geom_point(size=2,
                     position=position_dodge2(width=0.1))+
          labs(y=paste0("Probability of ",
                        ifelse(ytype[1]=='event',event_lab,nonevent_lab),
                        ifelse(is.null(time_units),'',
                               paste('',times, time_units, 'after baseline'))),
               x=ifelse(is.null(xvar_units),xlab,paste0(xlab,', ',xvar_units)),
               col=ifelse(is.null(flab),'',flab))+
          theme_Publication()
        
        if(color_palette=='viridis'){
          p=p+scale_color_viridis_d()
        } else {
          p=p+scale_color_brewer(palette=color_palette)
        }
        
      }
      
    }
    
  }
  
  if(is.null(fvar) | !separate_panels){
    print(p+scale_y_continuous(labels=scales::percent))
  } else if(separate_panels){
    print(p+scale_y_continuous(labels=scales::percent)+
            facet_wrap(~fvar))
  } 
  
}



  
# object=orsf
# xvar='bili'
# xlab="Participant bilirubin"
# xvar_units=NULL
# xvals=NULL
# nxpts=8
# ytype=c("event","nonevent"); ytype='nonevent'
# event_lab='death'
# nonevent_lab='survival'
# fvar='spiders'
# flab=NULL
# flvls=NULL
# time_units='years'
# xlvls=NULL
# sub_times=8
# separate_panels=F
# color_palette='Dark2'