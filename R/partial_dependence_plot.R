
#' Plot partial variable dependence using an oblique random survival forest
#' @param object an ORSF object (i.e. object returned from the ORSF function)
#' @param xvar a string giving the name of the x-axis variable
#' @param xvals a vector containing the values that partial dependence will be computed with.
#' @param times a vector of times to compute predicted survival probabilities
#' @param nonevent_lab the label that describes a non-event.
#' @param fvar (optional) a string indicating a variable to facet the plot with
#' @param flab the labels to be printed describing the facet variable. For a facet variable with k categories, flab should be a vector with k labels, given in the same order as the levels of the facet variable.
#' @param data a dataframe to be used for predictions and x-values in the plot
#' @param time_units the unit of time, e.g. days, since baseline.
#' @param xlab the label to be printed describing the x-axis variable
#' @param xvar_units the unit of measurement for the x-axis variable. For example, age is usually measured in years.
#' @return A printed summary of the oblique random survival forest.
# @export
#' @examples

# object=orsf
# xvar='bili'
# xvals=NULL
# ytype=c("event","nonevent")
# event_lab='death'
# nonevent_lab='survival'
# fvar='sex'
# flab=NULL
# time_units='years'
# xlab=xvar
# xvar_units=NULL
# xlvls=NULL
# sub_times=c(3,5)
# 
# if(is.factor(object$data[[xvar]])){
#   if(is.null(xvals)){
#     xvals=levels(object$data[[xvar]])
#   }
# } else {
#   if(is.null(xvals)){
#     #xvals=sort(unique(object$data[[xvar]]))
#     xvals=quantile(object$data[[xvar]],probs=seq(0.10,0.90,length.out=9))
#   }
# }
# 
# times=object$oob_times
# prds=object$oob_preds
# 
# if(!is.null(sub_times)){
#   prds=prds[,which(times%in%sub_times)]
#   times=times[times%in%sub_times]
# }
# 
# nx = length(xvals)
# nt = length(times)
# indx=1:nt
# 
# if(!is.null(fvar)){
#   gg_fvar = levels(object$data[[fvar]])
# } else {
#   gg_fvar = 0
# }
# 
# ggdat=expand.grid(x=xvals,time=times,fvar=gg_fvar)%>%
#   dplyr::arrange(x)
# 
# if(is.null(fvar))ggdat$fvar=NULL
# 
# ggdat$pred_mean=ggdat$pred_sd=0
# 
# if(is.null(fvar)){
#   for(i in 1:nx){
#     object$data[[xvar]]=xvals[i]
#     prds=predict(object,newdata=object$data,times=times)
#     ggdat[indx,'pred_mean']=apply(prds,2,mean)
#     ggdat[indx,'pred_sd']=apply(prds,2,sd)
#     indx=indx+nt
#   }
# } else {
#   for(i in 1:nx){
#     object$data[[xvar]]=xvals[i]
#     for(f in gg_fvar){
#       object$data[[fvar]]=f
#       object$data[[fvar]][1]=setdiff(gg_fvar,f)[1]
#       prds=predict(object,newdata=object$data,times=times)
#       ggdat[indx,'pred_mean']=apply(prds,2,mean)
#       ggdat[indx,'pred_sd']=apply(prds,2,sd)
#       indx=indx+nt
#     }
#   }
# }
# 
# 
# 
# 
# if(ytype[1]=='event'){
#   ggdat$pred_mean=1-ggdat$pred_mean
# }
# 
# ggdat <- mutate(ggdat,
#                 pred_se=pred_sd/sqrt(nrow(object$data)),
#                 ymin=pred_mean-1.96*pred_se,
#                 ymax=pred_mean+1.96*pred_se)
# 
# if(ntimes>1){
#   # Plot will display predictions from >1 time point
#   if(is.factor(object$data[[xvar]])){
#     # Categorical variable
#   } else {
#     # Continuous variable
#     color_label = paste0("Time since \nbaseline, ",time_units)
#     p=ggplot(ggdat,aes(x=x,y=pred_mean,col=factor(time)))+
#       scale_color_viridis_d()+
#       geom_line(size=1.2)+
#       labs(y=paste0("Probability of ",
#                     ifelse(ytype[1]=='event',event_lab,nonevent_lab)),
#            x=ifelse(is.null(xvar_units),xlab,paste0(xlab,', ',xvar_units)),
#            col=paste0("Time since \nbaseline",time_units))+
#       theme_Publication()+
#       theme(legend.position = 'right',
#             legend.direction = 'vertical')
#   }
# 
# } else if(ntimes==1){
#   # Plot will display predictions from 1 time point
#   if(is.factor(object$data[[xvar]])){
#     # Categorical variable
#   } else {
#     # Continuous variable
#     p=ggplot(ggdat,aes(x=x, y=pred_mean, ymin=ymin, ymax=ymax, group=time))+
#       geom_errorbar(col='grey60',width=0.1)+
#       geom_line(col='grey60')+
#       geom_point(col='red',size=2)+
#       labs(y=paste0("Probability of ",
#                     ifelse(ytype[1]=='event',event_lab,nonevent_lab)),
#            x=ifelse(is.null(xvar_units),xlab,paste0(xlab,', ',xvar_units)))+
#       theme_Publication()
#   }
# 
# }
# 
# if(is.null(fvar)){
#   p+scale_y_continuous(labels=scales::percent)
# } else {
#   p+scale_y_continuous(labels=scales::percent)+facet_wrap(~fvar)
# }

  
