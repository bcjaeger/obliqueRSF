
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
# data=na.omit(pbc)
# xvar='sex'
# xvals=NULL
# times=c(10)
# include.hist=TRUE
# nonevent_lab='survival'
# fvar=NULL
# flab=NULL
# time_units='years'
# xlab='Serum bilirubin'
# xvar_units='mg/ul'

# if(is.factor(data[[xvar]])){
#   
#   if(is.null(xvals)){
#     xvals=levels(data[[xvar]])
#   }
#   
# } else {
#  
#   if(is.null(xvals)){
#     xvals=sort(unique(data[[xvar]]))
#   }
#    
# }
# 
# nx = length(xvals)
# nt = length(times)
# indx=1:nt
# 
# ggdat=expand.grid(x=xvals,time=times)%>%dplyr::arrange(x)
# ggdat$pred_mean=ggdat$pred_sd=0
# 
# for(i in 1:nx){
#   
#   data[[xvar]]=xvals[i]
#   prds=predict(object,newdata=data,times=times)
#   ggdat[indx,'pred_mean']=apply(prds,2,mean)
#   ggdat[indx,'pred_sd']=apply(prds,2,sd)
#   indx=indx+nt
#   
# }
# 
# ggdat %>% 
#   mutate(pred_se=pred_sd/sqrt(nrow(data)),
#          ymin=pred_mean-1.96*pred_se,
#          ymax=pred_mean+1.96*pred_se) %>%
#   ggplot(aes(x=x, y=pred_mean, ymin=ymin, ymax=ymax, group=time))+
#   geom_errorbar(col='grey60',width=0.1)+
#   geom_line(col='grey60')+
#   geom_point(col='red',size=2)+
#   theme_Publication()
