
#' Oblique Random Survival Forests
#' 
#' Oblique random survival forest are ensembles for right-censored survival data that incorporate linear combinations of input variables into random survival forests (see Ishwaran et al., 2008 <doi:10.1214/08-AOAS169>). Regularized Cox proportional hazard models (see Simon et al., 2016 <doi:10.18637/jss.v039.i05>) identify optimal linear combinations of input variables in each recursive partitioning step while building survival trees (see Bou-hamad et al., 2011 <doi: 10.1214/09-SS047>).
#' 
#' @docType package
#' @author Byron C. Jaeger <bcjaeger@uab.edu>
#' @import Rcpp pec data.table stats missForest glmnet survival prodlim tidyr ggplot2
#' @importFrom purrr map reduce map_dbl
#' @importFrom dplyr arrange mutate
#' @importFrom rlang sym
#' @importFrom scales rescale
#' @importFrom ggthemes theme_foundation
#' @importFrom Rcpp evalCpp
#' @useDynLib obliqueRSF
#' @name obliqueRSF
NULL 