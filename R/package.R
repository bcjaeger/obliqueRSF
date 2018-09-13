
#' Oblique Random Survival Forests
#' 
#' Oblique random survival forest (ORSFs) are ensembles for right-censured survival data that use linear combinations of input variables to recursively partition a set of training data. Regularized Cox proportional hazard models identify optimal linear combinations of input variables in each recursive partitioning step while building survival trees.
#' 
#' @docType package
#' @author Byron C. Jaeger <bcjaeger@uab.edu>
#' @import Rcpp pec data.table stats missForest glmnet survival prodlim
#' @importFrom purrr map reduce map_dbl
#' @importFrom dplyr arrange 
#' @importFrom rlang sym
#' @importFrom Rcpp evalCpp
#' @useDynLib obliqueRSF
#' @name obliqueRSF
NULL 