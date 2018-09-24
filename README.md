
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

Oblique random survival forest (ORSFs) are ensembles for right-censured survival data that use linear combinations of input variables to recursively partition a set of training data. Regularized Cox proportional hazard models identify optimal linear combinations of input variables in each recursive partitioning step while building survival trees.

Installation
============

You can install obliqueRSF from github with:

``` r
# install.packages("devtools")
devtools::install_github("bcjaeger/obliqueRSF")
```

Usage
=====

The `ORSF` function is the center piece of the `obliqueRSF` package

``` r


data("pbc",package='survival')

# event is death
# censored at time of last contact or transplant
pbc$status[pbc$status>=1]=pbc$status[pbc$status>=1]-1

# format categorical variables as factors
pbc = pbc %>% 
  dplyr::select(-id)%>%
  mutate(trt=factor(trt),
         ascites=factor(ascites),
         hepato=factor(hepato),
         spiders=factor(spiders),
         edema=factor(edema),
         stage=factor(stage,ordered=TRUE))

set.seed(62689)
orsf=ORSF(data=pbc, ntree=100, alpha=0.50,verbose=FALSE)
#> 
#> performing imputation with missForest:
#>   missForest iteration 1 in progress...done!
#>   missForest iteration 2 in progress...done!
#>   missForest iteration 3 in progress...done!
#>   missForest iteration 4 in progress...done!
#>   missForest iteration 5 in progress...done!
#>   missForest iteration 6 in progress...done!
#>   missForest iteration 7 in progress...done!

print(orsf)
#> 
#> Oblique Random Survival Forest: ORSF(data = pbc, alpha = 0.5, ntree = 100, verbose = FALSE)
#> $concordance
#> [1] 0.7948444
#> 
#> $integrated_briscr
#> 
#> Integrated Brier score (crps):
#> 
#>           IBS[0;time=4795)
#> Reference            0.174
#> ORSF                 0.116
nboots=25
```

Performance
===========

`orsf` objects take a long time to grow, but they usually provide excellent predictions. Below we use the `pec` package to compare the integrated Brier scores over 25 replications of bootstrap cross validation using the complete cases in the `pbc` data set.

``` r

pbc_cmp=na.omit(pbc)

ntree=300
nsize=20

mdls=list(
  orsf=ORSF(data=pbc_cmp, ntree=ntree, use.cv = FALSE, 
            min_obs_to_split_node = nsize, verbose=FALSE),
  rsf =rfsrc(Surv(time,status)~.,data=pbc_cmp,
             ntree=ntree,nodesize=nsize),
  bws = pec::selectCox(Surv(time,status)~., data=pbc_cmp),
  cboost = pec::coxboost(Hist(time,status)~.,
                       data=pbc_cmp,penalty=1000)
)

# for reproducibility
set.seed(32989) 

bri_score = pec::pec(mdls, data=pbc_cmp, cens.model = 'cox',
                     formula = Surv(time,status)~age, 
                     splitMethod = 'BootCv', B=nboots)
#> Split sample loop (B=25)
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> Warning in fitter(X, Y, strata = Strata, offset = offset, weights =
#> weights, : Ran out of iterations and did not converge
#> 10
#> Warning in fitter(X, Y, strata = Strata, offset = offset, weights =
#> weights, : Ran out of iterations and did not converge
#> 20

print(bri_score)
#> 
#> Prediction error curves
#> 
#> Prediction models:
#> 
#> Reference      orsf       rsf       bws    cboost 
#> Reference      orsf       rsf       bws    cboost 
#> 
#> Right-censored response of a survival model
#> 
#> No.Observations: 276 
#> 
#> Pattern:
#>                 Freq
#>  event          111 
#>  right.censored 165 
#> 
#> IPCW: cox model
#> 
#> Method for estimating the prediction error:
#> 
#> Bootstrap cross-validation
#> 
#> Type: resampling
#> Bootstrap sample size:  276 
#> No. bootstrap samples:  25 
#> Sample size:  276 
#> 
#> Cumulative prediction error, aka Integrated Brier score  (IBS)
#>  aka Cumulative rank probability score
#> 
#> Range of integration: 0 and time=4509 :
#> 
#> 
#> Integrated Brier score (crps):
#> 
#>           IBS[0;time=4509)
#> Reference            0.189
#> orsf                 0.133
#> rsf                  0.134
#> bws                  0.147
#> cboost               0.134


cnc_index = pec::cindex(mdls, data=pbc_cmp, cens.model = 'cox',
                        formula = Surv(time,status)~age, 
                        splitMethod = 'BootCv', B=nboots)
#> 10
#> Warning in fitter(X, Y, strata = Strata, offset = offset, weights =
#> weights, : Ran out of iterations and did not converge
#> 20

print(cnc_index)
#> 
#> The c-index for right censored event times
#> 
#> Prediction models:
#> 
#>   orsf    rsf    bws cboost 
#>   orsf    rsf    bws cboost 
#> 
#> Right-censored response of a survival model
#> 
#> No.Observations: 276 
#> 
#> Pattern:
#>                 Freq
#>  event          111 
#>  right.censored 165 
#> 
#> Censoring model for IPCW: cox model  
#> 
#> Method for estimating the prediction error:
#> 
#> Bootstrap cross-validation
#> 
#> Type: resampling
#> Bootstrap sample size:  276 
#> No. bootstrap samples:  25 
#> Sample size:  276 
#> 
#> Estimated C-index in % at time=4191 
#> 
#>        AppCindex BootCvCindex
#> orsf        82.5         81.4
#> rsf         88.1         80.0
#> bws         75.5         76.2
#> cboost      79.3         79.7
#> 
#> AppCindex    : Apparent (training data) performance
#> BootCvCindex : Bootstrap crossvalidated performance
```
