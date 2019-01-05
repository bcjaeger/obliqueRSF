
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
# censor study participants at time of last contact or transplant
pbc$status[pbc$status>=1]=pbc$status[pbc$status>=1]-1

# format categorical variables as factors
# if we don't do this, missforest will not impute 0/1 variables 
# in the way that we would like it to.
pbc = pbc %>% 
  dplyr::select(-id)%>%
  mutate(trt=factor(trt),
         ascites=factor(ascites),
         hepato=factor(hepato),
         spiders=factor(spiders),
         edema=factor(edema),
         stage=factor(stage,ordered=TRUE),
         time=time/365.25) 

orsf=ORSF(data=pbc, # data to fit trees with
          ntree=200, # number of trees to fit
          eval_times=c(1:10), # when will predictions be made?
          # note: eval_times will be used to make figures
          verbose=F) # suppresses console output
#> 
#> performing imputation with missForest:
#>   missForest iteration 1 in progress...done!
#>   missForest iteration 2 in progress...done!
#>   missForest iteration 3 in progress...done!
#>   missForest iteration 4 in progress...done!
#>   missForest iteration 5 in progress...done!
```

The `vdplot` function allows you to quickly look at the patterns in the predictions from an ORSF.

``` r
# Variable dependence plot (vdplot)

# Survival probabilities for a continuous variable
# note the use of sub_times, which allows you to pick
# one or more times in the evaluation times of an ORSF object

vdplot(object=orsf, xvar='bili', xlab='Bilirubin levels', 
       xvar_units = 'mg/dl', sub_times = 5)
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](man/figures/vdplot1-1.png)

here is another application of `vdplot` with a continuous x-variable, but this time we will show predictions at three times: 1 year, 3 years, and 5 years since baseline.

``` r
vdplot(object=orsf, xvar='albumin', xlab='Serum albumin', 
       xvar_units = 'g/dl', sub_times = c(1,3,5))
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](man/figures/vdplot2-1.png)

The `vdplot` function also supports categorical x-variables. Setting the x-variable to `sex`, and the facet variable to `hepato`, we find a clear interaction between sex and the presence of an enlarged liver (i.e., `hepato=1`).

``` r
vdplot(object=orsf, xvar='sex', xlab=c("Sex"), xlvls=c("Male","Female"),
       fvar='hepato',flab=c("Normal size liver","Enlarged liver"))
```

![](man/figures/vdplot3-1.png)

Is this interaction something that is purely explained by sex, or is it a confounding effect from other variables that are not the same between the two sexes? We can address this using the partial dependence plot (`pdplot`) function:

``` r
pdplot(object=orsf, xvar='sex',xlab='Sex', xlvls=c("Male","Female"),
       fvar='hepato',flvls=c("Normal size liver","Enlarged liver"),
       sub_times=c(1,3,5,7,9))
```

![](man/figures/pdplot1-1.png)

Taking into account the effects of other variables in the data, the interaction between sex and hepato is attenuated.

Performance
===========

`orsf` objects take a long time to grow, but they usually provide excellent predictions. Below we use the `pec` package to compare the integrated Brier scores over 5 replications of bootstrap cross validation using the complete cases in the `pbc` data set.

``` r

# it takes a bit more code to incorporate imputation 
# into bootstrap cross-validation. However, results
# are similar using only the complete observations,
# so that's what we will use here

pbc_cmp=na.omit(pbc)

ntree=150
nsize=20

mdls=list(
  orsf=ORSF(data=pbc_cmp, ntree=ntree, use.cv = FALSE, 
            min_obs_to_split_node = nsize, verbose=FALSE),
  rsf =rfsrc(Surv(time,status)~.,data=pbc_cmp,
             ntree=ntree,nodesize=nsize),
  bws = pec::selectCox(Surv(time,status)~., data=pbc_cmp),
  cboost = pec::coxboost(Hist(time,status)~.,data=pbc_cmp)
)

# for reproducibility
set.seed(32989) 

bri_score = pec::pec(mdls, data=pbc_cmp, cens.model = 'cox',
                     formula = Surv(time,status)~1, 
                     splitMethod = 'BootCv', B=nboots)
#> No covariates  specified: Kaplan-Meier for censoring times used for weighting.
#> Split sample loop (B=5)
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> 1
#> 2
#> 3
#> 4
#> Warning in fitter(X, Y, strata = Strata, offset = offset, weights =
#> weights, : Ran out of iterations and did not converge
#> 5

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
#> IPCW: marginal model
#> 
#> Method for estimating the prediction error:
#> 
#> Bootstrap cross-validation
#> 
#> Type: resampling
#> Bootstrap sample size:  276 
#> No. bootstrap samples:  5 
#> Sample size:  276 
#> 
#> Cumulative prediction error, aka Integrated Brier score  (IBS)
#>  aka Cumulative rank probability score
#> 
#> Range of integration: 0 and time=12.4 :
#> 
#> 
#> Integrated Brier score (crps):
#> 
#>           IBS[0;time=12.4)
#> Reference            0.201
#> orsf                 0.138
#> rsf                  0.147
#> bws                  0.159
#> cboost               0.143


cnc_index = pec::cindex(mdls, data=pbc_cmp, cens.model = 'cox',
                        formula = Surv(time,status)~1, 
                        splitMethod = 'BootCv', B=nboots)
#> No covariates  specified: cens.model coerced to "marginal".
#> 1
#> 2
#> 3
#> 4
#> 5

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
#> Censoring model for IPCW: marginal model (Kaplan-Meier for censoring distribution) 
#> 
#> Method for estimating the prediction error:
#> 
#> Bootstrap cross-validation
#> 
#> Type: resampling
#> Bootstrap sample size:  276 
#> No. bootstrap samples:  5 
#> Sample size:  276 
#> 
#> Estimated C-index in % at time=11.5 
#> 
#>        AppCindex BootCvCindex
#> orsf        83.0         80.5
#> rsf         87.4         78.6
#> bws         75.1         72.5
#> cboost      78.8         77.2
#> 
#> AppCindex    : Apparent (training data) performance
#> BootCvCindex : Bootstrap crossvalidated performance
```
