
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
         stage=factor(stage,ordered=TRUE),
         time=time/365.25)

set.seed(62689)

orsf=ORSF(data=pbc, ntree=500, alpha=0.50,verbose=FALSE)
#> 
#> performing imputation with missForest:
#>   missForest iteration 1 in progress...done!
#>   missForest iteration 2 in progress...done!
#>   missForest iteration 3 in progress...done!
#>   missForest iteration 4 in progress...done!
#>   missForest iteration 5 in progress...done!
#>   missForest iteration 6 in progress...done!
#>   missForest iteration 7 in progress...done!
#>   missForest iteration 8 in progress...done!
#>   missForest iteration 9 in progress...done!
#>   missForest iteration 10 in progress...done!

prds=orsf$oob_preds
vdat=orsf$imputed_data
times=orsf$oob_times

# Create a categorical variable indicating stage 4
vdat$stg4=factor(ifelse(vdat$stage==4,"Stage 4", "Stage 3 or less"))

# there are very few values of ast above 300,
# but cutting at 300 makes the figure more interpretable
vdat$ast=pmin(vdat$ast,300)

# Variable dependence plot
vdplot(object=prds[,c(1:5)],
       xvar='ast',
       times=round(times)[1:5],
       data=vdat,
       fvar = 'stg4',
       flab = levels(vdat$stg4),
       time_units='years',
       xlab='Aspartate aminotransferase',
       xvar_units = '(U/ml)',
       nonevent_lab="survival",
       include.hist=TRUE)+
  theme(legend.position = 'bottom',
        legend.direction= 'horizontal')
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](README-usage-1.png)

``` r


nboots=50
```

Performance
===========

`orsf` objects take a long time to grow, but they usually provide excellent predictions. Below we use the `pec` package to compare the integrated Brier scores over 50 replications of bootstrap cross validation using the complete cases in the `pbc` data set.

``` r

# it takes a bit more code to incorporate imputation 
# into bootstrap cross-validation. However, results
# are similar using only the complete observations,
# so that's what we will use here

pbc_cmp=na.omit(pbc)

ntree=300
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
#> Split sample loop (B=50)
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> Warning in fitter(X, Y, strata = Strata, offset = offset, weights =
#> weights, : Ran out of iterations and did not converge
#> 10
#> Warning in fitter(X, Y, strata = Strata, offset = offset, weights =
#> weights, : Ran out of iterations and did not converge
#> 20
#> 30
#> 40
#> 50

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
#> No. bootstrap samples:  50 
#> Sample size:  276 
#> 
#> Cumulative prediction error, aka Integrated Brier score  (IBS)
#>  aka Cumulative rank probability score
#> 
#> Range of integration: 0 and time=11.7 :
#> 
#> 
#> Integrated Brier score (crps):
#> 
#>           IBS[0;time=11.7)
#> Reference            0.187
#> orsf                 0.129
#> rsf                  0.134
#> bws                  0.147
#> cboost               0.135


cnc_index = pec::cindex(mdls, data=pbc_cmp, cens.model = 'cox',
                        formula = Surv(time,status)~1, 
                        splitMethod = 'BootCv', B=nboots)
#> No covariates  specified: cens.model coerced to "marginal".
#> Warning in fitter(X, Y, strata = Strata, offset = offset, weights =
#> weights, : Ran out of iterations and did not converge
#> 10
#> Warning in fitter(X, Y, strata = Strata, offset = offset, weights =
#> weights, : Ran out of iterations and did not converge
#> 20
#> 30
#> 40
#> 50

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
#> No. bootstrap samples:  50 
#> Sample size:  276 
#> 
#> Estimated C-index in % at time=11.5 
#> 
#>        AppCindex BootCvCindex
#> orsf        82.7         82.3
#> rsf         87.5         80.4
#> bws         75.1         78.4
#> cboost      78.6         81.0
#> 
#> AppCindex    : Apparent (training data) performance
#> BootCvCindex : Bootstrap crossvalidated performance
```
