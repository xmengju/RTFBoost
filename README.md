RTFBoost: a package for (robust) tree-based functional boosting
algorithms
================
Xiaomeng Ju and Matias Salibian Barrera
2021-11-23

This repository contains `R` code implementing a robust tree-based
boosting algorithm for scalar-on-function regression.

## Install and load package

You can install the development version of the package in R using:

``` r
devtools::install_github("xmengju/RTFBoost", auth_token ='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')
```

Once installed you can load the package with:

``` r
library(RTFBoost)
```

## An example: fruit fly data

Below we illustrate the use of the package with the fruit fly dataset.
The original data is provided
[here](https://anson.ucdavis.edu/~mueller/data/data.html). The data set
consists of number of eggs laid daily for each of 1000 medflies until
time of death. We used a subset of the data including flies that lived
at least 30 days, to predict the number of living days. The number of
eggs laid in the first 30 days is treated as the functional predictor
with values evaluated at each day.

In the data set, `fruitfly$eggs` denotes the predictor, and
`fruitfly$lifetime` denotes the response. We plot the data, highlighting
10 randomly chosen predictor curves in red. We also plot the
distribution of the lifetime variable, measured in days.

``` r
data(fruitfly)
matplot(t(fruitfly$eggs), lty = 1, type = "l", ylab = "price", xlab = "days", 
        main = "number of eggs laid each day (10 days)", col='gray70', ylim=c(0, 150))
set.seed(123)
n <- nrow(fruitfly$eggs)
matplot(t(fruitfly$eggs[sample(n, 10), ]), type='l', col='red', add=TRUE, lty=1)
hist(fruitfly$lifetime, xlab = "days", 
     main = "Histogram of the lifetime variable")
```

<img src="README_files/figure-gfm/plot-1.png" width="45%" /><img src="README_files/figure-gfm/plot-2.png" width="45%" />

In order to train our predictor, we split the data set into a `training`
set (with 60% of the available data), a `validation` set and a `test`
set (both with 20% of the data). We first randomly select the
observations for each of these three sets:

``` r
n0 <- floor( 0.2 * n) 
set.seed(1)
idx_test <- sample(n, n0)
idx_train <- sample((1:n)[-idx_test], floor( 0.6 * n ) )
idx_val <- (1:n)[ -c(idx_test, idx_train) ] 
```

We now create the matrices of explanatory variables (`x`) and vectors of
responses (`y`) corresponding to this partition.
<!-- Note that `ytrain` and `yval` may contain outliers. -->

``` r
xtrain <- fruitfly$eggs[idx_train, ]
ytrain <- fruitfly$lifetime[idx_train ]
xval <- fruitfly$eggs[idx_val, ]
yval <- fruitfly$lifetime[idx_val ]
xtest <- fruitfly$eggs[idx_test, ]
ytest <- fruitfly$lifetime[idx_test ]
```

The `RTFBoost` function implements tree-based functional boosting with
four options:

-   `TFBoost(LS)`: TFBoost with the squared loss;
-   `TFBoost(LAD)`: TFBoost with the L1 loss;
-   `RTFBoost(LAD-M)`: two-stage robust TFBoost that involves a
    LAD-stage and an M-stage; and
-   `RTFBoost(RR)`: two-stage robust TFBoost that involves an S-stage
    and an M-stage.

These options correspond to setting `control$type = 'L2'`,
`control$type = 'LAD'`, `control$type = 'LAD-M'`, and
`control$type = ''RR'` respectively.

We now explain how to fit a `RTFBoost` estimator with different types
and compare it with the `fgam` estimator proposed in [McLean el
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924/) and
implemented in the `refund` package.

The `RTFBoost` functional allows using Type A or Type B trees as the
base learners. To specify the tree type, the user needs to set
`tree.type = A` or `tree.type = B` in `control$tree.control`. Below, we
will fit `RTFBoost` with the type B trees, which train much faster
compared to the type A trees.

The following are parameters required for our estimator

``` r
tree.type  <- "B" # type of the base learner
num.dir <- 200  # number of random directions for type B tree
gg <- 1:30  # specify the grid the functional predictor was evaluated on
tt <- c(0,30) # domain of the functional predictorn
niter <- 1000 # number of boosting iterations 
make.prediction <- TRUE # make predictions based on test data
shrinkage <- 0.05 # shrinkage parameter for boosting
nknot <- 3 # the number of interior knots for cubic B-spline basis
init.type <- "median" # intialize with the median of the training responses
```

The depth of the base learners in `RTFBoost` is set with the argument
`max.depth` in `control$tree.control`. We considered `max.depth` from 1
to 4, and chose the depth the minimizes the robust MSPE on the
validation set (ℐ<sub>val</sub>) at early stopping time. Denote
validation residuals as
*r*<sub>*i*</sub> = *F̂*(*x*<sub>*i*</sub>) − *y*<sub>*i*</sub>, *i* ∈ ℐ<sub>val</sub>
, the robust MSPE is defined as
*μ̂*<sub>*M*</sub>({*r*<sub>*i*</sub>, *i* ∈ ℐ<sub>val</sub>}) + *σ̂*<sub>*M*</sub><sup>2</sup>({*r*<sub>*i*</sub>, *i* ∈ ℐ<sub>val</sub>}),
where *μ*<sub>*M*</sub> is the M-location estimator and
*σ*<sub>*M*</sub> is the M-scale estimator. We use Tukey’s score
function and let asymptotic efficiency of *μ*<sub>*M*</sub> and
*σ*<sub>*M*</sub> be 95%. Below is the function to compute robust MSPE:

``` r
cal.rmspe <- function(r){
  tmp <-  RobStatTM::locScaleM(x=r, psi='bisquare', eff= 0.95)
  return( val.errors[dd] <- tmp$mu^2 + tmp$disper^2)
}
```

Below we train `RTFBoost` with `control$type = 'L2'` and select the
depth `max.depth`. The estimator is initialized at the median of all
training responses.

``` r
type <- "L2"
tree.depths <- 1:3
model.list.l2 <- vector('list', length(tree.depths))
val.errors <- rep(NA, length(tree.depths))  # vector of validation robust MSPE

for(dd in 1:length(tree.depths)) {
   model.list.l2[[dd]] <- 
   RTFBoost(x.train = xtrain, y.train = ytrain,  x.val = xval,  y.val = yval,
            x.test = xtest, y.test = ytest, grid = gg, t.range  = tt, 
            control = RTFBoost.control(make.prediction = make.prediction, niter = niter,                  
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = dd),
            shrinkage = shrinkage, nknot = nknot, type = type, init.type = init.type))
   val.errors[dd] <- cal.rmspe(model.list.l2[[dd]]$f.val - yval)
}

model.l2 <-  model.list.l2[[which.min(val.errors)]]
rm(model.list.l2)
```

A sanity check

``` r
all.equal(mean((model.l2$f.test - ytest)^2), tail(model.l2$err.test,1))
```

    ## [1] TRUE

Then we repeat the same procedure fitting `TFBoost(LAD)` and
`TFBoost(LAD-M)`.

``` r
type <- "LAD"
tree.depths <- 1:3
model.list.lad <- vector('list', length(tree.depths))
val.errors <- rep(NA, length(tree.depths))  # vector of validation robust MSPE

for(dd in 1:length(tree.depths)) {
   model.list.lad[[dd]] <- 
   RTFBoost(x.train = xtrain, y.train = ytrain,  x.val = xval,  y.val = yval,
            x.test = xtest, y.test = ytest, grid = gg, t.range  = tt, 
            control = RTFBoost.control(make.prediction = make.prediction, niter = niter,                  
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = dd),
            shrinkage = shrinkage, nknot = nknot, type = type, init.type = init.type))
   val.errors[dd] <- cal.rmspe(model.list.lad[[dd]]$f.val - yval)
}

model.lad <-  model.list.lad[[which.min(val.errors)]]
rm(model.list.lad)
```

``` r
type <- "LAD-M"
n.init <- 500  # maximum number of iterations in the LAD-stage
niter <- 1000  # maximum number of iterations in the M-stage
tree.depths <- 1:3
model.list.ladm <- vector('list', length(tree.depths))
val.errors <- rep(NA, length(tree.depths))  # vector of validation robust MSPE

for(dd in 1:length(tree.depths)) {
   model.list.ladm[[dd]] <- 
   RTFBoost(x.train = xtrain, y.train = ytrain,  x.val = xval,  y.val = yval,
            x.test = xtest, y.test = ytest, grid = gg, t.range  = tt, 
            control = RTFBoost.control(make.prediction = make.prediction, niter = niter, n.init = n.init,                
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = dd),
            shrinkage = shrinkage, nknot = nknot, type = type, init.type = init.type))
   val.errors[dd] <- cal.rmspe(model.list.ladm[[dd]]$f.val - yval)

}

model.ladm <-  model.list.ladm[[which.min(val.errors)]]
rm(model.list.ladm)
```

For `RTFBoost(RR)` we controlled the level of robustness of the S-stage
by setting `bb  = 0.5` in the S-stage, corresponding to a M-local
estimator achieving 50% breakdown point. To initialize `RTFBoost(RR)`,
we used a single-index Type A LAD tree. We set the minimum number of
observations per node of the initial tree is to 10 and selected the
maximum depth of the tree from 0, 1, 2, and 3 to minimize the robust
MSPE on the validation set. An initial tree with depth 0 corresponds to
initializing `RTFBoost(RR)` with the median of all training responses.
This step may take several minutes to run.

``` r
n.init <- 200  # maximum number of iterations in the S-stage
niter <- 1000  # maximum number of iterations in the M-stage
type <- "RR"
tree.depths <- 1:3
model.list.rr <- vector('list', length(tree.depths))
val.errors <- rep(NA, length(tree.depths))  # vector of validation robust MSPE

for(dd in 1:length(tree.depths)) {
   model.list.rr[[dd]] <- 
    RTFBoost.validation(x.train = xtrain,  y.train = ytrain,  x.val = xval,  y.val =  yval, x.test = xtest,  y.test = ytest, grid = gg, t.range = tt, 
                                max.depth.init.set = c(1,2,3), min.leafsize.init.set = c(10), control = RTFBoost.control(make.prediction = make.prediction, niter = niter, n.init = n.init,                
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = dd),
            shrinkage = shrinkage, nknot = nknot, type = type, init.type = "LADTree"), save.all = FALSE)
   val.errors[dd]  <- cal.rmspe(model.list.rr[[dd]]$model.best$f.val - yval)

}
```

    ## [1] 1 1 1
    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 1 1 1
    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 1 1 1
    ## [1] 1
    ## [1] 2
    ## [1] 3

``` r
model.rr <-  model.list.rr[[which.min(val.errors)]]$model.best
params <- model.list.rr[[which.min(val.errors)]]$params
print(params)
```

    ## [1] 0 0

`RTFBoost(RR)` initialized with the median of training responses.

``` r
n.init <- 200  # maximum number of iterations in the S-stage
niter <- 1000  # maximum number of iterations in the M-stage
type <- "RR"
tree.depths <- 1:3
model.list.rrm <- vector('list', length(tree.depths))
val.errors <- rep(NA, length(tree.depths))  # vector of validation robust MSPE

for(dd in 1:length(tree.depths)) {
   model.list.rrm[[dd]] <- 
    RTFBoost(x.train = xtrain,  y.train = ytrain,  x.val = xval,  y.val =  yval, x.test = xtest,  y.test = ytest, grid = gg, t.range = tt, 
                                control = RTFBoost.control(type = type, make.prediction = make.prediction, niter = niter, n.init = n.init,                
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = dd),
            shrinkage = shrinkage, nknot = nknot, init.type = "median"))
   val.errors[dd] <- cal.rmspe(model.list.rrm[[dd]]$f.val - yval)
}

model.rr.m <-  model.list.rrm[[which.min(val.errors)]]
```

We fit `fgam` which estimates the regression function using
tensor-product B-splines with roughness penalties. We set the marginal
basis to be cubic B-splines of dimension 15 and used REML to select the
parameter that penalizes the second order marginal differences.

``` r
library(refund)
nbasis_FGAM <- 15
xx <- xtrain
yy <- ytrain
model.fgam <- fgam(yy ~ af(xx, splinepars = list(k=c(nbasis_FGAM,nbasis_FGAM), 
                                                 m = list(c(2,2),c(2,2)))), 
                   method ="REML")
```

We compare test errors of different estimators:

``` r
err.l2 <-mean((model.l2$f.test- ytest)^2 )
err.lad <-mean( (model.lad$f.test- ytest)^2 )
err.ladm <-mean( (model.ladm$f.test- ytest)^2 )
err.rr <-mean( (model.rr$f.test- ytest)^2 )
err.rr.m <-mean( (model.rr.m$f.test- ytest)^2 )
err.fgam <-  mean( (predict(model.fgam,
                               newdata=list(xx =xtest),
                               type='response') - ytest)^2 )
print(c(err.l2, err.lad, err.ladm, err.rr, err.rr.m, err.fgam))
```

    ## [1] 150.3275 146.3308 147.3556 146.5420 146.5420 152.9956

``` r
err.l2 <-cal.rmspe((model.l2$f.test- ytest))
err.lad <-cal.rmspe( (model.lad$f.test- ytest))
err.ladm <-cal.rmspe( (model.ladm$f.test- ytest))
err.rr <-cal.rmspe( (model.rr$f.test- ytest))
err.rr.m <-cal.rmspe( (model.rr.m$f.test- ytest))

err.fgam <-  cal.rmspe( (predict(model.fgam,
                               newdata=list(xx =xtest),
                               type='response') - ytest))

print(c(err.l2, err.lad, err.ladm, err.rr, err.rr.m, err.fgam))
```

    ## [1] 114.59009  96.90441  96.75248  88.19130  88.19130 125.81592

For `RTFBoost`, we can also separate the process of training the
predictor and evaluating it on a test set. In this way we have the
flexibility to compute predictions on multiple test sets when needed. In
this case the user needs to use `save.tree = TRUE` in the argument
`control`. This option includes the tree objects at all iterations in
the returned object, which are needed for calculating predictions. The
following code illustrates this procedure running `TFBoost(L2)` with
`max.depth = 1`:

``` r
max.depth <- 1
type <- "L2"
model.l2.separate <-   RTFBoost(x.train = xtrain, y.train = ytrain,  x.val = xval,  y.val = yval,
            x.test = xtest, y.test = ytest, grid = gg, t.range  = tt, 
            control = RTFBoost.control(make.prediction = TRUE, niter = niter,                  
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = max.depth),
            shrinkage = shrinkage, nknot = nknot, type = type, init.type = "median", save.tree = TRUE))
predictions <- RTFBoost.predict(model.l2.separate, newx = xtest, grid = gg, t.range = tt)
```

A sanity check

``` r
all.equal(mean((predictions$f.new - ytest)^2), tail(model.l2.separate$err.test, 1))
```

    ## [1] TRUE

``` r
all.equal(mean((model.l2.separate$f.test - ytest)^2), tail(model.l2.separate$err.test, 1))
```

    ## [1] TRUE
