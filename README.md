RTFBoost: a package for (robust) tree-based functional boosting
algorithms
================
Xiaomeng Ju and Matias Salibian Barrera
2021-11-22

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
set.seed(123)
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
num.dir <- 20  # number of random directions for type B tree
gg <- 1:30  # specify the grid the functional predictor was evaluated on
tt <- c(0,30) # domain of the functional predictor
niter <- 1000 # number of boosting iterations 
make.prediction <- TRUE # make predictions based on test data
shrinkage <- 0.05 # shrinkage parameter for boosting
nknot <- 3 # the number of interior knots for cubic B-spline basis
```

The depth of the base learners in `RTFBoost` is set with the argument
`max.depth` in `control$tree.control`. We considered `max.depth` from 1
to 4, and chose the depth the minimizes the robust MSPE on the
validation set (ℐ<sub>val</sub>) at early stopping time. Denote
validation residuals as
*r*<sub>*i*</sub> = *F̂*(*x*<sub>*i*</sub>) − *y*<sub>*i*</sub>, *i* ∈ ℐ<sub>val</sub>
, the robust MSPE is defined as

$$
\\begin{aligned}
 AR(p): Y\_i &= c + \\epsilon\_i + \\phi\_i Y\_{i-1} \\dots \\\\
 Y\_{i} &= c + \\phi\_i Y\_{i-1} \\dots
\\end{aligned}
$$

$$ \\begin{align}
    \\mu\_M(r) + \\sigma\_M^2(r),
 \\end{align}
$$

where *μ*<sub>*M*</sub> is the M-location estimator and
*σ*<sub>*M*</sub> is the M-scale estimator. We use Tukey’s score
function and let asymptotic efficiency of *μ*<sub>*M*</sub> and
*σ*<sub>*M*</sub> be 95%.

Below we train `RTFBoost` with `control$type = 'L2'` and select the
depth `max.depth`. The estimator is initialized at the median of all
training responses.

``` r
type <- "L2"
tree.depths <- 1:4
model.list <- vector('list', length(tree.depths))
val.errors <- rep(NA, length(tree.depths))  # vector of validation robust MSPE

for(dd in 1:length(tree.depths)) {
   model.list[[dd]] <- 
   RTFBoost(x.train = xtrain, y.train = ytrain,  x.val = xval,  y.val = yval,
            x.test = xtest, y.test = ytest, grid = gg, t.range  = tt, 
            control = RTFBoost.control(make.prediction = make.prediction, niter = niter,                  
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = dd),
            shrinkage = shrinkage, nknot = nknot, type = type))
   tmp <-  RobStatTM::locScaleM(x=model.list[[dd]]$f.val - yval, psi='bisquare', eff=0.95)
   val.errors[dd] <- tmp$mu^2 + tmp$disper^2
}

model.l2 <-  model.list[[which.min(val.errors)]]
```

A sanity check

``` r
all.equal(mean((model.l2$f.test - ytest)^2), tail(model.l2$err.test,1))
```

Then we repeat the same procedure fitting `TFBoost(LAD)` and
`TFBoost(LAD-M)`.

``` r
type <- "LAD"
tree.depths <- 1:4
model.list <- vector('list', length(tree.depths))
val.errors <- rep(NA, length(tree.depths))  # vector of validation robust MSPE

for(dd in 1:length(tree.depths)) {
   model.list[[dd]] <- 
   RTFBoost(x.train = xtrain, y.train = ytrain,  x.val = xval,  y.val = yval,
            x.test = xtest, y.test = ytest, grid = gg, t.range  = tt, 
            control = RTFBoost.control(make.prediction = make.prediction, niter = niter,                  
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = dd),
            shrinkage = shrinkage, nknot = nknot, type = type))
   tmp <-  RobStatTM::locScaleM(x=model.list[[dd]]$f.val - yval, psi='bisquare', eff=0.95)
   val.errors[dd] <- tmp$mu^2 + tmp$disper^2
}

model.lad <-  model.list[[which.min(val.errors)]]
```

``` r
type <- "LAD-M"
tree.depths <- 1:4
model.list <- vector('list', length(tree.depths))
val.errors <- rep(NA, length(tree.depths))  # vector of validation robust MSPE

for(dd in 1:length(tree.depths)) {
   model.list[[dd]] <- 
   RTFBoost(x.train = xtrain, y.train = ytrain,  x.val = xval,  y.val = yval,
            x.test = xtest, y.test = ytest, grid = gg, t.range  = tt, 
            control = RTFBoost.control(make.prediction = make.prediction, niter = niter,                  
            tree.control =  TREE.control(tree.type = tree.type, num.dir = num.dir, max.depth = dd),
            shrinkage = shrinkage, nknot = nknot, type = type))
   tmp <-  RobStatTM::locScaleM(x=model.list[[dd]]$f.val - yval, psi='bisquare', eff=0.95)
   val.errors[dd] <- tmp$mu^2 + tmp$disper^2
}

model.lad <-  model.list[[which.min(val.errors)]]
```

For
