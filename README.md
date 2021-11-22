RTFBoost: a package for a tree-based functional boosting algorithm
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

## An example: Fruit fly data

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
`control$type = LAD`, `control$type = LAD-M`, and `control$type = RR`
respectively.

We now explain how to fit a `RTFBoost` estimator with different types
and compare it with the `fgam` estimator proposed in [McLean el
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924/) and
implemented in the `refund` package.

To specify the tree type, the user needs to set `tree_type = A` or
`tree_type = B` in `tree_control`. Below, we will fit `RTFBoost` with
the type B tree, which trains much faster compared to the type A tree.

The following are parameters required for our estimator

``` r
tree_type  <- "B" # type of the base learner
num_dir <- 20  # number of random directions for type B tree
gg <- 1:30  # specify the grid the functional predictor was evaluated on
tt <- c(0,24) # domain of the functional predictor
niter <- 1000 # number of boosting iterations 
make_prediction <- TRUE # make predictions based on test data
loss <-  "l2" # loss for the boosting algorithm ("l2", "lad", or one specified by user_func)
shrinkage <- 0.05 # shrinkage parameter for boosting
nknot <- 3 # the number of interior knots for cubic B-spline basis
```
