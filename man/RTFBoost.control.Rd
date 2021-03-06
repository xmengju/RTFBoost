% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RTFBoost.R
\name{RTFBoost.control}
\alias{RTFBoost.control}
\title{Tuning and control parameters for (robust) tree-based functional boosting algorithms}
\usage{
RTFBoost.control(
  make.prediction = TRUE,
  type = "RR",
  eff = 0.95,
  bb = 0.5,
  init.type = "median",
  init.tree.params = list(tree.type = "A", tree.nindex = 1, max.depth = 1, num.dir =
    20, min.leafsize = 10),
  loss.s2 = "tukey",
  tree.control = TREE.control(),
  shrinkage = 0.05,
  n.init = 100,
  niter = 100,
  nknot = 3,
  precision = 4,
  save.f = FALSE,
  trace = FALSE,
  save.tree = FALSE,
  error.type = "mse",
  trim.prop = NULL,
  trim.c = 3,
  stage.no = 0
)
}
\arguments{
\item{make.prediction}{a logical flag: if \code{TRUE} the algorithm make predictions using \code{x.test} (defaults to \code{TRUE}).}

\item{type}{a string that specifies the type of the RTFBoost algorithm. Valid options are
"L2", "LAD", "LAD-M", or "RR", corresponding to the algorithms of TFBoost(L2), TFBoost(LAD), TFBoost(LAD-M), and TFBoost(RR) respectively (defaults to "RR").}

\item{eff}{a numeric value  between 0 and 1 indicating the efficiency (measured in a linear model with Gaussian errors) of Tukey’s loss function used in the 2nd stage of TFBoost(LAD-M) and TFBoost(RR) (defaults to 0.95).}

\item{bb}{a numeric value indicating the breakdown point of the M-scale estimator used in the 1st stage of TFBoost(RR) (defaults to 0.5).}

\item{init.type}{a string indicating the initial estimator to be used.  Valid options are: "mean", "median" or "LADTree" (defaults to 'median')}

\item{init.tree.params}{a list of parameters  for the initial LADtree with components \code{tree.type}, \code{tree.nindex}, \code{max.depth}, \code{num.dir}, and \code{min.leafsize}, where}

\item{loss.s2}{a string indicating the loss function used in the M-step of TFBoost(LAD-M). Valid potions are "huber" and "tukey" (defaults to "tukey").
\code{tree.type} is a character ("A" or "B") denoting Type A or Type B tree (defaults to "A"),

\code{tree.nindex} is an integer specifying the number of indices required by the initial Type A tree (defaults to 1),

\code{num.dir} is an integer specifying the number of random directions required by the initial Type B tree (defaults to 20),

\code{max.depth} is an integer specifying the maximum depth for the initial Type A or Type B tree (defaults to 1).

\code{min.leafsize} is an integer specifying the minimum number of observations per leaf node for the initial Type A or Type B tree (defaults to 10).}

\item{tree.control}{a list of control parameters for the tree base learner of RTFBoost (defaults to \code{TREE.control()}, see \code{\link{TREE.control}}).}

\item{shrinkage}{a numeric shrinkage parameter in boosting (defaults to 0.05).}

\item{n.init}{an  integer specifying the number of iterations for the 1st stage of  TFBoost(RR) (defaults to 100).}

\item{niter}{an  integer specifying the number of iterations (for TFBoost(LAD-M) and TFBoost(RR)  \eqn{T_{1,max}} + \eqn{T_{2,max}}) (defaults to 200)}

\item{nknot}{an integer denoting the number of interior knots of the cubic B-spline basis (defaults to 3)}

\item{precision}{an integer denoting the number of significant digits to keep when using validation error to calculate early stopping time (defaults to 4)}

\item{save.f}{a logical flag:if \code{TRUE} save the function estimates at all iterations (defaults to \code{FALSE})}

\item{trace}{a logical flag:  if \code{TRUE} print the number of completed iterations for monitoring progress (defaults to \code{FALSE})}

\item{save.tree}{a logical flag:  if \code{TRUE} save the tree objects at all iterations, required when the user calls \code{RTFBoost.predict} with the returned object from \code{RTFBoost} (defaults to \code{FALSE})}

\item{error.type}{a  character  string indicating  the  type  of  error metrics to be evaluated on the test set. Valid options are:  "mse" (mean squared error), "aad" (average absolute deviation), and "tmse" (trimmed mean squared error)}

\item{trim.prop}{a numeric value used as the trimming proportion  if  "tmse"  is  used  as  the  performance  metric.
"tmse" calculates the mean-square error of residuals (r) of which |r| < quantile(|r|, 1-trim_prop) (e.g. trim_prop = 0.1 ignores 10\% of the data and calculates MSE of residuals whose absolute values are below 90\% quantile of |r|). If both \code{trim.prop} and \code{trim.c} are specified, \code{trim.c} will be used.}

\item{trim.c}{a numeric value used for trimming if "trmse" is used as the performance metric  (defaults  to  3). "tmse"  calculates  the  mean squared error  of  the residuals(r) between median(r) + trim_c mad(r) and median(r) - trim_c mad(r).  If both \code{trim.prop} and \code{trim.c} are specified, \code{trim.c} will be used.}

\item{stage.no}{a numeric value indicating which stages are executed, required by RTFBoost(LAD-M) and RTFBoost(RR): 0 runs both stages, 1 runs the first stage only, and 2 runs the second stage only   (defaults  to 0).}
}
\value{
A list of all input parameters
}
\description{
Allow users to set the tuning and control parameters for tree-based functional boosting algorithms, including TFBoost(LAD), TFBoost(LAD-M), and TFBoost(RR) for robust scalar-on-function regression, and TFBoost(L2) for L2 scalar-on-function regression.
}
\details{
Various tuning and control parameters for tree-based functional boosting algorithms implemented in the
function \link{RTFBoost}.
}
\author{
Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
}
