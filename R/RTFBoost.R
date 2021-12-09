#' Tuning and control parameters for (robust) tree-based functional boosting algorithms
#' 
#' Allow users to set the tuning and control parameters for tree-based functional boosting algorithms, including TFBoost(LAD), TFBoost(LAD-M), and TFBoost(RR) for robust scalar-on-function regression, and TFBoost(L2) for L2 scalar-on-function regression. 
#' 
#'
#' Various tuning and control parameters for tree-based functional boosting algorithms implemented in the
#' function \link{RTFBoost}. 
#' 
#' @param make.prediction a logical flag: if \code{TRUE} the algorithm make predictions using \code{x.test} (defaults to \code{TRUE}). 
#' @param type a string that specifies the type of the RTFBoost algorithm. Valid options are
#'  "L2", "LAD", "LAD-M", or "RR", corresponding to the algorithms of TFBoost(L2), TFBoost(LAD), TFBoost(LAD-M), and TFBoost(RR) respectively (defaults to "RR").
#' @param eff a numeric value  between 0 and 1 indicating the efficiency (measured in a linear model with Gaussian errors) of Tukey’s loss function used in the 2nd stage of TFBoost(LAD-M) and TFBoost(RR) (defaults to 0.95).
#' @param bb a numeric value indicating the breakdown point of the M-scale estimator used in the 1st stage of TFBoost(RR) (defaults to 0.5). 
#' @param init.type a string indicating the initial estimator to be used.  Valid options are: "mean", "median" or "LADTree" (defaults to 'median')
#' @param init.tree.params a list of parameters  for the initial LADtree with components \code{tree.type}, \code{tree.nindex}, \code{max.depth}, \code{num.dir}, and \code{min.leafsize}, where 
#' @param loss.s2 a string indicating the loss function used in the M-step of TFBoost(LAD-M). Valid potions are "huber" and "tukey" (defaults to "tukey"). 
#' \code{tree.type} is a character ("A" or "B") denoting Type A or Type B tree (defaults to "A"), 
#' 
#' \code{tree.nindex} is an integer specifying the number of indices required by the initial Type A tree (defaults to 1),  
#' 
#' \code{num.dir} is an integer specifying the number of random directions required by the initial Type B tree (defaults to 20),
#' 
#' \code{max.depth} is an integer specifying the maximum depth for the initial Type A or Type B tree (defaults to 1). 
#' 
#' \code{min.leafsize} is an integer specifying the minimum number of observations per leaf node for the initial Type A or Type B tree (defaults to 10).
#' @param tree.control a list of control parameters for the tree base learner of RTFBoost (defaults to \code{TREE.control()}, see \code{\link{TREE.control}}). 
#' @param shrinkage a numeric shrinkage parameter in boosting (defaults to 0.05). 
#' @param n.init an  integer specifying the number of iterations for the 1st stage of  TFBoost(RR) (defaults to 100). 
#' @param niter  an  integer specifying the number of iterations (for TFBoost(LAD-M) and TFBoost(RR)  \eqn{T_{1,max}} + \eqn{T_{2,max}}) (defaults to 200)
#' @param nknot an integer of the number of interior knots of the cubic B-spline basis (defaults to 3)
#' @param precision an integer denoting the number of significant digits to keep when using validation error to calculate early stopping time (defaults to 4)
#' @param nknot an integer denoting the number of interior knots of the cubic B-spline basis (defaults to 3)
#' @param save.f a logical flag:if \code{TRUE} save the function estimates at all iterations (defaults to \code{FALSE})
#' @param trace a logical flag:  if \code{TRUE} print the number of completed iterations for monitoring progress (defaults to \code{FALSE})
#' @param save.tree a logical flag:  if \code{TRUE} save the tree objects at all iterations, required when the user calls \code{RTFBoost.predict} with the returned object from \code{RTFBoost} (defaults to \code{FALSE})
#' @param error.type a  character  string indicating  the  type  of  error metrics to be evaluated on the test set. Valid options are:  "mse" (mean squared error), "aad" (average absolute deviation), and "tmse" (trimmed mean squared error)
#' @param trim.prop a numeric value used as the trimming proportion  if  "tmse"  is  used  as  the  performance  metric. 
#' "tmse" calculates the mean-square error of residuals (r) of which |r| < quantile(|r|, 1-trim_prop) (e.g. trim_prop = 0.1 ignores 10% of the data and calculates MSE of residuals whose absolute values are below 90% quantile of |r|). If both \code{trim.prop} and \code{trim.c} are specified, \code{trim.c} will be used.
#' @param trim.c a numeric value used for trimming if "trmse" is used as the performance metric  (defaults  to  3). "tmse"  calculates  the  mean squared error  of  the residuals(r) between median(r) + trim_c mad(r) and median(r) - trim_c mad(r).  If both \code{trim.prop} and \code{trim.c} are specified, \code{trim.c} will be used.
#' @return A list of all input parameters
#'
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' 
#' @export
#' 

RTFBoost.control <- function(make.prediction = TRUE, type = "RR", 
                             eff = 0.95, bb = 0.5, 
                             init.type = "median",
                             init.tree.params = list(tree.type = "A",  tree.nindex = 1, max.depth = 1, num.dir = 20, min.leafsize = 10), 
                             loss.s2 = "tukey", 
                             tree.control = TREE.control(), 
                             shrinkage  = 0.05,  n.init = 100, niter = 100, nknot = 3, precision = 4, 
                             save.f = FALSE, trace = FALSE, save.tree = FALSE, 
                             error.type = "mse", trim.prop = NULL, trim.c = 3){

  return(list(make.prediction =  make.prediction, type = type, eff = eff, bb = bb, 
              init.type = init.type, init.tree.params =  init.tree.params, loss.s2  = loss.s2 , tree.control = tree.control, 
              shrinkage  = shrinkage,  n.init = n.init, niter = niter, nknot = nknot, precision = precision,
              save.f = save.f , trace = trace, save.tree = save.tree, 
              error.type = error.type, trim.prop = trim.prop, trim.c = trim.c))
}


#' (Robust) tree-based functional boosting 
#' 
#' This function implements tree-based boosting algorithms for functional regression, including TFBoost(LAD), TFBoost(LAD-M), and TFBoost(RR) for robust scalar-on-function regression, and TFBoost(L2) for L2 scalar-on-function regression. 
#'
#' This function implements tree-based boosting algorithms for functional regression. 
#' The algorithms use functional multi-index regression trees as the base learners, which are constructed using functions available in the \code{rpart} and \code{lme4} packages (see https://arxiv.org/abs/2109.02989). 
#' 
#' @param x.train a matrix or data.frame containing observations of the functional predictor in the training data. 
#' @param z.train a  matrix or data.frame containing values of the scalar predictors in the training data (optional, defaults to \code{NULL}). 
#' @param y.train  a  vector containing values of the scalar response in the training data. 
#' @param x.val a matrix or data.frame containing observations of the functional predictor in the validation data. 
#' @param z.val a  matrix or data.frame containing values of the scalar predictors in the validation data (optional, required when \code{z.train} is provided, defaults to \code{NULL}). 
#' @param y.val  a  vector containing values of the scalar response in the validation data. 
#' @param x.test a matrix or data.frame containing observations of the functional predictor in the test data (required when \code{make.prediction} in control is \code{TRUE}).
#' @param z.test a  matrix or data.frame containing values of the scalar predictors in the test data (required when \code{make.prediction} in \code{control} is \code{TRUE}, and \code{z.train} and  \code{z.val}  are provided).
#' @param y.test  a  vector containing values of the scalar response in the test data.  (optional).
#' @param grid a vector of time points at which  \code{x.train} and \code{x.val} are evaluated. 
#' @param t.range a vector provided in the form of c(a, b), where a and b are numerical values specifying the boundaries of the domain of the functional predictor.  
#' @param control a named list of control parameters as returned by \link{RTFBoost.control}.
#' @param ...  optional additional arguments passed to \code{RTFBoost}, including components \code{user.loss} and \code{tree.init}.  
#' 
#' \code{user.loss} is a list that specifies a user-defined loss  to replace the L2 loss used by TFBoost(L2), including components 
#' \code{func} and \code{func.grad}, denoting the loss function and its derivative respectively; 
#' 
#' \code{tree.init} is a user-provided object that has the same form as the one returned from \link{TREE}. It is used as the tree to initialize the boosting estimator. 
#' 
#' @return A list with the following components:
#'
#' \item{loss.train}{a vector of training errors from all iterations}
#' \item{loss.val}{a vector of validation errors from all iterations}
#' \item{f.train}{a vector of predicted values for the training data at the early stopping iteration}
#' \item{f.val}{a vector of predicted values for the validation data at the early stopping iteration}
#' \item{f.test}{a vector of predicted values for the test data at the early stopping iteration (returned if \code{make.prediction = TRUE} in \code{control})}
#' \item{early.stop.s1}{the early stopping iteration for the 1st stage of TFBoost(LAD-M) and TFBoost(RR) (returned if \code{type = "RR"} or \code{type = "LAD-M"}  in \code{control})}
#' \item{early.stop}{the early stopping iteration (for TFBoost(LAD-M) and TFBoost(RR), equals to \code{n.init} + the early stopping time in the 2nd stage}
#' \item{err.train}{a vector of mean squared prediction errors for training data from all iterations}
#' \item{err.val}{a vector of mean squared prediction errors for validation data from all iterations}
#' \item{err.test}{a vector of mean squared prediction errors for test data  before and at the early stopping iteration (returned if \code{make.prediction = TRUE} in \code{control})}
#' \item{init.vals}{a constant (if \code{init.type = "mean"} or \code{init.type = "median"}  in \code{control}) or a tree object (if \code{init.type = "LADTree"} in \code{control}) to initialize the boosting algorithm}
#' \item{tree.objs}{a list of trees objects (returned if \code{save.tree =  TRUE} in \code{control}, one per iteration and returned from  \link{TREE})}
#' \item{alpha}{a vector of step sizes (one per iteration)}
#' \item{save.f.train}{a matrix of predicted values for the training data from all iterations (returned if \code{save.f = TRUE} in \code{control}, one column per iteration)}
#' \item{save.f.val}{a matrix of predicted values for the validation data from all iterations (returned if \code{save.f = TRUE} in \code{control}, one column per iteration)}
#' \item{save.f.test}{a matrix of predicted values for the test data before and at the early stopping iteration (returned if \code{save.f = TRUE} and \code{make.prediction = TRUE} in \code{control}, one column per iteration)}
#' \item{control}{\code{control} from the input arguments}
#'
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#'
#' @examples
#' 
#' \dontrun{
#' data(fruitfly)
#' n <- nrow(fruitfly$eggs) 
#' n0 <- floor( 0.2 * n) 
#' set.seed(123)
#' idx_test <- sample(n, n0)
#' idx_train <- sample((1:n)[-idx_test], floor( 0.6 * n ) )
#' idx_val <- (1:n)[ -c(idx_test, idx_train) ] 
#' xtrain <- fruitfly$eggs[idx_train, ]
#' ytrain <- fruitfly$lifetime[idx_train ]
#' xval <- fruitfly$eggs[idx_val, ]
#' yval <- fruitfly$lifetime[idx_val ]
#' xtest <- fruitfly$eggs[idx_test, ]
#' ytest <- fruitfly$lifetime[idx_test ]
#' gg <- 1:24  
#' tt <- c(0,24) 
#' niter <- 1000
#' my.control <- TFBoost.control(make_prediction = TRUE, 
#' tree_control = TREE.control(tree_type  = "B", d = 1, num_dir = 200), 
#' shrinkage = 0.05, nknot = 3, loss = "l2")
#'
#' model_TFBoost <- TFBoost(x_train = xtrain, y_train = ytrain,  x_val = xval,  y_val = yval, 
#'        x_test = xtest, y_test = ytest, grid = gg, t_range  = tt, niter = niter, 
#'        control = my.control)
#'}    
#' @export



RTFBoost <- function(x.train, z.train = NULL, y.train, x.val, z.val = NULL, y.val, x.test = NULL, z.test = NULL, y.test = NULL, grid, t.range, control = RTFBoost.control(), ...){
  
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
  }
  
  
  if("user.loss" %in% names(list(...))){
    user.loss <- list(...)$user.loss
  }
  
  if("tree.init" %in% names(list(...))){
    tree.init <- list(...)$tree.init
  }
  
  n.train <- nrow(x.train); n.val <- nrow(x.val)
  niter <- control$niter
  
  if(control$save.f){
    save.f.train <- matrix(NA, n.train, control$niter)
    save.f.val <- matrix(NA, n.val, control$niter)
  }
  
  if(missing(t.range)){
     t.range <- c(min(grid), max(grid))
  }
  
  if(length(z.train)!=0){
    if(!missing(z.train)){
      if(!is.matrix(z.train)){
        z.train = as.matrix(z.train)
      }
      if(!is.matrix(z.val)){
        z.val = as.matrix(z.val)
      }
    }
  }
  
  dd <- 4; p <- ncol(x.train)
  grid0 <- seq(t.range[1],t.range[2], 1/(10*(p-1))) # in case of not evenly spaced
  knot <- quantile(grid0, (1:control$nknot)/(control$nknot+1) )
  delta <- sort(c(rep(range(grid0), dd), knot)) # exterior knots
  B <- spline.des(delta, grid, dd)$design
  B <- compute.orthonormal(B, grid, t.range)
  
  train.predictors <- t(apply(x.train, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t.range)})}))
  val.predictors <- t(apply(x.val, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t.range)})}))
  
  if(missing(x.test)) {
    control$make.prediction <- FALSE 
  } 
  
  alpha <- rep(NA, niter)
  
  loss.train <- loss.val <- rep(NA, niter)
  err.train <- err.val <- rep(NA, niter)
  
  init.tmp <- init.boosting(control$type, control$init.type) 
  init.func <- init.tmp$init.func
  
  if(exists("user.loss")){
    func <- user.loss$func
    func.grad <- user.loss$func.grad
  }else{
    func <- init.tmp$func
    func.grad <- init.tmp$func.grad
  }
  
  
  if(control$type == "RR"|(control$type == c("LAD-M") & control$loss.s2 == "tukey")){
    cc.s <- as.numeric(RobStatTM::lmrobdet.control(bb=control$bb, family='bisquare')$tuning.chi)
    cc.m <-  as.numeric(RobStatTM::lmrobdet.control(efficiency= control$eff, family='bisquare')$tuning.psi)
    func.2 <- func.tukey  # function in the 2nd stage 
    func.grad.2 <- func.tukey.grad # grad in the 2nd stage 
  }
  
  
  if(control$type == c("LAD-M") & control$loss.s2 == "huber"){
    func.2 <- func.huber
    func.grad.2 <- func.huber.grad
    cc.m <- uniroot( function(e) (cal.efficiency(e, func.huber.grad,func.huber.grad.prime)-control$eff), lower=0.1, upper=5)$root
  }
  
  if(control$init.type!="LADTree"){
     init.vals <- function(newx,newz){rep(init.func(y.train), nrow(newx))}
  }else{
    if(is.null(z.train)){
      dat.tmp <- data.frame(x.train, y.train = y.train)
    }else{
      dat.tmp <- data.frame(x.train, y.train = y.train, z.train = z.train)
    }
    
    if(!exists("tree.init")){
      init.tree.params <- control$init.tree.params 
      tree.init <- TREE.init(x = train.predictors, y = y.train, z = z.train, random.seed = 0, 
                             tree.type = init.tree.params$tree.type, 
                             tree.nindex = init.tree.params$tree.nindex, 
                             max.depth = init.tree.params$max.depth, 
                             num.dir = init.tree.params$num.dir,
                             min.leafsize = init.tree.params$min.leafsize,
                             make.prediction = FALSE,  
                              init.nmulti = 3)
    }
    init.vals <- function(newx,newz){TREE.init.predict(model = tree.init, newx = newx, newz = newz)}
  }

  f.train <-  init.vals(train.predictors, z.train)
  f.val <-  init.vals(val.predictors, z.val)
    
  init.status <- 0; ss <- 1;  bb <- control$bb
  
  model <- list(tree.objs = list())
  # initial cc 
  if(control$type == "RR"){
    cc <- cc.s
  }else{  
    cc <- 1
  }

  for(i in 1: (control$niter)){
    
    
    if(control$save.f){
      save.f.train[,i] <- f.train
      save.f.val[,i] <- f.val
    }

    if(control$trace){
      
      if(i%%100 ==0 )
      print(paste(i,"th iteration out of", niter, "iterations"))
    }
    
    if(control$type == "RR" & init.status == 0) {
      ss <- cal.ss.rr(f.train, y.train,  cc = cc.s, bb)
    }

    u <- as.numeric(cal.neggrad(control$type, y.train, f.train, func.grad, init.status, ss, cc))
    model$tree.objs[[i]] <- TREE(x = train.predictors, y = u, z = z.train, random.seed = i, control = control$tree.control)
    model$tree.objs[[i]]$tree.model <- lean_rpart.fn(tree_fitted = model$tree.objs[[i]]$tree.model)

    h.train <- model$tree.objs[[i]]$pred
    h.val <- TREE.predict(model$tree.objs[[i]], newx = val.predictors, newz = z.val)

    alpha[i] <- cal.alpha(f.train, h.train, y.train, func, control$type, init.status, ss, bb, cc)
      
    f.train <-  f.train + control$shrinkage*alpha[i]*h.train
    f.val <-  f.val + control$shrinkage*alpha[i]*h.val
    
    # compute MSE error 
    err.train[i] <- mse(f.train - y.train)
    err.val[i] <- mse(f.val - y.val)
    
    # compute loss
    if(control$type == "RR"){
      if(init.status == 0){
        loss.val[i] <-cal.ss.rr(f.val, y.val,  cc, bb)
        loss.train[i] <-ss   
      }else{
        loss.train[i] <- mean(func((f.train - y.train)/ss, cc = cc))
        loss.val[i] <- mean(func((f.val - y.val)/ss, cc = cc))
      }
    }
    
    if(control$type %in% c("LAD","L2")){
      loss.train[i] <- mean(func(f.train - y.train))
      loss.val[i] <- mean(func((f.val - y.val)))
    }
    
    if(control$type == "LAD-M"){
      if(init.status == 0){
        loss.train[i] <- mean(func(f.train - y.train))
        loss.val[i] <- mean(func((f.val - y.val)))
      }else{
        loss.train[i] <- mean(func((f.train - y.train)/ss, cc = cc))
        loss.val[i] <- mean(func((f.val - y.val)/ss, cc = cc))
      }
    }
    
    if(i == 1){  # initialize
      
      if(control$type %in% c("RR","LAD-M")){
          when.init <- 1
      }
      
      early.stop <- 1
      f.train.early <- f.train
      f.val.early <- f.val
      
    }else{
      if(i <= control$n.init){
        
        if(round(loss.val[i], control$precision) < min(round(loss.val[1:(i-1)], control$precision))){
          
          if(control$type %in% c("RR","LAD-M")){
            when.init <- i
          }
          early.stop <- i
          f.train.early <- f.train
          f.val.early <- f.val
        }
      }else{
        #loss.val[n.init] is changed at the n.init iteration 
        if(control$type %in% c("LAD-M", "RR")){
          if(round(loss.val[i], control$precision) < min(round(loss.val[(control$n.init):(i-1)], control$precision))){
            early.stop <- i
            f.train.early  <- f.train
            f.val.early <- f.val
          }
        }else{
          if(round(loss.val[i], control$precision) < min(round(loss.val[1:(i-1)], control$precision))){
            early.stop <- i
            f.train.early  <- f.train
            f.val.early <- f.val
          }
        }
      }
    }
    
    
    if((control$type == "RR") && (i == control$n.init)){
      init.status <- 1
      f.train <- f.train.early  # reset the current one
      f.val <- f.val.early
      ss <-  RobStatTM::mscale(f.train - y.train, delta = bb)
      cc <- cc.m
      loss.val[i] <- mean(func((f.val - y.val)/ss, cc = cc))
    }
    
    if((control$type == "LAD-M") && (i == control$n.init)){
      init.status <- 1
      f.train <- f.train.early  # reset the current one
      f.val <- f.val.early
      ss <-  RobStatTM::mscale(f.train - y.train, delta = bb)
      cc <- cc.m
      func <- func.2
      func.grad <- func.grad.2
      loss.val[i] <- mean(func((f.val - y.val)/ss, cc = cc))
    }
    
    
  }
  
  
  f.train <- f.train.early
  f.val <- f.val.early
  

  model <- c(model, list(loss.train=loss.train, loss.val = loss.val, f.train =  f.train, f.val = f.val, early.stop = early.stop, err.train = err.train, 
                 err.val = err.val, init.vals = init.vals,  alpha = alpha, control = control))
  
  if(control$type %in% c("LAD-M","RR")){
    model$early.stop.s1 <- when.init
  }
  
  if(control$make.prediction){
    tmp_predict <- RTFBoost.predict(model, newx = x.test, newy = y.test, newz = z.test, grid = grid, t.range = t.range)
    model$f.test <- tmp_predict$f.new
    model$err.test <- tmp_predict$err.new
    if(control$save.f){
      model$save.f.test <- tmp_predict$save.f.new
      
    }
  }
  

  if(!control$save.tree){
    model$tree.objs <- NULL
  }

  if(control$save.f){
    model <- c(model, list(save.f.train = save.f.train, save.f.val = save.f.val))
  }

  return(model)
}
    
#' RTFBoost.predict
#'
#' A function to make predictions and calculate test error given an object returned by TFBoost and test data
#'
#' A function to make predictions and calculate test error given an object returned by TFBoost and test data
#'
#'@param model an object returned by \code{RTFBoost}
#'@param newx a matrix or data.frame containing observations of the functional predictor in the new data 
#'@param newy  a  vector containing values of the scalar response in the new data.  (optional,  defaults to \code{NULL})
#'@param newz a  matrix or data.frame containing values of the scalar predictors in the new data (required if \code{z.train} and  \code{z.val} were used to fit \code{model},  defaults to \code{NULL}).
#'@param grid  a vector of time points at which \code{x.test} are evaluated. 
#'@param t.range a vector provided in the form of c(a, b), where a and b are numerical values specifying the boundaries of the domain of the functional predictor. 
#'@return A list with with the following components:
#'
#' \item{f.new}{predicted values with model at the early stopping iteration using x_test (or x_test and z_test) as the predictors}
#' \item{err.new}{a vector of test errors before and at the early stopping iteration (returned if newy is provided)}
#' \item{save.f.new}{a matrix of test function estimates at all iterations (returned if save_f = TRUE in control)}
#'
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#' 
#' @examples
#' \dontrun{
#' data(GED)
#' n <- nrow(GED$price) 
#' n0 <- floor( 0.2 * n) 
#' set.seed(123)
#' idx_test <- sample(n, n0)
#' idx_train <- sample((1:n)[-idx_test], floor( 0.6 * n ) )
#' idx_val <- (1:n)[ -c(idx_test, idx_train) ] 
#' xtrain <- GED$price[idx_train, ]
#' ytrain <- GED$demand[idx_train ]
#' xval <- GED$price[idx_val, ]
#' yval <- GED$demand[idx_val ]
#' xtest <- GED$price[idx_test, ]
#' ytest <- GED$demand[idx_test ]
#' gg <- 1:24  
#' tt <- c(0,24) 
#' niter <- 1000
#' my.control <- TFBoost.control(make_prediction = TRUE, save_tree = TRUE,
#'                               tree_control = TREE.control(tree_type  = "B", d = 1, num_dir = 200), 
#'                              shrinkage = 0.05, nknot = 3, loss = "l2")
#' model_TFBoost <- TFBoost(x_train = xtrain, y_train = ytrain,  x_val = xval,  y_val = yval, 
#'                          grid = gg, t_range  = tt, niter = niter, 
#'                          control = my.control)

#' predictions<- TFBoost.predict(model_TFBoost, newx = xtest)

#'} 
#' @export
#' 

RTFBoost.predict <- function(model, newx, newz = NULL, newy = NULL, grid, t.range){

  control <- model$control
  
  if(control$save.f){
    save.f.new <- matrix(NA, nrow(newx), model$early.stop)
  }

  dd <- 4; p <- ncol(newx)
  grid0 <- seq(t.range[1],t.range[2], 1/(10*(p-1))) # in case of not evenly spaced
  knot <- quantile(grid0, (1:control$nknot)/(control$nknot+1) )
  delta <- sort(c(rep(range(grid0), dd), knot)) # exterior knots
  B <- spline.des(delta, grid, dd)$design
  B <- compute.orthonormal(B, grid, t.range)
  
  new.predictors <- t(apply(newx, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t.range)})}))
  f.new  <- model$init.vals(new.predictors, newz)
  err.new <- rep(NA, model$early.stop)
  
  for(i in 1: model$early.stop){
    
    if(control$save.f){
      save.f.new[,i] <- f.new
    }
  
    if(control$type %in% c("RR","LAD-M")){
      if( (i <= model$early.stop.s1) | (i >=  (control$n.init+1))){
        f.new <- f.new + control$shrinkage*model$alpha[i]* TREE.predict(model$tree.objs[[i]], newx =  new.predictors, newz = newz)
        if(!missing(newy)){
          err.new[i] <- switch(control$error.type,
                               "mse" = {
                                 mse( f.new- newy)
                               },
                               "aad" = {
                                 aad(f.new - newy)
                               },
                               "tmse" = {
                                 if(control$trim.c!=NULL){
                                   tmse(trim_prop = NULL, trim_c = control$trim.c, f.new - newy)
                                 }else{
                                   tmse(trim_prop = control$trim_prop, trim_c = NULL, f.new - newy)
                                 }
                               }
          )
        }
      }
    }else{
      f.new <- f.new + control$shrinkage*model$alpha[i]* TREE.predict(model$tree.objs[[i]], newx =  new.predictors, newz = newz)
      if(!missing(newy)){
        err.new[i] <- switch(control$error.type,
                             "mse" = {
                               mse( f.new- newy)
                             },
                             "aad" = {
                               aad(f.new - newy)
                             },
                             "tmse" = {
                               if(control$trim.c!=NULL){
                                 tmse(trim_prop = NULL, trim_c = control$trim.c, f.new - newy)
                               }else{
                                 tmse(trim_prop = control$trim_prop, trim_c = NULL, f.new - newy)
                               }
                             }
        )
      }
    }
    
   
    
  }
  
  res = list(f.new = f.new)
  if(control$save.f){
    res <- c(res,  list(save.f.new =  save.f.new))
  }
  
  if(!missing(newy)){
    res <- c(res, list(err.new = err.new))
  }
  return(res)
  
}

# transform a basis matrix to an orthonormal basis matrix
compute.orthonormal <- function(B, grid, t_range){
  
  d <- ncol(B)
  Phi_i <- B
  Psi <- matrix(NA, nrow = nrow(B), ncol(B))
  Psi[,1] <-   B[,1]/ sqrt(riemman(B[,1]*B[,1], grid, t_range))
  
  for(i in 2:d){
    for(j in 1:(i-1)){
      Phi_i[,i] <-   Phi_i[,i]  -  riemman(Phi_i[,i]*Psi[,j], grid, t_range) * Psi[,j]
    }
    Psi[,i] <-   Phi_i[,i]/ sqrt(riemman(Phi_i[,i]*Phi_i[,i], grid, t_range))
  }
  return(Psi)
}



#' Robust tree-based functional boosting with initialization parameters chosen on a validation set
#' 
#' A function to fit TFBoost(RR) (see also \code{\link{RTFBoost}}) where the initialization parameters are chosen
#' based on the performance on the validation set.
#'
#' This function runs the TFBoost(RR) algorithm (see \code{\link{RTFBoost}}) on different combinations of the
#' parameters for the initial fit, and chooses the optimal set based on the performance on the validation set.
#' 
#' @param x.train a matrix or data.frame containing observations of the functional predictor in the training data. 
#' @param z.train a  matrix or data.frame containing values of the scalar predictors in the training data (optional, defaults to \code{NULL}). 
#' @param y.train  a  vector containing values of the scalar response in the training data. 
#' @param x.val a matrix or data.frame containing observations of the functional predictor in the validation data. 
#' @param z.val a  matrix or data.frame containing values of the scalar predictors in the validation data (optional, required when \code{z.train} is provided, defaults to \code{NULL}). 
#' @param y.val  a  vector containing values of the scalar response in the validation data. 
#' @param x.test a matrix or data.frame containing observations of the functional predictor in the test data (required when \code{make.prediction} in control is \code{TRUE}).
#' @param z.test a  matrix or data.frame containing values of the scalar predictors in the test data (required when \code{make.prediction} in \code{control} is \code{TRUE}, and \code{z.train} and  \code{z.val}  are provided).
#' @param y.test  a  vector containing values of the scalar response in the test data.  (optional).
#' @param grid a vector of time points at which  \code{x.train} and \code{x.val} are evaluated. 
#' @param t.range a vector provided in the form of c(a, b), where a and b are numerical values specifying the boundaries of the domain of the functional predictor.  
#' @param control a named list of control parameters as returned by \link{RTFBoost.control}.
#' @param max.depth.init.set a vector of possible values of the maximum depth of the initial LADTree that the algorithm choses from (defaults to c(1,2,3,4)).
#' @param min.leafsize.init.set  a vector of possible values of the minimum observations per node of the initial LADTree that the algorithm choses from (defaults to c(10,20,30)).
#' @param save.all a logical flag: if \code{TRUE} the algorithm saves the predictions for median initialization and LADTree initialization with every combination of \code{max.depth.init.set} and \code{min.leaf.size.init.set} (defaults to \code{TRUE}). 
#' 
#' @return A list with the following components:
#'
#' \item{model.best}{an object returned by \link{RTFBoost} that is trained with selected initialization parameters}
#' \item{param}{a vector of selected initialization parameters (return (0,0) if selected initialization is the median of the training responses)}
#' \item{errs.val}{a vector of robust prediction errors for validation data from all iterations}
#' \item{errs.test}{a vector of robust prediction errors for test data from all iterations  before and at the early stopping iteration (returned if \code{make.prediction = TRUE} in \code{control})}
#' \item{pred.list.train}{a list of predicted values for the training data at the early stopping time  (returned if \code{save.all = TRUE}, one element per combination in the order of median initialization, and then \code{expand.grid(sort(min.leafsize.init.set,TRUE), max_depths= sort(max.depth.init.set))})}
#' \item{pred.list.val}{a list of predicted values for the validation data at the early stopping time  (returned if \code{save.all = TRUE}, one column per combination)}
#' \item{pred.list.test}{a list of predicted values for the test data at the early stopping time  (returned if \code{save.all = TRUE} and \code{make.prediction = TRUE} in \code{control}, one element per combination)}
#'
#'@export


RTFBoost.validation <- function(x.train, z.train = NULL, y.train,  x.val,  z.val = NULL, y.val, x.test, z.test = NULL, y.test, grid, t.range, 
                                max.depth.init.set = c(1,2,3,4), min.leafsize.init.set = c(10,20,30), control = RTFBoost.control(), save.all = TRUE){
  
  
  control.tmp <- control
  control.tmp$init.type <- "median"

  res_list_train <- list()
  res_list_val <- list()
  
  model_best <- RTFBoost(x.train = x.train, z.train = z.train, y.train = y.train,  
                         x.val = x.val, z.val = z.val, y.val = y.val,
                         x.test = x.test, z.test = z.test, y.test = y.test, 
                         grid = grid, t.range  = t.range,  control = control.tmp)
  
  if(save.all){
    res_list_train[[1]] <- model_best$f.train - y.train
    res_list_val[[1]] <- model_best$f.val - y.val
  }
  
  flagger_outlier <- which(abs(model_best$f.val - y.val)>3*mad(model_best$f.val - y.val))
  
  if(length(flagger_outlier) == length(y.val)){
    best_err <- Inf
  }else{
    if(length(flagger_outlier)>=1){
      best_err <- mean(abs(model_best$f.val[-flagger_outlier] - y.val[-flagger_outlier]))  
    }else{
      best_err <- mean(abs(model_best$f.val - y.val))
    }
  }
  
  params = c(0,0)
  
  # loss and error for all combinations 
  errs_val <-  rep(NA, 1+ length(min.leafsize.init.set)*length(max.depth.init.set))
  
  # test error for all combinations 
  if(control$make.prediction){
    errs_test <- rep(NA, 1+ length(min.leafsize.init.set)*length(max.depth.init.set))
  } 
  
  errs_val[1] <- best_err

  
  if(control$make.prediction){
    res_list_test <- list()
    res_list_test[[1]] <- model_best$f.test - y.test
    errs_test[1] <- as.numeric(model_best$err.test[model_best$early.stop])
  }
  
  deg <- 4; p <- ncol(x.train)
  grid0 <- seq(t.range[1],t.range[2], 1/(10*(p-1))) # in case of not evenly spaced
  knot <- quantile(grid0, (1:control$nknot)/(control$nknot+1) )
  delta <- sort(c(rep(range(grid0), deg), knot)) #exterior knots
  B<- spline.des(delta, grid, deg)$design
  B <- compute.orthonormal(B,grid, t.range)
  
  
  train.predictors <- t(apply(x.train, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t.range)})}))
  val.predictors <- t(apply(x.val, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t.range)})}))
  
  
  if(control$init.type == "LADTree") {
    
    model_pre_tree <- NA
    combs <- expand.grid(min_leafs= sort(min.leafsize.init.set,TRUE), max_depths= max.depth.init.set)
    j_tmp <- rep(1, nrow(combs)) #need to consider that tree or not
    
    tree_init_list <- list()
    for(j in 1:nrow(combs)) {
      min_leaf_size <- combs[j, 1]
      max_depths <- combs[j, 2]
      if(is.null(z.train)){
        dat_tmp <- data.frame(x.train, y.train = y.train)
      }else{
        dat_tmp <- data.frame(x.train, y.train = y.train, z.train = z.train)
      }
      
      tree_init_list[[j]] <- TREE.init(x = train.predictors, y = y.train, z = z.train, newx = val.predictors, newy = y.val, 
                                       newz = z.val, random.seed = 0, tree.type = control$init.tree.params$tree.type, 
                                       tree.nindex  = control$init.tree.params$tree.nindex, 
                                       max.depth = max_depths,  
                                       min.leafsize = min_leaf_size,
                                       num.dir = control$init.tree.params$num.dir, make.prediction = TRUE,
                                       init.nmulti = 3)
      
    }
    
    if(length(min.leafsize.init.set) >1){
      for(j in 1:length(max.depth.init.set)){
        for(k in 1:(length(min.leafsize.init.set)-1)){
          idx_jk <- which(combs[,1] == sort(min.leafsize.init.set, TRUE)[k] & combs[,2] == max.depth.init.set[j])
          idx_jk_plus <- which(combs[,1] == sort(min.leafsize.init.set, TRUE)[k+1] & combs[,2] == max.depth.init.set[j])
          equal_tmp<- all.equal(tree_init_list[[idx_jk]]$tree.model,tree_init_list[[idx_jk_plus]]$tree.model) == TRUE
          if(length(equal_tmp)==2 & sum(equal_tmp) == 0){
            j_tmp[ idx_jk_plus] <- 0
          }
        }
      }
    }
    
    if(length(max.depth.init.set) > 1){
      for(k in 1:length(min.leafsize.init.set)){
        for(j in 1:(length(max.depth.init.set)-1)){
          idx_kj <- which(combs[,1] == sort(min.leafsize.init.set, TRUE)[k] & combs[,2] == max.depth.init.set[j])
          idx_kj_plus<- which(combs[,1] == sort(min.leafsize.init.set, TRUE)[k] & combs[,2] == max.depth.init.set[j+1])
          equal_tmp<- all.equal(tree_init_list[[idx_kj]]$tree.model,tree_init_list[[idx_kj_plus]]$tree.model) == TRUE
          if(length(equal_tmp)==2 & sum(equal_tmp) == 0){
            j_tmp[idx_kj_plus] <- 0
          }
        }
      }
    }
    print(j_tmp)  # all the selected initial trees 
    
  
    for(j in 1:nrow(combs)) {
      
      print(j)
      if(j_tmp[j] == 1){
        
        min_leaf_size <- combs[j, 1]
        max_depths <- combs[j, 2]
        
        control.tmp$init.type <- "LADTree"
        control.tmp$init.tree.params$max.depth <- max_depths
        control.tmp$init.tree.params$min.leafsize  <- min_leaf_size
        
        
        model_tmp <- RTFBoost(x.train = x.train, y.train = y.train,  x.val = x.val,  y.val = y.val,
                              x.test = x.test, y.test = y.test, grid = grid, t.range  = t.range,
                              control = control.tmp, tree.init = tree_init_list[[j]])
        
        if(length(flagger_outlier) == length(y.val)){
          err_tmp <- Inf
        }else{
          if(length(flagger_outlier) >=1){
            err_tmp <- mean(abs(model_tmp$f.val[-flagger_outlier] - y.val[-flagger_outlier]))
          }else{
            err_tmp <- mean(abs(model_tmp$f.val - y.val))
          }
        }
        
        errs_val[j+1] <- err_tmp
        if(save.all){
          res_list_train[[j+1]] <- model_tmp$f.train - y.train
          res_list_val[[j+1]] <- model_tmp$f.val - y.val
          if(control$make.prediction){
            res_list_test[[j+1]] <- model_tmp$f.test - y.test
          }
        }
        
        if(control$make.prediction){
          errs_test[j+1] <- as.numeric(model_tmp$err.test[model_tmp$early.stop])
        }
        
        if(control$trace){
          print(paste("leaf size:", min_leaf_size, " depths:", max_depths, " err(val):", round(err_tmp,4), " best err(val) :", round(best_err,4) ,sep = ""))
        }
        
        if(err_tmp < best_err) {
          model_best <- model_tmp
          params <- combs[j, ]
          best_err <- err_tmp
          rm(model_tmp)
        }else{
          rm(model_tmp)
        }
      }
    }
  }
  
  # validation error 
  res <- list(model.best = model_best, params = params, errs.test = errs_test) 
           
  if(save.all){
    res$pred.list.train <- res_list_train
    res$pred.list.val <- res_list_val
    
    if(control$make.prediction){
      res$pred.list.test <- res_list_test
    }
  }  
  
  return(res)
}

