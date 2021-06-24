#' Tuning and control parameters for the tree-based functional boosting algorithm
#' 
#' Tuning and control parameters for the TFBoost algorithm.
#'
#' Various tuning and control parameters for the \code{\link{TFBoost}} algorithm implemented in the
#' function 
#' 
#' @param make_prediction logical indicating whether to make predictions using \code{x_test} (defaults to \code{TRUE})
#' @param tree_control control parameters for the tree learner of TFBoost (defaults to \code{TREE.control()}, see \code{\link{TREE.control}})
#' @param trim_prop trimming proportion if 'trmse' is used as the performance metric (numeric). 'trmse' calculates the root-mean-square error of residuals (r) of which |r| < quantile(|r|, 1-trim_prop)  (e.g. trim_prop = 0.1 ignores 10\% of the data and calculates RMSE of residuals whose absolute values are below 90\% quantile of |r|). If  both \code{trim_prop} and \code{trim_c} are specified, \code{trim_c} will be used.
#' @param trim_c the trimming constant if 'trmse' is used as the performance metric (numeric, defaults to 3). 'trmse' calculates the root-mean-square error of the residuals (r) between median(r) + trim_c mad(r) and median(r) - trim_c mad(r).  If  both \code{trim_prop} and \code{trim_c} are specified, \code{trim_c} will be used.
#' @param type  type of TFBoost (character, 'L2' or 'LAD' or 'RR', defaults to 'L2')
#' @param shrinkage shrinkage parameter in boosting (numeric, defaults to 0.05)
#' @param precision number of significant digits to keep when using validation error to calculate early stopping time (numeric, defaults to 4)
#' @param init_type type of initialization for TFBoost, at the mean or median of the training responses (character, 'mean' or 'median', defaults to 'mean')
#' @param nknot number of interior knots of the cubic B-spline basis (numeric, defaults to 3)
#' @param save_f logical indicating whether to save the function estimates at all iterations (defaults to \code{FALSE})
#' @param trace logical indicating whether to print the number of completed iterations for monitoring progress (defaults to \code{FALSE})
#' @param save_tree logical indicating whether to save the tree objects at all iterations, required when the user calls \code{TFBoost.predict} with the returned object from \code{TFBoost} (defaults to \code{FALSE})
#' @param error_type a character string (or vector of character strings) indicating the type of error metrics to be evaluated on the test set. Valid options are: "rmse" (root mean squared error), "aad" (average absolute deviation), and "trmse" (trimmed root mean squared error)
#' @return A list of all input parameters
#'
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' 
#' #' @examples
#'  \dontrun{
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
#' my.control <- TFBoost.control(make_prediction = TRUE, 
#' tree_control = TREE.control(tree_type  = "B", d = 1, num_dir = 200), 
#' shrinkage = 0.05, nknot = 3, loss = "l2")
#'
#' model_TFBoost <- TFBoost(x_train = xtrain, y_train = ytrain,  x_val = xval,  y_val = yval, 
#'        x_test = xtest, y_test = ytest, grid = gg, t_range  = tt, niter = niter, 
#'        control = my.control)
#'}        
#' @export
#' 

RTFBoost.control <- function(make_prediction = TRUE, eff_m = 0.95, bb = 0.5, trim_prop = NULL, trim_c = 3,
                             max_depth_init = 3, min_leaf_size_init = 10, 
                             tree_control = TREE.control(), type = "L2", shrinkage  = 0.05, precision = 4, 
                             init_type = "median", n_init = 20, niter = 100, nknot = 3, save_f = FALSE, 
                             trace = FALSE, save_tree = FALSE, error_type = "mse"){
  
  cc_s <- as.numeric(RobStatTM::lmrobdet.control(bb=.5, family='bisquare')$tuning.chi)
  cc_m <-  as.numeric(RobStatTM::lmrobdet.control(efficiency=eff_m, family='bisquare')$tuning.psi)
  
  return(list(make_prediction =  make_prediction,  cc_s = cc_s, cc_m = cc_m, bb = bb, trim_prop = trim_prop, 
              trim_c = trim_c , max_depth_init = max_depth_init , min_leaf_size_init = min_leaf_size_init ,
              tree_control = tree_control, type = type,  shrinkage = shrinkage, precision = precision, 
              init_type = init_type,  n_init = n_init, niter = niter, nknot = nknot, save_f = save_f, trace = trace, 
              save_tree = save_tree, error_type = error_type))
}

#' Tree-based functional boosting 
#' 
#' This function implements a tree-based boosting algorithm for functional regression,
#'
#' This function implements  a tree-based boosting algorithm for functional regression (TFBoost).
#' This function uses the functions available in the \code{rpart} package to construct functional multi-index regression trees.
#' 
#' @param x_train functional predictor matrix in the training data (matrix/dataframe)
#' @param z_train scalar predictor matrix in the training data (matrix/dataframe, optional)
#' @param y_train scalar response vector in the training data (vector/dataframe)
#' @param x_val functional predictor matrix in the validation data (matrix/dataframe)
#' @param z_val scalar predictor matrix in the validation data (matrix/dataframe, optonal)
#' @param y_val scalar response vectorin the validation data (vector/dataframe)
#' @param x_test functional predictor matrix for test data (matrix/dataframe, optional, required when \code{make_prediction} in control is \code{TRUE})
#' @param z_test scalar predictor matrix for test data (matrix/dataframe, optional, required when \code{make_prediction} in control is \code{TRUE})
#' @param y_test scalar response vector for test data (vector/dataframe, optional, required when \code{make_prediction} in control is \code{TRUE} and \code{z_train} and \code{z_val} are provided)
#' @param grid common grid that the \code{x_train}, \code{x_val}, and \code{x_test} are evaluated at (vector)
#' @param t_range domain of the functional predictor, provided in the form of c(left end of the domain, right end of the domain) (vector)
#' @param niter number of boosting iterations (numeric)
#' @param control a named list of control parameters, as returned by \code{\link{TFBoost.control}}
#' 
#' @return A list with the following components:
#'
#' \item{B}{predicted values with model at the early stopping iteration using x_test (or x_test and z_test) as the predictors}
#' \item{loss_train}{a vector of training errors for all iterations}
#' \item{loss_val}{a vector of validation errors for all iterations}
#' \item{f_train_t}{predicted values with model at the early stopping iteration on the training predictors}
#' \item{f_val_t}{predicted values with model at the early stopping iteration on the validation predictors}
#' \item{f_t_test}{predicted values with model at the early stopping iteration on the test predictors (returned if make_prediction = TRUE in control)}
#' \item{early_stop}{early stopping iteration}
#' \item{err_train}{a vector of training mean-squared-errors for all iterations}
#' \item{err_val}{a vector of validation mean-squared-errors for all iterations}
#' \item{err_test}{a vector of test mean-squared-errors before and at the early stopping iteration (returned if make_prediction = TRUE in control)}
#' \item{grid}{\code{grid} form the input arguments}
#' \item{t_range}{\code{t_range} from the input arguments}
#' \item{init_vals}{a constant to initialize the function estimates for training, validation, and test sets}
#' \item{tree_obj}{a list of fitted functional multi-index trees (one per iteration and returned from \link{TREE})}
#' \item{alpha}{a vector of boosting step sizes (one per iteration)}
#' \item{control}{\code{control} from the input arguments}
#' \item{save_f_train}{a matrix of training function estimates at all iterations (returned if save_f = TRUE in control)}
#' \item{save_f_val}{a matrix of validation function estimates at all iterations (returned if save_f = TRUE in control)}
#' \item{save_f_test}{a matrix of test function estimates before and at the early stopping iteration (returned if save_f = TRUE and make_prediction = TRUE in control)}
#'
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' 
#' @examples
#' 
#'  \dontrun{
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
#' my.control <- TFBoost.control(make_prediction = TRUE, 
#' tree_control = TREE.control(tree_type  = "B", d = 1, num_dir = 200), 
#' shrinkage = 0.05, nknot = 3, loss = "l2")
#'
#' model_TFBoost <- TFBoost(x_train = xtrain, y_train = ytrain,  x_val = xval,  y_val = yval, 
#'        x_test = xtest, y_test = ytest, grid = gg, t_range  = tt, niter = niter, 
#'        control = my.control)
#' }    
#' @export

RTFBoost <- function(x_train, z_train = NULL, y_train,  x_val,  z_val = NULL, y_val, x_test, z_test = NULL, y_test, grid, t_range,  tree_init = NULL, control = RTFBoost.control()){
  

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
  }
  
  ss <- NULL; cc <- NULL
  
  make_prediction <- control$make_prediction
  trim_prop <- control$trim_prop; trim_c <- control$trim_c 
  tree_control <- control$tree_control
  type <- control$type  
  trace <- control$trace
  nknot <- control$nknot
  save_f <- control$save_f
  save_tree <- control$save_tree
  n_init <- control$n_init
  niter <- control$niter 
  error_type <- control$error_type
  shrinkage <- control$shrinkage
  precision <- control$precision
  init_type <- control$init_type
  max_depth_init <- control$max_depth_init
  min_leaf_size_init <- control$min_leaf_size_init
  
  if(type == "RR"){
    cc <- control$cc_s
    cc_m <- control$cc_m
    bb <- control$bb
  }

  n_train <- nrow(x_train)
  n_val <- nrow(x_val)
  
  if(save_f){
    save_f_train <- matrix(NA, n_train, niter)
    save_f_val <- matrix(NA, n_val, niter )
  }
  
  if(missing(t_range)){
     t_range <- c(min(grid), max(grid))
  }
  
  if(!missing(z_train)){
    if(!is.matrix(z_train)){
      z_train = as.matrix(z_train, dimnames = list(NULL, names(z_train)))
    }
    if(!is.matrix(z_val)){
      z_val = as.matrix(z_val, dimnames = list(NULL, names(z_val)))
    }
  }
  

  dd <- 4; p <- ncol(x_train)
  grid0 <- seq(t_range[1],t_range[2], 1/(10*(p-1))) # in case of not evenly spaced
  knot <- quantile(grid0, (1:nknot)/(nknot+1) )
  delta <- sort(c(rep(range(grid0), dd), knot)) #exterior knots
  B<- spline.des(delta, grid, dd)$design
  B <- compute.orthonormal(B,grid, t_range)
  train_predictors <- t(apply(x_train, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))
  val_predictors <- t(apply(x_val, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))
  
  if(missing(x_test)) {
    make_prediction <- FALSE 
  } 
  
  tree.obj <- list()
  alpha <- rep(NA, niter)
  
  loss_train <- loss_val <- rep(NA, niter)
  err_train <- err_val <- rep(NA, niter)
  
  # initialize functions (loss and gradient)
  init_tmp <- init.boost(type) 
  func <- init_tmp$func
  func.grad <- init_tmp$func.grad
  # initialization
  if(init_type == "LADTree"){
    
      if(is.null(z_train)){
        dat_tmp <- data.frame(x_train, y_train = y_train)
      }else{
        dat_tmp <- data.frame(x_train, y_train = y_train, z_train = z_train)
      }
    
    if(is.null(tree_init)){
      tree_init <- TREE.init(x = train_predictors, y = y_train, z = z_train, newx = val_predictors, newy = y_val, 
                             newz = z_val, random.seed = 0, 
                             max_depth_init = max_depth_init,  
                             min_leaf_size_init = min_leaf_size_init,
                             num_dir = tree_control$num_dir, make_prediction = make_prediction)
    }

    f_train_early <- f_train_t <-  TREE.init.predict(tree_init, newx = train_predictors, newz = z_train)$pred
    f_val_early <- f_val_t <-  TREE.init.predict(tree_init, newx = val_predictors, newz = z_val)$pred
    model <- list( tree_init =  tree_init)
  }
  
  if(init_type == "median"){
    f_train_init <- median(y_train)
    f_train_early <-   f_train_t <- rep(f_train_init, length(y_train))
    f_val_early <- f_val_t <- rep(median(y_train), length(y_val))
    model <- list( f_train_init = f_train_init)
  }
  
  if(init_type == "mean"){
    f_train_init <- mean(y_train)
    f_train_early <-   f_train_t <- rep(f_train_init, length(y_train))
    f_val_early <- f_val_t <- rep(mean(y_train), length(y_val))
    model <- list(f_train_init = f_train_init)
  }
  
  init_status <- 0
  
  for(i in 1:niter){
    
    if(save_f){
      save_f_train[,i] <- f_train_t
      save_f_val[,i] <- f_val_t
    }
    
    if(type == "RR" & init_status == 0) {
      ss <- cal.ss.rr(f_train_t, y_train,  cc, bb)
    }
    
    if(trace){
      if(i%%100 ==0 )
      print(paste(i,"th iteration out of", niter, "iterations"))
    }
  
    u <- as.numeric(cal.neggrad(type, y_train,f_train_t, func.grad, init_status, ss, cc))
    tree.obj[[i]] <- TREE(x = train_predictors, y = u, z = z_train, newx = val_predictors, newz = z_val, random.seed = i, control = tree_control)
    h_train_t <-tree.obj[[i]]$pred_train
    h_val_t <- tree.obj[[i]]$pred_test

    alpha[i] <- cal.alpha(f_train_t = f_train_t, h_train = h_train_t, y_train = y_train, func = func, type = type,
                          init_status = init_status, ss = ss, cc = cc)
    f_train_t <-  f_train_t + shrinkage*alpha[i]*h_train_t
    f_val_t <-  f_val_t + shrinkage*alpha[i]*h_val_t
    #err_train[i] <- mse(f_train_t, y_train)
    #err_val[i] <- mse(f_val_t, y_val)
    err_train[i] <- mean(abs(f_train_t - y_train))
    err_val[i] <- mean(abs(f_val_t - y_val))
    
    if(type == "RR"){
      if(init_status == 0){
        loss_val[i] <-cal.ss.rr(f_val_t, y_val,  cc, bb)
        loss_train[i] <-ss   
      }else{
        loss_train[i] <- mean(func((f_train_t - y_train)/ss, cc = cc_m))
        loss_val[i] <- mean(func((f_val_t - y_val)/ss, cc = cc_m))
      }
    }else{
        loss_train[i] <- mean(func(f_train_t - y_train))
        loss_val[i] <- mean(func((f_val_t - y_val)))
    }

   if(i == 1){
      when_init <- 1
      early_stop <- 1
      f_train_early <- f_train_t
      f_val_early <- f_val_t
    }else{
      if(i<=n_init){
        if(round(loss_val[i], precision) < min(round(loss_val[1:(i-1)], precision))){
          when_init <- i
          early_stop <- i
          f_train_early <- f_train_t
          f_val_early <- f_val_t
        }
      }else{
        if(round(loss_val[i], precision) < min(round(loss_val[(n_init):(i-1)], precision))){
          early_stop <- i
          f_train_early  <- f_train_t
          f_val_early <- f_val_t
        }
      }
    }
    
    if((type == "RR") && (i == n_init)){
      init_status <- 1
      f_train_t <- f_train_early  # rest the current one
      f_val_t <- f_val_early
      ss <-  RobStatTM::mscale(f_train_t - y_train,  tuning.chi= cc, delta = bb)
      cc <- cc_m
      loss_val[i] <- mean(func((f_val_t - y_val)/ss, cc = cc_m))
    }
  }
  
  f_train_t <- f_train_early
  f_val_t <- f_val_early
  
  
  model <- c(model, list(B = B, loss_train=loss_train, loss_val = loss_val, f_train_t =  f_train_t, f_val_t = f_val_t, when_init = when_init, early_stop = early_stop, 
                err_train = err_train,  err_val = err_val, grid = grid,  t_range = t_range, tree.obj = tree.obj, 
                alpha = alpha, control = control))
  
  if(make_prediction){
    tmp_predict <- RTFBoost.predict(model, newx = x_test, newy = y_test, newz = z_test)
    model <- c(model, list(f_test_t = tmp_predict$pred, err_test = tmp_predict$err_test))
    if(save_f){
      model <- c(model, list(save_f_test = tmp_predict$save_f_test))
    }
  }
      
  if(save_f){
    model <- c(model, list(save_f_train = save_f_train, save_f_val = save_f_val))
  }
  if(!save_tree){
    model$tree.obj <- NULL
  }
  return(model)
}
    
#' TFBoost.predict
#'
#' A function to make predictions and calculate test error given an object returned by TFBoost and test data
#'
#' A function to make predictions and calculate test error given an object returned by TFBoost and test data
#'
#'@param model an object returned by TFBoost
#'@param newx functional predictor matrix for test data (matrix/dataframe)
#'@param newy scalar predictor matrix for test data (matrix/dataframe, optional, requires if model is trained with z_train and z_val)
#'@param newz scalar response vector for test data (vector/dataframe)
#'@return A list with with the following components:
#'
#' \item{f_t_test}{predicted values with model at the early stopping iteration using x_test (or x_test and z_test) as the predictors}
#' \item{err_test}{a vector of test errors before and at the early stopping iteration (returned if newy is provided)}
#' \item{save_f_test}{a matrix of test function estimates at all iterations (returned if save_f = TRUE in control)}
#'
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#' 
#' #' #' @examples
#'  \dontrun{
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

RTFBoost.predict <- function(model, newx, newy, newz = NULL){

  control <- model$control

  if(control$save_f){
    save_f_test <- matrix(NA, nrow(newx), model$early_stop)
  }

  test_predictors <- t(apply(newx, 1, function(xx){apply(model$B, 2, function(bb){riemman (xx*bb, model$grid, model$t_range)})}))
  
  if(control$init_type == "LADTree"){
    f_test_t  <- TREE.init.predict(model$tree_init, newx = test_predictors, newz = z_test)$pred
  }else{
    f_test_t  <- model$f_train_init
  }
  err_test <- data.frame(matrix(NA, nrow = model$early_stop, ncol = length(control$error_type)))
  colnames(err_test) <- control$error_type
  
  if((control$type == "RR") && (control$n_init < model$early_stop)){
    for(i in 1: model$when_init){
      if(control$save_f){
        save_f_test[,i] <- f_test_t
      }
      f_test_t <- f_test_t + control$shrinkage*model$alpha[i]* TREE.predict(model$tree.obj[[i]], newx =  test_predictors, newz = newz)$pred
      
      for(err_type in control$error_type){
        err_test[i,err_type] <- cal_error(control$trim_prop, control$trim_c, err_type, f_test_t, y_test)
      }

    }
    for(i in (control$n_init+1):model$early_stop){
      if(control$save_f){
        save_f_test[,i] <- f_test_t
      }
      f_test_t <- f_test_t + control$shrinkage*model$alpha[i]* TREE.predict(model$tree.obj[[i]], newx =  test_predictors, newz = newz)$pred
      
      if(!missing(newy)){
        for(err_type in control$error_type){
          err_test[i,err_type] <- cal_error(control$trim_prop, control$trim_c, err_type, f_test_t, y_test)
        }   
      }
    }
  }else{
    for(i in 1: model$early_stop){
      if(control$save_f){
        save_f_test[,i] <- f_test_t
      }
      f_test_t <- f_test_t + control$shrinkage*model$alpha[i]* TREE.predict(model$tree.obj[[i]], newx =  test_predictors, newz = newz)$pred
      
      if(!missing(newy)){
        for(err_type in control$error_type){
          err_test[i,err_type] <- cal_error(control$trim_prop, control$trim_c, err_type, f_test_t, y_test)
        }   
      }
    }
  }
  
  if((missing(newy)) & (control$save_f == FALSE)){
    return(f_test_t)
  }else{
    res <- list(pred = f_test_t)
    if(!missing(newy)){
      res <- c(res, list(err_test = err_test))
    }
    if(control$save_f){
      res <- c(res,  list(save_f_test =  save_f_test))
    }
    return(res)
  }
}


RTFBoost.validation <- function(functionx_train, z_train = NULL, y_train,  x_val,  z_val = NULL, y_val, x_test, z_test = NULL, y_test, grid, t_range, 
                                max_depth_init_set = c(1,2,3,4), min_leaf_size_init_set = c(10,20,30), control = RTFBoost.control()){

  control_tmp <- control

  control_tmp$init_type <- "median"
  model_best <- RTFBoost(x_train = x_train, y_train = y_train,  x_val = x_val,  y_val = y_val,
                         x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                         control = control_tmp)
  
  
  flagger_outlier <- which(abs(model_best$f_val_t - y_val)>3*mad(model_best$f_val_t - y_val))
  
  if(length(flagger_outlier)>=1){
    best_err <- mean(abs(model_best$f_val_t[-flagger_outlier] - y_val[-flagger_outlier]))  #test with tau-scale
  }else{
    best_err <- mean(abs(model_best$f_val_t - y_val))
  }
  
  params = c(0,0)
  errs_val <-  rep(NA, 1+ length(min_leaf_size_init_set)*length(max_depth_init_set))
  errs_test <- matrix(NA, 1+ length(min_leaf_size_init_set)*length(max_depth_init_set), length(control$error_type))
  errs_val[1] <- best_err
  
  if(control$make_prediction){
    errs_test[1,] <- as.numeric(model_best$err_test[control$niter,])
  }
  
  if(control$init_type == "LADTree") {
    model_pre_tree <- NA
    combs <- expand.grid(min_leafs= sort(min_leaf_size_init_set,TRUE), max_depths= max_depth_init_set)
    j_tmp <- rep(1, nrow(combs)) #need to consider that tree or not
    
    tree_init_list <- list()
    for(j in 1:nrow(combs)) {
      min_leaf_size <- combs[j, 1]
      max_depths <- combs[j, 2]
      if(is.null(z_train)){
        dat_tmp <- data.frame(x_train, y_train = y_train)
      }else{
        dat_tmp <- data.frame(x_train, y_train = y_train, z_train = z_train)
      }
      
      tree_init_list[[j]] <- TREE.init(x = train_predictors, y = y_train, z = z_train, newx = val_predictors, newy = y_val, 
                             newz = z_val, random.seed = 0, 
                             max_depth_init = max_depths,  
                             min_leaf_size_init = min_leaf_size,
                             num_dir = control$tree_control$num_dir, make_prediction = make_prediction)
     }
    
    for(j in 1:length(max_depth_init_set)){
      for(k in 1:(length(min_leaf_size_init_set)-1)){
        idx_jk <- which(combs[,1] == sort(min_leaf_size_init_set, TRUE)[k] & combs[,2] == max_depth_init_set[j])
        idx_jk_plus<- which(combs[,1] == sort(min_leaf_size_init_set, TRUE)[k+1] & combs[,2] == max_depth_init_set[j])
        equal_tmp<- all.equal(tree_init_list[[idx_jk]]$tree.model,tree_init_list[[idx_jk_plus]]$tree.model) == TRUE
        if(length(equal_tmp)==2 & sum(equal_tmp) == 0){
          j_tmp[ idx_jk_plus] <- 0
        }
      }
    }
    
    for(k in 1:length(min_leaf_size_init_set)){
      for(j in 1:(length(max_depth_init_set)-1)){
        idx_kj <- which(combs[,1] == sort(min_leaf_size_init_set, TRUE)[k] & combs[,2] == max_depth_init_set[j])
        idx_kj_plus<- which(combs[,1] == sort(min_leaf_size_init_set, TRUE)[k] & combs[,2] == max_depth_init_set[j+1])
        equal_tmp<- all.equal(tree_init_list[[idx_kj]]$tree.model,tree_init_list[[idx_kj_plus]]$tree.model) == TRUE
        if(length(equal_tmp)==2 & sum(equal_tmp) == 0){
          j_tmp[idx_kj_plus] <- 0
        }
      }
    }
    
    print(j_tmp)
    
    err_cvs <- matrix(NA, length(j_tmp))
    
    for(j in 1:nrow(combs)) {
      
      print(j)
      if(j_tmp[j] == 1){
        
        min_leaf_size <- combs[j, 1]
        max_depths <- combs[j, 2]
        control_tmp$max_depth_init <- max_depths
        control_tmp$min_leaf_size_init  <- min_leaf_size
        
        model_tmp <- RTFBoost(x_train = x_train, y_train = y_train,  x_val = x_val,  y_val = y_val,
                              x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, tree_init= tree_init_list[[j]],
                              control = control)
  
        if(length(flagger_outlier)>=1){
          err_tmp <- mean(abs(model_tmp$f_val_t[-flagger_outlier] - y_val[-flagger_outlier]))
        }else{
          err_tmp <- mean(abs(model_tmp$f_val_t - y_val))
        }
        
        errs_val[j+1] <- err_tmp
        
        if(control$make_prediction){
          errs_test[j+1,] <- as.numeric(model_tmp$err_test[control$niter,])
        }
        
        if(control$trace){
          print(paste("leaf size:", min_leaf_size, " depths:", max_depths, " err(val):", round(err_tmp,4), " best err(val) :", round(best_err,4) ,sep = ""))
        }
        err_cvs[j] <- err_tmp 
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
  
 
  model_best$err_cvs <- err_cvs
  model_best$params = params
  
  return(model_best)
}