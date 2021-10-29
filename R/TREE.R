#' Tuning and control parameters for the functional multi-index tree
#' 
#' Tuning and control parameters for the functional multi-index tree. 
#'
#' Various tuning and control parameters for the functional multi-index tree algorithm implemented in the
#' function \code{\link{TREE}}
#' 
#' @param tree.type  a character ("A" or "B") denoting Type A or Type B tree (defaults to "A").
#' @param tree.nindex an integer specifying the number of indices required by a Type A tree (defaults to 1).
#' @param max.depth an integer specifying the maximum depth for a Type A or Type B tree (defaults to 1).
#' @param num.dir an integer specifying the number of random directions required by the initial Type B tree (defaults to 20).
#' @param min.leafsize an integer specifying the minimum number of observations per leaf node for a Type A or Type B tree (defaults to 2).
#' @param nmulti  an integer specifying  the number of random initialization points to fit a Type A tree (defaults to 3). 
#' @param nscreen an integer specifying the number of random points from which to select the initial points to fit a Type A tree (defaults to 30). 
#' 
#' @return A list of all input parameters
#'
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' 
#' @export
#' 
TREE.control <- function(tree.type = "A",  tree.nindex = 1, max.depth = 1, num.dir = 20,  min.leafsize = 2,  nmulti= 3, nscreen = 30){
  
  return(list(tree.type = tree.type, tree.nindex = tree.nindex,  max.depth = max.depth, num.dir = num.dir,  min.leafsize = min.leafsize,  nmulti= nmulti, nscreen = nscreen))
}

#' Functional multi-index tree 
#' 
#' This function implements a algorithm for functional multi-index tree.
#'
#' This function implements a functional multi-index tree algorithm developed based on functions available in the \code{rpart} package. 
#' 
#' @param x a matrix or data.frame containing basis projections of the functional predictor in the training data. 
#' @param z a matrix or data.frame containing values of the scalar predictors in the training data (optional, defaults to \code{NULL}). 
#' @param y a vector containing values of the scalar response in the training data (optional, defaults to \code{NULL}).
#' @param random.seed an integer as the seed to generate random directions for the Type B tree. 
#' @param control a named list of control parameters, as returned by \code{\link{TREE.control}}
#' 
#' @return A list with the following components:
#' \item{beta.opt}{a matrix of index coefficients estimated by the Type A tree, returned if \code{control$tree.type = "A"}}
#' \item{betas.selected}{a matrix of index coefficients selected  by the Type B tree, returned if \code{control$tree.type = "B"}}
#' \item{pred}{a vector of predicted values for the training data}
#' \item{tree.model}{an rpart object of the fitted tree}
#' \item{control}{\code{control} from the input arguments}
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' 
#' @export

TREE <-function(x, y, z, random.seed, control = TREE.control()) {
  
  
  if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
  }
  
  p <- ncol(x)

  tree.nindex <- control$tree.nindex
  if(control$tree.type == "A"){
    
    fval.min <- Inf
    set.seed(random.seed)
    optim.fn <-  function(param){tree.loss(param, x, y, z, d = control$max.depth,  minbucket = control$min.leafsize, control$tree.nindex, type = "theta")}
    start_save <- screen.dir(control$nmulti, control$nscreen, p, tree.nindex , optim.fn, max_iter = 10)$init_points

    for(i in 1:control$nmulti) {
      
      param <- start_save[((i-1)*tree.nindex  +1):(i*tree.nindex ),] 
      optim.return <- tryCatch(Nelder_Mead(fn = optim.fn, par= as.numeric(param), lower = c(rep(-pi/2,tree.nindex ), rep(0, (p-2)*tree.nindex )), upper = c(rep(pi/2,tree.nindex ), rep(pi, (p-2)*tree.nindex)), control = list(maxfun = 1000)), warning = function(a){return(1)})
      
      while(class(optim.return) == "numeric"){
        print("resample!")
        param <- start_save[((i-1)*tree.nindex  +1):(i*tree.nindex),] <-  sample.dir(1, p, tree.nindex , type = "theta")
        optim.return <- tryCatch(Nelder_Mead(fn = optim.fn, par = as.numeric(param), lower = c(rep(-pi/2,tree.nindex ), rep(0, (p-2)*tree.nindex )), upper = c(rep(pi/2,tree.nindex ), rep(pi, (p-2)*tree.nindex )), control = list(maxfun = 1000)), warning = function(a){return(1)}) 
      }
      
      par.tmp <- matrix(optim.return$par, ncol = p-1)
      
      if(optim.return$fval < fval.min) {
          if(tree.nindex == 1){
            beta.opt <-  as.matrix(p_to_c(par.tmp))
          }else{
            beta.opt <- apply(par.tmp, 1, p_to_c) #ncol: number of index, nrow: p
          }
          fval.min <- optim.return$fval
        }
      }
    
    if(is.null(z)){
      index <- x %*% beta.opt
    }else{
      index <- cbind(x %*% beta.opt, z)
      colnames(index) <-  c(paste("X", 1:tree.nindex, sep = ""), colnames(z))
    }
    dat <- data.frame(index, y = y)
    
    tree.model <- rpart(y~., data = dat, control = rpart.control(maxdepth = control$max.depth, cp = 0, minbucket =  control$min.leafsize))
    model <- list(pred = predict(tree.model), beta.opt = beta.opt, control = control, tree.model = tree.model)

  }
  
  if(control$tree.type == "B"){
    
    set.seed(random.seed)
    betas <-  t(sample.dir(control$num.dir, p, num_index = 1, type = "beta"))
    
    if(is.null(z)){
      index <- x %*% betas
    }else{
      index <- cbind(x %*% betas, z)
      index <- matrix(index, ncol = control$num.dir + ncol(z))
      colnames(index) <-  c(paste("X", 1:control$num.dir, sep = ""), colnames(z))
    }
    dat <- data.frame(index, y = y)
    
    tree.model <- rpart(y~., data = dat, control = rpart.control(maxdepth = control$max.depth,  cp = 0, minbucket = control$min.leafsize))
    betas <-data.frame(betas)
    used.var <- setdiff(tree.model$frame$var, "<leaf>")
    used.var <- setdiff(used.var,colnames(z))
    
    if(length(used.var)==1){
      betas.selected <-data.frame(betas[, used.var])
      names(betas.selected) <- used.var 
    }else{
      betas.selected <-betas[, used.var]
    }
    
    model <- list(pred = predict(tree.model), betas.selected =  betas.selected, control = control, tree.model = tree.model)
  
  }
  return(model)
}

#' TREE.predict
#'
#' A function to make predictions given an object returned by TREE and test data
#'
#' A function to make predictions given an object returned by TREE and test data
#'
#'@param model an object returned by \code{TREE}
#'@param newx  a matrix or data.frame containing basis projections of the functional predictor in the new data. 
#'@param newz  a matrix or data.frame containing values of the scalar predictors in the new data (optional, defaults to \code{NULL}, required if \code{z} was used to fit \code{model}). 
#'@return 
#'
#' \item{pred}{a vector of predicted values using \code{newx} (or \code{newx} and \code{newz}) as the predictors}
#'
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#' 
#' @export
#' 
TREE.predict <- function(model, newx, newz){

    if(model$control$tree.type == "A"){
      if(is.null(newz)){
        index <- as.matrix(newx %*% model$beta.opt)
      }else{
        if(!is.matrix(newz)){
          newz = as.matrix(newz, dimnames = list(NULL, names(newz)))
        }
        index <- cbind(as.matrix(newx %*% model$beta.opt), newz)
        colnames(index)[1:ncol(model$beta.opt)] <- c(paste("X", 1:ncol(model$beta.opt), sep = ""))
      }
    }

    if(model$control$tree.type == "B"){
      betas <- data.frame(matrix(0, ncol(newx), model$control$num.dir))
      betas[, names(model$betas.selected)] <- model$betas.selected
      if(is.null(newz)){
        index <- newx %*% as.matrix(betas)
      }else{
        if(!is.matrix(newz)){
          newz = as.matrix(newz, dimnames = list(NULL, names(newz)))
        }
        index <- cbind(newx %*% as.matrix(betas), newz)
        colnames(index)[1:model$control$num.dir] <- c(paste("X", 1:model$control$num.dir, sep = ""))
      }
    }
   
    dat <- data.frame(index)
    pred <- predict(model$tree.model, newdata = dat)
    return(pred)
}


p_to_c <- function(theta){
  
  d <- length(theta)+1
  beta <- rep(NA, d)
  beta[1] <-cos(theta[1])
  tmp1 <- cumprod(sin(theta))
  
  if(d>2){
    for(i in 2:(d-1)){
      beta[i] <- cos(theta[i])*tmp1[i-1]
    }
  }
  beta[d] <- tmp1[d-1]
  
  return(beta)
}


c_to_p <- function(beta){
  
  d <- length(beta)
  theta <- rep(NA, d-1)
  tmp <- beta[d]
  
  for(i in (d-1):1){
    theta[i] <- atan(tmp/beta[i])
    if(theta[i] <0 & i> 1){
      theta[i] <-  theta[i] + pi
    }
    tmp <- tmp/sin(theta[i])
  }
  
  return(theta)
}



tree.loss<- function(param, x, y, z, d, minbucket, num_index, type = "theta") {
  
  param <- matrix(param, nrow = num_index) # every row is an index 
  
  if(type == "theta"){
    if(num_index == 1){
      beta <- p_to_c(param)
    }else{
      beta <- apply(param, 1,p_to_c) # apply to rows
    }
  }else{
    beta <- param
  }
  
  loss <- function(x, y, z,beta) {
    if(!is.null(z)){
      index <- cbind(x %*% beta, z)
    }else{
      index <- x %*% beta
    }
    dat <- data.frame(index, y = y)
    tree.model <- rpart(y~.,  data = dat, control = rpart.control(maxdepth = d,cp = 0, minbucket = minbucket))
    t.ret <- mean((y-predict(tree.model))^2)
    return(t.ret)
  }
  return(loss(x,y,z,beta))
}

sample.dir <- function(num_sample, p, num_index, type = "beta"){
  
  if(type == "theta"){
    tmppp <- matrix(runif(num_sample*(p-1)*num_index, 0,pi), nrow = num_sample*num_index)
    tmppp[,1] <- tmppp[,1] - pi/2
  }else{
    tmp <- matrix(rnorm(num_sample*p*num_index, 0,1), nrow = num_sample*num_index)
    tmpp <- apply(tmp, 1, function(x) {sqrt(sum(x^2))})
    tmppp <- apply(tmp, 2, function(x){x/tmpp})
    tmppp[,1] <- abs(tmppp[,1])
  }
  pars <- tmppp
  return(pars)
}

screen.dir <- function(nmulti, nscreen, p, num_index, optim.fn,  max_iter = 10){
  
  tmp <-  as.matrix(sample.dir(nscreen, p, num_index, type = "theta"))
  
  if(p == 2|| (num_index == 1 & nscreen == 1)){
    tmp <- t(tmp)
  }
  
  cal_nelder_mead <- function(i){
    par_val <- as.numeric(tmp[((i-1)*num_index +1):(i*num_index),])
    suppressWarnings(tmp.return <-Nelder_Mead(fn = optim.fn, par = par_val , lower = c(rep(-pi/2,num_index), rep(0, (p-2)*num_index)), upper = c(rep(pi/2,num_index), rep(pi, (p-2)*num_index)), control = list(maxfun = max_iter)))
    return(tmp.return) 
  }
  
  list_tmp <- lapply(1:nscreen, cal_nelder_mead)
  loss_sampled <- unlist(lapply(list_tmp, function(x){x$fval}))
  idx <- which(rank( loss_sampled, ties.method = "first") <= nmulti)
  idx_tmp  <- as.numeric(sapply(idx, function(x){((x-1)*num_index+1):(x*num_index)}))
  
  return(list(init_points = matrix(tmp[idx_tmp,], ncol = p-1)))
}


# set.seed(123)
# n <- 200; p <- 5
# x<- matrix(runif(n*p), ncol = p)
# y <- x[,1] + x[,2] + rnorm(n)
# train_idx <- sample(n,50)
# test_idx <- setdiff(1:n, train_idx)
# model <- TREE(x, y, newx = x[1:10,], nmulti = 3, random.seed = 42, control = TREE.control(d = 3, num_dir = 1000, method = "method3"))
# model$test_mse
