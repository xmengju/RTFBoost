tree.init.loss <- function(param, x, y, z, d, minbucket, num_index, type = "theta") {
  
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
    tree.model <-  rpart(y~., data = dat, control = rpart.control(maxdepth = d, minbucket = minbucket, xval = 0, cp = -Inf), method = alist)
    #t.ret <- mean((y-predict(tree.model))^2)
    t.ret <- mean(abs(y-predict(tree.model)))
    
    return(t.ret)
  }
  return(loss(x,y,z,beta))
}


TREE.init <- function(x, y, z, newx, newy, 
                      newz, random.seed = 0, 
                      tree.type, 
                      tree.nindex, 
                      max.depth, 
                      num.dir,
                      min.leafsize,
                      make.prediction,  
                      init.nmulti = 3){
  
  if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
  }
  
  set.seed(random.seed)
  p <- ncol(x)
  
  if(tree.type == "A"){
    
    fval.min <- Inf
    set.seed(random.seed)
    optim.fn <-  function(param){tree.init.loss(param, x, y, z, d = max.depth, minbucket = min.leafsize, num_index = tree.nindex, type = "theta")}
    start_save <- screen.dir(init.nmulti, nscreen = 30, p = p, num_index = tree.nindex, optim.fn, max_iter = 10)$init_points
    
    for(i in 1:init.nmulti) {
      
      param <- start_save[((i-1)*tree.nindex +1):(i*tree.nindex),] 
      optim.return <- tryCatch(Nelder_Mead(fn = optim.fn, par= as.numeric(param), lower = c(rep(-pi/2,tree.nindex), rep(0, (p-2)*tree.nindex)), upper = c(rep(pi/2,tree.nindex), rep(pi, (p-2)*tree.nindex)), control = list(maxfun = 1000)), warning = function(a){return(1)})
      
      while(class(optim.return) == "numeric"){
        print("resample!")
        param <- start_save[((i-1)*tree.nindex +1):(i*tree.nindex),] <-  sample.dir(1, p, tree.nindex, type = "theta")
        optim.return <- tryCatch(Nelder_Mead(fn = optim.fn, par = as.numeric(param), lower = c(rep(-pi/2,tree.nindex), rep(0, (p-2)*tree.nindex)), upper = c(rep(pi/2,tree.nindex), rep(pi, (p-2)*tree.nindex)), control = list(maxfun = 1000)), warning = function(a){return(1)}) 
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
    
    tree.model <- rpart(y~., data = dat,  control= rpart.control(maxdepth = max.depth, minbucket =  min.leafsize, xval = 0, cp = -Inf), method = alist)
    model <- list(beta.opt = beta.opt, tree.model = tree.model)
  }
  
  if(tree.type == "B"){
    
    betas <-  t(sample.dir(num.dir, p, num_index = 1, type = "beta"))
    
    if(is.null(z)){
      index <- x %*% betas
    }else{
      index <- cbind(x %*% betas, z)
      index <- matrix(index, ncol = num.dir + ncol(z))
      colnames(index) <-  c(paste("X", 1:num.dir, sep = ""), colnames(z))
    }
    
    dat <- data.frame(index, y = y)
    tree.model <- rpart(y~., data = dat, control = rpart.control(maxdepth = max.depth, minbucket = min.leafsize, xval = 0, cp = -Inf), method = alist)
    betas <-data.frame(betas)
    used.var <- setdiff(tree.model$frame$var, "<leaf>")
    used.var <- setdiff(used.var, colnames(z))
    
    if(length(used.var)==1){
      betas.selected <-data.frame(betas[, used.var])
      names(betas.selected) <- used.var
    }else{
      betas.selected <-betas[, used.var]
    }
    
    model <- list(betas.selected = betas.selected, num.dir = num.dir, tree.model = tree.model)
  }

  model$tree.type = tree.type 
  return(model)
}

TREE.init.predict <- function(model, newx, newz){

    if(model$tree.type == "A"){
      if(is.null(newz)){
        index <- as.matrix(newx %*% model$beta.opt)
      }else{
        index <- cbind(as.matrix(newx %*% model$beta.opt), newz)      
        colnames(index)[1:ncol(model$beta.opt)] <- c(paste("X", 1:ncol(model$beta.opt), sep = ""))

      }
    }
  
    if(model$tree.type == "B"){
      betas <- data.frame(matrix(0, ncol(newx), model$num.dir)) # default names
      betas[, names(model$betas.selected)] <- model$betas.selected
      
      if(is.null(newz)){
        index <- newx %*% as.matrix(betas)
      }else{
        if(!is.matrix(newz)){
          newz = as.matrix(newz)
        }
        index <- cbind(newx %*% as.matrix(betas), newz)
        colnames(index)[1:model$num.dir] <- c(paste("X", 1:model$num.dir, sep = ""))
        
      }
    }
  
  dat <- data.frame(index)
  pred <- predict(model$tree.model, newdata = dat)
  return(pred = pred)
}
