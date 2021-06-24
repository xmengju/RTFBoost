TREE.init <- function(x, y, z, newx, newy, newz, random.seed, max_depth_init,  min_leaf_size_init, num_dir, make_prediction){
  
  if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
  }
  
  set.seed(random.seed)
  p <- ncol(x)
  betas <-  t(sample.dir(num_dir, p, num_index = 1, type = "beta"))
  
  if(is.null(z)){
    index <- x %*% betas
  }else{
    index <- cbind(x %*% betas, z)
    index <- matrix(index, ncol = num_dir + ncol(z))
    colnames(index) <-  c(paste("X", 1:num_dir, sep = ""), colnames(z))
  }
  dat <- data.frame(index, y = y)
  tree.model <- rpart(y~., data = dat, control = rpart.control(maxdepth = max_depth_init, minbucket = min_leaf_size_init, xval = 0, cp = -Inf), method = alist)
  betas <-data.frame(betas)
  used.var <- setdiff(tree.model$frame$var, "<leaf>")
  used.var <- setdiff(used.var,colnames(z))
  
  if(length(used.var)==1){
    betas_selected <-data.frame(betas[, used.var])
    names(betas_selected) <- used.var 
  }else{
    betas_selected <-betas[, used.var]
  }
  
  model <- list(betas_selected =  betas_selected, pred_train = predict(tree.model), num_dir = num_dir, tree.model = tree.model)
  
  if(make_prediction){
    if(!missing(newx)){
      tmp <-  TREE.init.predict(model, newx = newx, newz = newz)
      model$pred_test <- tmp$pred
      if(!missing(newy)){
        model$test_mse <- mean((model$pred_test  - newy)^2)
      }     
    }
  }
  return(model)
}

TREE.init.predict <- function(model, newx, newz){

    betas <- data.frame(matrix(0, ncol(newx), model$num_dir))
    betas[, names(model$betas_selected)] <- model$betas_selected
    if(is.null(newz)){
      index <- newx %*% as.matrix(betas)
    }else{
      if(!is.matrix(newz)){
        newz = as.matrix(newz, dimnames = list(NULL, names(newz)))
      }
      index <- cbind(newx %*% as.matrix(betas), newz)
      colnames(index)[1:model$num_dir] <- c(paste("X", 1:model$num_dir, sep = ""))
  }

  
  dat <- data.frame(index)
  pred <- predict(model$tree.model, newdata = dat)
  return(list(pred = pred))
}
