
eigen_fun <- function(j){
  if(j == 1){
    fun <- function(t){
      rep(1,length(t))
    }
  }else{
    fun <- function(t){
      sqrt(2)*cos( (j-1)*pi*t)
    }
  }
  return(fun)
}

gen_x <- function(n, type = "clean"){
  
  t <- seq(0, 1, length.out = 100)
  X <- matrix(0, nrow = n, ncol =100)
  for(j in 1:50){
    X <- X + t(matrix(eigen_fun(j)(t))%*% t(matrix(rnorm(n, 0, 1/j^2))))
  }
  return(X)
}

gen_beta <- function(){
  t <- seq(0, 1, length.out = 100)
  beta <- rep(0, nrow = n, ncol =100)
  beta_coef <- c(0.3, 4*(-1)^{(2:50)+1}*(2:50)^{-2})
  for(j in 1:50){
    beta <- beta + eigen_fun(j)(t)*beta_coef[j]
  }
  return(beta)
}

generate_data <- function(n, type = "clean"){
  
  X <- gen_x(n, type = type)
  beta <- gen_beta()
  Y <- as.numeric(X%*%matrix(beta)) + rnorm(n)
  idx_train <- sample(1:n, round(0.6*n))
  idx_val <- sample(setdiff(1:n, idx_train), round(0.2*n))
  idx_test <- setdiff(1:n, c(idx_train, idx_val))
  x_train <- X[idx_train,];x_val <- X[idx_val,]; x_test <- X[idx_test,]
  y_train <- Y[idx_train]; y_val <- Y[idx_val]; y_test <- Y[idx_test]
  
  return(list(x = list(x_train = x_train, x_val= x_val, x_test = x_test), 
              y = list(y_train = y_train, y_val= y_val, y_test = y_test)))
}

conduct_methods <- function(x_train, y_train, x_val, y_val, x_test, y_test, d = 2, precision = 4, nknot = 3, n_init = 500,
                            shrinkage = 0.05; niter = 1000){
  
  control.tree.list <- set.control(d, precision,  shrinkage, nknot, n_init, niter)
  
  model.l2 <- RTFBoost(x_train = x_train, y_train = y_train,  x_val = x_val,  y_val = y_val,
                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                       control = control.tree.list[[1]])
  model.lad <- RTFBoost(x_train = x_train, y_train = y_train,  x_val = x_val,  y_val = y_val,
                        x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                        control = control.tree.list[[2]])
  
  control.rr <-  RTFBoost.control(make_prediction = TRUE, eff_m = 0.95, bb = 0.5, 
                                  tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                  type = "RR", shrinkage  = shrinkage, precision = precision, 
                                  init_type = "median", niter = niter, n_init = n_init,
                                  nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  model.rr <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val, y_val = y_val,
                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                       control = control.rr)
  
  return(list(model.l2= model.l2, model.lad = model.lad,  model.rr  =  model.rr ))
}

