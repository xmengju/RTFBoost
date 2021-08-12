### check the selection of optimal initial tree 
g_func_no = 3; type = "C2"; d <- 1; SNR = 5; nknot <- 3; seed <- 3

dir <- paste("Results/Results_0_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
file_name <- paste("g_func_no_", g_func_no, "_case_id_0_d_1",  "_type_",type, "_seed_", seed,".RData", sep = "")
load(paste(dir,"/",file_name,sep = ""))

control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no = g_func_no, SNR = 5)
dat <- dat.generate(seed = seed, control = control, type = type,direction = "symmetric")

x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
y_train <- dat$y$y_train;  y_val <- dat$y$y_val;  y_test <- dat$y$y_test

t_range <- c(0,1)
grid <- dat$tt

u <- 1:length(y_train)
range_beta <- 4:7
range_eta <- 20
norder  <- 4
model_RobustFPLM <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                                 norder = norder, fLoss = "lmrob", trace = TRUE)
tmp_RobustFPLM <- FPLMBsplines.predict(model_RobustFPLM, newx = x_test, newy = y_test)

dat2save$RFPLM$pred_test <- tmp_RobustFPLM$pred
dat2save$RFPLM$err_test <- tmp_RobustFPLM$error
library(rlist)
save(dat2save, file = paste(dir,"/",file_name,sep = ""))

# #load(paste(dir,"/",file_name,sep = ""))
# 
# 
# matplot(t(x_train),col = "grey", lty = 1, type = "l")
# matplot(t(x_train[dat$outliers_train,]),col = "red",add = TRUE, lty = 1, type = "l")
# 
# hist(y_train)
# 
# hist(dat$g$g_train)
# 
# d <- 4; shrinkage.l2 <- 0.05; precision <- 6; niter <- 1000; nknot<- 3
# 
# shrinkage.rr <- 0.05; n_init <- 500
# niter <- 1000
# control.rr <- RTFBoost.control(make_prediction = TRUE, eff_m= 0.95, bb = 0.5, 
#                                tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
#                                type = "RR", shrinkage  = shrinkage.rr, precision = precision, 
#                                init_type = "median", n_init = n_init, niter = niter, 
#                                nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
# 
# model.rr <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
#                      x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
#                      control = control.rr)
# 
# matplot(t(x_train), col = "grey", type = "l", lty = 1, ylim = c(-10,10))
# matplot(t(x_train[dat$outliers_train,]), col = "red", add = TRUE, type = "l", lty = 1)
# 
# matplot(t(x_test), col = "grey", type = "l", lty = 1, ylim = c(-10,10))
# matplot(t(x_test[c(139,398,479,891),]), col = "red", add = TRUE, type = "l", lty = 1)
# 
# 
# control.rr.tree <- RTFBoost.control(make_prediction = TRUE, eff_m= 0.95, bb = 0.5,
#                                     tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
#                                     type = "RR", shrinkage  = shrinkage.rr, precision = precision,
#                                     init_type = "LADTree", n_init = n_init, niter = niter,
#                                     max_depth_init = 4, min_leaf_size_init =  20, nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
# 
# 
# model.rr.tree <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
#                          x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
#                          control = control.rr.tree)
# 
# which(model.rr.tree$f_test_t - y_test>20)
# boxplot(model.rr.tree$f_test_t - y_test, ylim = c(-5,5))
# boxplot(model.rr$f_test_t - y_test, ylim = c(-5,5))
# 
# 
# flagger_outlier <- which(abs(model.rr$f_val_t - y_val)>3*mad(model.rr$f_val_t - y_val))
# 
# mean(abs(model.rr$f_val_t[-flagger_outlier] - y_val[-flagger_outlier]))
# mean(abs(model.rr.tree$f_val_t[-flagger_outlier] - y_val[-flagger_outlier]))
# 
# 
# mean((model.rr$f_val_t[-dat$outliers_val] - y_val[-dat$outliers_val])^2)
# mean((model.rr.tree$f_val_t[-dat$outliers_val] - y_val[-dat$outliers_val])^2)
# 
# 
# par(mfrow = c(1,2))
# hist(model.rr$f_val_t[-flagger_outlier] - y_val[-flagger_outlier])
# hist(model.rr.tree$f_val_t[-flagger_outlier] - y_val[-flagger_outlier])
# 
# 
# boxplot(model.rr$f_val_t - y_val, ylim = c(-20,20))
# boxplot(model.rr.tree$f_val_t - y_val)
# 
# hist(model.rr.tree$f_val_t - y_val)
