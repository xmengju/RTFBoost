#### 
rm(list = ls())
source("Code/Exp_A/dat_generate.R")

control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no = 1, SNR = 5)
seed <- 1
lim_range <- list(c(-10,15), c(-10,15), c(-5,15), c(-20,20), c(-5,20), c(-5,20))
names(lim_range) <- paste("C", 0:5, sep = "")
pdf("visualize_outliers.pdf", width = 7, height = 9)

for(type in paste("C", 0:5, sep = "")){
  par(mfrow = c(3,2))
  for(g_func_no in 1:5){
    control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no = g_func_no, SNR = 5)
    dat <- dat.generate(seed, control = control, type = type)
    if(g_func_no == 1){
      matplot(t(dat$x$x_train), lty = 1, type = "l", 
            main = paste("training data", type), ylab = "x(t)", xlab = "t",  col='gray70', ylim=lim_range[[type]])
      matplot(t(dat$x$x_train[dat$outliers_train,]), lty = 1, type = "l", 
            col='red', add = TRUE, ylim=c(-10, 10))
    }
  
    mycol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50")
    c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
    c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
    if(type == "C0"){
    hist(dat$y$y_train,  xlab = NULL, xlim = c(-20,20), 15,col =c1,  main = paste("Histogram of training y: r", g_func_no, sep = ""))
    }else{
    hist(dat$y$y_train[-dat$outliers_train],  xlab = NULL, xlim = c(-20,20),15, col =c1,  main= paste("Histogram of training y: r", g_func_no, sep = ""))
    hist(dat$y$y_train[dat$outliers_train], xlab = NULL,  col =c2, add = TRUE, 15)
    }
  }
  
  
}
dev.off()
