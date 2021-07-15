#### 
rm(list = ls())
source("Code/Exp_A/dat_generate.R")

control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no = 1, SNR = 5)
seed <- 1
lim_range <- list(c(-10,15), c(-5,15), c(-20,20), c(-5,20))
names(lim_range) <- paste("C", 1:4, sep = "")
pdf("visualize_outliers.pdf")
par(mfrow = c(2,2))
for(type in paste("C", 1:4, sep = "")){
  dat <- dat.generate(seed, control = control, type = type)
  matplot(t(dat$x$x_train), lty = 1, type = "l", 
          main = paste("training data", type), ylab = "x(t)", xlab = "t",  col='gray70', ylim=lim_range[[type]])
  matplot(t(dat$x$x_train[dat$outliers_train,]), lty = 1, type = "l", 
          col='red', add = TRUE, ylim=c(-10, 10))
}
dev.off()