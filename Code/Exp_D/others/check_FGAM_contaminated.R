## check the performance of FGAM

i <- 10
control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no = 1, SNR = 5)
dat <- dat.generate(seed = i, control = control, type = "C3", direction = "symmetric")

x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
y_train <- dat$y$y_train;  y_val <- dat$y$y_val;  y_test <- dat$y$y_test


matplot(t(x_train), col = "grey", type = "l", lty = 1)
matplot(t(x_train[dat$outliers_train,]), add = TRUE, col = "red", type = "l", lty = 1)


t_range <- c(0,1)
grid <- dat$tt

xx <- x_train
yy <- y_train
nbasises_FGAM <- 15
model_FGAM_1 <- fgam(yy ~ af(xx, splinepars=list(k=c(nbasises_FGAM,nbasises_FGAM),m=list(c(2,2),c(2,2)))), gamma = 1.2, method="REML")
pred_train <- predict(model_FGAM_1 ,newdata=list(xx =x_train),type='response')
pred_test <-  predict(model_FGAM_1 ,newdata=list(xx =x_test),type='response')
err_test <-  mean((pred_test- y_test)^2)
err_train <-   mean((pred_train-y_train)^2)
print(c(err_train, err_test))

i <- 50
dat <- dat.generate(seed = i, control = control, type = "C3", direction = "symmetric")

x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
y_train <- dat$y$y_train;  y_val <- dat$y$y_val;  y_test <- dat$y$y_test

t_range <- c(0,1)
grid <- dat$tt

xx <- x_train
yy <- y_train
nbasises_FGAM <- 15
model_FGAM_2 <- fgam(yy ~ af(xx, splinepars=list(k=c(nbasises_FGAM,nbasises_FGAM),m=list(c(2,2),c(2,2)))), gamma = 1.2, method="REML")
pred_train <- predict(model_FGAM_2 ,newdata=list(xx =x_train),type='response')
pred_test <-  predict(model_FGAM_2 ,newdata=list(xx =x_test),type='response')
err_test <-  mean((pred_test- y_test)^2)
err_train <-   mean((pred_train-y_train)^2)
print(c(err_train, err_test))
make_plot_fgam(model_FGAM_1, nlevels = 40)
make_plot_fgam(model_FGAM_2, nlevels = 40)




yyy <- seq(0,1,length.out = 200)
xxx <- seq(-10,10, length.out = 100)

zzz <- matrix(NA, length(xxx), length(yyy))
for(i in 1:length(yyy)){
  for(j in 1:length(xxx)){
    zzz[j,i] <- (sin(3*pi*yyy[i]/2) + sin(pi*yyy[i]/2))*xxx[j]
  }
}
contour(xxx,yyy, zzz, nlevels = 20)
lines( x_train[1,],seq(0,1,length.out = 100),col = "red")
lines( x_train[2,],seq(0,1,length.out = 100),col = "red")
lines( x_train[12,],seq(0,1,length.out = 100),col = "red")



par(mfrow = c(1,2))
make_plot_fgam(model_FGAM_1)
make_plot_fgam(model_FGAM_2, nlevels = 0)
lines( x_train[10,],seq(0,1,length.out = 100),col = "red")


matplot(t(x_train[c(1,3),]), type = "l", lty = 1)

make_plot_fgam <- function (object,  xval = NULL, tval = NULL, deriv2 = FALSE, 
          theta = 50, ticktype = "detailed",nlevels = 30,...) 
{
    tnames <- names(object$model)
    af.term <- tnames[grep(paste("", ".omat", sep = ""), 
                           tnames)[1]]
    af.term <- strsplit(af.term, ".omat")
    tnames <- object$fgam$labelmap
    afind <- grep(paste("te[(]", af.term, sep = ""), tnames)
    af.ind <- grep(paste("af[(]", af.term, sep = ""), names(object$fgam$ft))
   
    
    tname <- tnames[[afind]]
    basistype <- strsplit(tname, "[(]")[[1]][1]
    sstring <- paste(basistype, "[(]", af.term, "\\.omat,", af.term, 
                   "\\.tmat", "[)]:L\\.", af.term, sep = "")
    tind <- grepl(sstring, names(object$coef))

    temp <- list()
    temp[[paste("L.", af.term, sep = "")]] <- 1
   
    tvar <- tolower(af.term)
    
    if(object$fgam$ft[[af.ind]]$Qtransform){
      mtitle <- bquote(paste(hat(F),'(p,t),   p=',hat(G)[t],'(',.(tvar),')',sep=''))
    }else{
      mtitle <- bquote(paste(hat(F),'(',.(tvar),',t)',sep=''))
    }
    tvar <- paste('\n',tvar,sep='')
    vis.gam(object,view=c(paste(af.term,'.omat',sep=''),paste(af.term,'.tmat',sep='')),cond=temp,
            xlab=tvar,ylab='\nt', plot.type="contour", main=mtitle, nlevels = nlevels,...)
}

