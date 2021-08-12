
file_name <- "/Users/xmengju/Research/RTFBoost/Writing/Report_3/appendix.tex"

file.create(file_name)

for(g_func_no in 1:5){
  
  dat2plot_mean <- NULL
  dat2plot_sd <- NULL
  res_table <- NULL
  for(type in c("C0", "C1","C2","C3","C4", "C5")){
    res1 <- summarize_1_combine(g_func_no = g_func_no, type = type, summarize_type = "others")$res
    res0 <- summarize0(g_func_no = g_func_no, type = type)
    
    dat2plot_mean <- rbind(c(apply(res1, 2,mean), apply(res0, 2,mean)))
    dat2plot_sd <- rbind(c(apply(res1, 2,sd), apply(res0, 2,sd)))
    
    tmpp <-  round(dat2plot_mean,3)
    idx <- which(tmpp ==  sort(unique(tmpp))[1])
    idx <-  c(idx, which(tmpp == sort(unique(tmpp))[2]))
    
    tmppp <- paste(trimws(format(round(dat2plot_mean,3),nsmall= 3)), " (",trimws(format(round( dat2plot_sd ,3),nsmall= 3)),")", sep = "")
    tmppp[idx] <- paste("\\textbf{", trimws(format(round(dat2plot_mean[idx],3),nsmall= 3)), "} (",trimws(format(round(dat2plot_sd[idx],3),nsmall= 3)),")", sep = "")
    
    res_table <- rbind(res_table, tmppp)
  }
  
  colnames(res_table) <-  c("TFBoost(L2)","TFBoost(LAD)","TFBoost(LAD-M)","TFBoost(RR.5)","TFBoost(RR.2)","FPPR","FGAM","RFPLM", "MFLM","RFSIR")
  rownames(res_table)<- c("$C_0$","$C_1$","$C_2$","$C_3$","$C_4$", "$C_5$")
  caption_tmp <- paste("Summary statistics of test errors for data generated from $r_", g_func_no, "$; displayed in the form of mean (sd)", ".", sep = "" )
  
  table_tmp <- xtable(t(res_table), caption = caption_tmp,  table.placement="H", floating = getOption("xtable.floating", TRUE),
                      caption.placement = getOption("xtable.caption.placement", "bottom"))
  print(table_tmp, sanitize.text.function=function(x){x},  include.colnames=TRUE,
        append = TRUE, file = file_name, table.placement="H")
}


library(ggplot2)



make_figure <- function(g_func_no){
  
  y_lim <- list(c(0.1,0.4),  c(0.15,0.6),  c(0.25,0.8), 
                c(0.25,1.1),  c(0.45,1.5))
  tmp <- NULL
  for(type in c("C0","C1","C2","C3","C4","C5")){

    res2 <- summarize_1_combine(g_func_no = g_func_no, type = type, summarize_type = "all")
    res2 <- data.frame(res2$res, rep(type, nrow(res2$res)))
    
    colnames(res2) <- c("TFBoost(L2)", "TFBoost(LAD)", "TFBoost(LAD-M)","TFBoost(RR.5)", "TFBoost(RR.2)", "type")
    tmp  <- rbind(tmp, res2)
  }
  dat2plot <- melt(tmp, id.vars = "type", variable.name = "method")
  
  tmp <- NULL
  for(type in c("C0","C1","C2","C3","C4","C5")){
    res0 <- summarize0(g_func_no = g_func_no, type = type)
    res0 <- data.frame(res0, rep(type, nrow(res0)))
    colnames(res0) <- c( "FPPR","FGAM","RFPLM","MFLM", "RFSIR", "type")
    tmp  <- rbind(tmp, res0)
  }
  
  dat2plot <- rbind(dat2plot, melt(tmp, id.vars = "type", variable.name = "method"))
  dat2plot$value <- as.numeric(dat2plot$value)
  p <- ggplot(dat2plot, aes(fill=method, y=value)) + 
    geom_boxplot() + ylab("MSPE") + xlab(" ")  + ggtitle(paste("r", g_func_no, sep = ""))+
    facet_wrap(vars(type),  ncol = 3) + theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5), legend.position="bottom") +  scale_fill_brewer(palette="Spectral")
  
  return(p + ylim(y_lim[[g_func_no]][1],y_lim[[g_func_no]][2]))
  
}

library(reshape2)
library(ggplot2)
pp <- list()
for(g_func_no in 1:5){
  pp[[g_func_no]] <- make_figure(g_func_no = g_func_no)
}

pdf("make_figures_tmp.pdf", width = 15, height = 10)
for(i in 1:5){
  print(pp[[i]])
}
dev.off()
