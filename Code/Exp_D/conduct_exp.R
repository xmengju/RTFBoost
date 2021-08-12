rm(list = ls())
library(rpart)
library(splines)
library(MyFPLM)
library(rpart)
library(refund)
library(fda)
library(robustbase)
library(np)
library(RobStatTM)

source("Code/RTFBoost.R")
source("Code/TREE.init.R")
source("Code/TREE.R")
source("Code/user_rpart.R")
source("Code/utils.R")
source("Code/Exp_D/fun_exp.R")
source("Code/Exp_D/dat_generate.R")
source("Code/FSIM.R")
source("Code/RFSIR.R")

args <- commandArgs(trailingOnly = TRUE) 
nknot <-3; g_func_nos <- 1:5
ds <- c(1,2,3,4); SNRs <- c(5)
shrinkage <- 0.05 
types <- c("C0", "C1","C2","C3","C4","C5")
make_up  <- FALSE

validation_tree <- TRUE
competitor.control <- list(nbasises_FGAM = 15, niter_FPPR =15, nknots_FPPR = 3, k = 10, qs = 1:10)

if(length(args)== 0){
  g_func_nos <- c(1,2,3,4,5)
  SNRs <- 5
  ds <- c(1,2,3,4)
  test <- TRUE
  RNGkind(sample.kind = "default")
  exp_id <- 1  # row index in the conduct_sheet 
  case_id <- 1
  if(case_id == 1){
    n_init <- 5
    niter <- 10  
    ncol_num <- 6
  }
  if(case_id == 0){
    n_init <- 5 # doesn't matter
    niter <- 10  # doesn't matter
    ncol_num <- 10
  }
}else{
  test <- FALSE
  exp_id <- as.numeric(args[1]) 
  case_id <-  as.numeric(args[2]) 
  if(case_id == 1){
    n_init <- 1000
    niter <- 2000 
    ncol_num <- 6
  }
  if(case_id == 0){
    n_init <- 1000 # doesn't matter
    niter <- 2000  # doesn't matter
    ncol_num <- 10
  }
}

if(case_id == 0){
  methods <- c("FGAM","FPPR","RFSIR","MFLM","RFPLM")
}else{
  methods <- c("RTFBoost")
}


seeds <- 1:20

if(case_id == 0){
  ds <- 1
}

if(make_up == FALSE){
  all_exps <- expand.grid(g_func_nos, ds, types, seeds, SNRs)
  colnames(all_exps) <- c("g_func_no", "d","type","seed","SNR" )
}else{
  file_name <- paste("Code/Exp_C/make_up/case_", case_id,"_nknot_3.rds", sep = "")
  all_exps <- data.frame(readRDS(file_name))
}


if(nrow(all_exps)%%ncol_num!=0){
  if(nrow(all_exps)%%ncol_num!=0){
    conduct_sheet <- matrix(c(1:nrow(all_exps), rep(NA, ncol_num - nrow(all_exps)%%ncol_num)), ncol = ncol_num, byrow = TRUE)
  }else{
    conduct_sheet <- matrix(c(1:nrow(all_exps)), ncol = ncol_num, byrow = TRUE)
  }
}else{
  conduct_sheet <- matrix(1:nrow(all_exps), ncol = ncol_num, byrow = TRUE)
}


print(dim(conduct_sheet))

conduct.exp <- function(exp_id = 1, conduct_sheet){
  
  for(i in 1:ncol(conduct_sheet)) {
    
    print(c(i,"th trial"))
    if(!is.na(conduct_sheet[exp_id,i])){
      
      g_func_no <- all_exps[conduct_sheet[exp_id, i],]$g_func_no
      seed <- all_exps[conduct_sheet[exp_id, i],]$seed
      d <-  all_exps[conduct_sheet[exp_id, i],]$d
      type <-  all_exps[conduct_sheet[exp_id, i],]$type
      SNR <- all_exps[conduct_sheet[exp_id, i],]$SNR

      precision <- 6
   
      if(test){
        dir_name <- paste("Results_tmp/Results_", case_id, "_type_", type, "_g_",g_func_no,  "_SNR_", SNR, "_nknot", nknot,  sep = "")
      }else{
        dir_name <- paste("Results/Results_", case_id, "_type_", type, "_g_",g_func_no, "_SNR_", SNR, "_nknot", nknot,  sep = "")
      }
      
      if(dir.exists(dir_name) == FALSE){
        dir.create(dir_name)
      }
      
      start_time <- Sys.time()
      if(g_func_no == 5){
        control.tree.list <- set.control(d, precision, shrinkage, nknot, n_init = 1000, niter = 2000)
      }else{
      control.tree.list <- set.control(d, precision, shrinkage, nknot, n_init, niter)
    }
      dat2save <- tryCatch(do.exp(seed = seed, g_func_no = g_func_no, SNR = SNR, d = d, 
                                  control.tree.list = control.tree.list, methods = methods, type = type, 
                                  validation_tree = validation_tree, competitor.control =  competitor.control))
      
      end_time <- Sys.time()
      dat2save$time <- end_time - start_time
      print(dat2save$time)
      save(file = paste(dir_name, "/g_func_no_", g_func_no, "_case_id_",case_id, "_d_", d,"_type_", type,   "_seed_", seed, ".RData", sep = ""), dat2save)
      
    }
  }
}


conduct.exp(exp_id = exp_id, conduct_sheet)
