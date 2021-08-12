### make_up indices for seeds 1-20 RTFBoost
make_up_idx_case0 <- function(seeds, nknot = 3, g_func_nos = c(1:5)){
  
  all_exps <- NULL; SNR <- 5;nknot <- 3
  
  for(g_func_no in g_func_nos){
    for(type in c("C0","C1","C2","C3","C4","C5")){
      dir <- paste("Results/Results_0_type_", type,"_g_", g_func_no, "_SNR_", SNR, "_nknot", nknot, "/", sep = "")
      for(d in 1){
        pattern <-glob2rx(paste("g_func_no_",g_func_no,"_case_id_0_d_",d,  "*", sep = ""))
        file_list <- list.files(dir, pattern= pattern)
        for(i in seeds){
          file_name <- paste("g_func_no_",g_func_no,"_case_id_0_d_",d,"_type_", type, "_seed_",i, ".RData", sep = "")
          if(!file_name%in%file_list){
            all_exps <- rbind(all_exps, c(g_func_no, d, type, i, SNR))
          }
        }
      }
    }
  }
  colnames(all_exps) <-  c("g_func_no", "d","type","seed","SNR" )
  return(all_exps)
}


## checking missing values 
make_up_idx_case1 <- function(seeds, nknot = 3, g_func_nos = c(1:5)){
  
  all_exps <- NULL; SNR <- 5;nknot <- 3

  for(g_func_no in g_func_nos){
    for(type in c("C0","C1","C2","C3","C4","C5")){
      dir <- paste("Results/Results_1_type_", type,"_g_", g_func_no, "_SNR_", SNR, "_nknot", nknot, "/", sep = "")
      for(d in 1:4){
        pattern <-glob2rx(paste("g_func_no_",g_func_no,"_case_id_1_d_",d,  "*", sep = ""))
        file_list <- list.files(dir, pattern= pattern)
        for(i in seeds){
          file_name <- paste("g_func_no_",g_func_no,"_case_id_1_d_",d,"_type_", type, "_seed_",i, ".RData", sep = "")
          if(!file_name%in%file_list){
            all_exps <- rbind(all_exps, c(g_func_no, d, type, i, SNR))
          }
        }
      }
    }
  }
  colnames(all_exps) <-  c("g_func_no", "d","type","seed","SNR" )
  return(all_exps)
}


library(rlist)


# all exist 
seeds <- 1:20; nknot = 3
all_exps <- data.frame(make_up_idx_case0(seeds, nknot = 3, g_func_nos = c(1:5)))
all_exps$d <- as.numeric(all_exps$d)
all_exps$g_func_no <- as.numeric(all_exps$g_func_no)
all_exps$SNR <- as.numeric(all_exps$SNR)
all_exps$seed <- as.numeric(all_exps$seed)
file_name <- paste("Code/Exp_C/make_up/case_0_nknot_",nknot,".rds", sep = "")
list.save(all_exps,file = file_name)



# all exist 
seeds <- 1:20; nknot = 3
all_exps <- data.frame(make_up_idx_case1(seeds, nknot = 3, g_func_nos = c(1:5)))

all_exps$d <- as.numeric(all_exps$d)
all_exps$g_func_no <- as.numeric(all_exps$g_func_no)
all_exps$SNR <- as.numeric(all_exps$SNR)
all_exps$seed <- as.numeric(all_exps$seed)
file_name <- paste("Code/Exp_C/make_up/case_1_nknot_",nknot,".rds", sep = "")
list.save(all_exps,file = file_name)




