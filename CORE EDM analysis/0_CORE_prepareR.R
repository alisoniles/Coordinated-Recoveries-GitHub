#
# 0_CORE_prepareR.R: prepare CORE data set and load functions for analysis
#
# ver 1.0: initially written on 4/20/2021 by A.Iles

# load functions
source("0_CORE_functions.R")

# load data
load("Data/Rdata/norm_block_data.Rdata")  
data <- norm_block_data

# load the optimal E, tau and theta for the response processes
ETT <- read.csv("Data/csv_data/CORE_opt_E_tau_theta.csv")

#Create list of putative causal variables (forcing processes)
     aa <- colnames(data)
     bb <- data.frame(c(aa[11:length(aa)] ))
     colnames(bb) <- c("variable")
     cc <- str_split_fixed(bb$variable, "[.]", n=3)
FP_list <- cbind(bb,cc)
colnames(FP_list) <- c("FP", "cat", "subcat", "offset") 
FP_list$offset <- as.numeric(FP_list$offset)
rm(aa, bb, cc)

#