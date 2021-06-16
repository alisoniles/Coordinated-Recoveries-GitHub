#
# 0_CORE_functions.R: define functions
#
# ver 1.0: initially written on 04/15/2021 by A.Iles


# load dependent libraries
library(rEDM)
library(dplyr)
library(tseriesChaos)
library(multispatialCCM)
library(tidyr)
library(stringr)
library(zoo)

library(ggplot2)
library(ggforce)
library(gridExtra)
library(foreach)
library(parallel)
library(doParallel)
library(pracma) #detrending function

#Normalize each column in a data frame
normalize <- function(block)
{
  if(NCOL(block) > 1)
  {
    n <- NROW(block)
    means <- sapply(block, mean, na.rm = TRUE)
    sds <- sapply(block, sd, na.rm = TRUE)
    return((block - matrix(rep(means, each = n), nrow = n)) / 
             matrix(rep(sds, each = n), nrow = n))
  }
  else
    return((block - mean(block, na.rm = TRUE)) / sd(block, na.rm = TRUE))
}



#Function to remove to extra NAs, but leave one NA between each stock's data in an MPG
shape_stock_data <- function(D, min_stock_size = 10) {
  Dcc <- D[complete.cases(D[,3]),] 
  
  #get library of beginning and end points of each chunk of data
  lib <- create_lib(Dcc, min_stock_size)
  
  #insert NA's between each chunk of data
  x <- c(NA, NA, NA)
  for (r in 1:nrow(lib)){
    tmp <- Dcc[lib[r,1]:lib[r,2],]
    x <- rbind(x,tmp,c(NA, NA, NA))}
  
  data <- (x[2:(nrow(x)-1),])
  data <- as.data.frame(data)
  rownames(data) <- NULL 

  return(data)
}


#Create library list for EDM that identifies the beginning and end points of each chunk of data combined into the same embedding. Year must be in the first column of the time series. Complete cases only. Removes sections that are not continuous for at least 'cont_tp' time points.
create_lib <- function(D, min_cont_tp)
{
  lib <- matrix(NA, nrow = length(which(diff(D$year)!=1))+1, ncol = 2) 
  lib[,1] <- c(1, which(diff(D$year)!=1)+1)
  lib[,2] <- c(which(diff(D$year)!=1), nrow(D))
  
  minlib <- lib[,2]-lib[,1] #only include in the library the sections of data that are continuous for at least 'cont_tp' time points. 
  lib <- lib[minlib>=min_cont_tp,] 
  
  return(lib)
}



#' Randomization test for nonlinearity using S-maps and surrogate data
#' 
#' tests for nonlinearity using S-maps by comparing improvements in forecast skill (delta rho and delta mae) between 
#' linear and nonlinear models with a null distribution from surrogate data.
#' 
#' ts: the original time series
#' method: which algorithm to use to generate surrogate data
#' num_surr: the number of null surrogates to generate
#' E: the embedding dimension for s_map
#'  ... optional arguments to s_map

#ts <- Dt; method <-  "ebisuzaki"; num_surr <-  200; E <-  6; tau <-  1
test_nonlin <- function(ts, method = "ebisuzaki", num_surr = 200, 
                        theta = c(0, 0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9),
                        tau = tau, tp = 1, E = E, ...)
{
  
  compute_test_stat <- function(ts, ...)
  {
    results <- s_map(ts, tau = tau, tp = 1, E = E, ...)
    delta_rho <- max(as.numeric(results$rho))-results$rho[[1]]
    opt_theta <- results$theta[results$rho==max(as.numeric(results$rho))]
    return(list( 
      theta = as.numeric(results$theta), 
      rho = as.numeric(results$rho),
      opt_theta = opt_theta,
      delta_rho = delta_rho))
  }
  
  actual_stats <- compute_test_stat(ts, ...)
  delta_rho <- actual_stats$delta_rho; names(delta_rho) <- NULL
  theta <- actual_stats$theta
  rho <- actual_stats$rho
  opt_theta <- actual_stats$opt_theta
  
  surrogate_data <- make_surrogate_data(ts[,2], method = method, num_surr = num_surr) #generate surrogate data 
  null_stats <- data.frame(t(apply(surrogate_data, 2, compute_test_stat, ...)))
  null_stats_delta_rho <- lapply(null_stats, function(x) tail(as.numeric(unlist(x)), n=1))
  
  return(list(
    E = E,
    theta = theta,
    rho = rho,
    opt_theta = opt_theta,
    num_surr = num_surr,
    delta_rho = delta_rho, 
    delta_rho_p_value = (sum(null_stats_delta_rho>delta_rho) + 1) / (num_surr + 1)))
}   



#' Randomization test for nonlinearity using 'PredictNonlinear' wrapper function and surrogate data
#' 
#' tests for nonlinearity using S-maps by comparing improvements in forecast skill (delta rho and delta mae) between 
#' linear and nonlinear models with a null distribution from surrogate data.
#' 
#' D: the original time series with time in the first column
#' method: which algorithm to use to generate surrogate data, "ebisuzaki", "random_shuffle"
#' num_surr: the number of null surrogates to generate
#' E: the embedding dimension for s_map
#'  ... optional arguments to s_map

test_PredictNonlinear <- function(D, lib = lib, pred = lib, columns = columns, target = target, 
                        method = "ebisuzaki", num_surr = 200, 
                        theta = c(0, 0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9),
                        tau = 1, tp = 1, E = E, ...)
{
  actual_results <- PredictNonlinear(dataFrame = D, lib = lib, pred = lib, columns = colnames(D)[2], target = colnames(D)[2], E = E, showPlot = FALSE)
  delta_rho <- max(as.numeric(actual_results$rho))-actual_results$rho[[1]]
  opt_theta <- actual_results$Theta[actual_results$rho==max(as.numeric(actual_results$rho))]
  theta = as.numeric(actual_results$Theta)
  rho = as.numeric(actual_results$rho)
  
  surrogate_data <- make_surrogate_data(D[,2], method = method, num_surr = num_surr) #generate surrogate data 
  null_stats_delta_rho <- matrix(data=NA, nrow=num_surr, ncol=1)
  for(s in 1:num_surr){
    SD <- data.frame(cbind(D[,1], surrogate_data[,s]))
    colnames(SD) <- c("year", "surr")
    surr_results <- PredictNonlinear(dataFrame = SD, lib = lib, pred = lib, columns = "surr", target = "surr", E = E, showPlot = FALSE)
    null_stats_delta_rho[s] <- max(as.numeric(surr_results$rho))-surr_results$rho[[1]]
    }

return(list(
    E = E,
    theta = theta,
    rho = rho,
    opt_theta = opt_theta,
    num_surr = num_surr,
    delta_rho = delta_rho, 
    delta_rho_p_value = (sum(null_stats_delta_rho>delta_rho) + 1) / (num_surr + 1)))
}   



# perform normal multivariate S-map analysis to choose the theta with the optimal forecast skill
# block: SSR manifold with first column time and the target of the smap analysis in the second column
find_opt_theta  <- function(block, columns_names=columns_names, lib=t(c(1,nrow(block))), thetarange=seq(0,10,0.1), criterion="rmse") {
  # search optimal theta
  foreach(theta=thetarange, .combine=rbind, .packages=env_lib_ts) %do% {
    block_lnlp(block, 
               lib=lib,
               method="s-map",
               theta=theta, 
               first_column_time = TRUE,
               columns = columns_names) 
  } %>% slice(if_else(criterion=="rho",which.max(.[,criterion]),which.min(.[,criterion]))) #slice results based on max rho or min rmse
}




#Function to remove to extra NAs, but leave one NA between each stock's data in an MPG
#'data' with columns: year, forcing process and response process
CCM_shape_block_data <- function(data, min_cont_ts = 10) {

  data <- rbind(rep(NA, ncol(data)),data) #add an initial row of NA so the difference function works for the first row
  
  dataNA <- !complete.cases(data) # locate NA
  data[dataNA,] <- NA
  rownames(data) <- NULL
  
  #record the beginning and endpoints of each data chunk
  CC <- complete.cases(data[,1])  
  lib <- matrix(NA, nrow = length(which(diff(CC)==1)), ncol = 2)
  lib[,1] <- which(diff(CC)==1)+1
  lib[,2] <- which(diff(CC)==-1)
  colnames(lib) <- c("lower", "upper")
  
  #only include in the library the sections of data that are continuous for at least min_cont_ts time points. 
  minlib <- lib[,2]-lib[,1]+1
  #if there are no library chunks long enough, then return error
  if(sum(minlib>=min_cont_ts)<=1){
    print("Error - data chunks are too small for the given E and tau")
    return()
  }
  lib <- as.matrix(lib[(minlib>=min_cont_ts),])
  
  x <- rep(NA, ncol(data))
  for (r in 1:nrow(lib)){
    xtmp <- data[lib[r,1]:lib[r,2],]
    x <- rbind(x,rep(NA, ncol(data)),xtmp)}
  data <- (x[3:nrow(x),])
  data <- as.data.frame(data)
  rownames(data) <- NULL
  
  return(data)
}

#modified ccmtest function to test one-way interactions (we're only interested in the effects on salmon, not what salmon affect and most of the variables are exogenous anyway)
ccmtest_oneway <- function(CCM_boot_AcB) {
  #Tests for significant causal signal based on 95%
  #confidence intervals from bootstrapping.
  #Compares shortest library to longest
  pval_a_cause_b<-1-sum(CCM_boot_AcB$FULLinfo[1,]<
                          CCM_boot_AcB$FULLinfo[nrow(CCM_boot_AcB$FULLinfo),], na.rm=T)/
    ncol(CCM_boot_AcB$FULLinfo)
  return(pval_a_cause_b)
}





#Run the multispatial convergent cross mapping algorithm on A causing B and B causing A. 
# This function runs the multispatial CCM_boot and tests output from CCM_boot for significant causal signal using ccmtest to compare the 95% confidence intervals for estimated rho for the shortest and longest libraries calculated, and uses this to determine whether predictive power has significantly increased.  Reorganizes the output into a data frame for plotting. 
# Desired library lengths for which to compute CCM. Defaults to the maximum possible length ((tau * (E - 1) + (E + 1)):length(A) - E + 2) (though number of resulting predictions may be smaller because of gaps in the time series). Shortening this list (e.g., only predicting every nth element) will reduce run-time for the algorithm, but may also reduce ability to detect causal relations.
#A <- AB[,2]; B <- AB[,3]; An <- An; Bn <- Bn; EB <- EB; tau <- tau
CCM_boot_df_and_sig_test <- function(A, B, An, Bn, EB, tau)
{
  DL <- c((tau * (EB - 1) + (EB + 1)), (length(A) - EB + 2))
  
  CCM_BcA <- CCM_boot(B, A, EB, tau=tau, DesiredL=DL, iterations=500) # Does B "cause" A?
  CCM_A_B_sig_test<-ccmtest_oneway(CCM_BcA)
  
  #make data frame of results for output
  BcA_df <- as.data.frame(cbind(CCM_BcA$Lobs, CCM_BcA$rho, CCM_BcA$sdevrho, (CCM_BcA$rho-CCM_BcA$sdevrho), (CCM_BcA$rho+CCM_BcA$sdevrho)))
  colnames(BcA_df) <- c("Lobs", "rho", "sdevrho", "lower", "upper")
  
  BcA_df$model <- paste(Bn, " causes ", An, ", p = ", signif(CCM_A_B_sig_test,2) , sep="")
  BcA_df$response_var <- An
  BcA_df$causal_var <- Bn
  BcA_df$pval <- signif(CCM_A_B_sig_test,2)
  BcA_df$E <- EB
  BcA_df$tau <- tau
  
  #estimate slope of rho values as a function of library size for the largest 1/4 of library sizes
  D <- BcA_df[BcA_df$Lobs > (max(BcA_df$Lobs, na.rm=TRUE)*0.75),]
  BcA_df$slope <- lm(D$rho ~D$Lobs)$coefficients[2]
  
  return(BcA_df)
}




