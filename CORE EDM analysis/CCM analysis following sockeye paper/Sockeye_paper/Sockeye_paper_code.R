#### README
# The following code is intended for use with the R programming language. No conversion is necessary.
#
# The remaining dataset files will require conversion to .csv format. This can be done using Excel's "save as .csv" command. The expected datafiles are:
# "Dataset S1.xls" ---> "sockeye_ret_data.csv"
# "Dataset S2.xls" ---> "sockeye_data.csv"
# "Dataset S3.xls" ---> "env_data.csv"

#### LICENSE
# This software is Copyright © 2014 The Regents of the University of California. All Rights Reserved.
# 
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
# 
# Permission to make commercial use of this software may be obtained by contacting:
# 
# Technology Transfer Office
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# invent@ucsd.edu
# 
# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
# 
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

#---- packages used ----
library(rEDM)
library(rjags)
library(reshape2)
library(rgl)
library(ggplot2)
library(gridExtra)
library(xtable)

#---- function definitions ----
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

normalize_by_cycle_line <- function(ts)
{
     n <- length(ts)
     means <- rep.int(NA, times = 4)
     sds <- rep.int(NA, times = 4)
     mu <- rep.int(NA, times = n)
     sigma <- rep.int(NA, times = n)
     for(k in 1:4)
     {
          index <- seq(from = k, to = n, by = 4)
          means[k] <- mean(ts[index], na.rm = TRUE)
          sds[k] <- sd(ts[index], na.rm = TRUE)
          mu[index] <- means[k]
          sigma[index] <- sds[k]
     }
     ts <- (ts - mu) / sigma
     df <- data.frame(cbind(ts, mu, sigma))
     return(df)
}

compute_stats <- function(obs, pred)
{
     # computes performance metrics for how well predictions match observations
     # obs = vector of observations
     # pred = vector of prediction
     
     N = sum(is.finite(obs) & is.finite(pred))
     rho = cor(obs, pred, use = "pairwise.complete.obs")
     mae = mean(abs(obs-pred), na.rm = TRUE)
     return(data.frame(N = N, rho = rho, mae = mae))
}

preprocess_data <- function()
{
     preprocess_stock <- function(stock_df)
     {
          n <- NROW(stock_df)
          stock_df$rec45 <- stock_df$rec4 + stock_df$rec5
          stock_df$ret <- stock_df$rec4 + c(NA, stock_df$rec5[1:(n-1)]) # age-4 and age-5 fish (aligned to rec4)
          
          temp <- normalize_by_cycle_line(stock_df$rec45)
          stock_df$rec45_n <- temp$ts
          stock_df$rec45_mu <- temp$mu
          stock_df$rec45_sigma <- temp$sigma
          
          temp <- normalize_by_cycle_line(stock_df$rec4)
          stock_df$rec4_n <- temp$ts
          stock_df$rec4_mu <- temp$mu
          stock_df$rec4_sigma <- temp$sigma
          
          temp <- normalize_by_cycle_line(stock_df$rec5)
          stock_df$rec5_n <- temp$ts
          stock_df$rec5_mu <- temp$mu
          stock_df$rec5_sigma <- temp$sigma
          
          temp <- normalize_by_cycle_line(stock_df$eff)
          stock_df$eff_n <- temp$ts
          stock_df$eff_mu <- temp$mu
          stock_df$eff_sigma <- temp$sigma
          
          return(stock_df)
     }
     
     make_block <- function(stock_df, env_data)
     {
          discharge_names <- c("D_max", "D_apr", "D_may", "D_jun")
          temp_names <- c("ET_apr", "ET_may", "ET_jun", "PT_apr", "PT_may", "PT_jun", "PT_jul")
          pdo_names <- "PDO_win"
          discharge <- normalize(env_data[, discharge_names])
          temperature <- normalize(env_data[, temp_names])
          pdo <- normalize(env_data[, pdo_names])
          
          # line up environmental data
          # lag temperature and river discharge 2 years
          desired_years <- stock_df$yr + 2
          index_in_env_data <- match(desired_years, env_data$year)
          index_in_stock_df <- 1:length(desired_years)
          
          discharge_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = NCOL(discharge)))
          discharge_cols[index_in_stock_df,] <- discharge[index_in_env_data, ]
          stock_df[, discharge_names] <- discharge_cols
          
          temp_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = NCOL(temperature)))
          temp_cols[index_in_stock_df,] <- temperature[index_in_env_data, ]
          stock_df[, temp_names] <- temp_cols
          
          # lag PDO by 1 year (winter before smolt outmigration)
          desired_years <- stock_df$yr + 1
          index_in_env_data <- match(desired_years, env_data$year)
          pdo_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = 1))
          pdo_cols[index_in_stock_df,] <- pdo[index_in_env_data]
          stock_df[, pdo_names] <- pdo_cols
          
          return(stock_df)
     }
     
     data <- read.csv("sockeye_data.csv")
     
     # filter stocks we don't want
     stock_data <- split(data, data$stk)
     stock_data <- lapply(stock_data, preprocess_stock)
     
     # add env data
     env_data <- read.csv("env_data.csv")
     block_data <- lapply(stock_data, function(stock_df) { make_block(stock_df, env_data)})
     
     # save and return
     save(block_data, file = "block_data.Rdata")
     return()
}

compute_nonlinearity_aggregated <- function()
{
     load("block_data.Rdata")
     ret <- lapply(block_data, function(x) {
          temp <- x$ret
          temp <- temp[is.finite(temp)]
          return((temp - mean(temp)) / sd(temp))
     })
     x <- c()
     lib <- matrix(NA, nrow = 9, ncol = 2)
           last <- 0
           for(i in 1:9)
           {
                x <- c(x, ret[[i]])
                lib[i,] <- c(last+1, last + length(ret[[i]]))
                last <- lib[i,2]
           }
     simplex_output <- simplex(x, lib = lib, pred = lib, E = 1:6, exclusion_radius = 0, silent = TRUE)
     E <- simplex_output$E[which.max(simplex_output$rho)]
     smap_output <- s_map(x, lib = lib, pred = lib, E = E, exclusion_radius = 0, silent = TRUE)
     theta <- smap_output$theta[which.max(smap_output$rho)]
     
     save(simplex_output, E, smap_output, theta, file = "results_nonlinear_aggregated.Rdata")
     return()
}

test_nonlinearity_aggregated <- function(num_shuffles = 500)
{
     get_smap_stats <- function(x, lib, E = NULL)
     {
          if(is.null(E))
          {
               # compute E using simplex on recruits time series
               simplex_output <- simplex(x, E = 1:8, silent = TRUE)
               best_rho_E <- simplex_output$E[which.max(simplex_output$rho)]
               best_mae_E <- simplex_output$E[which.min(simplex_output$mae)]
               E <- min(best_rho_E, best_mae_E)
          }
          
          # compute theta using s-map and E 
          smap_output <- s_map(x, lib = lib, pred = lib, E = E, silent = TRUE)
          
          best_rho <- max(smap_output$rho)
          best_mae <- min(smap_output$mae)
          return(data.frame(delta_mae = best_mae - smap_output$mae[smap_output$theta == 0]))
     }
     
     load("block_data.Rdata")
     ret <- lapply(block_data, function(x) {
          temp <- x$ret
          temp <- temp[is.finite(temp)]
          return((temp - mean(temp)) / sd(temp))
     })
     x <- c()
     lib <- matrix(NA, nrow = 9, ncol = 2)
     last <- 0
     for(i in 1:9)
     {
          x <- c(x, ret[[i]])
          lib[i,] <- c(last+1, last + length(ret[[i]]))
          last <- lib[i,2]
     }
     E <- 4
     
     cat("calculating for actual data... ", sep = "")
     start_time <- proc.time()
     actual <- get_smap_stats(x, lib, E)
     delta_mae <- actual$delta_mae
     elapsed_time <- proc.time() - start_time
     cat("(", elapsed_time[3], " sec.)\n", sep = "")
     
     # null distribution
     cat("calculating for random shuffles... ", sep = "")
     start_time <- proc.time()
     null_dist <- do.call(rbind, lapply(1:num_shuffles, function(i) {
          x_shuffle <- c()
          for(i in 1:9)
          {
               n <- length(ret[[i]])
               x_shuffle <- c(x_shuffle, ret[[i]][sample(n, n)])
          }
          return(get_smap_stats(x_shuffle, lib, E))
     }))
     
     delta_mae_p = (sum(null_dist$delta_mae < delta_mae)+1) / num_shuffles
     elapsed_time <- proc.time() - start_time
     cat("(", elapsed_time[3], " sec.)\n", sep = "")
     
     save(delta_mae = delta_mae, delta_mae_p = delta_mae_p, 
          file = "test_nonlinear_aggregated.Rdata")
     return()
}

compute_nonlinearity_stock <- function()
{
     get_smap_stats <- function(x, E = NULL)
     {
          if(is.null(E))
          {
               # compute E using simplex on recruits time series
               simplex_output <- simplex(x, E = 1:8, silent = TRUE)
               best_rho_E <- simplex_output$E[which.max(simplex_output$rho)]
               best_mae_E <- simplex_output$E[which.min(simplex_output$mae)]
               E <- min(best_rho_E, best_mae_E)
          }
          
          # compute theta using s-map and E 
          smap_output <- s_map(x, E = E, silent = TRUE)
          
          best_rho <- max(smap_output$rho)
          best_mae <- min(smap_output$mae)
          return(data.frame(delta_mae = best_mae - smap_output$mae[smap_output$theta == 0]))
     }
     
     nonlinearity_for_stock <- function(stock_df, num_shuffles = 500, max_E = 8)
     {
          x <- stock_df$ret
          x <- x[is.finite(x)]
          n <- length(x)
          
          cat("calculating for actual data for ", as.character(stock_df$stk[1]), "... ", sep = "")
          start_time <- proc.time()
          simplex_output <- simplex(x, E = 1:8, silent = TRUE)
          best_rho_E <- simplex_output$E[which.max(simplex_output$rho)]
          best_mae_E <- simplex_output$E[which.min(simplex_output$mae)]
          E <- min(best_rho_E, best_mae_E)
          
          # compute theta using s-map and E 
          smap_output <- s_map(x, E = E, silent = TRUE)
          
          best_rho <- max(smap_output$rho)
          best_mae <- min(smap_output$mae)
          theta <- smap_output$theta[which.min(smap_output$mae)]
          delta_mae <- best_mae - smap_output$mae[smap_output$theta == 0]
          elapsed_time <- proc.time() - start_time
          cat("(", elapsed_time[3], " sec.)\n", sep = "")
          
          cat("calculating for random shuffles for ", as.character(stock_df$stk[1]), "... ", sep = "")
          start_time <- proc.time()
          null_dist <- do.call(rbind, lapply(1:num_shuffles, function(i) {
               x_shuffle <- x[sample(n, n)]
               return(get_smap_stats(x_shuffle, E))
          }))
          
          delta_mae_p = (sum(null_dist$delta_mae < delta_mae)+1) / num_shuffles
          elapsed_time <- proc.time() - start_time
          cat("(", elapsed_time[3], " sec.)\n", sep = "")
          
          return(list(simplex_output = simplex_output, 
                      smap_output = smap_output, 
                      E = E, 
                      theta = theta, 
                      delta_mae = delta_mae, 
                      delta_mae_p = delta_mae_p))
     }
     
     load("block_data.Rdata")
     nonlinearity_results <- lapply(block_data, nonlinearity_for_stock)
     saveRDS(nonlinearity_results, file = "results_nonlinearity_stock.RDS")
     return()
}

simple_EDM <- function()
{
     forecast <- function(stock_df)
     {  
          make_forecasts <- function(block, mu_4, sigma_4, mu_5, sigma_5)
          {
               rec4 <- block_lnlp_4(block, target_column = 2, columns = 1)
               rec5 <- block_lnlp_4(block, target_column = 3, columns = 1)
               
               rec4 <- rec4*sigma_4 + mu_4
               rec5 <- rec5*sigma_5 + mu_5
               return(rec4 + c(NA, rec5[1:(NROW(block)-1)]))
          }
          
          # set up recruits and spawners
          valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
          returns <- stock_df$ret[valid]
          years <- stock_df$yr[valid]
          spawners <- stock_df$eff_n[valid]
          recruits_4 <- stock_df$rec4_n[valid]
          mu_4 <- stock_df$rec4_mu[valid]
          sigma_4 <- stock_df$rec4_sigma[valid]
          recruits_5 <- stock_df$rec5_n[valid]
          mu_5 <- stock_df$rec5_mu[valid]
          sigma_5 <- stock_df$rec5_sigma[valid]
          
          # make block
          block <- data.frame(years = years, eff = spawners, 
                              rec4 = recruits_4, rec5 = recruits_5)
          
          if(length(returns) < 2) # check for enough data
               return(data.frame(year = NaN, obs = NaN, pred = NaN))
          
          forecasts <- make_forecasts(block, mu_4, sigma_4, mu_5, sigma_5)
          return(data.frame(year = years, obs = returns, pred = forecasts))
     }
     
     load("block_data.Rdata")
     
     # make forecasts for each stock
     results <- lapply(names(block_data), 
                       function(stk_name) {
                            cat("forecasting for ", stk_name, "... ", sep = "")
                            start_time <- proc.time()
                            output <- forecast(block_data[[stk_name]])
                            elapsed_time <- proc.time() - start_time
                            cat("(", elapsed_time[3], " sec.)\n", sep = "")
                            return(output)
                       })
     names(results) <- names(block_data)
     saveRDS(results, file = "results_simple_EDM.RDS")
     
     # compute stats
     stats <- do.call(rbind, lapply(results, function(stock_results) {
          compute_stats(stock_results$obs, stock_results$pred)
     }))
     stats$stk <- names(block_data)
     saveRDS(stats, file = "stats_simple_EDM.RDS")
     return()
}

multivariate_EDM <- function()
{
     forecast <- function(stock_df)
     {
          load("block_data.Rdata")
          env_names <- c("D_max", "D_apr", "D_may", "D_jun", 
                         "ET_apr", "ET_may", "ET_jun", 
                         "PT_apr", "PT_may", "PT_jun", "PT_jul", 
                         "PDO_win")
          
          # set up recruits and spawners
          valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
          years <- stock_df$yr[valid]
          returns <- stock_df$ret[valid]
          spawners <- stock_df$eff_n[valid]
          recruits_4 <- stock_df$rec4_n[valid]
          mu_4 <- stock_df$rec4_mu[valid]
          sigma_4 <- stock_df$rec4_sigma[valid]
          recruits_5 <- stock_df$rec5_n[valid]
          mu_5 <- stock_df$rec5_mu[valid]
          sigma_5 <- stock_df$rec5_sigma[valid]
          env <- normalize(stock_df[,env_names])
          
          # make block
          block <- data.frame(years = years, eff = spawners, 
                              rec4 = recruits_4, rec5 = recruits_5)
          block <- cbind(block, env[valid, ])
          
          if(length(returns) < 2) # check for enough data
               return(data.frame(year = NaN, obs = NaN, pred = NaN))
          
          columns <- list()
          for(E in 1:2)
          {
               columns <- c(columns, combn(env_names, E, simplify = FALSE))
          }
          columns <- lapply(columns, function(embedding) c("eff", embedding))
          columns <- c(columns, "eff")
          rec4_preds <- do.call(cbind, block_lnlp_4(block, target_column = 2, columns = columns))
          rec5_preds <- do.call(cbind, block_lnlp_4(block, target_column = 3, columns = columns))
          rec4_preds <- rec4_preds*sigma_4 + mu_4
          rec5_preds <- rec5_preds*sigma_5 + mu_5
          forecasts <- data.frame(rec4_preds + rbind(NA, rec5_preds[1:NROW(block)-1,]))
          names(forecasts) <- lapply(columns, function(v) paste(v, sep = "", collapse = ", "))
          output <- cbind(year = years, obs = returns, forecasts)
          
          return(output)
     }
     
     load("block_data.Rdata")
     
     # make forecasts for each stock
     results <- lapply(names(block_data), 
                       function(stk_name) {
                            cat("forecasting for ", stk_name, "... ", sep = "")
                            start_time <- proc.time()
                            output <- forecast(block_data[[stk_name]])
                            elapsed_time <- proc.time() - start_time
                            cat("(", elapsed_time[3], " sec.)\n", sep = "")
                            return(output)
                       })
     names(results) <- names(block_data)
     saveRDS(results, file = "results_multivariate_EDM.RDS")
     
     # compute stats
     stats <- lapply(names(block_data), function(stk_name) {
          output <- results[[stk_name]]
          stats <- do.call(rbind, lapply(3:NCOL(output), function(j) {
               compute_stats(output[,2], output[,j])
          }))
          stats$columns <- names(output)[3:NCOL(output)]
          stats$stk <- stk_name
          return(stats)        
     })
     
     stats <- lapply(stats, function(stats_df) {
          stats_df$E <- sapply(strsplit(stats_df$columns, ", "), length)
          with_eff_only <- subset(stats_df, E == 1)
          with_one_env_var <- subset(stats_df, E == 2)
          if(max(with_one_env_var$rho) <= with_eff_only$rho)
               return(subset(stats_df, E <= 2))
          best_env_var <- strsplit(with_one_env_var$columns[which.max(with_one_env_var$rho)], 
                                   ", ")[[1]][2]
          with_two_env_var <- subset(stats_df, E == 3)
          idx <- grep(best_env_var, with_two_env_var$columns)
          return(rbind(with_eff_only, with_one_env_var, with_two_env_var[idx,]))
     })
     
     saveRDS(stats, file = "stats_multivariate_EDM.RDS")
     return()
}




block_lnlp_4 <- function(block, target_column, columns, norm = FALSE)
{
     if(norm)
     {
          block[,columns] <- normalize(block[,columns])
     }
     
     lib_segments <- matrix(NA, nrow = 4, ncol = 2)
     segment_size <- NROW(block)/4
     start_index <- 1
     for(i in 1:4)
     {
          lib_segments[i,1] <- floor(start_index)
          end_index <- start_index - 1 + segment_size
          lib_segments[i,2] <- floor(end_index)
          start_index <- end_index+1
     }
     
     if(is.list(columns))
     {
          preds <- lapply(1:length(columns), function(x) {rep.int(NA, times = NROW(block))})
          for(i in 1:4)
          {
               pred_index <- lib_segments[i,1]:lib_segments[i,2]
               
               temp <- block_lnlp(block, lib = lib_segments[-i,], pred = lib_segments[i,], 
                                  target_column = target_column, tp = 0, 
                                  first_column_time = TRUE, 
                                  columns = columns, stats_only = FALSE)
               
               for(j in 1:length(columns))
                    preds[[j]][pred_index] <- temp[[j]]$model_output$pred[pred_index]
          }
     }
     else
     {
          preds <- rep.int(NA, times = NROW(block))
          for(i in 1:4)
          {
               pred_index <- lib_segments[i,1]:lib_segments[i,2]
               
               temp <- block_lnlp(block, lib = lib_segments[-i,], pred = lib_segments[i,], 
                                  target_column = target_column, tp = 0, 
                                  first_column_time = TRUE, 
                                  columns = columns, stats_only = FALSE)
               
               preds[pred_index] <- temp[[1]]$model_output$pred[pred_index]
          }
     }
     return(preds)
}

block_lnlp_4_v <- function(block, target_column, columns)
{
     lib_segments <- matrix(NA, nrow = 4, ncol = 2)
     segment_size <- NROW(block)/4
     start_index <- 1
     for(i in 1:4)
     {
          lib_segments[i,1] <- floor(start_index)
          end_index <- start_index - 1 + segment_size
          lib_segments[i,2] <- floor(end_index)
          start_index <- end_index+1
     }
     
     pred <- rep.int(NA, times = NROW(block))
     pred_var <- rep.int(NA, times = NROW(block))    
     for(i in 1:4)
     {
          pred_index <- lib_segments[i,1]:lib_segments[i,2]
          
          temp <- block_lnlp(block, lib = lib_segments[-i,], pred = lib_segments[i,], 
                             target_column = target_column, tp = 0, 
                             first_column_time = TRUE, 
                             columns = columns, stats_only = FALSE)
          
          pred[pred_index] <- temp[[1]]$model_output$pred[pred_index]
          pred_var[pred_index] <- temp[[1]]$model_output$pred_var[pred_index]
     }
     return(cbind(pred, pred_var))
}





plot_nonlinearity <- function()
{
     load("results_nonlinear_aggregated.Rdata")
     
     par(mfrow = c(1, 2), mar = c(3, 3, 0.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
     max_rho <- max(simplex_output$rho, na.rm = TRUE)
     min_rho <- min(simplex_output$rho, na.rm = TRUE)
     y_limits <- c(max(0, 1.1*min_rho - 0.1*max_rho), min(1.0, 1.1*max_rho - 0.1*min_rho))
     plot(simplex_output$E, pmax(0, simplex_output$rho), type = "l", 
          lwd = 1.5, ylim = y_limits, 
          xlab = "E", ylab = expression(rho))
     max_rho <- max(smap_output$rho, na.rm = TRUE)
     min_rho <- min(smap_output$rho, na.rm = TRUE)
     y_limits <- c(max(0, 1.1*min_rho - 0.1*max_rho), min(1.0, 1.1*max_rho - 0.1*min_rho))
     plot(smap_output$theta, 1 - smap_output$mae, type = "l", 
          lwd = 1.5, xlim = c(0, 4),  ylim = c(0.44, 0.49), 
          xlab = expression(theta), ylab = "1 - MAE")
     text(max(smap_output$theta), y_limits[2], 
          paste("E = ", E, sep = ""), adj = c(1, 1))
     return()
}


print_nonlinearity_table <- function()
{
     nonlinear_results <- readRDS("results_nonlinearity_stock.RDS")
     temp_table <- do.call(rbind, lapply(nonlinear_results, function(res) {
          return(data.frame(E = res$E, theta = res$theta, 
                            delta_mae = res$delta_mae, p_value = res$delta_mae_p))
     }))
     temp_table$stock <- rownames(temp_table)
     temp_table <- temp_table[order(temp_table$stock), 
                              c("stock", "E", "theta", "delta_mae", "p_value")]
     temp_table$"significantly nonlinear?" <- temp_table$p_value <= 0.05
     
     my_table <- xtable(temp_table, digits = 3)
     print(my_table, type = "html", file = "tables/Table_S1.html", include.rownames = FALSE)
     
     return()    
}


compute_ccm <- function()
{
     load("block_data.Rdata")
     env_names <- c("D_max", "D_apr", "D_may", "D_jun", 
                    "ET_apr", "ET_may", "ET_jun", 
                    "PT_apr", "PT_may", "PT_jun", "PT_jul", 
                    "PDO_win")
     
     ccm_table <- do.call(rbind, lapply(block_data, function(stock_df) {
          valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
          block <- stock_df[valid,]
          
          ccm_rhos <- do.call(cbind, lapply(env_names, function(env_var) {
               output <- block_lnlp(block, tp = 0, target_column = env_var, 
                                    columns = c("rec45_n", "eff_n"), silent = TRUE)
               return(output$rho)
          }))
          colnames(ccm_rhos) <- env_names
          ccm_rhos <- cbind(N = sum(valid), ccm_rhos)
          return(ccm_rhos)
     }))
     rownames(ccm_table) <- names(block_data)
     saveRDS(ccm_table, file = "results_ccm.RDS")
     return()
}

print_ccm_table <- function()
{
     ccm_table <- data.frame(readRDS("results_ccm.RDS"))
     ccm_table <- cbind("N" = ccm_table$N,
                        "95% p" = tanh(qnorm(0.95, sd = 1/sqrt(ccm_table$N - 3))), 
                        ccm_table[,2:NCOL(ccm_table)])
     my_table <- xtable(ccm_table, digits = 3)
     print(my_table, type = "html", file = "tables/Table_S3.html")
     return()
}

#---- initial set up ----
preprocess_data()

#---- nonlinear test ----
compute_nonlinearity_aggregated()
test_nonlinearity_aggregated()
compute_nonlinearity_stock()

#---- run EDM models ----
simple_EDM()
multivariate_EDM()

#---- run Ricker models ----
write_model_files()
standard_ricker()
extended_ricker()

# extract_results_for_best_models()

#---- produce figures ----
compute_seymour_ricker_params()
plot_seymour_ricker_halves() # figure 1a
compute_seymour_ricker_env_params()
plot_seymour_env_surface(plot_ricker = TRUE) # figure 1b
plot_seymour_env_surface(plot_ricker = FALSE) # figure 1c
plot_total_returns() # figure 2
plot_rho_comparison() # figure 4

#---- produce supplemental figures ----
plot_nonlinearity() # figure S2
plot_mae_comparison() # figure S3
compute_chilko_smolt_forecasts()
plot_chilko_smolt_model() # figure S4
plot_late_shuswap_CI() # figure S5

#---- produce tables ----
print_env_comparison_table() # table 1
print_nonlinearity_table() # table S1
print_comparison_table() # table S2
compute_ccm()
print_ccm_table() # table S3
print_EDM_env_models() # table S4






