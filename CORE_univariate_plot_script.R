
rm(list = ls())

par(mar = c(4, 4, 1, 1), mfrow = c(2, 3), mgp = c(2.5, 1, 0))  # set up margins for plotting
out <- lapply(names(simplex_out), function(var) {
  plot(simplex_out[[var]]$E, simplex_out[[var]]$rho, type = "l", xlab = "Embedding Dimension (E)", 
       ylab = "Forecast Skill (rho)", main = var)
})


par(mar = c(4, 4, 1, 1), mfrow = c(2, 3), mgp = c(2.5, 1, 0))  # set up margins for plotting
lapply(names(smap_out), function(var) {
  plot(smap_out[[var]]$theta, smap_out[[var]]$rho, type = "l", xlab = "Nonlinearity (theta)", 
       ylab = "Forecast Skill (rho)", main = var)
})

univ_cols <- c("rec_t", "rec_tm1", "rec_tm2", "rec_tm3")
univ_output <- block_lnlp(lagged_data, lib = composite_lib, pred = composite_pred, 
                          columns = univ_cols , target_column = 1, stats_only = FALSE, first_column_time = TRUE)

univ_model_output <- univ_output[[1]]$model_output 

univ_observed <- univ_model_output$obs
univ_predicted <- univ_model_output$pred

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), pty = "s")  # set up margins for plotting
plot_range <- range(c(univ_observed, univ_predicted), na.rm = TRUE)
plot(univ_observed, univ_predicted, xlim = plot_range, ylim = plot_range, xlab = "Observed", 
     ylab = "Predicted", asp = 1)
abline(a = 0, b = 1, lty = 2, col = "blue")

# multi_cols <- c("rec_t", "rec_tm1", "rec_tm2", "rec_tm3", "eff_t", "eff_tm1","eff_tm2","eff_tm3")
multi_output <- block_lnlp(lagged_data, lib = composite_lib, pred = composite_pred, columns = multi_cols , target_column = 1, stats_only = FALSE, first_column_time = TRUE)

multi_model_output <- multi_output[[1]]$model_output 

multi_observed <- multi_model_output$obs
multi_predicted <- multi_model_output$pred


par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), pty = "s")  # set up margins for plotting
plot_range <- range(c(univ_observed, univ_predicted), na.rm = TRUE)
plot(univ_observed, univ_predicted, xlim = plot_range, ylim = plot_range, xlab = "Observed", 
     ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "darkgrey", lwd = 2)
abline(lm(univ_predicted ~ univ_observed), col = "black", lty = 3, lwd = 2)

points(multi_observed, multi_predicted, pch = 2, col = "red")
abline(lm(multi_predicted ~ multi_observed), col = "red", lty = 3, lwd = 2)

legend("topleft", legend = c(paste("rho =",round(univ_output[[1]]$stats$rho,2)), 
                             paste("rho =", round(multi_output[[1]]$stats$rho, 2))), 
       lty = 3, lwd = 2, col = c("black", "red"), bty = "n")




