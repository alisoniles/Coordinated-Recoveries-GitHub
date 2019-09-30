

rm(list = ls())

library(rEDM)
data(tentmap_del)
str(tentmap_del)

#########
#########
#########

## Determine embedding dimension using Simplex Projection

## Split library and prediction sets of data
lib <- c(1, 100)
pred <- c(101, 500)

## use simplex algorithm to forecast points in pred set from nearest neighbors... 
#  ...and identify optimal embedding dimension
simplex_output <- simplex(tentmap_del, lib, pred)
str(simplex_output)

plot(rho ~ E, data = simplex_output, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")
## ... forecast skill peaks at two embeddings

simplex_output <- simplex(tentmap_del, lib, pred, E = 2, tp = 1:10)
plot(rho ~ tp, data = simplex_output, type = "l", xlab = "Time to Prediction (tp)", ylab = "Forecast Skill (rho)")
## and forecast skill decays rapidly with how far out we are predicting

## Use S-map with defined embedding dimension of 2
## ... and identify level of non-linearity
smap_output <- s_map(tentmap_del, lib, pred, E = 2)
str(smap_output)
plot(rho ~ theta, data = smap_output, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")
## forecast skill continues to improve with higher theta 
## in other words, the more we weight local over distant and use fewer points to define function

## Add stochastic noise so there is tradeoff
## Between allowing the S-map to be more and more wiggly 
## and the number of data points that are used to determine the function.
tentmap_noisy <- tentmap_del + rnorm(length(tentmap_del), sd = sd(tentmap_del) * 0.2)

smap_output <- s_map(tentmap_noisy, lib, pred, E = 2)
plot(rho ~ theta, data = smap_output, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")
## Now we see that an intermediate theta (some wiggliness) maximizes predictive skill

#########
#########
#########

## Multivariate Models

data(block_3sp)
str(block_3sp)

lib <- c(1, 100)
pred <- c(101, 200)

cols <- c(1, 2, 4)
target <- 1

block_lnlp_output <- block_lnlp(block_3sp, lib = lib, pred = pred, 
                                columns = cols,  target_column = target, 
                                first_column_time = TRUE)
str(block_lnlp_output)

block_lnlp_output <- block_lnlp(block_3sp, lib = lib, pred = pred, 
                                columns = c("x_t", "x_t-1", "y_t"),
                                target_column = "x_t", stats_only = FALSE)

str(block_lnlp_output)
model_output <- block_lnlp_output[[1]]$model_output # KI moved the [[1]] before the $

observed <- model_output$obs
predicted <- model_output$pred
time <- model_output$time

plot_range <- range(c(observed, predicted), na.rm = TRUE)
plot(observed, predicted, xlim = plot_range, ylim = plot_range, xlab = "Observed", 
     ylab = "Predicted", asp = 1)
abline(a = 0, b = 1, lty = 2, col = "blue")

plot(time, observed, type = "l", xlim = c(101,150))
lines(time, predicted, col="blue")

############### ############### ############### ############### ############### 
## Real Data ##
### Tilman ####
############### ############### ############### ############### ############### 

rm(list = ls(all.names = TRUE))

data(e120_biodiversity)

normalize <- function(x, ...) {
  (x - mean(x, ...))/sd(x, ...)
}

# separate time column from data
vars <- c("AbvBioAnnProd", "noh020tot", "invrichness", "SummerPrecip.mm.")
composite_ts <- e120_biodiversity[, vars]

# normalize each time series within a plot
data_by_plot <- split(composite_ts, e120_biodiversity$Plot)
normalized_data <- lapply(data_by_plot, function(df) sapply(df, normalize))
composite_ts <- cbind(Year = e120_biodiversity$Year, data.frame(do.call(rbind, 
                                                                        normalized_data)))
# KI Note: Now we have the year and 4 variables with no lagged variables

# To prevent lagged vectors from being constructed that span separate plots, 
# we need to create an appropriate index variable that identifies different segments.

segments_end <- cumsum(sapply(data_by_plot, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)

# Choose random segments for prediction
set.seed(2312)
rndlib <- sample(1:NROW(segments), floor(NROW(segments) * 0.75))
composite_lib <- segments[rndlib, ]
composite_pred <- segments[-rndlib, ]

# Now we have normalized the data, split data into segments by plot, 
# and created a series of chunks for training and prediction...

# Because the time series for precipitation does not vary among replicates, 
# we also need to construct separate variables for analyzing precipitation dynamics:
precip_ts <- unique(e120_biodiversity[, c("Year", "SummerPrecip.mm.")])
precip_ts <- precip_ts[order(precip_ts$Year), ]

#################
## Run simplex ##
## and S-map ####
#################

simplex_out <- lapply(names(composite_ts)[2:4], function(var) {
simplex(composite_ts[, c("Year", var)], E = 2:4, lib = composite_lib, pred = composite_pred)
})
simplex_out[[length(simplex_out) + 1]] <- simplex(precip_ts, E = 2:5)
names(simplex_out) <- names(composite_ts)[-1]

par(mar = c(4, 4, 1, 1), mfrow = c(2, 3), mgp = c(2.5, 1, 0))  # set up margins for plotting
out <- lapply(names(simplex_out), function(var) {
  plot(simplex_out[[var]]$E, simplex_out[[var]]$rho, type = "l", xlab = "Embedding Dimension (E)", 
       ylab = "Forecast Skill (rho)", main = var)
})

best_E <- sapply(simplex_out, function(df) {df$E[which.max(df$rho)]})
best_E

smap_out <- lapply(names(composite_ts)[2:4], function(var) {
  s_map(composite_ts[, c("Year", var)], E = best_E[var], lib = composite_lib, 
        pred = composite_pred)
})
smap_out[[length(smap_out) + 1]] <- s_map(precip_ts, E = best_E[length(smap_out) + 
                                                                  1])
names(smap_out) <- names(simplex_out)

par(mar = c(4, 4, 1, 1), mfrow = c(2, 2), mgp = c(2.5, 1, 0))  # set up margins for plotting
lapply(names(smap_out), function(var) {
  plot(smap_out[[var]]$theta, smap_out[[var]]$rho, type = "l", xlab = "Nonlinearity (theta)", 
       ylab = "Forecast Skill (rho)", main = var)
})



################## 
## Multivariate ##
##### Models #####
################## 

data_by_plot <- split(composite_ts, e120_biodiversity$Plot)

block_data <- do.call(rbind, lapply(data_by_plot, function(df) {
  n <- NROW(df)
  temp <- data.frame(Year = df$Year)
  temp$AB_tm <- df$AbvBioAnnProd
  temp$AB_tm1 <- c(NA, temp$AB_tm[-n])
  temp$AB_tm2 <- c(NA, temp$AB_tm1[-n])
  temp$AB_tm3 <- c(NA, temp$AB_tm2[-n])
  
  temp$NO_tm <- df$noh020tot
  temp$NO_tm1 <- c(NA, temp$NO_tm[-n])
  temp$NO_tm2 <- c(NA, temp$NO_tm1[-n])
  temp$NO_tm3 <- c(NA, temp$NO_tm2[-n])
  
  temp$IV_tm <- df$invrichness
  temp$IV_tm1 <- c(NA, temp$IV_tm[-n])
  temp$IV_tm2 <- c(NA, temp$IV_tm1[-n])
  temp$IV_tm3 <- c(NA, temp$IV_tm2[-n])
  
  temp$PR_tm <- df$SummerPrecip.mm
  temp$PR_tm1 <- c(NA, temp$PR_tm[-n])
  temp$PR_tm2 <- c(NA, temp$PR_tm1[-n])
  temp$PR_tm3 <- c(NA, temp$PR_tm2[-n])
  
  return(temp)
}))
head(block_data[, 1:5], 20)

AB_columns <- c("AB_tm", "AB_tm1", "AB_tm2")
AB_output <- block_lnlp(block_data, lib = composite_lib, pred = composite_pred, 
                        columns = AB_columns, target_column = 1, stats_only = FALSE, first_column_time = TRUE)

model_output <- AB_output[[1]]$model_output 

observed <- model_output$obs
predicted <- model_output$pred
time <- model_output$time

plot_range <- range(c(observed, predicted), na.rm = TRUE)
plot(observed, predicted, xlim = plot_range, ylim = plot_range, xlab = "Observed", 
     ylab = "Predicted", asp = 1)
abline(a = 0, b = 1, lty = 2, col = "blue")

ABNO_columns <- c("AB_tm", "AB_tm1", "AB_tm2", "NO_tm", "NO_tm1", "NO_tm2")
ABNO_output <- block_lnlp(block_data, lib = composite_lib, pred = composite_pred, 
                          columns = ABNO_columns, target_column = 1, stats_only = FALSE, first_column_time = TRUE)

observed_AB <- AB_output[[1]]$model_output$obs
predicted_AB <- AB_output[[1]]$model_output$pred

observed_ABNO <- ABNO_output[[1]]$model_output$obs
predicted_ABNO <- ABNO_output[[1]]$model_output$pred

dev.off()

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0), pty = "s")  # set up margins for plotting
plot_range <- range(c(observed_AB, predicted_AB), na.rm = TRUE)
plot(observed_AB, predicted_AB, xlim = plot_range, ylim = plot_range, xlab = "Observed", 
     ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "darkgrey", lwd = 2)
abline(lm(predicted_AB ~ observed_AB), col = "black", lty = 3, lwd = 2)

points(observed_ABNO, predicted_ABNO, pch = 2, col = "red")
abline(lm(predicted_ABNO ~ observed_ABNO), col = "red", lty = 3, lwd = 2)

legend("topleft", legend = c(paste("rho =", round(AB_output[[1]]$stats$rho, 
                                                  2)), paste("rho =", round(ABNO_output[[1]]$stats$rho, 2))), lty = 3, lwd = 2, 
       col = c("black", "red"), bty = "n")
