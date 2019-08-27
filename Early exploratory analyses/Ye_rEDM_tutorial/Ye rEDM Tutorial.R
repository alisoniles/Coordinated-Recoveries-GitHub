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

block_lnlp_output <- block_lnlp(block_3sp, lib = lib, pred = pred, 
                                columns = c("x_t", "x_t-1", "y_t"),
                                target_column = "x_t", stats_only = FALSE)

model_output <- block_lnlp_output$model_output[[1]]

observed <- model_output$obs
predicted <- model_output$pred

plot_range <- range(c(observed, predicted), na.rm = TRUE)
plot(observed, predicted, xlim = plot_range, ylim = plot_range, xlab = "Observed", 
     ylab = "Predicted", asp = 1)
abline(a = 0, b = 1, lty = 2, col = "blue")
