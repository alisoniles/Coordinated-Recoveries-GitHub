library(rEDM)
data(tentmap_del)
str(tentmap_del)

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

smap_output <- s_map(tentmap_del, lib, pred, E = 2)
str(smap_output)
plot(rho ~ theta, data = smap_output, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

tentmap_noisy <- tentmap_del + rnorm(length(tentmap_del), sd = sd(tentmap_del) * 0.2)

## Add stochastic noise so there is tradeoff between non-linearity (theta) and forecast skill
smap_output <- s_map(tentmap_noisy, lib, pred, E = 2)
plot(rho ~ theta, data = smap_output, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

#########
# 
data(block_3sp)
str(block_3sp)


