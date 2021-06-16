# Function for calculate absolute and relative causal strengths
# of a pair of time series
#
# Syntax: causal_matrix = causal_decomposition(s1,s2,rnoise,nensemble)
# Parameters:
# s1: time series 1
# s2: time series 2
# rnoise: the level of noise used in ensemble empirical mode decomposition
#         the rnoise value represents the fraction of standard deviation of
#         time series (e.g., 0.1)
# nensemble: the number of ensemble to average the results of noise-assisted
#            empirical mode decomposition.
#            The error of IMFs equals to rnoise/sqrt(nensemble)
#
# Output:
# The output of causal_decomposition function contains a four by n matrix,
# where n is the number of IMFs decomposed from the data. The first column
# in the matrix indicates the relative causal strength from first time series
# to the second time series and vise versa in the second column. The third
# column indicates the absolute causal strength from the first time series
# to the second time series, and vice versa in the fourth column.
#
# Example: load ecosystem_data;
#          causal_matrix = causal_decomposition(DIDI,PARA,0.35,1000);
#
# Ver 1.0: Albert C. Yang, MD, PhD 7/11/2018
#
# Referece

causal_decomposition <- function (s1,s2,rnoise,nensemble){

# normalize time series
s1 <- s1-mean(s1)
s1 <- s1/std(s1)
s2 <- s2-mean(s2)
s2 <- s2/std(s2)


# Decompose input data to Intrinsic Mode Functions (IMFs) with the 
# Complete Ensemble Empirical Mode Decomposition with Adaptive Noise (CEEMDAN) algorithm, a variant of EEMD.
imfs1 <- eemd(s1, noise_strength = rnoise, ensemble_size = nensemble)
imfs2 <- eemd(s2, noise_strength = rnoise, ensemble_size = nensemble)

# determine imf size (less one because of the residual column)
imfsize <- size(imfs1)
imfsize <- imfsize[2]-1 

# calculate phase coherence between paired IMFs
ps <- phasefcimf(imfs1,imfs2,imfsize)

# calculate variance of IMFs
v1 <- nvar(imfs1[,1:imfsize])
v2 <- nvar(imfs2[,1:imfsize])

# initiate phase vectors
p12 <- matrix(nrow=imfsize,ncol=1)
p21 <- matrix(nrow=imfsize,ncol=1)

# compute causal strength for each paired IMFs
for (i in 1:imfsize){

   # remove an IMF and do the redecomposition
   s1r <- s1-imfs1[,i]
   s2r <- s2-imfs2[,i]
   imfs1r <- eemd(s1r, noise_strength = rnoise, ensemble_size = nensemble)
   imfs2r <- eemd(s2r, noise_strength = rnoise, ensemble_size = nensemble)

   # recalculate phase coherence between paired IMFs
   ps12 <- phasefcimf(imfs1,imfs2r,imfsize)
   ps21 <- phasefcimf(imfs1r,imfs2,imfsize)

   # calculate absolute causal strength using variance-weighted Euclidian distance
   # between the phase coherence of the original IMFs and redecomposed IMFs
   p12[i] <- sqrt(sum(v1*v2*(ps12-ps)^2)/sum(v1*v2))
   p21[i] <- sqrt(sum(v1*v2*(ps21-ps)^2)/sum(v1*v2))

}

causal_matrix <- matrix(nrow=imfsize,ncol=4)
# calculate relative causal strengths from absolute causal strengths (p12 and p21)
# the output is 4 by n matrix (n: number of IMFs)
# the first column indicates relative causal strength from time series 1 to 2
# the second column indicates relative causal strength from time series 2 to 1
# the third column indicates absolute causal strength from time series 1 to 2
# the fourth column indicates absolute causal strength from time series 2 to 1
for (i in 1:imfsize){

    if (p12[i]<0.05 & p21[i]<0.05){
        alpha <- 1
    } else {
        alpha <- 0
    }

    causal_matrix[i,] <- rbind((alpha+p12[i])/(2*alpha+p12[i]+p21[i]), (alpha+p21[i])/(2*alpha+p12[i]+p21[i]), p12[i], p21[i])
}

return (causal_matrix)
}
