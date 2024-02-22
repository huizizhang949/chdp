library(pbmcapply)

# ------------- single-cell clustering with Gaussian kernel to incorporate latent time as a covariate -------
# prepare the datasets, and latent time
G=10
y1 <- sim1.data[1:120,1:G]; z1 <- sim1.data[1:120,G+1]; t1 <- sim1.data[1:120,G+2]
y2 <- sim1.data[121:240,1:G]; z2 <- sim1.data[121:240,G+1]; t2 <- sim1.data[121:240,G+2]
p_j1 <- sim1.data[1:120,(G+3):(G+4)]; p_j2 <- sim1.data[121:240,(G+3):(G+4)]

# -------- clustering estimate -----------
# run consensus clustering, with 100 chains and 5000 iterations
Width <- 2
Depth <- 2000

Width <- 100
Depth <- 5000

# for reproducibility
set.seed(14663)
seeds <- sample(1:1e8, Width)
# parallel computing on 8 cores
consensus_result <- pbmclapply(1:Width, function(i) {

  set.seed(seeds[i])
  b_inits <- runif(2,-5,5)
  alpha_inits <- runif(1,0,10)
  alpha_0_inits <- runif(1,0,10)

  gkernelHDP_mcmc(Y = list(t(y1), t(y2)), t = list(t1, t2), J = 4, niter = Depth, burn_in = 0, thinning = 1,
                 empirical = TRUE, empirical_z = TRUE,
                 b_initial = b_inits, alpha_initial = alpha_inits, alpha_0_initial = alpha_0_inits,
                 beta.mean = 0.6)
}, mc.cores = 8)

# choose W and D based on the plot (the first element is the baseline)
Ds <- c(1,100,500,1000,1500,2000)
Ws <- c(1,2)

Ds <- c(1, 100, 500, seq(1000, 5000, by=1000))
Ws <- c(1,10,20,50,100)

psm_list <- psm_list_consensus(Ws, Ds, consensus_result = consensus_result)
# plot mean absolute difference to determine suitable width and depth
plot_consensus(Ws, Ds, psm_list)

# after choosing suitable W and D, compute optimal clustering
# below assume the best width and depth are the largest possible values
opt <- opt_cl_consensus(W = Width, D = Depth, consensus_result = consensus_result)
# optimal clustering
opt$opt_cl
# posterior similarity matrix based on the chosen W and D
opt$psm

# plot psm
# without distinguishing datasets
plotpsm(psm = opt$psm)
# distinguish datasets
plot_psm(psm.tot = opt$psm, size = rep(120,2))

# compare with the truth, note that the cluster labels can be different while the partition is the same
# e.g. cluster 1 in opt_cl may correspond to cluster 2 in the truth
table(c(z1,z2), unlist(opt$opt_cl))
mclust::adjustedRandIndex(c(z1,z2), unlist(opt$opt_cl))

# summary of cluster sizes
cluster_size(opt_cl = opt$opt_cl, data_names = c('data1', 'data2'))

# ----------- for inference of the other parameters, fix the optimal clustering and run a post-processing step ---------

set.seed(3)
post_result <- gkernelHDP_mcmc(Y = list(t(y1), t(y2)), t = list(t1, t2), J = length(unique(unlist(opt$opt_cl))),
                               niter = 1000, burn_in = 0, thinning = 1,
                               empirical = TRUE, Z_fix = opt$opt_cl, beta.mean = 0.6)

# plot posterior samples for time-dependent probabilities
# the function also accepts truth if available. remember to match the labels in optimal clustering
# to the true labels when providing the true probabilities
plot_probx(post_result = post_result, opt_cl = opt$opt_cl, x = list(t1, t2), mfrow = c(1,2),
           data_names = c('data1', 'data2'), xlab = 't', thinning = 5,
           truth = list(p_j1[,c(2,1)], p_j2[,c(2,1)]))




