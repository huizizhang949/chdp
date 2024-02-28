# ------------- single-cell clustering with Gaussian kernel to incorporate latent time as a covariate -------
# prepare the datasets, and latent time
G=10
y1 <- sim1.data[1:120,1:G]; z1 <- sim1.data[1:120,G+1]; t1 <- sim1.data[1:120,G+2]
y2 <- sim1.data[121:240,1:G]; z2 <- sim1.data[121:240,G+1]; t2 <- sim1.data[121:240,G+2]
p_j1 <- sim1.data[1:120,(G+3):(G+4)]; p_j2 <- sim1.data[121:240,(G+3):(G+4)]

load(file='man/consensus_result.RData')
load(file='man/post_result.RData')
load(file='man/opt.RData')
# -------- clustering estimate -----------
# an example of running consensus clustering, with 10 chains and 500 iterations
Width <- 10
Depth <- 500

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
# save(consensus_result, file='man/consensus_result.RData')


# choose W and D based on the plot (the first element is the baseline)
# specify candidate values
Ds <- c(1,seq(50,500,by=50))
Ws <- c(1,seq(2,10,by=2))

# a list of posterior similarity matrices for different combination W and D
psm_list <- psm_list_consensus(Ws = Ws, Ds = Ds, consensus_result = consensus_result)
# plot mean absolute difference to determine suitable width and depth
plot_consensus(Ws = Ws, Ds = Ds, psm_list = psm_list)

# after choosing suitable W and D, compute optimal clustering
# below assume the best width and depth are the largest possible values
opt <- opt_cl_consensus(W = Width, D = Depth, consensus_result = consensus_result)
# save(opt, file='man/opt.RData')

# optimal clustering
opt_cl <- opt$opt_cl
# posterior similarity matrix based on the chosen W and D
opt$psm

# plot psm
# without distinguishing datasets
mcclust.ext::plotpsm(psm = opt$psm)
# distinguish datasets
plot_psm(psm.tot = opt$psm, size = rep(120,2))

# compare with the truth, note that the cluster labels can be different while the partition is the same
# e.g. cluster 1 in opt_cl may correspond to cluster 2 in the truth
table(c(z1,z2), unlist(opt_cl))
mclust::adjustedRandIndex(c(z1,z2), unlist(opt_cl))

# summary of cluster sizes
cluster_size(opt_cl = opt_cl, data_names = c('data 1', 'data 2'))

# ----------- for inference of the other parameters, fix the optimal clustering and run a post-processing step ---------

# post-processing: fix Z to the optimal one from consensus clustering
set.seed(3)
post_result <- gkernelHDP_mcmc(Y = list(t(y1), t(y2)), t = list(t1, t2),
                               niter = 2000, burn_in = 1500, thinning = 1,
                               empirical = TRUE, Z_fix = opt_cl, beta.mean = 0.6)
# save(post_result,file='man/post_result.RData')

# plot posterior samples for time-dependent probabilities
# the function also accepts truth if available. remember to match the labels in optimal clustering
# to the true labels when providing the true probabilities
plot_c_prob(post_result = post_result, opt_cl = opt_cl, t = list(t1, t2), mfrow = c(1,2),
            data_names = c('data 1', 'data 2'), xlab = 't', thinning = 5,
            truth = list(p_j1[,c(2,1)], p_j2[,c(2,1)]))

# compute posterior mean of the probability (PP) of belonging to each cluster, for every cell
# in the mean time, perform PCA after concatenating the datasets (along gene axis)
# plot first principal component (PC1) against time, colored by PP.
PP_mean_result <- plot_pp(post_result = post_result, Y = list(t(y1), t(y2)), t = list(t1, t2),
                          data_names = c('data1', 'data2'), nrow = 1, mc.cores = 4, xlab='t')
# PP
PP_mean <- PP_mean_result$PP_mean
# ggplot objects
# first dataset
print(PP_mean_result$gglist[[1]])
# second dataset
print(PP_mean_result$gglist[[2]])

# compute the posterior mean and highest posterior density (HPD) interval for mean latent count
# as a function of covariate, also return the results of posterior mean and HPD.
mean_latent_counts <- plot_mean_latent_count(post_result = post_result, t = list(t1, t2), gene_ix = c(1,2),
                       prob = 0.99, xlab='t', mc.cores = 4)
# to draw plot manually
d=1;g=1
# mean
plot(t1[order(t1)], mean_latent_counts[[d]][order(t1),g,1],type='l', xlab='t', ylab='Mean latent counts',
     ylim = range(mean_latent_counts[[d]][,g,]))
# lower bound
lines(t1[order(t1)], mean_latent_counts[[d]][order(t1),g,2],col='red',lty=2)
# upper bound
lines(t1[order(t1)], mean_latent_counts[[d]][order(t1),g,3],col='red',lty=2)


# compute posterior mean of the latent count
Y_latent <- latent_count(post_result = post_result, Y = list(t(y1), t(y2)), opt_cl = opt_cl)

# visualize observed and latent counts on 2D using t-sne
plot_tsne(Y = list(t(y1), t(y2)), Y_latent = Y_latent, opt_cl = opt_cl, color_pal = c16)

