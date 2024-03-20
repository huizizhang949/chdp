# ------------- single-cell clustering with Gaussian kernel to incorporate latent time as a covariate -------
# prepare the datasets, and latent time
G=100; C=c(200,300); J=3
Y1 <- sim1.data[1:C[1],1:G]; z1 <- sim1.data[1:C[1],G+1]; t1 <- sim1.data[1:C[1],G+2]
Y2 <- sim1.data[(C[1]+1):sum(C),1:G]; z2 <- sim1.data[(C[1]+1):sum(C),G+1]; t2 <- sim1.data[(C[1]+1):sum(C),G+2]
p_j1 <- sim1.data[1:C[1],-(1:(G+2))]; p_j2 <- sim1.data[(C[1]+1):sum(C),-(1:(G+2))]

load(file='man/opt.RData')
load(file='man/psm_list.RData')
load(file='man/post_result.RData')
# -------- clustering estimate -----------
# an example of running consensus clustering, with 10 chains and 500 iterations
Width <- 100
Depth <- 200

# for reproducibility
set.seed(14663)
seeds <- sample(1:1e8, Width)
# parallel computing on 8 cores
consensus_result <- pblapply(1:Width, function(i) {

  set.seed(seeds[i])
  b_inits <- runif(2,-5,5)
  alpha_inits <- runif(1,0,10)
  alpha_0_inits <- runif(1,0,10)

  gkernelHDP_mcmc(Y = list(t(Y1), t(Y2)), t = list(t1, t2), J = 7, niter = Depth, burn_in = 0, thinning = 1,
                 empirical = TRUE, empirical_z = TRUE,
                 b_initial = b_inits, alpha_initial = alpha_inits, alpha_0_initial = alpha_0_inits,
                 beta.mean = 0.6)
}, cl=8)
# save(consensus_result, file='man/consensus_result.RData')

# choose W and D based on the plot (the first element is the baseline)
# specify candidate values
Ds <- c(1,seq(20,200,by=10))
Ws <- c(1,seq(10,100,by=5))

# a list of posterior similarity matrices for different combination W and D
psm_list <- psm_list_consensus(Ws = Ws, Ds = Ds, consensus_result = consensus_result)
# plot mean absolute difference to determine suitable width and depth
pdf(file='demo/figures/plot_consensus.pdf',w=10,h=8)
plot_consensus(Ws = Ws, Ds = Ds, psm_list = psm_list)
dev.off()

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
pdf('demo/figures/psm_all_nb.pdf',w=6,h=6)
mcclust.ext::plotpsm(psm = opt$psm)
dev.off()
# distinguish datasets
pdf('demo/figures/psm_single_nb.pdf',w=6,h=6)
plot_psm(psm.tot = opt$psm, size = C)
dev.off()

# compare with the truth, note that the cluster labels can be different while the partition is the same
# e.g. cluster 1 in opt_cl may correspond to cluster 2 in the truth
table(c(z1,z2), unlist(opt_cl))
mclust::adjustedRandIndex(c(z1,z2), unlist(opt_cl))

# summary of cluster sizes
pdf('demo/figures/cluster_size_nb.pdf',w=10,h=6)
cluster_size(opt_cl = opt_cl, data_names = c('data 1', 'data 2'))
dev.off()
# ----------- for inference of the other parameters, fix the optimal clustering and run a post-processing step ---------

# post-processing: fix Z to the optimal one from consensus clustering
set.seed(3)
post_result <- gkernelHDP_mcmc(Y = list(t(Y1), t(Y2)), t = list(t1, t2),
                               niter = 7000, burn_in = 5000, thinning = 1,
                               empirical = TRUE, Z_fix = opt_cl, beta.mean = 0.6)
# 14.661 mins
# save(post_result,file='man/post_result.RData')

# plot posterior samples for time-dependent probabilities
# the function also accepts truth if available. remember to match the labels in optimal clustering
# to the true labels when providing the true probabilities
# save individual figure by hand, w=8, h=6
plot_c_prob(post_result = post_result, opt_cl = opt_cl, t = list(t1, t2), mfrow = c(1,2),
            data_names = c('data 1', 'data 2'), xlab = 't', thinning = 5, plot.empty.cluster = FALSE,
            truth = list(p_j1[,c(2,1,3)], p_j2[,c(2,1,3)]))

# compute posterior mean of the probability (PP) of belonging to each cluster, for every cell
# in the mean time, perform PCA after concatenating the datasets (along gene axis)
# plot first principal component (PC1) against time, with observations colored by PP
# save individual figure by hand, w=12, h=6
PP_mean_result <- plot_pp(post_result = post_result, Y = list(t(Y1), t(Y2)), t = list(t1, t2),
                          data_names = c('data 1', 'data 2'), nrow = 1, mc.cores = 4, xlab='t')
# PP
PP_mean <- PP_mean_result$PP_mean
# ggplot objects
# first dataset
print(PP_mean_result$gglist[[1]])
# second dataset
print(PP_mean_result$gglist[[2]])

# compute and plot the posterior mean and highest posterior density (HPD) interval for mean latent count
# as a function of covariate, also return the results of posterior mean and HPD.
# save individual figure by hand, w=12, h=6
mean_latent_counts <- plot_mean_latent_count(post_result = post_result, t = list(t1, t2), gene_ix = c(1,2),
                       prob = 0.99, xlab='t')
# to draw plot manually
d=1;g=1
# mean
plot(t1[order(t1)], mean_latent_counts[[d]][order(t1),g,1],type='l', xlab='t', ylab='Mean latent counts',
     ylim = range(mean_latent_counts[[d]][,g,]))
# lower bound
lines(t1[order(t1)], mean_latent_counts[[d]][order(t1),g,2],col='red',lty=2)
# upper bound
lines(t1[order(t1)], mean_latent_counts[[d]][order(t1),g,3],col='red',lty=2)


# compute posterior mean of the latent count, conditional on allocations
Y_latent <- latent_count(post_result = post_result, Y = list(t(Y1), t(Y2)), opt_cl = opt_cl, mc.cores = 2)

# visualize observed and latent counts on 2D using t-sne
pdf('demo/figures/tsne_nb.pdf',w=10,h=6)
set.seed(1)
plot_tsne(Y = list(t(Y1), t(Y2)), Y_latent = Y_latent, opt_cl = opt_cl, color_pal = c15)
dev.off()

# ------ marker genes ------------
# ---- globally differentially expressed (DE) and dispersed (DD) marker genes ------
global_result <- global_marker_genes(post_result = post_result, threshold = c(2.5,2.5),
                                     alpha_m = 0.3, alpha_d = 0.3, mc.cores = 2)

# plot tail probabilities against mean absolute log-fold change (LFC) for all genes, in terms of
# mean expression and dispersion parameter. Also summarize the overlaps between global DE and DD genes
pdf('demo/figures/global_nb.pdf',w=12,h=6)
plot_global_marker_genes(global_output = global_result)
dev.off()
# show four heatmaps for estimated mean and dispersion parameters
# first two show estimates for all genes across clusters, with genes ordered by decreasing tail probabilities
# from top to bottom. A black horizontal line separates global and non-global genes.
# last two plots show marker genes only, where the estimates are normalized such that
# the mean of the estimated parameters across clusters is zero.

# save individual figure by hand, w=8, h=8
global_marker_genes_heatmaps(global_output = global_result, post_result = post_result)

# plot heatmaps of observed counts, with genes ordered by decreasing tail probabilities (DE) from top to bottom,
# a red dashed line separates global and non-globle DE genes. Cells are separated
# by clusters (yellow vertical solid line) and also separated by datasets within each cluster (yellow vertical dashed line)
pdf('demo/figures/obs_heatmap.pdf',w=12,h=6)
observed_counts_heatmap(Y = list(t(Y1),t(Y2)), opt_cl = opt_cl, global_output = global_result)
dev.off()

# heatmap for estimated latent counts, in the same fashion as observed counts
pdf('demo/figures/latent_heatmap.pdf',w=12,h=6)
latent_counts_heatmap(Y_latent = Y_latent, opt_cl = opt_cl, global_output = global_result)
dev.off()
# ---- locally differentially expressed (DE) and dispersed (DD) marker genes ------
local_result <- local_marker_genes(post_result = post_result, threshold = c(1.2,1.2),
                                   alpha_m = 0.3, alpha_d = 0.3, mc.cores = 2)

# Return four ggplot objects
# 1. Tail probabilities against mean absolute LFC based on mean expression for each cluster, and the threshold to decide local DE genes.
# 2. A bar-chart showing the number of local DE genes for each cluster.
# 3. Tail probabilities against mean absolute LFC based on dispersion for each cluster, and the threshold to decide local DD genes.
# 4. A bar-chart showing the number of local DD genes for each cluster.
ggs_local <- plot_local_marker_genes(local_output = local_result, nrow = 1)
# plot
gridExtra::grid.arrange(grobs=ggs_local,nrow=2)

# Return two lists of ggplot objects:
# mu.ggs: Each item is the heatmap of posterior mean of log mean expression,
# for local DE genes of the cluster. Genes are ordered by decreasing tail probabilities from top to bottom.
# phi.ggs: Each item is the heatmap of posterior mean of log dispersion,
# for local DD genes of the cluster. Genes are ordered by decreasing tail probabilities from top to bottom.
heatmaps_ggs_local <- local_marker_genes_heatmaps(local_output = local_result, post_result = post_result)
# plot
gridExtra::grid.arrange(grobs=heatmaps_ggs_local$mu.ggs,nrow=1)
gridExtra::grid.arrange(grobs=heatmaps_ggs_local$phi.ggs,nrow=1)


# ----------- posterior predictive checks ---------
# use mixed predictive distribution, where dispersion is generated from its prior

# ----- a single replicate --------------------
# generate a single replicate, compare the relationship between statistics between replicate and observed data
# compare differences in statistics between replicate and observed data
pdf('demo/figures/ppc_single_nb.pdf',w=12,h=6)
set.seed(2)
plot_ppc_single(post_result = post_result, Y = list(t(Y1), t(Y2)), opt_cl = opt_cl, data_names = c('data 1', 'data 2'))
dev.off()
# ------- multiple replicates ---------------
# generate multiple replicates, and compute key statistics (gene-wise) for replicated and observed data.
set.seed(3)
ppc_multiple_df <- ppc_multiple(post_result = post_result, Y = list(t(Y1), t(Y2)),
                                opt_cl = opt_cl, number_rep = 200, mc.cores = 2)
# plot the results from multiple replicates. compare density plots for each statistic between
# replicated (grey) and observed data (red)
pdf('demo/figures/ppc_multiple_nb.pdf',w=12,h=6)
plot_ppc_multiple(ppc_multiple_df = ppc_multiple_df, data_names = c('data 1', 'data 2'))
dev.off()

# ----- posterior predictive p-values ----------
# compute posterior predictive p-values for three discrepancy measures (Chi-squared statistic, Freeman-Tukey statistic
# and dropout probabilities)
set.seed(4)
ppp_mixed_result <- ppp_mixed(post_result = post_result, Y = list(t(Y1), t(Y2)),
                              opt_cl = opt_cl, number_rep = 200, mc.cores = 2)

# plot histograms of p-values for each discrepancy measure
ppp_hist(ppp_output = ppp_mixed_result)

# compute p-values condition on optimal clustering
set.seed(5)
ppp_mixed_cl_result <- ppp_mixed_cl(post_result = post_result, Y = list(t(Y1), t(Y2)),
                                  opt_cl = opt_cl, number_rep = 200, mc.cores = 2)

# pdf(file='~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Mac/Project/code/rpkg/ppp_cl_dropout.pdf',w=12,h=6)
# plot histograms for each discrepancy measure, for each cluster and each dataset
ppp_hist_cl(ppp_output = ppp_mixed_cl_result, nrow = 1, data_names = c('data 1', 'data 2'))
# dev.off()






