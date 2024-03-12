
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chdp

<!-- badges: start -->
<!-- badges: end -->

The chdp package is developed for implementing covariate-dependent
hierarchical Dirichlet process, which allows for clustering across
related datasets and incorporating the information of external
covariate. The covariate can be flexibly included through different
kernel functions. Two different applications are included in the
package, one for single-cell clustering with a Gaussian kernel, and the
other for clustering time-series data based on vector autoregression
(VAR) with a periodic kernel.

## Installation

You can install the development version of chdp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("huizizhang949/chdp")
```

# Example 1 - clustering single-cell data

First, prepare the datasets including the external covariate - latent
time:

``` r
library(chdp)
# G - dimension (number of genes)
# C - number of data points in each dataset
# J - number of clusters
# Y - matrix: rows are observations
# z - data allocation
# t - covariate (latent time)
# p_j - covariate-dependent time probabilities
G=100; C=c(200,300); J=3
# first dataset
Y1 <- sim1.data[1:C[1],1:G]; z1 <- sim1.data[1:C[1],G+1]; 
t1 <- sim1.data[1:C[1],G+2]; p_j1 <- sim1.data[1:C[1],-(1:(G+2))]
# second dataset
Y2 <- sim1.data[(C[1]+1):sum(C),1:G]; z2 <- sim1.data[(C[1]+1):sum(C),G+1]; 
t2 <- sim1.data[(C[1]+1):sum(C),G+2]; p_j2 <- sim1.data[(C[1]+1):sum(C),-(1:(G+2))]
```

## Consensus clustering

To implement consensus clustering, first choose the width (number of
chains) and depth (length of each chain). Below is an example of running
consensus clustering, with 100 chains and 200 iterations on 8 cores

``` r
Width <- 100
Depth <- 200
```

``` r
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
```

Specify a range of candidate values (the first element is the baseline):

``` r
Ds <- c(1,seq(20,200,by=10))
Ws <- c(1,seq(10,100,by=5))
```

Compute posterior similarity matrices (PSM) for different combination W
and D

``` r
psm_list <- psm_list_consensus(Ws = Ws, Ds = Ds, consensus_result = consensus_result)
```

Choose a suitable value based on mean absolute difference (MAD) between
PSMs

``` r
plot_consensus(Ws = Ws, Ds = Ds, psm_list = psm_list)
```

<p align="center">
<img src="demo/figures/plot_consensus.png" width="60%" style="display: block; margin: auto;" />
</p>

After choosing suitable W and D, compute optimal clustering

``` r
# below assume the best width and depth are the largest possible values
opt <- opt_cl_consensus(W = Width, D = Depth, consensus_result = consensus_result)
# optimal clustering
opt_cl <- opt$opt_cl
# posterior similarity matrix based on the chosen W and D
opt$psm
```

To understand the uncertainty in data partition

``` r
# plot psm
# without distinguishing datasets
mcclust.ext::plotpsm(psm = opt$psm)
# distinguish datasets
plot_psm(psm.tot = opt$psm, size = C)
```

<p align="center">

<div class="figure" style="text-align: center">

<img src="demo/figures/psm_all_nb.png" alt=" " width="49%" height="20%" /><img src="demo/figures/psm_single_nb.png" alt=" " width="49%" height="20%" />
<p class="caption">
</p>

</div>

</p>

And compare compare with the truth, if known

``` r
# note that the cluster labels can be different while the partition is the same
table(c(z1,z2), unlist(opt_cl))
mclust::adjustedRandIndex(c(z1,z2), unlist(opt_cl))
```

Blow shows summary of cluster sizes

``` r
cluster_size(opt_cl = opt_cl, data_names = c('data 1', 'data 2'))
```

<p align="center">
<img src="demo/figures/cluster_size_nb.png" width="60%" style="display: block; margin: auto;" />
</p>

## Post-processing

For inference of the other parameters, fix the optimal clustering and
run a post-processing step

``` r
# post-processing: fix Z to the optimal one from consensus clustering
set.seed(3)
post_result <- gkernelHDP_mcmc(Y = list(t(Y1), t(Y2)), t = list(t1, t2),
                               niter = 7000, burn_in = 5000, thinning = 1,
                               empirical = TRUE, Z_fix = opt_cl, beta.mean = 0.6)
```

### Covariate-related results and latent counts

Check posterior samples for time-dependent probabilities

``` r
# the function also accepts truth if available. remember to match the labels in optimal clustering
# to the true labels when providing the true probabilities
plot_c_prob(post_result = post_result, opt_cl = opt_cl, t = list(t1, t2), mfrow = c(1,2),
            data_names = c('data 1', 'data 2'), xlab = 't', thinning = 5, plot.empty.cluster = FALSE,
            truth = list(p_j1[,c(2,1,3)], p_j2[,c(2,1,3)]))
```

<p align="center">
<img src="demo/figures/plot_c_prob_nb1.png" width="49%" height="20%" style="display: block; margin: auto;" /><img src="demo/figures/plot_c_prob_nb2.png" width="49%" height="20%" style="display: block; margin: auto;" />
</p>

Below we compute posterior mean of the probability (PP) of belonging to
each cluster, for every cell. In the mean time, perform PCA after
concatenating the datasets (along gene axis).

``` r
# plot first principal component (PC1) against time, with observations colored by PP
PP_mean_result <- plot_pp(post_result = post_result, Y = list(t(Y1), t(Y2)), t = list(t1, t2),
                          data_names = c('data 1', 'data 2'), nrow = 1, mc.cores = 4, xlab='t')
```

<p align="center">
<img src="demo/figures/plot_pp_mean_nb1.png" width="49%" height="20%" style="display: block; margin: auto;" /><img src="demo/figures/plot_pp_mean_nb2.png" width="49%" height="20%" style="display: block; margin: auto;" />
</p>

And further we can compute posterior mean and highest posterior density
(HPD) interval for mean latent count as a function of covariate, without
conditioning on allocations.

``` r
# the functionalso return the results of posterior mean and HPD.
mean_latent_counts <- plot_mean_latent_count(post_result = post_result, t = list(t1, t2), gene_ix = c(1,2),
                       prob = 0.99, xlab='t')
```

<p align="center">
<img src="demo/figures/mean_latent_count_nb1.png" width="49%" height="20%" style="display: block; margin: auto;" /><img src="demo/figures/mean_latent_count_nb2.png" width="49%" height="20%" style="display: block; margin: auto;" />
</p>

Compute posterior mean of the latent count, conditional on allocations

``` r
Y_latent <- latent_count(post_result = post_result, Y = list(t(Y1), t(Y2)), opt_cl = opt_cl, mc.cores = 2)
```

Visualize observed and latent counts on 2D using t-sne

<p align="center">
<img src="demo/figures/tsne_nb.png" width="60%" style="display: block; margin: auto;" />
</p>

### Marker genes

To identify global marker genes

``` r
# ---- globally differentially expressed (DE) and dispersed (DD) marker genes ------
global_result <- global_marker_genes(post_result = post_result, threshold = c(2.5,2.5),
                                     alpha_m = 0.3, alpha_d = 0.3, mc.cores = 2)
```

Get a summary of global marker genes

``` r
# plot tail probabilities against mean absolute log-fold change (LFC) for all genes, in terms of
# mean expression and dispersion parameter. Also summarize the overlaps between global DE and DD genes
plot_global_marker_genes(global_output = global_result)
```

<p align="center">
<img src="demo/figures/global_nb.png" width="60%" style="display: block; margin: auto;" />
</p>

Below shows four heatmaps for estimated (log) mean and dispersion
parameters

``` r
# first two show estimates for all genes across clusters, with genes ordered by decreasing tail probabilities
# from top to bottom. A black horizontal line separates global and non-global genes.
# last two plots show marker genes only, where the estimates are normalized such that
# the mean of the estimated parameters across clusters is zero.
global_marker_genes_heatmaps(global_output = global_result, post_result = post_result)
```

<p align="center">

<div class="figure" style="text-align: center">

<img src="demo/figures/DE_global.png" alt=" " width="40%" height="20%" /><img src="demo/figures/DD_global.png" alt=" " width="40%" height="20%" /><img src="demo/figures/DE_relative.png" alt=" " width="40%" height="20%" /><img src="demo/figures/DD_relative.png" alt=" " width="40%" height="20%" />
<p class="caption">
</p>

</div>

</p>

Given marker genes, we plot heatmaps of observed counts, with genes
ordered by decreasing tail probabilities (DE) from top to bottom

``` r
# a red dashed line separates global and non-globle DE genes. Cells are separated
# by clusters (yellow vertical solid line) and also separated by datasets within each cluster (yellow vertical dashed line)
observed_counts_heatmap(Y = list(t(Y1),t(Y2)), opt_cl = opt_cl, global_output = global_result)
```

<p align="center">
<img src="demo/figures/obs_heatmap.png" width="60%" style="display: block; margin: auto;" />
</p>

And do the same for estimated latent counts

``` r
latent_counts_heatmap(Y_latent = Y_latent, opt_cl = opt_cl, global_output = global_result)
```

<p align="center">
<img src="demo/figures/latent_heatmap.png" width="60%" style="display: block; margin: auto;" />
</p>

As for local marker genes,

``` r
# ---- locally differentially expressed (DE) and dispersed (DD) marker genes ------
local_result <- local_marker_genes(post_result = post_result, threshold = c(1.2,1.2),
                                   alpha_m = 0.3, alpha_d = 0.3, mc.cores = 2)
```

Below functions produce similar plots as for the global marker genes

``` r
# Return four ggplot objects
# 1. Tail probabilities against mean absolute LFC based on mean expression for each cluster, and the threshold to decide local DE genes.
# 2. A bar-chart showing the number of local DE genes for each cluster.
# 3. Tail probabilities against mean absolute LFC based on dispersion for each cluster, and the threshold to decide local DD genes.
# 4. A bar-chart showing the number of local DD genes for each cluster.
ggs_local <- plot_local_marker_genes(local_output = local_result, nrow = 1)
# plot
gridExtra::grid.arrange(grobs=ggs_local,nrow=2)
```

``` r
# Return two lists of ggplot objects:
# mu.ggs: Each item is the heatmap of posterior mean of log mean expression,
# for local DE genes of the cluster. Genes are ordered by decreasing tail probabilities from top to bottom.
# phi.ggs: Each item is the heatmap of posterior mean of log dispersion,
# for local DD genes of the cluster. Genes are ordered by decreasing tail probabilities from top to bottom.
heatmaps_ggs_local <- local_marker_genes_heatmaps(local_output = local_result, post_result = post_result)
# plot
gridExtra::grid.arrange(grobs=heatmaps_ggs_local$mu.ggs,nrow=1)
gridExtra::grid.arrange(grobs=heatmaps_ggs_local$phi.ggs,nrow=1)
```

## Posterior predictive check

We apply mixed predictive distribution, where dispersion is generated
from its prior

``` r
# ----- a single replicate --------------------
# generate a single replicate
# 1. compare the relationship between statistics between replicate (black) and observed data (red)
# 2. compare differences in statistics between replicate and observed data
set.seed(2)
plot_ppc_single(post_result = post_result, Y = list(t(Y1), t(Y2)), opt_cl = opt_cl, data_names = c('data 1', 'data 2'))
```

<p align="center">
<img src="demo/figures/ppc_single_nb1.png" width="60%" style="display: block; margin: auto;" /><img src="demo/figures/ppc_single_nb2.png" width="60%" style="display: block; margin: auto;" />
</p>

``` r
# ------- multiple replicates ---------------
# generate multiple replicates, and compute key statistics (gene-wise) for replicated and observed data.
set.seed(3)
ppc_multiple_df <- ppc_multiple(post_result = post_result, Y = list(t(Y1), t(Y2)),
                                opt_cl = opt_cl, number_rep = 200, mc.cores = 2)
# plot the results from multiple replicates. compare density plots for each statistic between
# replicated (grey) and observed data (red)
plot_ppc_multiple(ppc_multiple_df = ppc_multiple_df, data_names = c('data 1', 'data 2'))
```

<p align="center">
<img src="demo/figures/ppc_multiple_nb1.png" width="60%" style="display: block; margin: auto;" /><img src="demo/figures/ppc_multiple_nb2.png" width="60%" style="display: block; margin: auto;" />
</p>

Finally, posterior predictive p-values are also available from the
following functions

``` r
# ----- posterior predictive p-values ----------
# compute posterior predictive p-values for three discrepancy measures 
# based on Chi-squared statistic, Freeman-Tukey statistic and dropout probabilities
set.seed(4)
ppp_mixed_result <- ppp_mixed(post_result = post_result, Y = list(t(Y1), t(Y2)),
                              opt_cl = opt_cl, number_rep = 200, mc.cores = 2)

# plot histograms of p-values for each discrepancy measure
ppp_hist(ppp_output = ppp_mixed_result)

# compute p-values condition on optimal clustering
set.seed(5)
ppp_mixed_cl_result <- ppp_mixed_cl(post_result = post_result, Y = list(t(Y1), t(Y2)),
                                  opt_cl = opt_cl, number_rep = 200, mc.cores = 2)

# plot histograms for each discrepancy measure, for each cluster and each dataset
ppp_hist_cl(ppp_output = ppp_mixed_cl_result, nrow = 1, data_names = c('data 1', 'data 2'))
```

# Example 2 - clustering time-series data with VAR likelihood

Prepare the datasets including the external covariate - time:

``` r
rm(list=ls())
# G - dimension (number of features)
# C - number of data points in each dataset
# J - number of clusters
# Y - matrix: rows are observations
# z - data allocation
# t - covariate (latent time)
# p_j - covariate-dependent time probabilities
G=2; C=rep(151,2); J=3
# first dataset
Y1 <- sim2.data[1:C[1],1:G]; z1 <- sim2.data[1:C[1],G+1]; 
t1 <- sim2.data[1:C[1],G+2]; p_j1 <- sim2.data[1:C[1],-(1:(G+2))]
# second dataset
Y2 <- sim2.data[(C[1]+1):sum(C),1:G]; z2 <- sim2.data[(C[1]+1):sum(C),G+1]; 
t2 <- sim2.data[(C[1]+1):sum(C),G+2]; p_j2 <- sim2.data[(C[1]+1):sum(C),-(1:(G+2))]

# the first observation will not be clustered in each dataset
C=C-1
```

## Clustering estimate

Similar to above, we first estimate clustering and check the uncertainty
in data partition

``` r
set.seed(3901)
test1 <- pkernelHDP_mcmc(Y=list(Y1,Y2), t=list(t1,t2), niter=5000, J=6, burn_in = 3000, thinning = 1,
                         empirical_z=TRUE, target_accept=0.234, mu_h = -2)

# compute optimal clustering, here the argument pkernel_output should contain a list
# of outputs from pkernelHDP_mcmc.
opt <- opt_cl_var(pkernel_output = list(test1))

# optimal clustering
opt_cl <- opt$opt_cl
# posterior similarity matrix based on the chosen W and D
opt$psm

# plot psm
# without distinguishing datasets
mcclust.ext::plotpsm(psm = opt$psm)

# distinguish datasets, note that the first observation in each dataset is not clustered
plot_psm(psm.tot = opt$psm, size = C)
```

<p align="center">

<div class="figure" style="text-align: center">

<img src="demo/figures/psm_all_var.png" alt=" " width="49%" height="20%" /><img src="demo/figures/psm_single_var.png" alt=" " width="49%" height="20%" />
<p class="caption">
</p>

</div>

</p>

And compare compare with the truth, if known

``` r
# note that the cluster labels can be different while the partition is the same
# and first observation in each dataset is not clustered
table(c(z1[-1],z2[-1]), unlist(opt_cl))
mclust::adjustedRandIndex(c(z1[-1],z2[-1]), unlist(opt_cl))
```

Blow shows summary of cluster sizes

``` r
cluster_size(opt_cl = opt_cl, data_names = c('data 1', 'data 2'))
```

<p align="center">
<img src="demo/figures/cluster_size_var.png" width="60%" style="display: block; margin: auto;" />
</p>

Specifically, for time-series data, we visualize the cluster by plotting
each feature against time and coloring observations by clusters

``` r
plot_var_cluster(Y=list(Y1,Y2), t=list(t1,t2), opt_cl = opt_cl, color_pal = c25)
```

<p align="center">
<img src="demo/figures/plot_cluster_var.png" width="60%" style="display: block; margin: auto;" />
</p>

## Post-processing

For inference of the other parameters, fix the optimal clustering and
run a post-processing step

``` r
# post-processing: fix Z to the optimal one from the above step
set.seed(666)
post_result <- pkernelHDP_mcmc(Y=list(Y1,Y2), t=list(t1,t2), niter=16000, burn_in = 12000, thinning = 2,
                              Z_fix = opt_cl, target_accept=0.234, mu_h = -2)
```

### Covariate-related results

Check posterior samples for time-dependent probabilities, which shows
periodicity as expected.

``` r
# the function also accepts truth if available. remember to match the labels in optimal clustering
# to the true labels when providing the true probabilities

# remove the time for the first observation (not clustered)
plot_c_prob(post_result = post_result, opt_cl = opt_cl, t = list(t1[-1], t2[-1]), mfrow = c(1,2),
            data_names = c('data 1', 'data 2'), xlab = 't', thinning = 5, plot.empty.cluster = FALSE,
            truth = list(p_j1, p_j2))
```

<p align="center">
<img src="demo/figures/plot_c_prob_var1.png" width="49%" height="20%" style="display: block; margin: auto;" /><img src="demo/figures/plot_c_prob_var2.png" width="49%" height="20%" style="display: block; margin: auto;" />
</p>

### Conditional mean

We examine the mean of each observation in the normal likelihood, which
relies on the past observation.

``` r
# ------ conditional mean ------
# calculates the mean conditional on the past observation as well as credible intervals
conditional_mean_result <- conditional_mean(Y=list(Y1,Y2), t=list(t1,t2), opt_cl = opt_cl, post_result = post_result,
                                            prob = c(0.005,0.995), mc.cores = 2)

# plots the posterior mean and credible intervals for the conditional mean, for the specified dataset and cluster
plot_conditional_mean(conditional_mean_output = conditional_mean_result, data = 1, cluster = 1)
```

<p align="center">
<img src="demo/figures/plot_conditional_mean.png" width="60%" style="display: block; margin: auto;" />
</p>

### Prediction

Predict future observations (future 10 observations with the same time
increment as the observed data). Note that uncertainty increases as
predictions are made further away.

``` r
# predict future observations, and computes credible intervals.
set.seed(2)
pred_trend_result <- pred_trend(post_result = post_result, Y=list(Y1,Y2),
                                t_pred = 1+(diff(t1)[1])*(1:10), prob = c(0.005,0.995), mc.cores = 2)

# plot posterior predictive means and associated credible intervals, in addition to observed data (black), for each feature.
plot_pred_trend(pred_trend_output = pred_trend_result, Y=list(Y1,Y2), t_obs = list(t1,t2))
```

<p align="center">

<div class="figure" style="text-align: center">

<img src="demo/figures/plot_pred_trend1.png" alt=" " width="49%" height="20%" /><img src="demo/figures/plot_pred_trend2.png" alt=" " width="49%" height="20%" />
<p class="caption">
</p>

</div>

</p>

To calculate time-dependent probabilities of belonging to each cluster
for future time points

``` r
pred_c_prob_result <- pred_c_prob(post_result = post_result, t_pred = seq(1,2,by=diff(t1)[1]), mc.cores = 2)

# plot posterior samples of time-dependent probabilities for future time points and observed time points,
# for the specified dataset and cluster
plot_pred_c_prob(pred_c_prob_output = pred_c_prob_result, t_obs = t1, data = 1, cluster = 1, thinning = 5, color_pal = c25)
```

<p align="center">
<img src="demo/figures/plot_pred_c_prob.png" width="40%" style="display: block; margin: auto;" />
</p>

## Posterior predictive check

For multiple replicated datasets, we compare the mean and credible
intervals of the replicates with the observed data

``` r
# generate multiple replicated datasets for posterior predictive check, and computes
# the mean and credible intervals
set.seed(2)
ppc_var_result <- ppc_var(post_result = post_result, Y=list(Y1,Y2), t = list(t1,t2),
                          opt_cl = opt_cl, number_rep = 200, prob = c(0.005,0.995), mc.cores = 2)

# plot observed data in black, and the posterior mean of the replicated datasets (red line) and credible intervals (red area)
plot_ppc_var(ppc_var_output = ppc_var_result, Y=list(Y1,Y2))
```

<p align="center">

<div class="figure" style="text-align: center">

<img src="demo/figures/plot_ppc_var1.png" alt=" " width="49%" height="20%" /><img src="demo/figures/plot_ppc_var2.png" alt=" " width="49%" height="20%" />
<p class="caption">
</p>

</div>

</p>
