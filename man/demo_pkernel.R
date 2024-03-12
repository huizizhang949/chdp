# -------------- clustering with VAR likelihood and time as external covariate ------------
# prepare the datasets, and time
G=2; C=rep(151,2); J=3
Y1 <- sim2.data[1:C[1],1:G]; z1 <- sim2.data[1:C[1],G+1]; t1 <- sim2.data[1:C[1],G+2]
Y2 <- sim2.data[(C[1]+1):sum(C),1:G]; z2 <- sim2.data[(C[1]+1):sum(C),G+1]; t2 <- sim2.data[(C[1]+1):sum(C),G+2]
p_j1 <- sim2.data[1:C[1],-(1:(G+2))]; p_j2 <- sim2.data[(C[1]+1):sum(C),-(1:(G+2))]

# the first observation will not be clustered in each dataset
C=C-1

load(file='man/test1.RData')
load(file='man/opt_var.RData') #opt
load(file='man/post_result_var.RData') #opt
# -------- clustering estimate -----------
set.seed(493)
test1 <- pkernelHDP_mcmc(Y=list(Y1,Y2), t=list(t1,t2), niter=2000, J=4, burn_in = 1000, thinning = 1,
                         empirical_z=TRUE, target_accept=0.234)
# 32.834 secs
# save(test1,file='man/test1.RData')

# compute optimal clustering, here the argument pkernel_output should contain a list
# of outputs from pkernelHDP_mcmc.
opt <- opt_cl_var(pkernel_output = list(test1))
# save(opt,file='man/opt_var.RData')

# optimal clustering
opt_cl <- opt$opt_cl
# posterior similarity matrix based on the chosen W and D
opt$psm

# plot psm
# without distinguishing datasets
mcclust.ext::plotpsm(psm = opt$psm)
# distinguish datasets, note that the first observation in each dataset is not clustered
plot_psm(psm.tot = opt$psm, size = C)

# compare with the truth, note that the cluster labels can be different while the partition is the same
# e.g. cluster 1 in opt_cl may correspond to cluster 2 in the truth
table(c(z1[-1],z2[-1]), unlist(opt_cl))
mclust::adjustedRandIndex(c(z1[-1],z2[-1]), unlist(opt_cl))

# summary of cluster sizes
cluster_size(opt_cl = opt_cl, data_names = c('data 1', 'data 2'))

# for each dataset, plot each feature against time and color observations by time
plot_var_cluster(Y=list(Y1,Y2), t=list(t1,t2), opt_cl = opt_cl, color_pal = c25)

# ----------- for inference of the other parameters, fix the optimal clustering and run a post-processing step ---------

# post-processing: fix Z to the optimal one from the above step
set.seed(666)
post_result <- pkernelHDP_mcmc(Y=list(Y1,Y2), t=list(t1,t2), niter=10000, burn_in = 6000, thinning = 2,
                              Z_fix = opt_cl, target_accept=0.234)

# 47.131 secs
save(post_result,file='man/post_result_var.RData')

# plot posterior samples for time-dependent probabilities
# the function also accepts truth if available. remember to match the labels in optimal clustering
# to the true labels when providing the true probabilities

# remove the time for the first observation (not clustered)
plot_c_prob(post_result = post_result, opt_cl = opt_cl, t = list(t1[-1], t2[-1]), mfrow = c(1,2),
            data_names = c('data 1', 'data 2'), xlab = 't', thinning = 5, plot.empty.cluster = TRUE,
            truth = list(p_j1[,c(2,1)], p_j2[,c(2,1)]))


# ------ conditional mean ------
# calculates the mean conditional on the past observation as well as credible intervals
conditional_mean_result <- conditional_mean(Y=list(Y1,Y2), t=list(t1,t2), opt_cl = opt_cl, post_result = post_result,
                                            prob = c(0.005,0.995), mc.cores = 2)

# plots the posterior mean and credible intervals for the conditional mean
plot_conditional_mean(conditional_mean_output = conditional_mean_result, data = 1, cluster = 1)



# ----- prediction ----------
# predict future observations, and computes credible intervals.
set.seed(2)
pred_trend_result <- pred_trend(post_result = post_result, Y=list(Y1,Y2),
                                t_pred = seq(1,1.5,by=diff(t1)[1]), prob = c(0.005,0.995), mc.cores = 2)

# plot posterior predictive means and associated credible intervals, in addition to observed data (black), for each feature.
plot_pred_trend(pred_trend_output = pred_trend_result, Y=list(Y1,Y2), t_obs = list(t1,t2))


# calculate time-dependent probabilities of belonging to each cluster for future time points
pred_c_prob_result <- pred_c_prob(post_result = post_result, t_pred = seq(1,1.5,by=diff(t1)[1]), mc.cores = 2)

# plot posterior samples of time-dependent probabilities for future time points and observed time points,
# for the specified dataset and cluster
plot_pred_c_prob(pred_c_prob_output = pred_c_prob_result, t_obs = t1, data = 1, cluster = 1, thinning = 5, color_pal = c25)


# ------ posterior predictive check --------
# generate multiple replicated datasets for posterior predictive check, and computes
# the mean and credible intervals
set.seed(2)
ppc_var_result <- ppc_var(post_result = post_result, Y=list(Y1,Y2), t = list(t1,t2),
                          opt_cl = opt_cl, number_rep = 200, prob = c(0.005,0.995), mc.cores = 2)

# plot observed data in black, and the posterior mean of the replicated datasets (red line) and credible intervals (red area)
plot_ppc_var(ppc_var_output = ppc_var_result, Y=list(Y1,Y2))






