#' Compute and plot posterior mean of covariate-dependent mean latent count
#'
#' @description
#' This function computes posterior mean and highest posterior density interval of
#' mean latent count, and plots it versus covariate
#'
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#' @param t a list of covariates for each dataset.
#' @param gene_ix index of genes to plot.
#' @param prob the target probability content of the interval.
#' @param xlab (optional) string. The label for x-axis in the plot.
#' @param data_names optional. The labels for each dataset in the plot.
#' @param gene_names optional. The labels for each gene in the plot
#' @param nrow parameter used to cut the plot window into subpanels.
#'
#' @return
#' a lineplot showing posterior mean of mean latent count in black, and highest
#' posterior density interval in red dashed line. Also return a list of arrays for each dataset.
#' Each array is of the form \eqn{a[c,g,i]} where \eqn{c}, \eqn{g}, \eqn{i} correspond to cell, gene and
#' statistics (mean, lower bound, upper bound), respectively.
#' @export
#'
#' @examples
#' plot_mean_latent_count(post_result = post_result, t = list(t1, t2), gene_ix = c(1,2),
#'                        prob = 0.99, xlab='t')
plot_mean_latent_count <- function(post_result, t, gene_ix, prob=0.95,
                                   xlab=NULL, data_names=NULL, gene_names=NULL, nrow=1){

  C <- sapply(t, length)
  L <- post_result$output_index
  D <- length(t)

  if(is.null(data_names)){
    data_names <- paste0('data ',1:D)
  }
  if(is.null(gene_names)){
    gene_names <- paste0('gene ', gene_ix)
  }

  if(is.null(xlab)){
    xlab <- 'x'
  }

  post_list <- lapply(1:D, function(d) {

    # a list of n_cell * n_gene matrix of mean latent count
    temp <- lapply(1:L, function(i) {

      post_result$P_C_J_D_output[[d]][[i]] %*% post_result$mu_star_1_J_output[[i]][,gene_ix]

    })

    # turn into an array
    arr <- simplify2array(temp)

    # compute mean, and quantiles
    mean <- apply(arr,1:2, mean)
    lower <- apply(arr,1:2,function(x) {
      coda::HPDinterval(coda::mcmc(x), prob = prob)[1]
    })
    upper <- apply(arr,1:2,function(x) {
      coda::HPDinterval(coda::mcmc(x), prob = prob)[2]
    })


    # turn result into an cell*gene*3 array, last dimension is for 3 statistics,
    return(simplify2array(list(mean,lower,upper)))
  })

  # plot for each dataset
  for(d in 1:D){
    par(mfrow=c(nrow, floor(length(gene_ix)/nrow)),cex.axis=1.2,cex.lab=1.2)
    for (idx in 1:length(gene_ix)) {
      o <- order(t[[d]])
      plot(t[[d]][o],post_list[[d]][o,idx,1],type='l',main=paste0(data_names[d],': ',gene_names[idx]),
           xlab=xlab, ylab = 'Mean latent count',lwd=2,ylim=range(post_list[[d]][,idx,]))
      lines(t[[d]][o],post_list[[d]][o,idx,2],col='red',lty=2) #lower
      lines(t[[d]][o],post_list[[d]][o,idx,3],col='red',lty=2) #upper
    }

  }

  return(post_list)

}



#' Compute posterior mean of latent counts
#'
#' @description
#' The function computes the posterior mean of latent counts, based on allocations.
#'
#' @importFrom pbapply pblapply
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#' @param Y a list of matrices. Each matrix is a gene-by-cell matrix of mRNA counts corresponding to a dataset.
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param mc.cores number of cores used for parallel computing.
#'
#' @return a list of matrices, each is a gene-by-cell matrix of posterior mean of latent counts for a dataset.
#' @export
#'
#' @examples
#' latent_count(post_result = post_result, Y = list(t(y1), t(y2)), opt_cl = opt_cl)
latent_count <- function(post_result, Y, opt_cl, mc.cores=8){

  D <- length(Y)
  C <- sapply(Y, ncol)
  G <- nrow(Y[[1]])
  J <- length(unique(unlist(opt_cl)))


  # posterior sample
  mu_sample <- post_result$mu_star_1_J_output
  phi_sample <- post_result$phi_star_1_J_output
  beta_sample <- post_result$Beta_output
  Z <- opt_cl
  L <- post_result$output_index


  # latent Y
  Y_latent <- lapply(1:D, function(d) {

    print(paste('Compute for data',d))

    temp <- pblapply(1:L, function(i) {
      matrix(as.vector(Y[[d]])*(as.vector(t(mu_sample[[i]][Z[[d]],]))+as.vector(t(phi_sample[[i]][Z[[d]],])))/
               (as.vector(t(mu_sample[[i]][Z[[d]],]))*rep(beta_sample[[i]][[d]],each=G)+as.vector(t(phi_sample[[i]][Z[[d]],]))) +

               as.vector(t(mu_sample[[i]][Z[[d]],]))*(as.vector(t(phi_sample[[i]][Z[[d]],]))*(1-rep(beta_sample[[i]][[d]],each=G)))/
               (as.vector(t(mu_sample[[i]][Z[[d]],]))*rep(beta_sample[[i]][[d]],each=G)+as.vector(t(phi_sample[[i]][Z[[d]],]))),

             nrow = G,
             ncol = C[d]
      )
    }, cl = mc.cores)

    return(Reduce('+', temp)/L)
  })

  return(Y_latent)
}

#' Plot observed and estimated latent counts on 2D using t-sne
#'
#' @param Y a list of matrices. Each matrix is a gene-by-cell matrix of mRNA counts corresponding to a dataset.
#' @param Y_latent a list of matrices, each is a gene-by-cell matrix of posterior mean of latent counts for a dataset.
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param data_names optional. The labels for each dataset in the plot.
#' @param color_pal optional. A vector of color names to map to clusters.
#'
#' @return two plots showing t-sne representation of the observed counts (first) and latent counts (second),
#' where observations are colored by optimal clustering, and shaped by different datasets.
#' @export
#'
#' @examples
#' plot_tsne(Y = list(t(y1), t(y2)), Y_latent = Y_latent, opt_cl = opt_cl, color_pal = c15)
plot_tsne <- function(Y, Y_latent, opt_cl, data_names=NULL, color_pal=NULL){

  # combine datasets
  Y_all <- t(do.call(cbind, Y))
  Y_latent_all <- t(do.call(cbind, Y_latent))
  C <- sapply(Y, ncol)
  D <- length(Y)

  if(is.null(data_names)){
    data_names <- paste0('data ',1:D)
  }

  # observed data
  tsne_obs <- Rtsne::Rtsne(Y_all, check_duplicates = FALSE)
  df_obs <- data.frame(tsne1=tsne_obs$Y[,1],tsne2=tsne_obs$Y[,2],
                       cluster=as.factor(unlist(opt_cl)),dataset=rep(data_names,times=C))

  g1 <- ggplot(data = df_obs, aes(x=tsne1,y=tsne2,col=cluster,pch=dataset))+
    geom_point(size=0.5)+
    theme_bw()+
    # scale_color_manual(values=c16)+
    # scale_shape_manual(values=c(1, 3))+
    labs(title='Observed counts')

  # latent counts
  tsne_latent <- Rtsne::Rtsne(Y_latent_all, check_duplicates = FALSE)
  df_latent <- data.frame(tsne1=tsne_latent$Y[,1],tsne2=tsne_latent$Y[,2],
                          cluster=as.factor(unlist(opt_cl)),dataset=rep(data_names,times=C))

  g2 <- ggplot(data = df_latent, aes(x=tsne1,y=tsne2,col=cluster,pch=dataset))+
    geom_point(size=0.5)+
    theme_bw()+
    # scale_color_manual(values=c16)+
    # scale_shape_manual(values=c(1, 3))+
    labs(title='Latent counts')

  if(!is.null(color_pal)){
    g1 <- g1+scale_color_manual(values=color_pal)
    g2 <- g2+scale_color_manual(values=color_pal)
  }

  print(ggarrange(g1, g2, nrow=1, common.legend = TRUE, legend='bottom'))
}















