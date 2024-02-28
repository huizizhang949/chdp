#' Compute posterior probability of each cluster and plot
#'
#' @description
#' The function computes the posterior mean of the probability (PP) of belonging to each cluster, for every cell.
#' Perform principal component analysis (PCA) on the combined data (concatenated along features) and plot first PC
#' against the covariate.
#'
#' @importFrom pbmcapply pbmclapply
#' @param Y a list of matrices. Each matrix is a gene-by-cell matrix of mRNA counts corresponding to a dataset.
#' @param t a list of covariates for each dataset.
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#' @param data_names the labels for each dataset in the plot.
#' @param nrow parameter used to cut the plot window into subpanels.
#' @param mc.cores number of cores used for parallel computing.
#' @param size parameter used to control the point size in the plot.
#' @param alpha parameter used to control the transparancy of the points in the plot.
#' @param xlab string. The label for x-axis in the plot.
#'
#' @return
#' one plot for each datases that contains subpanels for clusters, showing PC1 against covariate with
#' observations colored by PP.
#' \item{PP_mean}{A list of matrices. Each is a cell-by-cluster matrix of posterior probabilities, for one dataset.}
#' \item{gglist}{A list of items for each dataset. Each item has ggplot objects saving the plots.}
#'
#' @export
#'
#' @examples
#' plot_pp(post_result = post_result, Y= list(t(y1), t(y2)), t = list(t1, t2),
#'         data_names = c('data1', 'data2'), nrow = 1, mc.cores = 4, xlab='t')
plot_pp <- function(post_result, Y, t, data_names=NULL, nrow=1, mc.cores=8, size=1, alpha=1, xlab=NULL){

  Y_all <- t(do.call(cbind, Y))
  D <- length(Y)
  C <- sapply(Y, ncol)
  J <- ncol(post_result$P_C_J_D_output[[1]][[1]])
  L <- post_result$output_index

  ind <- 1:L

  if(is.null(xlab)){
    xlab <- 'x'
  }
  if(is.null(data_names)){
    data_names <- paste0('data ',1:D)
  }
  # PCA
  pc <- prcomp(Y_all,scale. = T)

  # separate by datasets
  C_cum <- c(0,cumsum(C))
  pcs <- lapply(1:D, function(d) {

    temp <- pc$x[(C_cum[d]+1):C_cum[d+1],]

    return(temp)
  })

  # compute posterior mean of the proabability
  PP_mean <- lapply(1:D, function(d) {

    # for each sample
    PP_mcmc <- pbmclapply(ind, function(i) {

      # a n_cell * J matrix
      loop.result <- vapply(1:C[d], function(cc) {

        ## Set log probability
        LP <- vapply(1:J, function(j) {

          sum(dnbinom(Y[[d]][,cc],mu=post_result$mu_star_1_J_output[[i]][j,]*post_result$Beta_output[[i]][[d]][cc],
                      size = post_result$phi_star_1_J_output[[i]][j,], log = TRUE)) +
            log(post_result$Q_J_D_output[[i]][j,d]) - (t[[d]][cc]-post_result$t_star_J_D_output[[i]][j,d])^2/2/post_result$sigma_star_2_J_D_output[[i]][j,d]

        }, FUN.VALUE = numeric(1))

        # Compute the normalizing constant
        nc <- -max(LP)
        P <- exp(LP+nc)/sum(exp(LP+nc)) # LP is a vector of length J
        return(P)
      }, FUN.VALUE = numeric(J))

      return(t(loop.result))
    }, mc.cores = mc.cores)

    # average over all samples
    return(Reduce('+',PP_mcmc)/L)
  })

  # plot for each data
  my_palette <- colorRampPalette(c("grey", 'blue','red'))

  ggs <- lapply(1:D, function(d) {

   f <- lapply(1:J, function(j) {
      ggplot(data = data.frame(t=t[[d]],PC1=pcs[[d]][,1],PP=PP_mean[[d]][,j]),aes(x=t,y=PC1,col=PP))+
        geom_point(size=size,alpha=alpha)+
        scale_colour_gradientn(colours = my_palette(n=400), limits=c(0, 1))+
        theme_bw()+
        labs(title=paste('HET: Cluster',j), x=xlab)
    })

  })

  for(d in 1:D){
    print(ggarrange(plotlist=ggs[[d]],nrow=nrow,common.legend = T,legend = 'bottom'))
  }

  return(list(PP_mean=PP_mean, gglist=ggs))
}

