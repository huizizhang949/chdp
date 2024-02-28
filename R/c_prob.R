#' Plot posterior samples for covariate-dependent probabilities
#'
#' @description
#' This function plots covariate-dependent probabilities against the covariate
#' based on the MCMC samples from the post-processing step.
#'
#'
#' @param post_result output from \code{gkernelHDP_mcmc} or \code{pkernelHDP_mcmc} where clustering is fixed to the optimal one.
#' @param opt_cl optimal clustering used for the post-processing step.
#' @param t a list of vectors. Each vector denotes the covariate in one dataset.
#' @param mfrow parameter used to cut the plot window into subpanels
#' @param data_names the labels for each dataset in the plot.
#' @param xlab string. The label for x-axis in the plot.
#' @param thinning if provided, apply thinning when plotting MCMC samples.
#' @param truth a list of matrices. Each matrix stores the covariate-dependent probabilities for one dataset,
#' where rows correspond to observations and columns correspond to clusters.
#' @param plot.empty.cluster if \code{FALSE}, only plot for occupied clusters.
#'
#' @return plots showing MCMC samples of covariate-dependent probabilities for every cluster in each dataset.
#' @export
#'
#' @examples
#' plot_c_prob(post_result = post_result, opt_cl = opt$opt_cl, t = list(t1, t2),
#'             mfrow = c(1,2), data_names = c('data1', 'data2'),
#'             xlab = 't', thinning = 5, truth = list(p_j1, p_j2))
plot_c_prob <- function(post_result, opt_cl, t, mfrow=c(2,2), data_names=NULL, xlab=NULL, thinning=NULL, truth=NULL,
                       plot.empty.cluster=FALSE){

  # number of clusters
  J <- ncol(post_result$P_C_J_D_output[[1]][[1]])
  # number of datasets
  D <- length(post_result$P_C_J_D_output)
  # number of MCMC samples
  L <- post_result$output_index

  if(is.null(xlab)){
    xlab <- 'x'
  }
  if(is.null(data_names)){
    data_names <- paste0('data ',1:D)
  }
  # apply thinning if provided
  if(is.null(thinning)){
    ix <- 1:L
  }else{
    ix <- seq(thinning, L, by=thinning)
  }

  for(d in 1:D){
    par(mfrow=mfrow)
    if(plot.empty.cluster){
      Js <- 1:J
    }else{
      Js <- as.numeric(names(table(opt_cl[[d]])))
    }

    for (j in Js) {
      o <- order(t[[d]])

      plot(sort(t[[d]]),sort(t[[d]]),ylim=c(0,1),type='n',main=paste0(data_names[d],': cluster ',j),
           xlab=xlab, ylab='probability')

      for (i in ix) {
        lines(sort(t[[d]]),post_result$P_C_J_D_output[[d]][[i]][o,j],col='grey')
      }
      if(!is.null(truth)){
        lines(sort(t[[d]]),truth[[d]][o,j],col='red')
      }
    }
  }

}

