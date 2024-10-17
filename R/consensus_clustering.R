#' Compute posterior similarity matrix from consensus clustering
#'
#' @description
#' For each combination of width and depth, compute the corresponding posterior similarity matrix.
#'
#' @param Ws a vector of candidate values for widths.
#' @param Ds a vector of candidate values for depths.
#' @param consensus_result a list of outputs from \code{gkernelHDP_mcmc}.
#' @param mc.cores number of cores used for parallel computing.
#'
#' @return a list of list of posterior similarity matrices. The outer level corresponds to different widths and the inner
#' level corresponds to different depths.
#' @export
#'
#' @examples
#' psm_list_consensus(Ws = Ws, Ds = Ds, consensus_result = consensus_result)
psm_list_consensus <- function(Ws, Ds, consensus_result, mc.cores = 8){

  # ---- a list of length = max(Ws), each contains the Z_output (a list of length = niter, within each there two lists of allocations for two datasets) ----
  mcmc_z_result <- lapply(1:max(Ws), function(w) {
    consensus_result[[w]]$Z_output
  })

  # ---- a list of length = no. of different Ds, each is a matrix of max(Ws)*n_obs ----
  # depend on values in Ds and max(Ws)
  z_mat_list_by_d <- lapply(Ds, function(d) {
    t(sapply(mcmc_z_result, function(z_output) {
      z_result <- unlist(z_output[[d]])
    }))
  })

  # ---- a list of list, the first level is of length = n_Ws, the second level is of length = n_Ds ----
  # within the list of list is a posterior similarity matrix based on one combination of d and w
  # depend on the values in Ds and Ws

  psm_list <- NULL
  for (i in 1:length(Ws)) {
    psm_list[[i]] <- pbapply::pblapply(1:length(Ds), cl=mc.cores, function(j) {
      if(Ws[i]==1) {
        mat <- mcclust::comp.psm(matrix(z_mat_list_by_d[[j]][1,],nrow=1))
      }else{
        mat <- mcclust::comp.psm(z_mat_list_by_d[[j]][1:Ws[i],])
      }
      return(mat)
    })
  }

  return(psm_list)

}


#' Plot mean absolute difference in posterior similarity matrices from consensus clustering
#'
#' @description
#' Compute mean absolute difference in posterior similarity matrices between difference values
#' of widths (or depths) and draw line plots.
#'
#' @import ggplot2
#'
#' @param Ws a vector of candidate values for widths.
#' @param Ds a vector of candidate values for depths.
#' @param psm_list output from \code{psm_list_consensus}.
#'
#' @return two line plots. The first plot shows mean absolute difference between different depths, the second for different widths.
#' @export
#'
#' @examples plot_consensus(Ws = Ws, Ds = Ds, psm_list = psm_list)
plot_consensus <- function(Ws, Ds, psm_list){

  # prepare the matrix for mean absolute difference of the consensus matrix computed from the d_{j} iteration of w_i chains
  # to that from the d_{j-1} iteration of w_i chains
  mean_abs_mat_d <- data.frame(tidyr::expand_grid(Ws,Ds,MAD=NA))
  mean_abs_mat_d$MAD <- apply(mean_abs_mat_d, 1, function(vec) {
    i <- which(Ws==vec[1])
    j <- which(Ds==vec[2])
    # the first d (smallest one) is not compared
    if(j==1) {
      return(NA)
    }else{
      mean(sqrt((psm_list[[i]][[j]]-psm_list[[i]][[j-1]])^2))
    }
  })

  mean_abs_mat_d <- mean_abs_mat_d[!is.na(mean_abs_mat_d$MAD),]

  # compare successive w
  mean_abs_mat_w <- data.frame(tidyr::expand_grid(Ws,Ds,MAD=NA))
  mean_abs_mat_w$MAD <- apply(mean_abs_mat_w, 1, function(vec) {
    i <- which(Ws==vec[1])
    j <- which(Ds==vec[2])
    # the first w (smallest one) is not compared
    if(i==1) {
      return(NA)
    }else{
      mean(sqrt((psm_list[[i]][[j]]-psm_list[[i-1]][[j]])^2))
    }
  })

  mean_abs_mat_w <- mean_abs_mat_w[!is.na(mean_abs_mat_w$MAD),]

  g1 <- ggplot(data = mean_abs_mat_d,aes(x=Ds,y=MAD,col=factor(Ws)))+
    geom_line()+
    scale_colour_discrete(name = "W")+
    labs(x='D')+
    theme_bw()+theme(legend.position='bottom')
  g2 <- ggplot(data = mean_abs_mat_w,aes(x=Ws,y=MAD,col=factor(Ds)))+
    geom_line()+
    scale_colour_discrete(name = "D")+
    labs(x='W')+
    theme_bw()+theme(legend.position='bottom')

  print(gridExtra::grid.arrange(g1,g2,ncol=2))
  return(list(d=mean_abs_mat_d,w=mean_abs_mat_w))
}


#' Compute optimal clustering from consensus clustering results
#'
#' @description
#' Provide the optimal clustering and posterior similarity matrix.
#'
#' @param Width a suitable value for width.
#' @param Depth a suitable value for depth.
#' @param consensus_result a list of outputs from \code{gkernelHDP_mcmc}.
#'
#' @return a list containing the following:
#' \item{opt_cl}{a list of vectors. Each vector stores the optimal clustering for one dataset. }
#' \item{psm}{posterior similarity matrix based on the provided width and depth.}
#' @export
#'
#' @examples
#' opt_cl_consensus(W = 100, D = 1000, consensus_result = consensus_result)
opt_cl_consensus <- function(Width, Depth, consensus_result){

  # size of each data
  C <- sapply(consensus_result[[1]]$Z_output[[1]],length)
  # number of datasets
  D <- length(consensus_result[[1]]$Z_output[[1]])

  n_cum <- c(0,cumsum(C))

  # MCMC samples from the given width and depth
  mcmc_z_result <- t(sapply(1:max(Width), function(w) {
    unlist(consensus_result[[w]]$Z_output[[Depth]])
  }))

  # psm
  psm <- mcclust::comp.psm(mcmc_z_result)

  # optimal clustering
  VI_test <- mcclust.ext::minVI(psm,mcmc_z_result,method=('all'),include.greedy=FALSE)
  VI_test <- VI_test$cl[1,]

  # get the optimal clustering for each dataset
  opt_cl <- NULL
  for(d in 1:D){
    opt_cl[[d]] <- VI_test[(n_cum[d]+1):n_cum[d+1]]
  }

  return(list(opt_cl=opt_cl, psm=psm))
}









