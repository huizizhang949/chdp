#' Compute optimal clustering for model based on VAR likelihood
#'
#' @description
#' Provide the optimal clustering and posterior similarity matrix.
#'
#' @param pkernel_output a list of outputs from \code{pkernelHDP_mcmc} with unknown allocations.
#'
#' @return a list containing the following:
#' \item{opt_cl}{a list of vectors. Each vector stores the optimal clustering for one dataset. }
#' \item{psm}{posterior similarity matrix combining samples from all chains.}
#' @export
#'
#' @examples
#' opt_cl_var(pkernel_output = pkernel_output)
opt_cl_var <- function(pkernel_output){

  # number of chains
  W <- length(pkernel_output)
  # size of each data
  C <- sapply(pkernel_output[[1]]$Z_output[[1]],length)
  # number of datasets
  D <- length(pkernel_output[[1]]$Z_output[[1]])

  n_cum <- c(0,cumsum(C))

  # MCMC samples of Z from all chains
  mcmc_z_result <- lapply(1:W, function(w) {

    temp <- pkernel_output[[w]]$Z_output
    # row corresponds to samples, columns are observations
    z_mat <- lapply(temp, function(l) {
      unlist(l)
    })

    z_mat <- do.call(rbind, z_mat)
    return(z_mat)

  })

  mcmc_z_result <- do.call(rbind,mcmc_z_result)


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
