#' Posterior predictive checks based on multiple replicates
#'
#' @description
#' The function generates multiple replicated datasets and compute four statistics for the replicates and observed data.
#'
#' @importFrom pbapply pblapply
#'
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#' @param Y a list of matrices. Each matrix is a gene-by-cell matrix of mRNA counts corresponding to a dataset.
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param number_rep number of replicates to generate.
#' @param mc.cores number of cores used for parallel computing.
#'
#' @return The output contains the following components:
#' \item{rep_Y_statistics}{statistics for replicated datasets.}
#' \item{Y_statistics}{statistics for the observed dataset.}
#' \item{number_rep}{number of replicated datasets.}
#' @export
#'
#' @examples
#' ppc_multiple_df <- ppc_multiple(post_result = post_result, Y = list(t(y1), t(y2)),
#'                                 opt_cl = opt_cl, number_rep = 10, mc.cores = 2)
ppc_multiple <- function(post_result, Y, opt_cl, number_rep, mc.cores=8){


  G <- nrow(Y[[1]])
  D <- length(Y)
  C <- sapply(Y, ncol)
  J <- ncol(post_result$P_C_J_D_output[[1]][[1]])
  L <- post_result$output_index
  Z <- opt_cl

  # samples
  mu_star_1_J_output <- post_result$mu_star_1_J_output
  alpha_phi_2_output <- post_result$alpha_phi_2_output
  Beta_output <- post_result$Beta_output
  b_output <- post_result$b_output


  # Generate some random index of length number_rep
  index <- sample(1:L, number_rep, replace = FALSE)

  seeds <- sample(1:1e8, number_rep)
  ##----------------------- Compute discrepancy measures ----------------------------------

  rep_Y_statistics <- pblapply(1:number_rep,
                                 FUN = function(i){

                                   set.seed(seeds[i])

                                   ii <- index[i]

                                   # Set theta
                                   theta <- list('mu' = mu_star_1_J_output[[ii]],
                                                 'b' = b_output[[ii]],
                                                 'alpha_phi_2' = alpha_phi_2_output[ii],
                                                 'Beta' = Beta_output[[ii]])

                                   # Simulate phi
                                   if(length(theta$b) == 3){

                                     phi <- matrix(rlnorm(n = J*G,
                                                          meanlog = theta$b[1]+theta$b[2]*log(theta$mu)+theta$b[3]*log(theta$mu)^2,
                                                          sdlog = sqrt(theta$alpha_phi_2)),
                                                   nrow = J,
                                                   ncol = G)
                                   }else{

                                     phi <- matrix(rlnorm(n = J*G,
                                                          meanlog = theta$b[1]+theta$b[2]*log(theta$mu),
                                                          sdlog = sqrt(theta$alpha_phi_2)),
                                                   nrow = J,
                                                   ncol = G)
                                   }

                                   # Replicated Y
                                   Y_rep <- lapply(1:D,
                                                   function(d){

                                                     matrix(rnbinom(n = G*C[d],
                                                                    mu = as.vector(t(theta$mu[Z[[d]],]))*rep(theta$Beta[[d]], each = G),
                                                                    size = as.vector(t(phi[Z[[d]],]))),

                                                            nrow = G,
                                                            ncol = C[d])
                                                   })

                                   # Statistics for replicated Y
                                   Y_rep_statistics_t <- lapply(1:D, function(d){

                                     data.frame(dataset = d,
                                                gene = 1:G,
                                                mean.log.shifted.counts = apply(Y_rep[[d]], 1, function(x) mean(log(x+1))),
                                                sd.log.shifted.counts = apply(Y_rep[[d]], 1, function(x) sd(log(x+1))),
                                                log.mean.counts = apply(Y_rep[[d]], 1, function(x) log(mean(x))),
                                                dropout.probability = apply(Y_rep[[d]], 1, function(x) length(which(x == 0))/C[d]))
                                   })

                                   # Combine all data
                                   Y_rep_statistics_t <- do.call(rbind,
                                                                 Y_rep_statistics_t)
                                   Y_rep_statistics_t$ind <- i

                                   return(Y_rep_statistics_t)
                                 }, cl = mc.cores)


  # Combine all replicated data
  rep_Y_statistics <- do.call(rbind,
                              rep_Y_statistics)

  ##-- For observed Y
  Y_statistics <- lapply(1:D, function(d){

    rel.df.d <- data.frame(dataset = d,
                           gene = 1:G,
                           mean.log.shifted.counts = apply(Y[[d]],
                                                           1,
                                                           function(x) mean(log(x+1))),

                           sd.log.shifted.counts = apply(Y[[d]],
                                                         1,
                                                         function(x) sd(log(x+1))),

                           log.mean.counts = apply(Y[[d]],
                                                   1,
                                                   function(x) log(mean(x))),

                           dropout.probability = apply(Y[[d]],
                                                       1,
                                                       function(x) length(which(x == 0))/C[d]))
  })

  Y_statistics <- do.call(rbind,
                          Y_statistics)
  Y_statistics$ind <- 0

  ##-- Return final output
  return(list('rep_Y_statistics' = rep_Y_statistics,
              'Y_statistics' = Y_statistics,
              'number_rep' = number_rep))
}


#' Plot the results for posterior predictive checks from multiple replicates
#'
#' @description
#' The function compares the kernel density estimates of three statistics: mean and mean and standard
#' deviation of log shifted counts and dropout probabilities, between the observed (red) and replicated (grey) data.
#'
#' @param ppc_multiple_df output from \code{ppc_multiple}.
#' @param data_names optional. The labels for each dataset in the plot.
#'
#' @return two plots for each dataset, comparing the density plots for each statistics between the
#' observed and replicated datasets.
#' @export
#'
#' @examples
#' plot_ppc_multiple(ppc_multiple_df = ppc_multiple_df, data_names = c('data 1', 'data 2'))
plot_ppc_multiple <- function(ppc_multiple_df, data_names=NULL){

  # Statistics
  Y_statistics <- ppc_multiple_df$Y_statistics
  rep_Y_statistics <- ppc_multiple_df$rep_Y_statistics


  Y_statistics$ind <- 'observed data'
  rep_Y_statistics$ind <- paste('replicated data', rep_Y_statistics$ind)

  D <- max(Y_statistics$dataset)

  if(is.null(data_names)){
    data_names <- paste0('data ',1:D)
  }

  # Kernel density plot
  df <- rbind(rep_Y_statistics, Y_statistics)

  colnames(df)[7] <- 'source'

  for(d in 1:D){

    df_d <- dplyr::filter(df, dataset == d)

    plot1 <- ggplot2::ggplot()+
      geom_density(mapping = aes(x = mean.log.shifted.counts,
                                 colour = source),

                   data = df_d,
                   size = 1.2)+
      theme_bw()+
      xlab('mean of log shifted counts')+
      scale_color_manual(values=c(rep('grey',ppc_multiple_df$number_rep), 'red'))+
      theme(legend.position = "none")+
      ggtitle(data_names[d])


    plot2 <- ggplot()+
      geom_density(mapping = aes(x = sd.log.shifted.counts,
                                 colour = source),

                   data = df_d,
                   size = 1.2)+
      theme_bw()+
      xlab('standard deviaion of log shifted counts')+
      scale_color_manual(values=c(rep('grey',ppc_multiple_df$number_rep), 'red'))+
      theme(legend.position = "none")+
      ggtitle(data_names[d])


    plot3 <- ggplot()+
      geom_density(mapping = aes(x = dropout.probability,
                                 colour = source),

                   data = df_d,
                   size = 1.2)+
      theme_bw()+
      xlab('dropout probabilities')+
      scale_color_manual(values=c(rep('grey',ppc_multiple_df$number_rep), 'red'))+
      theme(legend.position = "none")+
      ggtitle(data_names[d])

    gridExtra::grid.arrange(grobs=list(plot1,plot2,plot3), nrow = 1)
  }


}
