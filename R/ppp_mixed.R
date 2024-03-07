Discrepancy_measures <- function(Y,
                                 theta){

  ##-- Dimensions
  D <- length(Y)
  J <- nrow(theta$mu)
  G <- ncol(theta$mu)

  df <- lapply(1:D,function(d) {
    ##-- Indicator function to distinguish between zero and non-zero counts
    decision1 <- ifelse(Y[[d]]==0, 1, 0)
    ##-- Expectations
    Expectations <- t(apply(theta$mu[theta$Z[[d]],1:G], 2, function(x) x*theta$Beta[[d]]))

    ##-- D1
    D1 <- rowSums((Y[[d]]-Expectations)^2/Expectations)
    D2 <- rowSums((sqrt(Y[[d]])-sqrt(Expectations))^2)

    ##-- Probability
    probabilty <- dnbinom(x = 0,
                          mu = Expectations,
                          size = t(theta$phi[theta$Z[[d]],1:G]))

    ##-- D3
    D3 <- rowSums(abs(decision1-probabilty))

    ##-- Return Discrepancy measures
    return(data.frame('dataset' = d,
                      'gene' = 1:G,
                      'D1' = D1,
                      'D2' = D2,
                      'D3' = D3))
  })

  return(do.call(rbind,df))
}


#' Compute posterior predictive p-values based on mixed predictive distribution
#'
#' @description
#' This function computes posterior predictive p-values based on multiple replicates.
#' Three discrepancy measures are calculated based on Chi-squared statistic, Freeman-Tukey statistic
#' and dropout probabilities, for each gene in each dataset.
#'
#' @importFrom magrittr %>%
#'
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#' @param Y a list of matrices. Each matrix is a gene-by-cell matrix of mRNA counts corresponding to a dataset.
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param number_rep number of replicates to generate.
#' @param mc.cores number of cores used for parallel computing.
#'
#' @return a data frame of calculated p-values for each discrepancy measure.
#' @export
#'
#' @examples
#' ppp_mixed_result <- ppp_mixed(post_result = post_result, Y = list(t(y1), t(y2)),
#'                               opt_cl = opt_cl, number_rep = 200, mc.cores = 2)
ppp_mixed <- function(post_result, Y, opt_cl, number_rep, mc.cores = 8){

  ##-- Input data
  C <- sapply(opt_cl, function(x) length(x))
  J <- length(unique(unlist(opt_cl)))
  G <- nrow(Y[[1]])
  D <- length(Y)
  L <- post_result$output_index

  ##-- Generate some random index of length number_rep
  index <- sample(1:L, number_rep, replace = FALSE)
  seeds <- sample(1:1e8, number_rep)

  # Samples
  beta.sample <- lapply(index,
                        function(t) post_result$Beta_output[[t]])

  mu.sample <- lapply(index,
                      function(t) post_result$mu_star_1_J_output[[t]])

  alpha.phi2.sample <- sapply(index,
                              function(t) post_result$alpha_phi2_output[t])

  b.sample <- lapply(index,
                     function(t) post_result$b_output[[t]])

  # sample of phi

  if(length(b.sample[[1]]) == 3){

    phi.sample <- lapply(1:number_rep,
                         function(i){

                           matrix(rlnorm(n = J*G,
                                         meanlog = b.sample[[i]][1]+b.sample[[i]][2]*log(mu.sample[[i]])+b.sample[[i]][3]*log(mu.sample[[i]])^2,
                                         sdlog = sqrt(alpha.phi2.sample[i])),
                                  nrow = J,
                                  ncol = G)
                         })

  }else{

    phi.sample <- lapply(1:number_rep,
                         function(i){

                           matrix(rlnorm(n = J*G,
                                         meanlog = b.sample[[i]][1]+b.sample[[i]][2]*log(mu.sample[[i]]),
                                         sdlog = sqrt(alpha.phi2.sample[i])),
                                  nrow = J,
                                  ncol = G)

                         })
  }

  print('Generate replicates and compute discrepancy measures')
  # generate replicated Y and compute discrepancy measures
  y.rep.Dis <- pblapply(1:number_rep,
                  FUN = function(i){

                    set.seed(seeds[i])

                    dat <- lapply(1:D,
                                  function(d){

                                    matrix(rnbinom(n = G*C[d],
                                                   mu = as.vector(t(mu.sample[[i]][opt_cl[[d]],]))*rep(beta.sample[[i]][[d]], each = G),
                                                   size = as.vector(t(phi.sample[[i]][opt_cl[[d]],]))),

                                           nrow = G,
                                           ncol = C[d])
                                  })


                    df <- Discrepancy_measures(Y = dat,
                                               theta = list('Z' = opt_cl,
                                                            'mu' = mu.sample[[i]],
                                                            'phi' = phi.sample[[i]],
                                                            'Beta' = beta.sample[[i]])
                    )

                    return(df)

                  }, cl = mc.cores)

  y.rep.Dis <- do.call(rbind, y.rep.Dis)

  print('Compute discrepancy measures for observed data')

  # Discrepancy for the observed data
  y.Dis <- pblapply(1:number_rep,
                  FUN = function(i){

                    df <- Discrepancy_measures(Y = Y,
                                               theta = list('Z' = opt_cl,
                                                            'mu' = mu.sample[[i]],
                                                            'phi' = phi.sample[[i]],
                                                            'Beta' = beta.sample[[i]])
                    )

                    return(df)

                  }, cl = mc.cores)

  y.Dis <- do.call(rbind, y.Dis)


  #--------------------------------- Compare statistics -----------------------------
  print('Compute p-values')

  p_value <- NULL

  for(d in 1:D){

    p_value[[d]] <- lapply(1:G,
                           function(g){


                             # For y rep
                             p_value_rep_dg <- y.rep.Dis %>%

                               dplyr::filter(dataset == d,
                                      gene == g)

                             # For y obs
                             p_value_obs_dg <- y.Dis %>%

                               dplyr::filter(dataset == d,
                                      gene == g)

                             # Output a data frame
                             return(data.frame(dataset = d,
                                               gene = g,
                                               Dis.1 = mean(p_value_rep_dg$D1 >= p_value_obs_dg$D1),
                                               Dis.2 = mean(p_value_rep_dg$D2 >= p_value_obs_dg$D2),
                                               Dis.3 = mean(p_value_rep_dg$D3 >= p_value_obs_dg$D3)))

                           })

    p_value[[d]] <- do.call(rbind,
                            p_value[[d]])

  }

  return(do.call(rbind,
                 p_value))
}


#' Histograms of posterior predictive p-values from multiple replicates
#'
#' @description
#' This function plots histograms of p-values for each discrepancy measure.
#'
#' @importFrom graphics hist
#'
#' @param ppp_output output from \code{ppp_mixed}.
#'
#' @return three histograms of p-values for discrepancy measures based on Chi-squared statistic, Freeman-Tukey statistic
#' and dropout probabilities, respectively.
#' @export
#'
#' @examples
#' ppp_hist(ppp_output = ppp_mixed_result)
ppp_hist <- function(ppp_output){

  ##-- Distribution of p-values for posterior predictive checks

  par(mfrow=c(1,3))

  ##-- Chi-square
  hist(ppp_output$Dis.1,

       main="Chi-Square",
       xlab="p-values for each gene",
       xlim=c(0,1),
       freq=FALSE,
       cex.lab = 1.5)

  ##-- Freeman-Tukey
  hist(ppp_output$Dis.2,

       main="Freeman-Tukey",
       xlab="p-values for each gene",
       xlim=c(0,1),
       freq=FALSE,
       cex.lab = 1.5)

  ##-- Dropout probabilities
  hist(ppp_output$Dis.3,

       main="Dropout probabilties",
       xlab="p-values for each gene",
       xlim=c(0,1),
       freq=FALSE,
       cex.lab = 1.5)


}

