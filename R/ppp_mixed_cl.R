Discrepancy_measures_cl <- function(Y,
                                    theta){

  ##-- Dimensions
  D <- length(Y)
  J <- nrow(theta$mu)
  G <- ncol(theta$mu)

  df <- lapply(1:D,function(d) {

    # only consider non-empty clusters
    Js <- sort(unique(theta$Z[[d]]),decreasing = FALSE)

    temp <- lapply(Js, function(j){

      sub_ind <- which(theta$Z[[d]]==j)
      Y_sub <- Y[[d]][,sub_ind]

      ##-- Indicator function to distinguish between zero and non-zero counts
      decision1 <- ifelse(Y_sub==0, 1, 0)

      ##-- Expectations
      Expectations <- t(apply(theta$mu[theta$Z[[d]][sub_ind],1:G], 2, function(x) x*theta$Beta[[d]][sub_ind]))

      ##-- D1
      D1 <- rowSums((Y_sub-Expectations)^2/Expectations)
      D2 <- rowSums((sqrt(Y_sub)-sqrt(Expectations))^2)

      ##-- Probability
      probabilty <- dnbinom(x = 0,
                            mu = Expectations,
                            size = t(theta$phi[theta$Z[[d]][sub_ind],1:G]))

      ##-- D3
      D3 <- rowSums(abs(decision1-probabilty))

      ##-- Return Discrepancy measures
      return(data.frame('dataset' = d,
                        'cluster' = j,
                        'gene' = 1:G,
                        'D1' = D1,
                        'D2' = D2,
                        'D3' = D3))
    })

    return(do.call(rbind,temp))

  })

  return(do.call(rbind,df))
}


#' Compute posterior predictive p-values based on mixed predictive distribution and conditioned
#' on optimal clustering
#'
#' @description
#' This function computes posterior predictive p-values based on multiple replicates and optimal clustering.
#' Three discrepancy measures are calculated for each cluster,
#' based on Chi-squared statistic, Freeman-Tukey statistic
#' and dropout probabilities.
#'
#'
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#' @param Y a list of matrices. Each matrix is a gene-by-cell matrix of mRNA counts corresponding to a dataset.
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param number_rep number of replicates to generate.
#' @param mc.cores number of cores used for parallel computing.
#'
#' @return a data frame of calculated p-values for each discrepancy measure, conditioned on
#' clusters and datasets.
#' @export
#'
#' @examples
#' ppp_mixed_cl_result <- ppp_mixed_cl(post_result = post_result, Y = list(t(y1), t(y2)),
#'                                     opt_cl = opt_cl, number_rep = 200, mc.cores = 2)

ppp_mixed_cl <- function(post_result, Y, opt_cl, number_rep, mc.cores = 8){

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

  alpha.phi_2.sample <- sapply(index,
                              function(t) post_result$alpha_phi_2_output[t])

  b.sample <- lapply(index,
                     function(t) post_result$b_output[[t]])

  # sample of phi

  if(length(b.sample[[1]]) == 3){

    phi.sample <- lapply(1:number_rep,
                         function(i){

                           matrix(rlnorm(n = J*G,
                                         meanlog = b.sample[[i]][1]+b.sample[[i]][2]*log(mu.sample[[i]])+b.sample[[i]][3]*log(mu.sample[[i]])^2,
                                         sdlog = sqrt(alpha.phi_2.sample[i])),
                                  nrow = J,
                                  ncol = G)
                         })

  }else{

    phi.sample <- lapply(1:number_rep,
                         function(i){

                           matrix(rlnorm(n = J*G,
                                         meanlog = b.sample[[i]][1]+b.sample[[i]][2]*log(mu.sample[[i]]),
                                         sdlog = sqrt(alpha.phi_2.sample[i])),
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


                          df <- Discrepancy_measures_cl(Y = dat,
                                                     theta = list('Z' = opt_cl,
                                                                  'mu' = mu.sample[[i]],
                                                                  'phi' = phi.sample[[i]],
                                                                  'Beta' = beta.sample[[i]])
                          )

                          return(df)

                        }, cl=mc.cores)

  y.rep.Dis <- do.call(rbind, y.rep.Dis)

  print('Compute discrepancy measures for observed data')

  # Discrepancy for the observed data
  y.Dis <- pblapply(1:number_rep,
                    FUN = function(i){

                      df <- Discrepancy_measures_cl(Y = Y,
                                                 theta = list('Z' = opt_cl,
                                                              'mu' = mu.sample[[i]],
                                                              'phi' = phi.sample[[i]],
                                                              'Beta' = beta.sample[[i]])
                      )

                      return(df)

                    }, cl=mc.cores)

  y.Dis <- do.call(rbind, y.Dis)

  y.rep.Dis <- data.table::data.table(y.rep.Dis)
  y.Dis <- data.table::data.table(y.Dis)
  #--------------------------------- Compare statistics -----------------------------
  print('Compute p-values')

  p_value <- NULL

  for(d in 1:D){

    Js <- sort(unique(opt_cl[[d]]),decreasing = FALSE)

    p_value[[d]] <- lapply(1:G,
                           function(g){

                            temp <- lapply(Js, function(j){
                              # For y rep
                              p_value_rep_dgj <- y.rep.Dis[y.rep.Dis$dataset == d & y.rep.Dis$gene == g & y.rep.Dis$cluster==j,]

                              # For y obs
                              p_value_obs_dgj <- y.Dis[y.Dis$dataset == d & y.Dis$gene == g & y.Dis$cluster==j,]

                              # Output a data frame
                              return(data.frame(dataset = d,
                                                gene = g,
                                                cluster = j,
                                                Dis.1 = mean(p_value_rep_dgj$D1 >= p_value_obs_dgj$D1),
                                                Dis.2 = mean(p_value_rep_dgj$D2 >= p_value_obs_dgj$D2),
                                                Dis.3 = mean(p_value_rep_dgj$D3 >= p_value_obs_dgj$D3)))

                            })

                            return(do.call(rbind, temp))

                           })

    p_value[[d]] <- do.call(rbind,
                            p_value[[d]])

  }

  return(do.call(rbind,
                 p_value))
}


#' Histograms of posterior predictive p-values from multiple replicates, conditioned on optimal clustering
#'
#' @description
#' This function plots histograms of p-values for each discrepancy measure and each cluster.
#'
#' @importFrom graphics hist
#'
#' @param ppp_output output from \code{ppp_mixed_cl}.
#' @param nrow parameter used to cut the plot window into subpanels.
#' @param data_names optional. The labels for each dataset in the plot.
#'
#' @return three plots for each discrepancy measures, containing subpanels for each cluster.
#' @export
#'
#' @examples
#' ppp_hist_cl(ppp_output = ppp_mix_cl_result, nrow = 1, data_names = c('data 1', 'data 2'))
ppp_hist_cl <- function(ppp_output, nrow=1, data_names=NULL){

  D <- max(ppp_output$dataset)

  if(is.null(data_names)){
    data_names <- paste0('data ',1:D)
  }


  ##-- Distribution of p-values for posterior predictive checks
  ##-- Chi-square
  for(d in 1:D){

    # only consider non-empty clusters
    Js <- sort(unique(ppp_output$cluster[ppp_output$dataset==d]),decreasing = FALSE)

    par(mfrow=c(nrow,floor(length(Js/nrow))))

    for(j in Js){
      hist(ppp_output$Dis.1[ppp_output$cluster==j & ppp_output$dataset==d],

           main=paste0("Chi-Square\n", data_names[d],": cluster ", j),
           xlab="p-values for each gene",
           xlim=c(0,1),
           freq=FALSE,
           cex.lab = 1.5)

    }
  }

  ##-- Freeman-Tukey
  for(d in 1:D){

    # only consider non-empty clusters
    Js <- sort(unique(ppp_output$cluster[ppp_output$dataset==d]),decreasing = FALSE)

    par(mfrow=c(nrow,floor(length(Js/nrow))))

    for(j in Js){
      hist(ppp_output$Dis.2[ppp_output$cluster==j & ppp_output$dataset==d],

           main=paste0("Freeman-Tukey\n", data_names[d],": cluster ", j),
           xlab="p-values for each gene",
           xlim=c(0,1),
           freq=FALSE,
           cex.lab = 1.5)

    }
  }

  ##---Dropout probabilities
  for(d in 1:D){

    # only consider non-empty clusters
    Js <- sort(unique(ppp_output$cluster[ppp_output$dataset==d]),decreasing = FALSE)

    par(mfrow=c(nrow,floor(length(Js/nrow))))

    for(j in Js){
      hist(ppp_output$Dis.3[ppp_output$cluster==j & ppp_output$dataset==d],

           main=paste0("Dropout probabilities\n", data_names[d],": cluster ", j),
           xlab="p-values for each gene",
           xlim=c(0,1),
           freq=FALSE,
           cex.lab = 1.5)

    }
  }

}

