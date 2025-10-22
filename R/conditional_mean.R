#' Compute conditional mean in the VAR likelihood
#'
#' @description
#' This function calculates the mean conditional on the past observation as well as credible intervals
#'
#' @importFrom stats quantile
#' @param Y a list of two matrices for two datasets. The columns correspond to features.
#' @param t a list of two vectors. Each vector is the external covariate (time) for individual dataset.
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param post_result output from \code{pkernelHDP_mcmc} for the post-processing step.
#' @param prob probability corresponding to the credible interval.
#' @param mc.cores number of cores used for parallel computing.
#'
#' @return a list of three components:
#' \item{Y_df}{a dataframe combining observations, optimal clustering and time across datasets.}
#' \item{Y_conditional_mean}{a dataframe of posterior mean of the conditional mean.}
#' \item{quantiles}{a list of dataframes of credible intervals for the conditional mean. Each dataframe corresponds to one feature.}
#' @export
#'
#' @examples
#' conditional_mean_result <- conditional_mean(Y=list(Y1,Y2), t=list(t1,t2), opt_cl = opt_cl,
#'     post_result = post_result, prob = 0.95, mc.cores = 2)

conditional_mean <- function(Y, t, opt_cl, post_result, prob=0.95, mc.cores=8){

  D <- length(Y); G <- ncol(Y[[1]])
  L <- post_result$output_index
  Z <- unlist(opt_cl)

  # design matrix
  X <- lapply(1:D, function(d) {

    temp <- as.matrix(cbind(rep(1, nrow(Y[[d]])-1), Y[[d]][-nrow(Y[[d]]),]))

  })

  X <- do.call(rbind,X)

  # for each dataset, save the data, t, optimal cluster, dataset indicator
  Y_df <- lapply(1:D, function(d) {

    df_d <- data.frame(Y[[d]][-1,]); df_d$cluster <- as.factor(opt_cl[[d]]); df_d$t <- t[[d]][-1]
    colnames(df_d)[1:G] <- paste0('y',1:G)
    df_d$dataset <- d

    return(df_d)
  })

  Y_df <- do.call(rbind,Y_df)

  # all samples for conditional means
  Y_mean_list <- pbapply::pblapply(1:L, cl=mc.cores,function(i) {

    L.s <- post_result$L_1_J_output[[i]]
    # a matrix of n_obs * G of Y_mean from sample i.
    mat <- t(vapply(1:nrow(Y_df),function(cc) {

      val <- matrix(X[cc,]%*%L.s[[Z[cc]]],nrow=1)

      return(val)

    },FUN.VALUE = numeric(G)))

  })

  # posterior mean
  Y_mean <- as.data.frame(Reduce('+',Y_mean_list)/L)
  colnames(Y_mean) <- paste0('y',1:G)
  Y_mean[,(G+1):(G+3)] <- Y_df[,c('cluster','t','dataset')]

  # an array of dimension n_obs * G * n_sample
  arr <- array(unlist(Y_mean_list), c(nrow(Y_df),G,L))
  quant <- apply(arr,1:2,function(x) as.numeric(coda::HPDinterval(coda::mcmc(x),prob=prob)))
  # turn into a list of length = G
  quant_list <- lapply(1:G, function(p) {
    df_quant <- as.data.frame(t(quant[,,p]))
    colnames(df_quant) <- c('lower','upper')
    df_quant[,3:5] <- Y_df[,c('cluster','t','dataset')]

    return(df_quant)
  })

  return(list(Y_df=Y_df, Y_conditional_mean=Y_mean, quantiles=quant_list))
}

#' Plot conditional mean
#'
#' @description
#' This function plots the posterior mean and credible intervals for the conditional mean.
#'
#' @param conditional_mean_output output from \code{conditional_mean}.
#' @param data numeric value to indicate which dataset to visualize.
#' @param cluster numeric value to indicate which cluster to visualize.
#'
#' @return plots showing the posterior means (red) and credible intervals (black bar) against time, with
#' observations in grey line.
#' @export
#'
#' @examples
#' plot_conditional_mean(conditional_mean_output = conditional_mean_result, data = 1, cluster = 1)
plot_conditional_mean <- function(conditional_mean_output, data, cluster){

  Y_df <- conditional_mean_output$Y_df
  Y_mean <- conditional_mean_output$Y_conditional_mean
  quant_list <- conditional_mean_output$quantiles

  G <- length(quant_list)

  d=data;j=cluster
  ix <- which(Y_df$cl==j & Y_df$dataset==d)
  ixd <- which(Y_df$dataset==d)
  ggs <- lapply(1:G,function(p) {

    temp <- quant_list[[p]][ix,]

    ggplot(data=Y_df[ixd,])+
      geom_line(aes(x=t,y=Y_df[ixd,p]),col='grey')+
      geom_errorbar(data=temp,aes(x=t,ymin=lower,ymax=upper), width=0.01,linewidth=0.5)+
      geom_point(data=Y_mean[ix,],aes(x=t,y=Y_mean[ix,p]),col='red',size=0.3)+
      theme_bw()+
      labs(y=paste0('y',p))

  })

  gridExtra::grid.arrange(grobs=ggs,nrow=G)
}

