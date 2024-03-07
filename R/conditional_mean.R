#' Title
#'
#' @param Y
#' @param t
#' @param opt_cl
#' @param post_result
#' @param prob
#'
#' @return
#' @export
#'
#' @examples
conditional_mean <- function(Y, t, opt_cl, post_result, prob=c(0.005,0.995), mc.cores=8){

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
  quant <- apply(arr,1:2,quantile,prob=prob)
  # turn into a list of length = G
  quant_list <- lapply(1:G, function(p) {
    df_quant <- as.data.frame(t(quant[,,p]))
    colnames(df_quant) <- c('lower','upper')
    df_quant[,3:5] <- Y_df[,c('cluster','t','dataset')]

    return(df_quant)
  })

  return(list(Y_df=Y_df, Y_conditional_mean=Y_mean, quantiles=quant_list))
}

#' Title
#'
#' @param conditional_mean_result
#' @param data
#' @param cluster
#'
#' @return
#' @export
#'
#' @examples
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

