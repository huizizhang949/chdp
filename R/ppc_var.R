#' Title
#'
#' @param post_result
#' @param Y
#' @param t
#' @param opt_cl
#' @param number_rep
#' @param prob
#' @param mc.cores
#'
#' @return
#' @export
#'
#' @examples
ppc_var <- function(post_result, Y, t, opt_cl, number_rep, prob=c(0.005,0.995), mc.cores=8){

  L <- post_result$output_index
  G <- ncol(Y[[1]])
  Z <- unlist(opt_cl)
  D <- length(Y)

  # design matrix
  X <- lapply(1:D, function(d) {

    temp <- as.matrix(cbind(rep(1, nrow(Y[[d]])-1), Y[[d]][-nrow(Y[[d]]),]))

  })

  X <- do.call(rbind,X)

  # for each dataset, save the data, t, optimal cluster, dataset indicator
  Y_df <- lapply(1:D, function(d) {

    df_d <- data.frame(Y[[d]][-1,]); df_d$t <- t[[d]][-1]
    colnames(df_d)[1:G] <- paste0('y',1:G)
    df_d$dataset <- d

    return(df_d)
  })

  Y_df <- do.call(rbind,Y_df)

  # generate y_{rep,t} conditional on y_{obs,t-1}
  index <- sample(1:L,size=number_rep)
  # simulate the data Y.rep.list from number_rep MCMC samples
  Y.rep.list <- pbapply::pblapply(index, cl=mc.cores, function(i) {

    L.s <- post_result$L_1_J_output[[i]]
    Sigma.s <- post_result$Sigma_1_J_output[[i]]
    # a matrix of n_obs * G of simulated Y from sample i
    mat <- t(vapply(1:nrow(Y_df),function(cc) {

      mu <- matrix(X[cc,]%*%L.s[[Z[cc]]],nrow=1)

      y <- mvnfast::rmvn(1, mu=mu, sigma = Sigma.s[[Z[cc]]])

      return(y)

    },FUN.VALUE = numeric(G)))

    mat <- as.data.frame(mat)
    colnames(mat) <- paste0('y',1:G)
    mat[,(G+1):(G+2)] <- Y_df[,c('t','dataset')]

    return(mat)
  })

  # ----------- plot CI and posterior mean for Y.rep -----------
  Y.rep.mean <- as.data.frame(Reduce('+',Y.rep.list)/number_rep)
  colnames(Y.rep.mean)[1:G] <- paste0('y',1:G)
  Y.rep.mean[,(G+1):(G+2)] <- Y_df[,c('t','dataset')]

  # an array of dimension n_obs * G * number_rep
  arr <- array(unlist(lapply(Y.rep.list,function(x) x[,1:G])), c(nrow(Y_df),G,number_rep))
  quant.rep <- apply(arr,1:2,quantile,prob=prob)
  # turn into a list of length=G
  quant.rep.list <- lapply(1:G, function(p) {
    df <- as.data.frame(t(quant.rep[,,p]))
    colnames(df) <- c('lower','upper')
    df[,3:4] <- Y_df[,c('t','dataset')]

    return(df)
  })

  return(list(Y.rep.mean=Y.rep.mean, quantiles=quant.rep.list))
}


#' Title
#'
#' @param ppc_var_output
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
plot_ppc_var <- function(ppc_var_output, Y){

  D <- length(Y)
  G <- ncol(Y[[1]])
  Y.rep.mean <- ppc_var_output$Y.rep.mean
  quant.rep.list <- ppc_var_output$quantiles
  t <- Y.rep.mean$t

  Y_df <- lapply(1:D, function(d) {

    df_d <- data.frame(Y[[d]][-1,]);
    colnames(df_d)[1:G] <- paste0('y',1:G)
    df_d$dataset <- d

    return(df_d)
  })

  Y_df <- do.call(rbind,Y_df); Y_df$t <- t

  for(d in 1:D) {

    ix <- which(Y_df$dataset==d)

    ggs_d <- lapply(1:G,function(p) {

      temp <- quant.rep.list[[p]][ix,]

      ggplot(data=Y_df[ix,])+
        geom_point(aes(x=t,y=Y_df[ix,p]),col='black',size=0.3)+
        geom_line(data=Y.rep.mean[ix,],aes(x=t,y=Y.rep.mean[ix,p]),col='red',lwd=0.3)+
        geom_ribbon(data=temp,aes(x=t,ymin=lower,ymax=upper),fill='red',alpha=0.3)+
        theme_bw()+
        labs(y=paste0('y',p))

    })

    gridExtra::grid.arrange(grobs=ggs_d,nrow=G)

  }
}

