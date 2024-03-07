#' Title
#'
#' @param post_result
#' @param Y
#' @param t_pred
#' @param prob
#' @param mc.cores
#'
#' @return
#' @export
#'
#' @examples
pred_trend <- function(post_result, Y, t_pred, prob=c(0.005,0.995), mc.cores=8){

  L <- post_result$output_index
  D <- length(Y)
  J <- ncol(post_result$P_C_J_D_output[[1]][[1]])
  G <- ncol(Y[[1]])

  # kernel functiom
  pkernel <- function(t,mu,lambda,sigma_2){

    val1 <- -2/sigma_2
    val2 <- sin((t-mu)/lambda)

    return(exp(val1*val2^2))
  }

  Y.pred.list <- pbapply::pblapply(1:L, cl=mc.cores, function(i) {

    # kernel params
    mu.s <- post_result$mu_J_D_output[[i]]
    lambda.s <- post_result$lambda_J_D_output[[i]]
    sigma_2.s <- post_result$sigma_2_J_D_output[[i]]
    q.s <- post_result$Q_J_D_output[[i]]

    # p(t) at future t_pred, a list of length=D
    p_d <- lapply(1:D,function(d) {

      # a n_future_obs x J matrix
      t(sapply(t_pred, function(ts) {
        q.s[,d]*pkernel(ts,mu.s[,d],lambda.s[,d],sigma_2.s[,d])/(sum(q.s[,d]*pkernel(ts,mu.s[,d],lambda.s[,d],sigma_2.s[,d])))
      }))

    })

    # cluster membership for each new obs
    Z_d <- lapply(1:D,function(d) {
      sapply(1:length(t_pred),function(cc) sample(1:J,1,prob=p_d[[d]][cc,]))
    })

    # cluster-specific params
    L.s <- post_result$L_1_J_output[[i]]
    Sigma.s <- post_result$Sigma_1_J_output[[i]]

    y.pred <- lapply(1:D, function(d) {
      mat <- matrix(NA,nrow=length(t_pred),ncol=G)
      # first future observation is based on the last observed value
      x <- matrix(c(1, as.numeric(Y[[d]][nrow(Y[[d]]),])),nrow=1)
      mat[1,] <- mvtnorm::rmvnorm(1, x%*%L.s[[Z_d[[d]][1]]],sigma = Sigma.s[[Z_d[[d]][1]]])

      return(mat)
    })
    for(cc in 2:length(t_pred)){

      preds <- lapply(1:D, function(d) {
        x <- matrix(c(1,y.pred[[d]][cc-1,]),nrow=1)
        val <- mvtnorm::rmvnorm(1,x%*%L.s[[Z_d[[d]][cc]]],sigma = Sigma.s[[Z_d[[d]][cc]]])

        return(val)
      })

      for(d in 1:D){
        y.pred[[d]][cc,] <- preds[[d]]
      }
    }

    y.pred <- do.call(rbind,y.pred)
    y.pred <- as.data.frame(y.pred)
    colnames(y.pred) <- paste0('y',1:G)
    y.pred$t <- rep(t_pred,D)
    y.pred$dataset <- rep(1:D,each=length(t_pred))

    return(y.pred)
  })

  # posterior mean of predictions
  Y.pred.mean <- as.data.frame(Reduce('+',Y.pred.list)/L)
  colnames(Y.pred.mean)[1:G] <- paste0('y',1:G)
  Y.pred.mean[,(G+1):(G+2)] <- Y.pred.list[[1]][,c('t','dataset')]

  # an array of dimension n_future_obs * G * n_sample
  arr <- array(unlist(lapply(Y.pred.list,function(x) x[,1:G])), c(D*length(t_pred),G,L))
  quant.pred <- apply(arr,1:2,quantile,prob=prob)
  # turn into a list of length=G
  quant.pred.list <- lapply(1:G, function(p) {
    df <- as.data.frame(t(quant.pred[,,p]))
    colnames(df) <- c('lower','upper')
    df[,3:4] <- Y.pred.mean[,c('t','dataset')]

    return(df)
  })

  return(list(t_pred=t_pred, Y.pred.mean=Y.pred.mean, quantiles=quant.pred.list))

}


#' Title
#'
#' @param pred_trend_result
#' @param Y
#' @param t_obs
#'
#' @return
#' @export
#'
#' @examples
plot_pred_trend <- function(pred_trend_output, Y, t_obs){

  t_pred <- pred_trend_output$t_pred
  Y.pred.mean <- pred_trend_output$Y.pred.mean
  quant.pred.list <- pred_trend_output$quantiles
  D <- length(Y)
  G <- ncol(Y[[1]])

  # for each dataset, save the data, t, dataset indicator
  Y_df <- lapply(1:D, function(d) {

    df_d <- data.frame(Y[[d]][-1,]); df_d$t <- t_obs[[d]][-1]
    colnames(df_d)[1:G] <- paste0('y',1:G)
    df_d$dataset <- d

    return(df_d)
  })

  Y_df <- do.call(rbind,Y_df)

  for(d in 1:D) {

    ix <- which(Y.pred.mean$dataset==d)

    ggs_d <- lapply(1:G,function(p) {

      temp <- quant.pred.list[[p]][ix,]

      ggplot(data=Y_df[Y_df$d==d,])+
        geom_point(aes(x=t,y=Y_df[Y_df$d==d,p]),col='black',size=0.3)+
        geom_line(data=Y.pred.mean[ix,],aes(x=t,y=Y.pred.mean[ix,p]),col='red',lwd=0.3)+
        geom_ribbon(data=temp,aes(x=t,ymin=lower,ymax=upper),fill='red',alpha=0.3)+
        theme_bw()+
        labs(y=paste0('y',p))

    })

    gridExtra::grid.arrange(grobs=ggs_d,nrow=G)

  }

}




#' Title
#'
#' @param post_result
#' @param t_pred
#' @param mc.cores
#'
#' @return
#' @export
#'
#' @examples
pred_c_prob <- function(post_result, t_pred, mc.cores=8){

  L <- post_result$output_index
  D <- length(post_result$P_C_J_D_output)

  # kernel function
  pkernel <- function(t,mu,lambda,sigma_2){

    val1 <- -2/sigma_2
    val2 <- sin((t-mu)/lambda)

    return(exp(val1*val2^2))
  }

  # a list of length=D, each list corresponds to p(t) for each dataset, within each list:
  # a list of length=n_sample, each is a matrix of n_future_time_points * J
  pt_predict <- lapply(1:D,function(d) {

    print(paste('Compute for data',d))

    temp <- pbapply::pblapply(1:L, cl=mc.cores,function(i) {

      mu.s <- post_result$mu_J_D_output[[i]]
      lambda.s <- post_result$lambda_J_D_output[[i]]
      sigma_2.s <- post_result$sigma_2_J_D_output[[i]]
      q.s <- post_result$Q_J_D_output[[i]]

      val <- lapply(t_pred, function(ts) {
        q.s[,d]*pkernel(ts,mu.s[,d],lambda.s[,d],sigma_2.s[,d])/(sum(q.s[,d]*pkernel(ts,mu.s[,d],lambda.s[,d],sigma_2.s[,d])))
      })

      val <- do.call(rbind,val)

      return(val)
    })

    return(temp)

  })

  return(list(t_pred=t_pred,pt_predict=pt_predict,pt_obs=post_result$P_C_J_D_output))
}


#' Title
#'
#' @param pred_c_prob_result
#' @param data
#' @param t_obs
#' @param cluster
#' @param thinning
#' @param color_pal
#'
#' @return
#' @export
#'
#' @examples
plot_pred_c_prob <- function(pred_c_prob_output, data, t_obs, cluster, thinning=NULL, color_pal=NULL){

  pt_predict <- pred_c_prob_output$pt_predict
  t_pred <- pred_c_prob_output$t_pred
  L <- length(pt_predict[[1]])
  pt_obs <- pred_c_prob_output$pt_obs

  d=data;j=cluster
  if(is.null(thinning)){
    inds <- seq(thinning, L, by=thinning)
  }else{
    inds <- 1:L
  }

  plot(1,1,xlim=c(0,max(t_pred)),ylim=c(0,1),type='n',xlab='t',ylab='p(t)',main=paste0('Dataset ', d,': cluster ',j))
  # prediction
  for(i in inds){
    lines(t_pred,pt_predict[[d]][[i]][,j],col='grey')
  }
  # observed
  for (i in inds) {
    lines(t_obs[-1],pt_obs[[d]][[i]][,j],col=ifelse(is.null(color_pal),'red',color_pal[j]))
  }

}


