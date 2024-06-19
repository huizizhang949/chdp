#' Find local marker genes
#'
#' @description
#' Identify locally differentially expressed and dispersed genes.
#'
#' @importFrom pbapply pblapply
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#' @param threshold a vector of length 2, corresponding to the threshold in the calculation of tail probabilities.
#' @param alpha_m optional. The threshold to classify locally differentially expressed (DE) genes. If not provided, the value will be
#' chosen to achieve an expected false discovery rate (EFDR) of 0.05.
#' @param alpha_d optional. The threshold to classify locally differentially dispersed (DD) genes. If not provided, the value will be
#' chosen to achieve an expected false discovery rate (EFDR) of 0.05.
#' @param mc.cores number of cores used for parallel computing.
#'
#' @return The output contains the following items:
#' \item{mu.output}{a dataframe containing mean absolute log-fold change, tail probabilities and classification of local DE genes for each cluster.}
#' \item{phi.output}{a dataframe containing mean absolute log-fold change, tail probabilities and classification of local DD genes for eachc luster.}
#' \item{alpha_m}{threshold to classify local DE genes.}
#' \item{alpha_d}{threshold to classify local DD genes.}
#' @export
#'
#' @examples
#' local_result <- local_marker_genes(post_result = post_result,
#'   threshold = c(1.2,1.2), alpha_m = 0.3, alpha_d = 0.3, mc.cores = 2)
local_marker_genes <- function(post_result, threshold=c(2.5,2.5), alpha_m=NULL, alpha_d=NULL,
                               mc.cores=8){

  J <- ncol(post_result$P_C_J_D_output[[1]][[1]])
  G <- ncol(post_result$mu_star_1_J_output[[1]])
  L <- post_result$output_index

  # threholds
  tau <- threshold[1]; omega <- threshold[2]
  # mcmc ssamples
  mu_mcmc <- post_result$mu_star_1_J_output
  phi_mcmc <- post_result$phi_star_1_J_output

  pairs <- combn(1:J,2)

  # compute tail probabilities for each gene, for every pair
  loop_local_output <- lapply(1:G, function(g) {

    temp <- apply(pairs, 2, function(vec) {
      j <- vec[1]
      j1 <- vec[2]

      sample_value <- lapply(1:L,function(iter) {

        lfc.mu <- abs(log(mu_mcmc[[iter]][j,g])-log(mu_mcmc[[iter]][j1,g]))
        lfc.phi <- abs(log(phi_mcmc[[iter]][j,g])-log(phi_mcmc[[iter]][j1,g]))
        return(c(lfc.mu, lfc.phi))
      })

      sample_value <- do.call(rbind,sample_value)

      df.mu <- data.frame(gene=g, j=j, j1=j1, mean.lfc=mean(sample_value[,1]), tail.prob=mean(sample_value[,1]>tau))
      df.phi <- data.frame(gene=g, j=j, j1=j1, mean.lfc=mean(sample_value[,2]), tail.prob=mean(sample_value[,2]>omega))

      return(list(df.mu=df.mu, df.phi=df.phi))
    })

    # for one gene, all pairwise comparisons and tail_probabilities and mean log-fold change
    df.mu <- do.call(rbind,lapply(temp,function(l) l$df.mu))
    df.phi <- do.call(rbind,lapply(temp,function(l) l$df.phi))

    # for this gene: for each j, find the minimum tail probabilities, for mu and phi
    temp_j <- lapply(1:J, function(j) {
      # find the pairs that contain cluster j
      ind <- c(1:ncol(pairs))[apply(pairs, 2, function(vec) any(vec==j))]
      df.mu.sub <- df.mu[ind,]
      df.phi.sub <- df.phi[ind,]

      # find the minimum value
      min_val_mu <- min(df.mu.sub$tail.prob)
      min_val_phi <- min(df.phi.sub$tail.prob)

      # find the index for the minimum value, there may be multiple minima
      min_ind_mu <- which(df.mu.sub$tail.prob==min_val_mu)
      min_ind_phi <- which(df.phi.sub$tail.prob==min_val_phi)

      # compute average lfc and store minimum tail probabilities
      mu.output <- data.frame(gene=g, cluster=j, mean.lfc=mean(df.mu.sub$mean.lfc[min_ind_mu]), tail.prob=min_val_mu)
      phi.output <- data.frame(gene=g, cluster=j, mean.lfc=mean(df.phi.sub$mean.lfc[min_ind_phi]), tail.prob=min_val_phi)


      return(list(mu.output=mu.output,phi.output=phi.output))
    })

    # summary for this gene
    mu.output <- do.call(rbind,lapply(temp_j,function(l) l$mu.output))
    phi.output <- do.call(rbind,lapply(temp_j,function(l) l$phi.output))


    return(list(mu.output=mu.output, phi.output=phi.output))
  })


  # outputs for all genes
  mu.output <- do.call(rbind,lapply(loop_local_output, function(l) l$mu.output))
  phi.output <- do.call(rbind,lapply(loop_local_output, function(l) l$phi.output))

  # function to compute expected false discovery rate (efdr)
  efdr <- function(P_g_star=NULL,alpha_m=NULL,L_g_star=NULL,alpha_d=NULL){
    if(is.null(L_g_star)&is.null(alpha_d)){
      vec_star <- P_g_star
      alpha <- alpha_m
    }else{
      vec_star <- L_g_star
      alpha <- alpha_d
    }
    value <- sum((1-vec_star)*ifelse(vec_star>alpha,1,0))/(sum(1-vec_star))
    return(value)
  }

  # choose suitable alpha_m and alpha_d
  if(is.null(alpha_m)){
    alphas <- seq(0,1,length.out = 1000)
    efdr_local_mu <- sapply(alphas, function(val) {
      efdr(P_g_star = mu.output$tail.prob,alpha_m = val)
    })
    alpha_m <- alphas[which.min(abs(efdr_local_mu-0.05))]
  }

  if(is.null(alpha_d)){
    alphas <- seq(0,1,length.out = 1000)
    efdr_local_phi <- sapply(alphas, function(val) {
      efdr(L_g_star = phi.output$tail.prob,alpha_d = val)
    })
    alpha_d <- alphas[which.min(abs(efdr_local_phi-0.05))]
  }

  # classify genes
  mu.output$DE <- ifelse(mu.output$tail.prob>alpha_m,'Yes','No')
  phi.output$DD <- ifelse(phi.output$tail.prob>alpha_d,'Yes','No')

  return(list(mu.output=mu.output,phi.output=phi.output, alpha_m=alpha_m,alpha_d=alpha_d))
}





#' Plot the results for local marker genes
#'
#' @description
#' Plot mean absolute log-fold change (LFC), and the number of local genes for each cluster.
#'
#' @param local_output output from \code{local_marker_genes}.
#' @param nrow number of rows for displaying the plots of tail probabilities against mean absolute LFC
#'
#' @return the output contains four ggplot objects: 1. Tail probabilities against mean absolute LFC based on mean expression for each cluster,
#' and the threshold to decide local DE genes. 2. A bar-chart showing the number of local DE genes for each cluster.
#' 3. Tail probabilities against mean absolute LFC based on dispersion for each cluster,
#' and the threshold to decide local DD genes. 4. A bar-chart showing the number of local DD genes for each cluster.
#' @export
#'
#' @examples
#' ggs_local <- plot_local_marker_genes(local_output = local_result)
#' gridExtra::grid.arrange(grobs=ggs_local,nrow=2)
plot_local_marker_genes <- function(local_output, nrow=2){

  J <- max(local_output$mu.output$cluster)
  # ---- mean expression -----
  # plot tail probabilities against mean absolute LFC
  p1 <- ggplot(data=local_output$mu.output,aes(x=mean.lfc,y=tail.prob,col=DE))+
    geom_point(show.legend = F,size=0.5)+
    facet_wrap(~cluster,nrow=nrow)+
    theme_bw()+
    scale_color_manual(values = c('Yes'='red','No'='grey'))+
    geom_hline(yintercept = local_output$alpha_m,linetype='dashed')+
    labs(x='Mean absolute log-fold change',y='Tail probability',title = 'Mean expression')

  # count how many local DE genes for each cluster
  df_mu <- subset(local_output$mu.output,DE=='Yes')
  df_mu_summary <- data.frame(cluster=as.factor(1:J),size=sapply(1:J,function(j) sum(df_mu$cluster==j)))

  p2 <- ggplot(data=df_mu_summary,aes(x=cluster,y=size))+
    geom_bar(stat='identity',fill='steelblue')+
    theme_bw()+
    labs(y='Number of marker genes',title='')

  # ---- dispersion -------
  p3 <- ggplot(data=local_output$phi.output,aes(x=mean.lfc,y=tail.prob,col=DD))+
    geom_point(show.legend = F,size=0.5)+
    facet_wrap(~cluster,nrow=nrow)+
    theme_bw()+
    scale_color_manual(values = c('Yes'='red','No'='grey'))+
    geom_hline(yintercept = local_output$alpha_d,linetype='dashed')+
    labs(x='Mean absolute log-fold change',y='Tail probability',title='Dispersion')

  # count how many local DE genes for each cluster
  df_phi <- subset(local_output$phi.output,DD=='Yes')
  df_phi_summary <- data.frame(cluster=as.factor(1:J),size=sapply(1:J,function(j) sum(df_phi$cluster==j)))

  p4 <- ggplot(data=df_phi_summary,aes(x=cluster,y=size))+
    geom_bar(stat='identity',fill='steelblue')+
    theme_bw()+
    labs(y='Number of marker genes',title='')


  return(list(p1,p2,p3,p4))


}




#' Plot heatmaps for estimated mean expression and dispersion across clusters, for local marker genes
#'
#' @description
#' The function shows the heatmaps for estimated mean and dispersion parameters of local marker genes, for each cluster.
#'
#' @param local_output output from \code{global_marker_genes}.
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#'
#' @return The output contains two components:
#' \item{mu.ggs}{a list of ggplot objects. Each object is the heatmap of posterior mean of log mean expression,
#' for local DE genes of the cluster. Genes are ordered by decreasing tail probabilities from top to bottom.}
#' \item{phi.ggs}{a list of ggplot objects. Each object is the heatmap of posterior mean of log dispersion,
#' for local DD genes of the cluster. Genes are ordered by decreasing tail probabilities from top to bottom}
#' @export
#'
#' @examples
#' heatmaps_ggs_local <- local_marker_genes_heatmaps(local_output = local_result, post_result = post_result)
#' gridExtra::grid.arrange(grobs=heatmaps_ggs_local$mu.ggs,nrow=1)
#' gridExtra::grid.arrange(grobs=heatmaps_ggs_local$phi.ggs,nrow=1)
local_marker_genes_heatmaps <- function(local_output, post_result){

  J <- ncol(post_result$P_C_J_D_output[[1]][[1]])
  L <- post_result$output_index
  mu.sample <- post_result$mu_star_1_J_output
  phi.sample <- post_result$phi_star_1_J_output

  # threshold
  alpha_m <- local_output$alpha_m
  alpha_d <- local_output$alpha_d

  # ------ posterior mean of log(mu or phi) ----------
  # DE: compute posterior mean of logarithm of mu
  logmu.sample <- rapply(mu.sample, log, how = 'replace')
  # G x J
  mean_logmu <- t(Reduce('+',logmu.sample)/L)

  # plot
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  hp1 <- lapply(1:J, function(j) {

    # summaries for cluster j, and local DE genes only
    df <- dplyr::filter(local_output$mu.output, cluster==j & tail.prob > alpha_m)
    # order by tail probabilities
    df <- df[order(df$tail.prob, decreasing = TRUE),]

    heat_map <- pheatmap(mean_logmu[df$gene,],
                         scale = "none",cluster_rows  = FALSE,cluster_cols = FALSE,
                         angle_col = 90,fontsize = 8,border_color=NA,color = myPalette(400),
                         main = paste('cluster',j),silent = T)

    return(ggplotify::as.ggplot(heat_map))
  })


  # DD: compute posterior mean of logarithm of phi
  logphi.samples <- rapply(phi.sample, log, how = 'replace')
  # G x J
  mean_logphi <- t(Reduce('+',logphi.samples)/L)

  # plot
  hp2 <- lapply(1:J, function(j) {

    # summaries for cluster j, and local DE genes only
    df <- dplyr::filter(local_output$phi.output, cluster==j & tail.prob > alpha_d)
    # order by tail probabilities
    df <- df[order(df$tail.prob, decreasing = TRUE),]

    heat_map <- pheatmap(mean_logphi[df$gene,],
                         scale = "none",cluster_rows  = FALSE,cluster_cols = FALSE,
                         angle_col = 90,fontsize = 8,border_color=NA,color = myPalette(400),
                         main = paste('cluster',j),silent = T)

    return(ggplotify::as.ggplot(heat_map))
  })

  return(list(mu.ggs=hp1, phi.ggs=hp2))
}








