#' Find global marker genes
#'
#' @description
#' Identify globally differentially expressed and dispersed genes。
#'
#' @importFrom pbapply pblapply
#' @importFrom utils combn
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#' @param threshold a vector of length 2, corresponding to the threshold in the calculation of tail probabilities.
#' @param alpha_m optional. The threshold to classify globally differentially expressed (DE) genes. If not provided, the value will be
#' chosen to achieve an expected false discovery rate (EFDR) of 0.05.
#' @param alpha_d optional. The threshold to classify globally differentially dispersed (DD) genes. If not provided, the value will be
#' chosen to achieve an expected false discovery rate (EFDR) of 0.05.
#' @param mc.cores number of cores used for parallel computing.
#'
#' @return The output contains the following items:
#' \item{mu.output}{a dataframe containing mean absolute log-fold change, tail probabilities and classification of global DE genes.}
#' \item{phi.output}{a dataframe containing mean absolute log-fold change, tail probabilities and classification of global DD genes.}
#' \item{alpha_m}{threshold to classify global DE genes.}
#' \item{alpha_d}{threshold to classify global DD genes.}
#' @export
#'
#' @examples
#' global_result <- global_marker_genes(post_result = post_result,
#'   threshold = c(2.5,2.5), alpha_m = 0.3, alpha_d = 0.3, mc.cores = 2)
global_marker_genes <- function(post_result, threshold=c(2.5,2.5), alpha_m=NULL, alpha_d=NULL,
                                mc.cores=8){


  J <- ncol(post_result$P_C_J_D_output[[1]][[1]])
  G <- ncol(post_result$mu_star_1_J_output[[1]])
  L <- post_result$output_index

  # threholds
  tau <- threshold[1]; omega <- threshold[2]
  # mcmc samples
  mu_mcmc <- post_result$mu_star_1_J_output
  phi_mcmc <- post_result$phi_star_1_J_output

  # all combination of j, j'
  pairs <- combn(1:J,2)

  loop_global_output <- pblapply(1:G, cl=mc.cores, function(g) {

    # for each pair of clusters, compute mean log-fold change and tail probability
    temp <- lapply(1:ncol(pairs), function(i) {
      j <- pairs[1,i]
      j1 <- pairs[2,i]

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

    # find the maximum value (may have multiple maxima)
    ind <- which(df.mu$tail.prob==max(df.mu$tail.prob))
    mu.output <- data.frame(gene=g, mean.lfc=mean(df.mu$mean.lfc[ind]), tail.prob=max(df.mu$tail.prob))
    ind <- which(df.phi$tail.prob==max(df.phi$tail.prob))
    phi.output <- data.frame(gene=g, mean.lfc=mean(df.phi$mean.lfc[ind]), tail.prob=max(df.phi$tail.prob))

    return(list(mu.output=mu.output, phi.output=phi.output))
  })

  # outputs for all genes
  mu.output <- do.call(rbind,lapply(loop_global_output, function(l) l$mu.output))
  phi.output <- do.call(rbind,lapply(loop_global_output, function(l) l$phi.output))


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
    efdr_global_mu <- sapply(alphas, function(val) {
      efdr(P_g_star = mu.output$tail.prob,alpha_m = val)
    })
    alpha_m <- alphas[which.min(abs(efdr_global_mu-0.05))]
  }

  if(is.null(alpha_d)){
    alphas <- seq(0,1,length.out = 1000)
    efdr_global_phi <- sapply(alphas, function(val) {
      efdr(L_g_star = phi.output$tail.prob,alpha_d = val)
    })
    alpha_d <- alphas[which.min(abs(efdr_global_phi-0.05))]
  }

  # classify genes
  mu.output$DE <- ifelse(mu.output$tail.prob>alpha_m,'Yes','No')
  phi.output$DD <- ifelse(phi.output$tail.prob>alpha_d,'Yes','No')

  return(list(mu.output=mu.output,phi.output=phi.output, alpha_m=alpha_m,alpha_d=alpha_d))
}

#' Plot the results for global marker genes
#'
#' @description
#' Plot mean absolute log-fold change (LFC), as well as the overlap between global DE and DD genes.
#'
#'
#' @param global_output output from \code{global_marker_genes}.
#'
#' @return the output contains three ggplot objects:: 1. Tail probabilities against mean absolute LFC based on mean expression, and the threshold value to decide global DE genes.
#' 2. Tail probabilities against mean absolute LFC based on dispersion, and the threshold value to decide global DD genes.
#' 3. A bar-chart summarizing the overlap between DE and DD genes.
#' @export
#'
#' @examples
#' plot_global_marker_genes(global_output = global_result)
plot_global_marker_genes <- function(global_output){

  p1 <- ggplot(data = global_output$mu.output,aes(x=mean.lfc,y=tail.prob,col=DE))+
    geom_point(show.legend = F,size=0.5)+
    theme_bw()+
    scale_color_manual(values = c('Yes'='red','No'='grey'))+
    geom_hline(yintercept = global_output$alpha_m,linetype='dashed')+
    labs(x='Mean absolute log-fold change',y='Tail probability',title = 'Global DE')

  p2 <- ggplot(data = global_output$phi.output,aes(x=mean.lfc,y=tail.prob,col=DD))+
    geom_point(show.legend = F, size=0.5)+
    theme_bw()+
    scale_color_manual(values = c('Yes'='red','No'='grey'))+
    geom_hline(yintercept = global_output$alpha_d,linetype='dashed')+
    labs(x='Mean absolute log-fold change',y='Tail probability',title = 'Global DD')

  df <- data.frame(DE=ifelse(global_output$mu.output$DE=='Yes','DE','non_DE'),
                   DD=ifelse(global_output$phi.output$DD=='Yes','DD','non_DD'))

  p3 <- ggplot(data = df, aes(x=DE,fill=DD))+
    geom_bar(position = 'dodge')+
    theme_bw()+
    scale_fill_manual(values = c('DD'='red','non_DD'='grey'))+
    labs(x='Differentially expressed or not',y='Number of genes',
         fill='Classify phi', title='Classfication of genes')+
    theme(legend.position = 'right')

  return(list(p1,p2,p3))

}


#' Plot heatmaps for estimated mean expression and dispersion across clusters, for global marker genes
#'
#' @description
#' The function shows the heatmaps for estimated mean and dispersion parameters for global marker genes.
#'
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @param global_output output from \code{global_marker_genes}.
#' @param post_result output from \code{gkernelHDP_mcmc} for the post-processing step.
#'
#' @return Four plots. First two plots show the heatmaps of posterior mean of log mean and dispersion parameters, respectively,
#' where genes are ordered by decreasing tail probabilities from top to bottom with a black horizontal line to distinguish
#' global and non-global marker genes. The last two plots show the heatmaps for global marker genes only. For each gene, the estimates
#' are normalized such that the mean of the estimated parameters across clusters is zero.
#' @export
#'
#' @examples
#' global_marker_genes_heatmaps(global_output = global_result, post_result = post_result)
global_marker_genes_heatmaps <- function(global_output, post_result){

  J <- ncol(post_result$P_C_J_D_output[[1]][[1]])
  L <- post_result$output_index
  mu.sample <- post_result$mu_star_1_J_output
  phi.sample <- post_result$phi_star_1_J_output

  # threshold
  alpha_m <- global_output$alpha_m
  alpha_d <- global_output$alpha_d

  # ------ posterior mean of log(mu or phi) ----------
  # DE: compute posterior mean of logarithm of mu
  logmu.sample <- rapply(mu.sample, log, how = 'replace')
  # G x J
  mean_logmu <- t(Reduce('+',logmu.sample)/L)

  # ordered by decreasing tail probability
  mean_logmu <- mean_logmu[order(global_output$mu.output$tail.prob,decreasing = TRUE),]

  # plot
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  hp1 <- pheatmap(mean_logmu, scale = "none", cluster_rows  = FALSE,cluster_cols = FALSE,
                  angle_col = 90,fontsize = 8,border_color=NA,color = myPalette(400),
                  main = 'DE',silent = T)
  grid::grid.newpage()
  grid::grid.draw(hp1$gtable)
  grid::downViewport("matrix.4-3-4-3")
  grid::grid.lines(x=c(0,1), y=1-mean(global_output$mu.output$tail.prob>alpha_m), gp=grid::gpar(col="black", lwd=2,lty=2))
  grid::popViewport()

  # DD: compute posterior mean of logarithm of phi
  logphi.samples <- rapply(phi.sample, log, how = 'replace')
  # G x J
  mean_logphi <- t(Reduce('+',logphi.samples)/L)

  # ordered by decreasing tail probability
  mean_logphi <- mean_logphi[order(global_output$phi.output$tail.prob,decreasing = TRUE),]

  # plot
  hp2 <- pheatmap(mean_logphi, scale = "none", cluster_rows  = FALSE,cluster_cols = FALSE,
                  angle_col = 90,fontsize = 8,border_color=NA,color = myPalette(400),
                  main = 'DD',silent = T)
  grid::grid.newpage()
  grid::grid.draw(hp2$gtable)
  grid::downViewport("matrix.4-3-4-3")
  grid::grid.lines(x=c(0,1), y=1-mean(global_output$phi.output$tail.prob>alpha_d), gp=grid::gpar(col="black", lwd=2,lty=2))
  grid::popViewport()


  # ----- relative values for global genes only -----
  relative_logmu <- mean_logmu-rowMeans(mean_logmu)

  hp3 <- pheatmap(relative_logmu[1:sum(global_output$mu.output$tail.prob>alpha_m),],scale = "none",cluster_rows  = FALSE,
                  cluster_cols = FALSE,angle_col = 90,fontsize = 8,border_color=NA,color = myPalette(400),
                  main = 'DE (relative)',silent = T)
  grid::grid.newpage()
  grid::grid.draw(hp3$gtable)

  relative_logphi <- mean_logphi-rowMeans(mean_logphi)
  hp4 <- pheatmap(relative_logphi[1:sum(global_output$phi.output$tail.prob>alpha_d),],scale = "none",cluster_rows  = FALSE,
                  cluster_cols = FALSE,angle_col = 90,fontsize = 8,border_color=NA,color = myPalette(400),
                  main = 'DD (relative)',silent = T)
  grid::grid.newpage()
  grid::grid.draw(hp4$gtable)





}


















