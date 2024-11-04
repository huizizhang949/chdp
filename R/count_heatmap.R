#' Heatmaps for observed and latent counts
#'
#' @description
#' Heatmap for observed and latent counts after ordering cells by clusters and by datasets, and genes by
#' tail probabilities. The counts are on the log scale after adding a pseudo-count of 1.
#'
#' @param Y a list of matrices. Each matrix is a gene-by-cell matrix of mRNA counts corresponding to a dataset.
#' @param Y_latent a list of matrices. Each matrix is a gene-by-cell matrix of latent mRNA counts corresponding to a dataset.
#' Can use the direct output from \code{latent_count}.
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param global_output output from \code{global_marker_genes}.
#'
#' @return A heatmap with rows for genes and columns for cells. Genes are ordered by decreasing tail probabilities (DE)
#' from top to bottom, and a red dashed line separates global and non-globle DE genes. Cells are separated
#' by clusters (yellow vertical solid line) and also separated by datasets within each cluster (yellow vertical dashed line).
#' @export
#'
#' @examples
#' observed_counts_heatmap(Y = list(t(Y1),t(Y2)), opt_cl = opt_cl, global_output = global_result)
#' latent_counts_heatmap(Y_latent = Y_latent, opt_cl = opt_cl, global_output = global_result)
#' @name counts_heatmap
observed_counts_heatmap <- function(Y, opt_cl, global_output){

  D <- length(Y)
  Y_all <- do.call(cbind, Y)
  J <- length(unique(unlist(opt_cl)))
  C <- sapply(opt_cl, length)
  n_cum <- c(0,cumsum(C))
  alpha_m <- global_output$alpha_m

  # index of cells in each cluster
  ind_by_cl_data <- unlist(lapply(1:J, function(j) {

    temp <- lapply(1:D, function(d) {
      ind <- which(opt_cl[[d]]==j)+n_cum[d]
      return(ind)
    })

    return(unlist(temp))
  }))

  # order matrix for Y_all, genes ordered by decreasing tail probabilities, cells ordered by clusters
  Y_all <- as.matrix(Y_all)[order(global_output$mu.output$tail.prob,decreasing = TRUE),ind_by_cl_data]

  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

  # vertical solid line to separate clusters
  cl_size <- as.numeric(table(unlist(opt_cl)))
  cl_size_cum <- cumsum(cl_size)
  cum_prop <- cl_size_cum/ncol(Y_all)

  # vertical dashed line to separate datasets
  cl_size_by_d <- unlist(lapply(1:J, function(j) {

    temp <- lapply(1:D, function(d) {
      sum(opt_cl[[d]]==j)
    })

    return(unlist(temp))
  }))

  # remove those indices separateing clusters
  cl_size_by_d_cum <- cumsum(cl_size_by_d)
  cl_size_by_d_cum <- setdiff(cl_size_by_d_cum, cl_size_cum)
  cum_prop_d <- cl_size_by_d_cum/ncol(Y_all)

  heatmap_obs <- pheatmap(log(Y_all+1),scale = "none",cluster_rows  = FALSE,
                          cluster_cols = FALSE,
                          angle_col = 90,fontsize = 8,border_color=NA,color = myPalette(n=400),
                          main = 'Observed counts',silent = T)

  grid::grid.newpage()
  grid::grid.draw(heatmap_obs$gtable)
  grid::downViewport("matrix.4-3-4-3")
  # vertical solid line to separate clusters
  for (i in 1:length(cum_prop[-J])) {
    grid::grid.lines(x=cum_prop[i],y=c(0,1), gp=grid::gpar(col="yellow", lwd=2))
  }
  # vertical dashed line to separate datasets
  for (i in 1:length(cum_prop_d)) {
    grid::grid.lines(x=cum_prop_d[i],y=c(0,1), gp=grid::gpar(col="yellow", lwd=2,lty=2))
  }
  grid::grid.lines(x=c(0,1),y=1-mean(global_output$mu.output$tail.prob>alpha_m), gp=grid::gpar(col="red", lwd=2,lty=2))
  grid::popViewport()

}

#' @rdname counts_heatmap
#' @export
latent_counts_heatmap <- function(Y_latent, opt_cl, global_output){

  D <- length(Y_latent)
  Y_latent_all <- do.call(cbind, Y_latent)
  J <- length(unique(unlist(opt_cl)))
  C <- sapply(opt_cl, length)
  n_cum <- c(0,cumsum(C))
  alpha_m <- global_output$alpha_m

  # index of cells in each cluster
  ind_by_latent_cl_data <- unlist(lapply(1:J, function(j) {

    temp <- lapply(1:D, function(d) {
      ind <- which(opt_cl[[d]]==j)+n_cum[d]
      return(ind)
    })

    return(unlist(temp))
  }))

  # order matrix for Y_latent_all, genes ordered bY_latent decreasing tail probabilities, cells ordered bY_latent clusters
  Y_latent_all <- as.matrix(Y_latent_all)[order(global_output$mu.output$tail.prob,decreasing = TRUE),ind_by_latent_cl_data]

  my_latentPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

  # vertical solid line to separate clusters
  cl_size <- as.numeric(table(unlist(opt_cl)))
  cl_size_cum <- cumsum(cl_size)
  cum_prop <- cl_size_cum/ncol(Y_latent_all)

  # vertical dashed line to separate datasets
  cl_size_by_latent_d <- unlist(lapply(1:J, function(j) {

    temp <- lapply(1:D, function(d) {
      sum(opt_cl[[d]]==j)
    })

    return(unlist(temp))
  }))

  # remove those indices separating clusters
  cl_size_by_latent_d_cum <- cumsum(cl_size_by_latent_d)
  cl_size_by_latent_d_cum <- setdiff(cl_size_by_latent_d_cum, cl_size_cum)
  cum_prop_d <- cl_size_by_latent_d_cum/ncol(Y_latent_all)

  heatmap_latent <- pheatmap(log(Y_latent_all+1),scale = "none",cluster_rows  = FALSE,
                          cluster_cols = FALSE,
                          angle_col = 90,fontsize = 8,border_color=NA,color = my_latentPalette(n=400),
                          main = 'Latent counts',silent = T)

  grid::grid.newpage()
  grid::grid.draw(heatmap_latent$gtable)
  grid::downViewport("matrix.4-3-4-3")
  # vertical solid line to separate clusters
  for (i in 1:length(cum_prop[-J])) {
    grid::grid.lines(x=cum_prop[i],y=c(0,1), gp=grid::gpar(col="yellow", lwd=2))
  }
  # vertical dashed line to separate datasets
  for (i in 1:length(cum_prop_d)) {
    grid::grid.lines(x=cum_prop_d[i],y=c(0,1), gp=grid::gpar(col="yellow", lwd=2,lty=2))
  }
  grid::grid.lines(x=c(0,1),y=1-mean(global_output$mu.output$tail.prob>alpha_m), gp=grid::gpar(col="red", lwd=2,lty=2))
  grid::popViewport()

}
