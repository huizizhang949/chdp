#' Plot a heatmap of posterior similarity matrix across related groups
#'
#' @description
#' Produce a heatmap of posterior similarity matrix showing both within-dataset and between-dataset similarity,
#' where data points are ordered by hierarchical clustering.
#'
#' @importFrom stats hclust as.dist
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @param psm.tot a posterior similarity matrix for all observations across datasets. Can be obtained from
#' MCMC samples of clusterings by calling \code{comp.psm} in package \code{mcclust.ext}. Observations from
#' the same dataset should be grouped together.
#' @param size a vector of the size of each dataset. The order of the dataset should be the same as \code{psm.tot}.
#'
#' @return a heatmap of posterior similarity matrix, with datasets separated by black solid lines. The main
#' diagonal blocks correspond to within-dataset posterior similarity matrix.
#' @export
#'
#' @examples plot_psm(psm.tot = psm, size = c(120, 120))
plot_psm <- function(psm.tot,size){
  # size: the size of each dataset
  if(any(psm.tot !=t(psm.tot)) | any(psm.tot >1) | any(psm.tot < 0) | sum(diag(psm.tot)) != nrow(psm.tot) ){
    stop("psm.tot must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}

  my_palette <- colorRampPalette(c("white", 'yellow','red'))(n = 400)

  # total number of observations
  n_all <- sum(size)

  # proportion of each dataset
  n_cum <- cumsum(size)

  n_prop <- n_cum/n_all

  n_cum <- c(0,n_cum)

  D <- length(size)

  hc <- lapply(1:D, function(d) {
    hclust(as.dist(1-psm.tot[(n_cum[d]+1):n_cum[d+1],(n_cum[d]+1):n_cum[d+1]]))
  })

  ind <- unlist(sapply(1:D,function(d) hc[[d]]$order+n_cum[d]))

  x <- pheatmap(psm.tot[ind,ind],scale = "none",cluster_rows  = FALSE,
                cluster_cols = FALSE,
                angle_col = 90,fontsize = 8,color = my_palette,border_color=NA,silent=FALSE)
  grid::downViewport("matrix.4-3-4-3")
  for (i in n_prop[-length(n_prop)]) {
    grid::grid.lines(x=c(0,1),y=i, gp=grid::gpar(col="black", lwd=2))
    grid::grid.lines(x=i,y=c(0,1), gp=grid::gpar(col="black", lwd=2))
  }
  grid::popViewport()
}
