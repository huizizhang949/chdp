#' Plot each feature against time with for time-series data and show clustering
#'
#' @description
#' This function shows the line plots of each feature against time, with observations colored by clustering, for each dataset.
#'
#'
#' @param Y a list of two matrices for two datasets. The columns correspond to features.
#' @param t a list of two vectors. Each vector is the external covariate (time) for individual dataset.
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param data_names optional. The labels for each dataset in the plot.
#' @param color_pal optional. A vector of color names to map to clusters.
#' @param nrow_legend number of rows in legend.
#'
#' @return Line plots for each dataset, showing each feature against time, and clusters.
#' @export
#'
#' @examples
#' plot_var_cluster(Y=list(Y1,Y2), t=list(t1,t2), opt_cl = opt_cl, color_pal = c25)
plot_var_cluster <- function(Y, t, opt_cl, data_names=NULL, color_pal=NULL, nrow_legend=1){

  D <- length(Y); G <- ncol(Y[[1]])
  if(is.null(data_names)){
    data_names <- paste0('data ',1:D)
  }

  # for each dataset, save the data, t, optimal cluster, dataset indicator
  Y_df <- lapply(1:D, function(d) {

    df_d <- data.frame(Y[[d]][-1,]); df_d$cluster <- as.factor(opt_cl[[d]]); df_d$t <- t[[d]][-1]
    colnames(df_d)[1:G] <- paste0('y',1:G)
    df_d$dataset <- data_names[d]

    return(df_d)
  })

  Y_df <- do.call(rbind,Y_df)

  # turn into a long format
  Y_df_long <- reshape2::melt(Y_df, id.vars = c('t','cluster','dataset'), variable.name = "y")

  gg <- ggplot(data=Y_df_long)+
    geom_line(aes(x=t,y=value),col='grey',linetype=1,linewidth=0.3)+
    geom_point(aes(x=t,y=value,shape=cluster,col=cluster),size=0.6)+
    # scale_color_manual(values=color_pal)+
    # scale_shape_manual(values=seq(1,23))+
    theme_bw()+
    labs(color='cluster',shape='cluster',y='')+
    theme(legend.position = 'bottom')+
    guides(color = guide_legend(nrow = nrow_legend))+
    facet_grid(y ~ dataset)

  if(!is.null(color_pal)){
    gg <- gg+scale_color_manual(values=color_pal)
  }

  print(gg)
}
