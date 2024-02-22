#' Produce a summary of cluster size
#'
#' @description
#' The function draws two bar plots to summarize the size of cluster.
#'
#' @param opt_cl a list of vectors. Each vector stores the optimal clustering for one dataset.
#' @param data_names the labels for each dataset in the plot.
#'
#' @return two bar plots. The first one shows the size of each cluster. The second one shows the proportion of each dataset
#' in every cluster, and compares that to the overall proportion (black solid line).
#' @export
#'
#' @examples cluster_size(opt_cl = opt_cl, data_names = c('data1', 'data2'))
cluster_size <- function(opt_cl, data_names){

  # number of identified clusters
  J <- length(unique(unlist(opt$opt_cl)))
  # cluster size
  cl_size <- table(unlist(opt$opt_cl))
  # plot
  p1 <- ggplot(data = data.frame(cluster=factor(1:J),size=as.numeric(unname(cl_size))),aes(x=cluster,y=size))+
    geom_bar(stat = 'identity',fill="steelblue")+
    geom_text(aes(label=size), vjust=-0.3, size=3.5)+
    theme_bw()

  # cluster size within each dataset
  cl_size1 <- table(opt_cl[[1]])
  cl_size2 <- table(opt_cl[[2]])

  # proportion of each dataset in every cluster
  df <- data.frame(cluster=rep(factor(1:J),2),
                   proportion=c(as.numeric(unname(cl_size1))/unname(cl_size),as.numeric(unname(cl_size2))/unname(cl_size)),
                   dataset=rep(data_names,each=J))
  p2 <- ggplot(data = df, aes(x=cluster,y=proportion,fill=dataset))+
    geom_bar(stat = 'identity')+
    geom_hline(yintercept = length(opt_cl[[2]])/length(unlist(opt_cl)))+
    theme_bw()

  gridExtra::grid.arrange(p1,p2,ncol=2)

}
