#' Simulated data from negative binomial likelihood with a Gaussian kernel
#'
#' @description
#' Two simulated datasets generated from a mixture of three binomial distributions with an overall mean of 0.6
#' for capture efficiency. The Gaussian kernel is applied to introduce dependence on the external covariate (latent time).
#'
#' @format A matrix with 500 rows and 105 columns, with rows for observations (cells):
#' \describe{
#' \enumerate{
#' \item First 200 rows are from the first dataset, with the rest from the second dataset.
#' \item First 100 columns correspond to genes, where first 70 genes are global DE and DD genes.
#' \item The last five columns denote cell allocations, latent time, and time-dependent probabilities
#' of belonging to three clusters.
#' }
#' }
"sim1.data"



#' Simulated data from VAR likelihood with a periodic kernel
#'
#' @description
#' Two simulated datasets from a mixture of two normal distributions with vector autoregression,
#' for \eqn{G=2} features. Both clusters are shared across two datasets.
#' The periodic kernel is applied to introduce dependence on the external covariate (time).
#'
#' @format A matrix with 300 rows and 6 columns, with rows for observations:
#' \describe{
#' \enumerate{
#' \item First 150 rows are from the first dataset, with the rest from the second dataset.
#' \item First 2 columns correspond to features.
#' \item The last four columns denote subject allocations, time, and time-dependent probabilities
#' of belonging to two clusters.
#' }
#' }
"sim2.data"

#' A vector of names for 15 distinguished colors
"c15"

#' A vector of names for 25 distinguished colors
"c25"


