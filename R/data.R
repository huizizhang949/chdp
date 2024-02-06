#' Simulated data from negative binomial likelihood with a Gaussian kernel
#'
#' @description
#' Two simulated datasets from a mixture of two binomial distributions, for \eqn{G=10} genes.
#' Both clusters are shared across two datasets. The Gaussian kernel is applied to introduce
#' dependence on the external covariate (latent time).
#'
#' @format A matrix with 240 rows and 14 columns, with rows for observations (cells):
#' \describe{
#' \enumerate{
#' \item First 120 rows are from the first dataset, with the rest from the second dataset.
#' \item First 10 columns correspond to genes.
#' \item The last four columns denote cell allocations, latent time, and time-dependent probabilities
#' of belonging to two clusters.
#' }
#' }
"sim1.data"



#' Simulated data from normal likelihood (vector autoregression) with a periodic kernel
#'
#' @description
#' Two simulated datasets from a mixture of two normal distributions, for \eqn{G=2} features.
#' Both clusters are shared across two datasets. The periodic kernel is applied to introduce
#' dependence on the external covariate (time).
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
