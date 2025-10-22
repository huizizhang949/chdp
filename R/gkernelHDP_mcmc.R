#' Perform MCMC for covariate-dependent hierarchical Dirichlet process to cluster across two single-cell RNA-sequencing datasets
#'
#' @description
#' Implement Gibbs sampling and adaptive Metropolis-Hastings algorithm to draw each parameter sequentially.
#' A Gaussian kernel is applied to introduce dependence on the external covariate, i.e., latent time, and the
#' algorithm provides posterior samples for all parameters including cluster allocations. The function is also
#' used in the post-processing step to fix allocations to the optimal clustering, in order to infer cluster-specific parameters.
#'
#' @importFrom stats coef deviance df.residual dnbinom kmeans lm plnorm pnorm qlnorm rbinom rgamma rlnorm rnorm runif var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import irlba
#'
#' @param Y a list of two matrices. Each is a gene-by-cell \eqn{(G \times C_d)} matrix of mRNA counts for individual dataset \eqn{d}.
#' Both datasets should share the same set of genes.
#' @param t a list of two vectors. Each vector is the external covariate (latent time) for individual dataset.
#' @param niter integer. Total number of MCMC iterations.
#' @param J integer. Truncation level, i.e., the number of clusters to be found. Only required when estimating \code{Z}.
#' @param burn_in the length of burn-in period during MCMC.
#' @param thinning the thinning applied after burn-in period.
#' @param empirical logical. If TRUE, empirical values are used to set up the prior for \eqn{\alpha_{\phi}^2} and \eqn{\mathbf{b}}.
#' @param empirical_z if \code{Z_fix} is not provided and \code{empirical_z=TRUE}, t-sne is first performed to reduce
#' the combined datasets into two dimensions, and k-means is applied to find \eqn{J} clusters on
#' the lower-dimensional embeddings. If \code{empirical_z=FALSE}, clusters are randomly initialized.
#' @param Z_fix optional. Should be provided when \code{empirical_z} and \code{J} are not available. A list of vectors,
#' where each vector stores the allocations in \eqn{1,\ldots,J} in one dataset. Typically used in the
#' post-processing step where allocations are fixed to the optimal clustering to infer component-specific parameters.
#' @param BB_SIZE optional. Default to TRUE. See \code{bayNorm} for details.
#' @param baynorm_beta,baynorm_tot a list of bayNorm estimates of capture efficiencies for each data (will be used as initial values),
#' bayNorm results obtained from applying bayNorm to combined datasets.
#' @param alpha_initial,alpha_0_initial initial values for concentration parameters \eqn{\alpha,\alpha_0}. Default to 1.
#' @param quadratic If \code{TRUE}, a quadratic relationship is assumed between \eqn{\phi^*_{j,g}} and \eqn{\mu^*_{j,g}}.
#' @param MH.variance additional variance added to the empirical covariance (variance) in adaptive Metropolis-Hastings.
#' @param mu_r,sigma_r mean and standard deviation of the hyper-prior (normal) for \eqn{r_j}.
#' @param eta_1,eta_2 shape and scale of the hyper-prior (inverse-gamma) for \eqn{s^2}.
#' @param mu_h,sigma_h mean and standard deviation of the hyper-prior (normal) for \eqn{h_j}.
#' @param kappa_1,kappa_2 shape and scale of the hyper-prior (inverse-gamma) for \eqn{m^2}.
#' @param beta.mean an estimate of global mean capture efficiency across cells, default to be 0.06. Used to set up prior parameters for capture efficiency.
#' @param alpha_mu_2 optional. Prior variance for \eqn{\mu^*_{j,g}}. If not provided, empirical values will be used.
#' @param partial_pca if \code{empirical_z=TRUE}, whether truncated PCA should be used to calculate principal components (requires the irlba package).
#' See \code{Rtsne}.
#'
#' @usage gkernelHDP_mcmc(Y, t, J = NULL, niter, burn_in = 1000, thinning = 1,
#'     empirical = TRUE, empirical_z = TRUE, Z_fix = NULL,
#'     quadratic = FALSE, MH.variance = 0.01, BB_SIZE = TRUE,
#'     baynorm.beta = baynorm.beta, baynorm_tot = baynorm_tot,
#'     alpha_initial = 1, alpha_0_initial = 1,
#'     mu_r = 0.5, sigma_r = 0.5, eta_1 = 5, eta_2 = 1, mu_h = -5,
#'     sigma_h = 0.5, kappa_1 = 5, kappa_2 = 1,
#'     beta.mean = 0.06, alpha_mu_2 = NULL,
#'     partial_pca = FALSE)
#'
#' @return \code{gkernelHDP_mcmc} returns a list containing the following components:
#' \item{output_index}{total number of saved MCMC samples, taking into account of burn-in and thinning.}
#' \item{b_output}{a list of vectors for sampled \eqn{\mathbf{b}}. The list is of length \code{output_index}.}
#' \item{alpha_phi_2_output}{a vector for sampled \eqn{\alpha_{\phi}^2}. The vector is of length \code{output_index}.}
#' \item{Z_output}{If \code{Z_fix} is not provided, a list of length \code{output_index}. Each element of
#' the list is a list of two vectors saving the allocations in individual dataset. Otherwise is \code{NULL}.}
#' \item{P_C_J_D_output}{a list of length two, each corresponding to a single dataset. Each element of the list is
#' a list of length \code{output_index}, saving the samples for time-dependent probabilities in a \eqn{C_d \times J} matrix .}
#' \item{P_output}{a list of length \code{output_index}. Each element of the list is a vector of length \eqn{J} for component probabilities \eqn{\mathbf{p}^J}.}
#' \item{alpha_output}{a vector of sampled concentration parameter \eqn{\alpha}.}
#' \item{alpha_0_output}{samples for concentration parameter \eqn{\alpha_0}. Similar to \code{alpha_output}.}
#' \item{mu_star_1_J_output}{a list of length \code{output_index}. Each element of the list is a \eqn{J \times G} matrix for mean expression \eqn{\mu_{j,g}^*}.}
#' \item{phi_star_1_J_output}{samples for dispersion parameter \eqn{\phi_{j,g}^*}. Similar to \code{mu_star_1_J_output}.}
#' \item{Beta_output}{samples for capture efficiencies. Similar to \code{Z_output}.}
#' \item{Q_J_D_output}{a list of length \code{output_index}. Each element of the list is a \eqn{J \times D} matrix for \eqn{q_{j,d}}.}
#' \item{Xi_C_D_output}{samples for latent variables. Similar to \code{Z_output}.}
#' \item{t_star_J_D_output}{samples for kernel parameters (mean). Similar to \code{Q_J_D_output}.}
#' \item{sigma_star_2_J_D_output}{samples for kernel parameters (variance). Similar to \code{Q_J_D_output}.}
#' \item{U_C_J_D_output}{samples for latent variables. Similar to \code{P_C_J_D_output}.}
#' \item{r_J_output}{a list of length \code{output_index}. Each element of the list is a vector of length \eqn{J} for \eqn{r_j}.}
#' \item{s_2_output}{a vector of sampled \eqn{s^2}.}
#' \item{h_J_output}{samples for \eqn{h_j}. Similar to \code{r_J_output}.}
#' \item{m_2_output}{samples for \eqn{m^2}. Similar to \code{s_2_output}.}
#' \item{acceptance_count_avg}{average acceptance probabilities over iterations for parameters drawn using adaptive Metroplis-Hastings.}
#' \item{alpha_mu_2}{prior variance for \eqn{\mu^*_{j,g}}.}
#' \item{v_1, v_2}{prior parameters for \eqn{\alpha_{\phi}^2}.}
#' \item{m_b}{prior mean for \eqn{\mathbf{b}}.}
#' \item{a_d_beta,b_d_beta}{prior parameters for \eqn{\beta_{c,d}.}}
#'
#' @export
#'
#' @examples
#' gkernelHDP_mcmc(Y = list(t(y1), t(y2)), t = list(t1, t2), J = 4, niter = 1000,
#'                 burn_in = 0, thinning = 1, empirical = TRUE, empirical_z = TRUE)
gkernelHDP_mcmc <- function(Y, t, J = NULL, niter, burn_in = 1000, thinning = 1,
                            empirical = TRUE, empirical_z = TRUE, Z_fix = NULL,
                            quadratic = FALSE, MH.variance = 0.01, BB_SIZE = TRUE,
                            baynorm.beta, baynorm_tot,
                            alpha_initial = 1, alpha_0_initial = 1,
                            mu_r = 0.5, sigma_r = 0.5, eta_1 = 5, eta_2 = 1, mu_h = -5,
                            sigma_h = 0.5, kappa_1 = 5, kappa_2 = 1,
                            beta.mean = 0.06, alpha_mu_2 = NULL,
                            partial_pca = FALSE){

  # Show time
  start_time <- Sys.time()
  print(paste('Start Initialization:', start_time))
  cat('\n')

  # Sample size after burn-in and thinning
  # Apply burn-in first, then apply thinning to the rest samples
  # output_size <- floor((niter-burn_in)/thinning)

  # Define Gaussian kernel (Radial basis function)
  rbf <- function(t_c_d, t_j_d, sigma_j_d_2) {exp(-(t_c_d-t_j_d)^2/(2*sigma_j_d_2))}

  # number of dataset, dimension of the data
  D <- length(Y); G <- nrow(Y[[1]])

  # Data size
  C_d <- sapply(Y, ncol)

  if(is.null(J) & is.null(Z_fix)){
    stop('At least one of J or Z_fix should be provided!')
  }

  if(is.null(J)){
    J <- length(unique(unlist(Z_fix)))
  }

  #----------------------- bayNorm estimates -------------

  # Estimate of mu
  baynorm.mu.estimate <- baynorm_tot$PRIORS$MME_prior[,1]

  # Estimate of phi
  if(!BB_SIZE){
    baynorm.phi.estimate <- baynorm_tot$PRIORS$MME_prior[,2]
  }else{
    baynorm.phi.estimate <- baynorm_tot$PRIORS$MME_SIZE_adjust
  }

  #------------------------ Step 1: Prepare for outputs -----------------
  b_output <- NULL
  alpha_phi_2_output <- c()
  Z_output <- NULL

  P_C_J_D_output <- rep(list(NULL),D)

  P_output <- NULL
  alpha_output <- c()
  alpha_0_output <- c()
  mu_star_1_J_output <- NULL
  phi_star_1_J_output <- NULL
  Beta_output <- NULL


  Q_J_D_output <- NULL
  Xi_C_D_output <- NULL
  t_star_J_D_output <- NULL
  sigma_star_2_J_D_output <- NULL

  U_C_J_D_output <- rep(list(NULL),D)

  r_J_output <- NULL
  s_2_output <- c()
  h_J_output <- NULL
  m_2_output <- c()

  sd_P_output <- c()
  #----------------------- Step 2: Initial values in MCMC-----------

  # ------ Initial Z --------
  # Use K-means for initialization of Z, t_star_J_D, sigma_star_2_J_D, r_J, s^2, h_j, m^2
  # after dimension reduction (t-sne)
  if(is.null(Z_fix)){
    if(empirical_z==TRUE){
      Y_all <- t(do.call(cbind, Y))
      tsne_results <- Rtsne::Rtsne(Y_all, perplexity=30, check_duplicates = FALSE, partial_pca = partial_pca)
      km_cluster <- kmeans(tsne_results$Y,centers=J)$cluster

      # Allocation variables Z
      Z_initial <- NULL
      Z_initial[[1]] <- km_cluster[1:C_d[1]]
      Z_initial[[2]] <- km_cluster[(C_d[1]+1):(C_d[2]+C_d[1])]

    }else{
      Z_initial <- lapply(1:D, function(d) sample(1:J, size = C_d[[d]],replace = TRUE,prob = rep(1/J,J)))
    }

  }else{
    Z_initial <- Z_fix
  }

  # ---------- Initial kernel and hyper-paramerters ---------
  # Compute the mean of t within each cluster in each dataset to initialize t_star_J_D
  # There may be empty clusters
  # Fill the non-empty values into the initialization matrix, impute NA values with dataset-specific mean
  t_star_J_D_list <- lapply(1:D, function(d) {

    val <- rep(NA, J)
    # Mean of t within each cluster of each dataset
    t_star_J <- tapply(t[[d]], Z_initial[[d]], mean)
    val[as.numeric(names(t_star_J))] <- t_star_J
    val[is.na(val)] <- mean(t[[d]])

    return(val)
  })

  # matrix of [J, D]
  t_star_J_D_initial <- do.call(cbind,t_star_J_D_list)

  # Initials for r_J as an empirical mean from t_star_J_D_initial
  r_J_initial <- apply(t_star_J_D_initial, 1, mean)

  # Initials for s^2
  s_2_initial <- mean(apply(t_star_J_D_initial, 1, var))

  # Initials for sigma_star_2_J_D
  sigma_star_2_J_D_list <- lapply(1:D, function(d) {

    val <- rep(NA, J)
    # Mean of t within each cluster of each dataset
    sigma_star_J <- tapply(t[[d]], Z_initial[[d]], var)
    val[as.numeric(names(sigma_star_J))] <- sigma_star_J
    val[is.na(val)] <- var(t[[d]])

    return(val)
  })

  # matrix of [J, D]
  sigma_star_2_J_D_initial <- do.call(cbind,sigma_star_2_J_D_list)

  # Initials for h_J from log(sigma_star_2_J_D)
  h_J_initial <- apply(log(sigma_star_2_J_D_initial), 1, mean)

  # Initials for m^2
  m_2_initial <- mean(apply(log(sigma_star_2_J_D_initial), 1, var))

  # ----------- Initial P and q ---------------
  # Component probabilities p_j, add 1 to avoid zero p if cluster is empty
  P_initial <- sapply(1:J, function(j) mean(unlist(Z_initial)==j))
  if(any(P_initial==0)){
    # Add 1 to avoid zero p if cluster is empty
    P_initial <- sapply(1:J, function(j) sum(unlist(Z_initial)==j)+1)/(sum(C_d)+J)
  }

  # Dataset-specific vector q_j_d
  Q_J_D_initial <- do.call(cbind, lapply(1:D, function(d) {
    val <- sapply(1:J,function(j) sum(Z_initial[[d]]==j))+1
    return(val)
  }))

  # ----------- Initial mu, phi ----------

  mu.estimate <- ifelse(baynorm.mu.estimate == 0, 0.01, baynorm.mu.estimate)
  phi.estimate <- ifelse(baynorm.phi.estimate == 0, 0.01, baynorm.phi.estimate)

  #cluster-specific parameters in NB likelihood
  mu_star_1_J_initial <- t(matrix(mu.estimate, nrow=G, ncol=J))
  phi_star_1_J_initial <- t(matrix(phi.estimate, nrow=G, ncol=J))

  # ----------- Initial b ---------
  if(quadratic){
    b_initial <- c(0,1,0)
  }else{
    b_initial <- c(0,1)
  }

  # --------- Prior: alpha_mu_2, v_1, v_2 and m_b. and Initial alpha_phi_2 -------------
  # alpha_mu_2
  if(is.null(alpha_mu_2)){
    # prior log(mu) ~ N(0, alpha_mu^2), so alpha_mu^2 = mean(log(mu)^2)
    log.mu.finite.index <- which(is.finite(log(baynorm.mu.estimate)))
    alpha_mu_2 <- mean(c(log(baynorm.mu.estimate[log.mu.finite.index]))^2)
  }

  # m_b, v1 and v2: based on non-zero mu
  x.1 <- log(baynorm.mu.estimate[is.finite(log(baynorm.mu.estimate))])
  y.1 <- log(baynorm.phi.estimate[is.finite(log(baynorm.mu.estimate))])

  if(quadratic==FALSE){
    lm.1 <- lm(y.1 ~ x.1)
  }else{
    lm.1 <- lm(y.1 ~ x.1+I(x.1^2))
  }

  # estimate of alpha_phi_2
  rse.lm.1.squared <- deviance(lm.1)/df.residual(lm.1)
  # variance of alpha_phi^2 is 1, mean is rse.lm.2.squared in the inverse gamma prior for alpha_phi^2
  variance <- 1
  v_1_empirical <- rse.lm.1.squared^2/variance + 2
  v_2_empirical <- (v_1_empirical - 1)*rse.lm.1.squared

  if(empirical == TRUE){
    v_1 <- v_1_empirical
    v_2 <- v_2_empirical
  }else{
    v_1 <- 2; v_2 <- 1
  }

  alpha_phi_2_initial <- rse.lm.1.squared

  if(empirical == TRUE){
    m_b <- as.numeric(coef(lm.1))
  }else{
    if(quadratic==TRUE){
      m_b <- c(-1,2,0)
    }else{
      m_b <- c(-1,2)
    }
  }

  # --------- Prior: a_d_beta and b_d_beta and Initials for beta -------------
  Beta_initial <- baynorm.beta
  loop.result <- lapply(1:D, function(d){

    ##-- Initial mean and variance
    baynorm_mean_capeff <- beta.mean
    baynorm_var_capeff <- var(baynorm.beta[[d]]) + 0.01

    ##-- Initial a_beta and b_beta
    a_beta <- ((1-baynorm_mean_capeff)/baynorm_var_capeff - 1/baynorm_mean_capeff)*baynorm_mean_capeff^2
    b_beta <- a_beta*(1/baynorm_mean_capeff - 1)

    while((baynorm_var_capeff >= baynorm_mean_capeff*(1-baynorm_mean_capeff)) | (a_beta < 1) | (b_beta < 1)){

      ##-- While the above condition is true, half the variance
      baynorm_var_capeff <- baynorm_var_capeff/2

      ##-- Compute the new a_beta and b_beta
      a_beta <- ((1-baynorm_mean_capeff)/baynorm_var_capeff - 1/baynorm_mean_capeff)*baynorm_mean_capeff^2
      b_beta <- a_beta*(1/baynorm_mean_capeff - 1)

    }

    return(c('a_beta' = a_beta,'b_beta' = b_beta))
  })

  a_d_beta <- unlist(lapply(loop.result, function(l) l['a_beta']))
  b_d_beta <- unlist(lapply(loop.result, function(l) l['b_beta']))

  if(any(c(a_d_beta,b_d_beta)<0)) {
    stop('Prior setting for beta inappropriate!!')
  }

  #-------------- Average acceptance probabilities ---------------------

  # Acceptance probability
  acceptance_count_avg <- data.frame(P_accept = rep(0,niter), alpha_accept = rep(0,niter),
                                     alpha_0_accept = rep(0,niter), unique_accept = rep(0,niter),
                                     Beta_accept = rep(0,niter), sigma_accept = rep(0, niter))

  #----------------------- Step 3: Set the initial_values as new values ----------------------------
  t_star_J_D_new <- t_star_J_D_initial
  sigma_star_2_J_D_new <- sigma_star_2_J_D_initial
  r_J_new <- r_J_initial
  s_2_new <- s_2_initial
  h_J_new <- h_J_initial
  m_2_new <- m_2_initial
  Q_J_D_new <- Q_J_D_initial

  Xi_C_D_new <- Xi_C_D_update_gkernel(Q_J_D = Q_J_D_new, C_d = C_d, t = t, t_star_J_D = t_star_J_D_new,
                                      sigma_star_2_J_D = sigma_star_2_J_D_new)

  U_C_J_D_new <- U_C_J_D_update_gkernel(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C_d = C_d, t = t,
                                        t_star_J_D = t_star_J_D_new, sigma_star_2_J_D = sigma_star_2_J_D_new,
                                        rbf = rbf)

  b_new <- b_initial
  alpha_phi_2_new <- alpha_phi_2_initial
  Z_new <- Z_initial
  P_new <- P_initial
  alpha_new <- alpha_initial
  alpha_0_new <- alpha_0_initial
  mu_star_1_J_new <- mu_star_1_J_initial
  phi_star_1_J_new <- phi_star_1_J_initial
  Beta_new <- Beta_initial

  # Total count of acceptance
  P_count <- 0; alpha_count <- 0; alpha_0_count <- 0; unique_count <- 0; Beta_count <- 0
  sigma_star_2_count <- 0

  #----------------------- Step 4: Prepare for the adaptive covariance update -----------------------------

  # 0) For sigma_star_2_J_D, matrix
  # X = -log(1/sigma_star - 1/upper)
  mean_X_sigma_star_2_new <- log(sigma_star_2_J_D_new)
  M_2_sigma_star_2_new <- matrix(0, nrow = J, ncol = D)
  variance_sigma_star_2_new <- matrix(0, nrow = J, ncol = D)

  for (d in 1:D) {
    for (j in 1:J) {
      ind <- c(1:C_d[d])[-log(U_C_J_D_new[[d]][,j])<Xi_C_D_new[[d]]*Q_J_D_new[j,d]]
      if(length(ind)!=0) {
        uppers <- -(t[[d]]-t_star_J_D_new[j,d])^2/2/(log(-log(U_C_J_D_new[[d]][,j]))-log(Xi_C_D_new[[d]])-log(Q_J_D_new[j,d]))
        upper <- min(uppers[ind])
        mean_X_sigma_star_2_new[j,d] <- -log(1/sigma_star_2_J_D_new[j,d] - 1/upper)
      }
    }
  }

  # 1) For Component probabilities
  sd_P_new <- 0.001
  mean_X_component_new <- log(matrix(P_initial[1:(J-1)]/P_initial[J], nrow = 1)) # 1x(J-1)
  tilde_s_component_new <- t(mean_X_component_new)%*%mean_X_component_new
  # At the 1st iteration, the covariance based on the initial values are 0
  covariance_component_new <- matrix(0, nrow = J-1, ncol = J-1)

  # 2) For alpha
  mean_X_alpha_new <- log(alpha_new)
  M_2_alpha_new <- 0
  variance_alpha_new <- 0

  # 3) For alpha_0
  mean_X_alpha_0_new <- log(alpha_0_new)
  M_2_alpha_0_new <- 0
  variance_alpha_0_new <- 0

  # 4) Unique parameters
  # JxG 2x2 matrices
  covariance_unique_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
  tilde_s_unique_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
  mean_X_unique_new <- rep(list(rep(list(matrix(0,nrow=1,ncol=2)),G)),J)
  for(j in 1:J){
    mean_X_unique_new[[j]] <- lapply(1:G, function(g) {
      matrix(c(log(mu_star_1_J_new[j,g]),log(phi_star_1_J_new[j,g])),nrow=1)
    })
    tilde_s_unique_new[[j]] <- lapply(1:G, function(g) {
      t(mean_X_unique_new[[j]][[g]])%*%mean_X_unique_new[[j]][[g]]
    })
  }

  # 5) Capture efficiency
  mean_X_capture_new <- lapply(1:D, function(d) log(Beta_new[[d]]/(1-Beta_new[[d]])))
  M_2_capture_new <- lapply(1:D, function(d) rep(0, C_d[d]))
  variance_capture_new <- lapply(1:D, function(d) rep(0, C_d[d]))

  # Index for each saved MCMC sample
  output_index <- 0

  start_time_mcmc <- Sys.time()
  print(paste('Start MCMC:', start_time_mcmc))
  cat('\n')

  pb <- txtProgressBar(min = 1, max = niter+1, style = 3)

  #----------------------- Step 5: Updates -----------------------------
  # Iteration starts with iter_num = 2
  for(iter in 2:(niter+1)){

    setTxtProgressBar(pb, iter)
    # Starting value of the output index = 1
    # If the current iteration is greater than the burn in and divisible by the thinning index
    criterion <- (iter-1 > burn_in & (iter-1-burn_in)%%thinning == 0)

    if(criterion){
      output_index <- output_index + 1
      update <- TRUE
    }else{
      update <- FALSE
    }


    # 0) --- Update latent/auxiliary variables Xi_C_D ----
    Xi_C_D_new <- Xi_C_D_update_gkernel(Q_J_D = Q_J_D_new, C_d = C_d, t = t, t_star_J_D = t_star_J_D_new,
                                        sigma_star_2_J_D = sigma_star_2_J_D_new)

    # 1) --- Update dataset-specific vector q_j_d ----
    Q_J_D_new <- Q_J_D_update_gkernel(Z = Z_new, alpha = alpha_new, P = P_new, Xi_C_D = Xi_C_D_new,
                                      t = t, t_star_J_D = t_star_J_D_new,
                                      sigma_star_2_J_D = sigma_star_2_J_D_new)

    # 2) ---- Update the allocation variable, if Z_fix is not provided -----
    if(is.null(Z_fix)){
      Z_new <- allocation_variables_update_gkernel(Y = Y, t = t, mu_star_1_J = mu_star_1_J_new, phi_star_1_J = phi_star_1_J_new,
                                                   Beta = Beta_new, Q_J_D = Q_J_D_new, t_star_J_D = t_star_J_D_new,
                                                   sigma_star_2_J_D = sigma_star_2_J_D_new)
    }


    # 3) ---- Update kernel parameters -----
    # 3-1) latent variable U_C_J_D
    U_C_J_D_new <- U_C_J_D_update_gkernel(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C_d = C_d, t = t,
                                          t_star_J_D = t_star_J_D_new, sigma_star_2_J_D = sigma_star_2_J_D_new,
                                          rbf = rbf)

    # 3_2) t_star_J_D
    t_star_J_D_new <- t_star_J_D_update(r_J = r_J_new, s_2 = s_2_new, Z = Z_new, t = t, C_d = C_d,
                                        sigma_star_2_J_D = sigma_star_2_J_D_new, U_C_J_D = U_C_J_D_new,
                                        Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new)

    # Update U again
    U_C_J_D_new <- U_C_J_D_update_gkernel(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C_d = C_d, t = t,
                                          t_star_J_D = t_star_J_D_new, sigma_star_2_J_D = sigma_star_2_J_D_new,
                                          rbf = rbf)

    # 3_3) sigma_star_2_J_D
    sigma_star_output <- sigma_star_2_J_D_update(sigma_star_2_J_D_old = sigma_star_2_J_D_new,
                                                 h_J = h_J_new, m_2 = m_2_new, Z = Z_new, t = t, C_d = C_d,
                                                 t_star_J_D = t_star_J_D_new, U_C_J_D = U_C_J_D_new,
                                                 Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new,
                                                 X_mean = mean_X_sigma_star_2_new, M_2 = M_2_sigma_star_2_new,
                                                 variance = variance_sigma_star_2_new, iter_num = iter,
                                                 MH.variance = MH.variance)


    sigma_star_2_J_D_new <- sigma_star_output$sigma_star_2_J_D_new
    mean_X_sigma_star_2_new <- sigma_star_output$X_mean_new
    M_2_sigma_star_2_new <- sigma_star_output$M_2_new
    variance_sigma_star_2_new <- sigma_star_output$variance_new
    sigma_star_2_count <- sigma_star_2_count + sigma_star_output$accept
    acceptance_count_avg$sigma_accept[iter-1] <- sigma_star_2_count/((iter-1)*J*D)


    # ---- Calculate P_C_J_D -------
    P_C_J_D_new <- lapply(1:D, function(d) {
      temp <- t(sapply(t[[d]],function(time) {
        # log-scale
        LP <- log(Q_J_D_new[,d])-(time-t_star_J_D_new[,d])^2/2/sigma_star_2_J_D_new[,d]
        nc <- -max(LP)
        P <- exp(LP+nc)/sum(exp(LP+nc))
        return(P)
      }))
    })


    # 4) ----- Update hyper parameters in priors for t_star, sigma_star_2 ------
    # 4-1) r_J
    r_J_new <- r_J_update_gkernel(t_star_J_D = t_star_J_D_new, mu_r = mu_r, sigma_r = sigma_r, s_2 = s_2_new)

    # 4-2) s_2
    s_2_new <- s_2_update_gkernel(t_star_J_D = t_star_J_D_new, eta_1 = eta_1, eta_2 = eta_2, r_J = r_J_new)

    # 4-3) h_J
    h_J_new <- h_J_update_gkernel(sigma_star_2_J_D = sigma_star_2_J_D_new, mu_h = mu_h, sigma_h = sigma_h, m_2 = m_2_new)

    # 4_4) m_2
    m_2_new <- m_2_update_gkernel(sigma_star_2_J_D = sigma_star_2_J_D_new, kappa_1 = kappa_1, kappa_2 = kappa_2, h_J = h_J_new)


    # 5) ----- Update the component probabilities P ----------
    component_output <- component_probabilities_update(P = P_new, Q_J_D = Q_J_D_new, alpha_0 = alpha_0_new,
                                                       alpha = alpha_new, covariance = covariance_component_new,
                                                       mean_x = mean_X_component_new,
                                                       tilde_s = tilde_s_component_new,
                                                       iter_num = iter, sd_P = sd_P_new,
                                                       MH.variance = MH.variance)
    P_new <- component_output$P_new
    tilde_s_component_new <- component_output$tilde_s_new
    mean_X_component_new <- component_output$mean_x_new
    covariance_component_new <- component_output$covariance_new
    sd_P_new <- component_output$sd_P_new
    P_count <- P_count + component_output$accept
    acceptance_count_avg$P_accept[iter-1] <- P_count/(iter-1)

    sd_P_output[iter-1] <- sd_P_new

    # 6) ------- Update alpha -----------
    alpha_output_sim <- alpha_update(Q_J_D = Q_J_D_new, P = P_new, alpha = alpha_new,
                                    X_mean = mean_X_alpha_new, M_2 = M_2_alpha_new,
                                    variance = variance_alpha_new, iter_num = iter,
                                    MH.variance = MH.variance)

    alpha_new <- alpha_output_sim$alpha_new
    mean_X_alpha_new <- alpha_output_sim$X_mean_new
    M_2_alpha_new <- alpha_output_sim$M_2_new
    variance_alpha_new <- alpha_output_sim$variance_new
    alpha_count <- alpha_count + alpha_output_sim$accept
    acceptance_count_avg$alpha_accept[iter-1] <- alpha_count/(iter-1)

    # 7) -------- Update alpha_0 ----------
    alpha_0_output_sim <- alpha_0_update(P = P_new, alpha_0 = alpha_0_new,
                                         X_mean = mean_X_alpha_0_new,
                                         M_2 = M_2_alpha_0_new,
                                         variance = variance_alpha_0_new, iter_num = iter,
                                         MH.variance = MH.variance)

    alpha_0_new <- alpha_0_output_sim$alpha_0_new
    mean_X_alpha_0_new <- alpha_0_output_sim$X_mean_new
    M_2_alpha_0_new <- alpha_0_output_sim$M_2_new
    variance_alpha_0_new <- alpha_0_output_sim$variance_new
    alpha_0_count <- alpha_0_count + alpha_0_output_sim$accept
    acceptance_count_avg$alpha_0_accept[iter-1] <- alpha_0_count/(iter-1)

    # 8) ------ Update mean_dispersion -------
    mean_dispersion_output <- mean_dispersion(mu_star_1_J = mu_star_1_J_new,
                                              phi_star_1_J = phi_star_1_J_new,
                                              v_1 = v_1, v_2 = v_2, m_b = m_b,
                                              quadratic = quadratic)

    alpha_phi_2_new <- mean_dispersion_output$alpha_phi_2
    b_new <- mean_dispersion_output$b

    # 9) -------- Update unique parameters --------
    unique_output_sim <- unique_parameters_update(mu_star_1_J = mu_star_1_J_new,
                                                  phi_star_1_J = phi_star_1_J_new,
                                                  mean_X_mu_phi = mean_X_unique_new,
                                                  tilde_s_mu_phi = tilde_s_unique_new,
                                                  J = J, G = G, Z = Z_new, b = b_new,
                                                  alpha_phi_2 = alpha_phi_2_new,
                                                  Beta = Beta_new, alpha_mu_2 = alpha_mu_2,
                                                  covariance = covariance_unique_new, Y = Y,
                                                  iter_num = iter, quadratic = quadratic,
                                                  MH.variance = MH.variance)

    mu_star_1_J_new <- unique_output_sim$mu_star_1_J_new
    phi_star_1_J_new <- unique_output_sim$phi_star_1_J_new
    tilde_s_unique_new <- unique_output_sim$tilde_s_mu_phi_new
    mean_X_unique_new <- unique_output_sim$mean_X_mu_phi_new
    covariance_unique_new <- unique_output_sim$covariance_new
    unique_count <- unique_count + unique_output_sim$accept_count
    acceptance_count_avg$unique_accept[iter-1] <- unique_count/((iter-1)*J*G)

    # 10) ------ Update capture efficiency -----------
    capture_output_sim <- capture_efficiencies_update(Beta = Beta_new, Y = Y, Z = Z_new,
                                                      mu_star_1_J = mu_star_1_J_new,
                                                      phi_star_1_J = phi_star_1_J_new,
                                                      a_d_beta = a_d_beta,
                                                      b_d_beta = b_d_beta, iter_num = iter,
                                                      M_2 = M_2_capture_new,
                                                      mean_X = mean_X_capture_new,
                                                      variance = variance_capture_new,
                                                      MH.variance = MH.variance)

    Beta_new <- capture_output_sim$Beta_new
    mean_X_capture_new <- capture_output_sim$mean_X_new
    M_2_capture_new <- capture_output_sim$M_2_new
    variance_capture_new <- capture_output_sim$variance_new
    Beta_count <- Beta_count+capture_output_sim$accept_count
    acceptance_count_avg$Beta_accept[iter-1] <- Beta_count/((iter-1)*(sum(C_d)))


    #-------------------------- Step 6: Save simulated values ------------------------
    if(update == TRUE){

      alpha_phi_2_output[output_index] <- alpha_phi_2_new
      b_output[[output_index]] <- as.vector(b_new)

      if(is.null(Z_fix)){
        Z_output[[output_index]] <- Z_new
      }

      P_C_J_D_output[[1]][[output_index]] <- P_C_J_D_new[[1]]
      P_C_J_D_output[[2]][[output_index]] <- P_C_J_D_new[[2]]

      P_output[[output_index]] <- P_new
      alpha_output[output_index] <- alpha_new
      alpha_0_output[output_index] <- alpha_0_new
      mu_star_1_J_output[[output_index]] <- mu_star_1_J_new
      phi_star_1_J_output[[output_index]] <- phi_star_1_J_new
      Beta_output[[output_index]] <- Beta_new

      Q_J_D_output[[output_index]] <- Q_J_D_new
      Xi_C_D_output[[output_index]] <- Xi_C_D_new

      t_star_J_D_output[[output_index]] <- t_star_J_D_new
      sigma_star_2_J_D_output[[output_index]] <- sigma_star_2_J_D_new

      U_C_J_D_output[[1]][[output_index]] <- U_C_J_D_new[[1]]
      U_C_J_D_output[[2]][[output_index]] <- U_C_J_D_new[[2]]

      r_J_output[[output_index]] <- as.vector(unname(r_J_new))
      s_2_output[output_index] <- s_2_new

      h_J_output[[output_index]] <- as.vector(unname(h_J_new))
      m_2_output[output_index] <- m_2_new

    }

  }

  ## Return the list
  my_list <- list('b_output' = b_output, 'alpha_phi_2_output' = alpha_phi_2_output, 'Z_output' = Z_output,
                  'P_C_J_D_output' = P_C_J_D_output,  'P_output' = P_output, 'alpha_output' = alpha_output,
                  'alpha_0_output' = alpha_0_output, 'mu_star_1_J_output' = mu_star_1_J_output,
                  'phi_star_1_J_output' = phi_star_1_J_output, 'Beta_output' = Beta_output,
                  'Q_J_D_output' = Q_J_D_output, 'Xi_C_D_output' = Xi_C_D_output, 't_star_J_D_output' = t_star_J_D_output,
                  'sigma_star_2_J_D_output' = sigma_star_2_J_D_output, 'U_C_J_D_output' = U_C_J_D_output,
                  'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output,
                  'acceptance_count_avg' = acceptance_count_avg,
                  'output_index' = output_index,
                  'alpha_mu_2' = alpha_mu_2, 'v_1' = v_1, 'v_2' = v_2, 'm_b' = m_b,
                  'a_d_beta' = a_d_beta, 'b_d_beta' = b_d_beta)

  close(pb)

  cat('\n')
  end_time <- Sys.time()
  print(paste('End:',end_time))
  cat('\n')
  diff_time <- difftime(end_time,start_time_mcmc)
  print(paste('MCMC running time:', round(diff_time, digits = 3),units(diff_time)))
  return(my_list)

}

