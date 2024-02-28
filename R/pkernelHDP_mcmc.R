#' Perform MCMC for covariate-dependent hierarchical Dirichlet process to cluster across time-series data
#'
#' @description
#' Implement Gibbs sampling and adaptive Metropolis-Hastings algorithm to draw each parameter sequentially.
#' A periodic kernel is applied to introduce dependence on the external covariate, i.e., time, and the
#' algorithm provides posterior samples for all parameters including cluster allocations. The function is also
#' used in the post-processing step to fix allocations to the optimal clustering, in order to infer cluster-specific parameters.
#'
#' @importFrom stats pgamma qgamma rexp
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom SciViews ln
#'
#' @param Y a list of two matrices of dimension \eqn{C_d \times G}.
#' The rows in each matrix correspond to observations and columns for features.
#' Both datasets should share the same set of features.
#' @param t a list of two vectors. Each vector is the external covariate (time) for individual dataset.
#' @param niter integer. Total number of MCMC iterations.
#' @param J integer. Truncation level, i.e., the number of clusters to be found. Only required when estimating \code{Z}.
#' @param burn_in the length of burn-in period during MCMC.
#' @param thinning the thinning applied after burn-in period.
#' @param empirical_z optional. Should be provided when \code{Z_fix} is not.
#' If \code{Z_fix} is not provided and \code{empirical_z=TRUE}, a simple Gaussian mixture model
#' is fitted using package \code{mclust}. If \code{empirical_z=FALSE}, clusters are randomly initialized.
#' @param Z_fix optional. Should be provided when \code{empirical_z} and \code{J} are not available. A list of vectors,
#' where each vector stores the allocations in one dataset. Typically used in the
#' post-processing step where allocations are fixed to the optimal clustering to infer cluster-specific parameters.
#' @param alpha_initial,alpha_0_initial initial values for concentration parameters \eqn{\alpha,\alpha_0}.
#' @param MH.variance additional variance added to the empirical covariance (variance) in adaptive Metropolis-Hastings.
#' @param target_accept the targeted average acceptance rate in adaptive Metropolis-Hastings (Algorithm 5).
#' @param mu_r,sigma_r mean and standard deviation of the hyper-prior (normal) for \eqn{r_j}.
#' @param eta_1,eta_2 shape and scale of the hyper-prior (inverse-gamma) for \eqn{s^2}.
#' @param mu_h,sigma_h mean and standard deviation of the hyper-prior (log-normal) for \eqn{h_j}.
#' @param kappa_1,kappa_2 shape and scale of the hyper-prior (inverse-gamma) for \eqn{m^2}.
#' @param auto.save logical. Whether intermediate results should be saved during MCMC sampling.
#' @param partial.save.name if \code{auto.save=TRUE}, the path and name to save the intermediate results.
#' @param save_frequency if \code{auto.save=TRUE}, save the intermediate results every \code{save_frequency} iterations.
#' @param save_ind optional. Used for consensus clustering, to save results at iteration index given by \code{save_ind} before burn-in.
#'
#' @usage pkernelHDP_mcmc(Y, t, J = NULL, niter, burn_in = 1000, thinning = 1,
#'     empirical_z = NULL, Z_fix = NULL, alpha_initial = 1, alpha_0_initial = 1,
#'     MH.variance = 0.01, target_accept = 0.234,
#'     mu_r = -2, sigma_r = 0.5, eta_1 = 5, eta_2 = 1, mu_h = -1,
#'     sigma_h = 0.5, kappa_1 = 26, kappa_2 = 1, auto.save = FALSE,
#'     partial.save.name = NULL, save_frequency = 100, save_ind = NULL)
#'
#' @return \code{pkernelHDP_mcmc} returns a list containing the following components:
#' \item{output_index}{total number of saved MCMC samples, taking into account of burn-in and thinning.}
#' \item{Z_output}{If \code{Z_fix} is not provided, a list of length \code{output_index}. Each element of
#' the list is a list of two vectors saving the allocations in individual dataset. Otherwise is \code{NULL}.}
#' \item{P_C_J_D_output}{a list of length two, each corresponding to a single dataset. Each element of the list is
#' a list of length \code{output_index}, saving the samples for time-dependent probabilities in a \eqn{C_d \times J} matrix .}
#' \item{P_output}{a list of length \code{output_index}. Each element of the list is a vector of length \eqn{J} for component probabilities \eqn{\mathbf{p}^J}.}
#' \item{alpha_output}{a vector of sampled concentration parameter \eqn{\alpha}.}
#' \item{alpha_0_output}{samples for concentration parameter \eqn{\alpha_0}. Similar to \code{alpha_output}.}
#' \item{L_1_J_output}{a list of length \code{output_index}.
#' Each element of the list is a list of \eqn{J} coefficient matrices \eqn{L_j^*}.}
#' \item{Sigma_1_J_output}{samples for covariance matrix \eqn{\Sigma_j^*}. Similar to \code{L_1_J_output}.}
#' \item{Q_J_D_output}{a list of length \code{output_index}. Each element of the list is a \eqn{J \times D} matrix for \eqn{q_{j,d}}.}
#' \item{Xi_C_D_output}{samples for latent variables. Similar to \code{Z_output}.}
#' \item{mu_J_D_output}{samples for kernel parameters (\eqn{\arg\max value}). Similar to \code{Q_J_D_output}.}
#' \item{lambda_J_D_output}{samples for kernel parameters (period). Similar to \code{Q_J_D_output}.}
#' \item{sigma_2_J_D_output}{samples for kernel parameters (variance). Similar to \code{Q_J_D_output}.}
#' \item{U_C_J_D_output}{samples for latent variables. Similar to \code{P_C_J_D_output}.}
#' \item{r_J_output}{a list of length \code{output_index}. Each element of the list is a vector of length \eqn{J} for \eqn{r_j}.}
#' \item{s_2_output}{a vector of sampled \eqn{s^2}.}
#' \item{h_J_output}{samples for \eqn{h_j}. Similar to \code{r_J_output}.}
#' \item{m_2_output}{samples for \eqn{m^2}. Similar to \code{s_2_output}.}
#' \item{acceptance_count_avg}{average acceptance probabilities over iterations for parameters drawn using adaptive Metroplis-Hastings.}
#' \item{L0, V0}{prior parameters for \eqn{L_j^*}.}
#' \item{Phi0, omega0}{prior parameters for \eqn{\Sigma_j^*}.}
#'
#' @export
pkernelHDP_mcmc <- function(Y, t, J = NULL, niter, burn_in = 1000, thinning = 1,
                            empirical_z = NULL, Z_fix = NULL,
                            alpha_initial = 1, alpha_0_initial = 1,
                            MH.variance = 0.01, target_accept = 0.234,
                            mu_r = -2, sigma_r = 0.5, eta_1 = 5, eta_2 = 1, mu_h = -1,
                            sigma_h = 0.5, kappa_1 = 26, kappa_2 = 1,
                            auto.save = FALSE, partial.save.name = NULL, save_frequency = 100,
                            save_ind = NULL){

  # Show time
  start_time <- Sys.time()
  print(paste('Start Initialization:', start_time))
  cat('\n')

  # Sample size after burn-in and thinning
  # apply burn-in first, then apply thinning to the rest samples
  # output_size <- floor((niter-burn_in)/thinning)

  # Define periodic kernel
  pkernel <- function(t_c_d, mu_j_d, lambda_j_d, sigma_2_j_d) {
    val1 <- -2/sigma_2_j_d
    val2 <- sin((t_c_d-mu_j_d)/lambda_j_d)

    return(exp(val1*val2^2))
  }

  # Number of datasets, and number of features
  D <- length(Y); G <- ncol(Y[[1]])

  # Check the length matches
  for(d in 1:D){
    if(length(t[[d]])!=nrow(Y[[d]])){
      stop(paste('length of t does not match nrow of Y for data',d))
    }
  }

  # First observation is not going to be clustered
  C_d <- rep(0,2); C_d[1] <- nrow(Y[[1]])-1; C_d[2] <- nrow(Y[[2]])-1

  # Prepare covariate matrix C_d[d] * (n_feature+1)
  X <- NULL
  for(d in 1:D){
    # First column is the intercept
    X[[d]] <- as.matrix(cbind(rep(1,C_d[d]),Y[[d]][-nrow(Y[[d]]),]))
  }

  # Remove first observation from Y and t
  Y[[1]] <- Y[[1]][-1,]; Y[[2]] <- Y[[2]][-1,]
  t[[1]] <- t[[1]][-1]; t[[2]] <- t[[2]][-1]

  if(is.null(J) & is.null(Z_fix)){
    stop('At least one of J or Z_fix should be provided!')
  }

  if(is.null(J)){
    J <- length(unique(unlist(Z_fix)))
  }

  #------------------------ Step 1: Prepare for outputs -----------------
  Z_output <- NULL

  P_C_J_D_output <- list(list(NULL),list(NULL))

  P_output <- NULL
  alpha_output <- c()
  alpha_0_output <- c()

  # Cluster-specific
  L_1_J_output <- NULL
  Sigma_1_J_output <- NULL

  Q_J_D_output <- NULL
  Xi_C_D_output <- NULL

  # Kernel parameters
  mu_J_D_output <- NULL
  sigma_2_J_D_output <- NULL
  lambda_J_D_output <- NULL

  U_C_J_D_output <- list(list(NULL),list(NULL))

  r_J_output <- NULL
  s_2_output <- c()
  h_J_output <- NULL
  m_2_output <- c()

  sd_mu_output <- NULL
  sd_lambda_output <- NULL
  sd_P_output <- c()
  sd_h_output <- NULL
  sd_m_2_output <- c()

  #----------------------- Step 2: Initial values in MCMC-----------

  # Use mclust for initialization of Z, lambda_J_D, sigma_2_J_D, r_J, s^2, h_j, m^2
  if(is.null(Z_fix)){
    if(empirical_z==TRUE){
      Y_all <- rbind(Y[[1]],Y[[2]])
      mc_cluster <- mclust::Mclust(Y_all, G=J, modelNames = 'VVV')$classification
      Z_initial <- NULL
      Z_initial[[1]] <- mc_cluster[1:C_d[1]]
      Z_initial[[2]] <- mc_cluster[(C_d[1]+1):(C_d[2]+C_d[1])]
    }else{
      Z_initial <- NULL
      Z_initial[[1]] <- sample(1:J, size = C_d[[1]],replace = TRUE,prob = rep(1/J,J))
      Z_initial[[2]] <- sample(1:J, size = C_d[[2]],replace = TRUE,prob = rep(1/J,J))
    }
  }else{
    Z_initial <- Z_fix
  }

  # Assume a common lambda within each dataset; period T=pi*lambda; t \in (0,1)
  lambda_J_D_initial <- matrix(NA, nrow = J, ncol = D)
  # Sample the number of complete periods on (0,1), different for each dataset
  n_period <- sample(2:5,size=2,replace = FALSE)
  for(d in 1:D){
    for(j in 1:J){
      # T = 1/n_period
      lambda_J_D_initial[j,d] <- 1/n_period[d]/pi
    }
  }

  # Initials for r_J as an empirical mean from log(lambda_J_D_initial)
  r_J_initial <- apply(ln(lambda_J_D_initial), 1, mean)

  # Initials for s^2
  s_2_initial <- mean(apply(ln(lambda_J_D_initial), 1, var))

  # Initialize mu_J_D from the prior: Unif(-pi*lambda/2, pi*lambda/2)
  mu_J_D_initial <- matrix(NA, nrow = J, ncol = D)

  for(d in 1:D){
    for(j in 1:J){
      mu_J_D_initial[j,d] <- runif(1, min=-pi*lambda_J_D_initial[j,d]/2, max=pi*lambda_J_D_initial[j,d]/2)
    }
  }

  # Initials for sigma_2_J_D
  sigma_2_J_D_initial <- matrix(NA, nrow = J, ncol = D)

  # Variance of sin((t-mu)/lambda) within each cluster of each dataset
  # For empty clusters, return var of sin((t-mu)/lambda) in the whole dataset
  temp1 <- unlist(lapply(1:J, function(j) {
    ind <- which(Z_initial[[1]]==j)

    if(length(ind)>1){
      val <- sin((t[[1]][ind]-mu_J_D_initial[j,1])/lambda_J_D_initial[j,1])
      return(var(val))
    }else{
      lambda_empty <- 1/n_period[1] / pi
      mu_empty <- runif(1,min=-pi*lambda_empty/2,max=pi*lambda_empty/2)
      val <- sin((t[[1]]-mu_empty)/lambda_empty)

      return(var(val))
    }
  }))

  temp2 <- unlist(lapply(1:J, function(j) {
    ind <- which(Z_initial[[2]]==j)

    if(length(ind)>1){
      val <- sin((t[[2]][ind]-mu_J_D_initial[j,2])/lambda_J_D_initial[j,2])
      return(var(val))
    }else{
      lambda_empty <- 1/n_period[2] / pi
      mu_empty <- runif(1,min=-pi*lambda_empty/2,max=pi*lambda_empty/2)
      val <- sin((t[[2]]-mu_empty)/lambda_empty)

      return(var(val))
    }
  }))

  # -2/sigma_2=-1/(2*var), so sigma_2=4*var
  sigma_2_J_D_initial[,1] <- 4*temp1
  sigma_2_J_D_initial[,2] <- 4*temp2

  # Initials for h_J from sigma_2_J_D
  h_J_initial <- apply(sigma_2_J_D_initial, 1, mean)

  # Initials for m^2
  m_2_initial <- mean(apply(sigma_2_J_D_initial, 1, var))

  # Component probabilities p_j, add 1 to avoid zero p if cluster is empty
  P_initial <- sapply(1:J, function(j) mean(unlist(Z_initial)==j))
  if(any(P_initial==0)){
    # add 1 to avoid zero p if cluster is empty
    P_initial <- sapply(1:J, function(j) sum(unlist(Z_initial)==j)+1)/(sum(C_d)+J)
  }

  # Dataset-specific vector q_j_d, and compute component probabilities p_j_d
  Q_J_D_initial <- matrix(NA,nrow = J, ncol = D)
  Q_J_D_initial[,1] <- sapply(1:J,function(j) sum(Z_initial[[1]]==j))+1
  Q_J_D_initial[,2] <- sapply(1:J,function(j) sum(Z_initial[[2]]==j))+1


  P_C_J_D_initial <- NULL
  P_C_J_D_initial[[1]] <- t(sapply(t[[1]],function(time) {
    Q_J_D_initial[,1]*pkernel(time,mu_J_D_initial[,1],lambda_J_D_initial[,1],sigma_2_J_D_initial[,1])/(
      sum(Q_J_D_initial[,1]*pkernel(time,mu_J_D_initial[,1],lambda_J_D_initial[,1],sigma_2_J_D_initial[,1])))
  }))
  P_C_J_D_initial[[2]] <- t(sapply(t[[2]],function(time) {
    Q_J_D_initial[,2]*pkernel(time,mu_J_D_initial[,2],lambda_J_D_initial[,2],sigma_2_J_D_initial[,2])/(
      sum(Q_J_D_initial[,2]*pkernel(time,mu_J_D_initial[,2],lambda_J_D_initial[,2],sigma_2_J_D_initial[,2])))
  }))


  # Set up prior for L_j based on MLE from the whole dataset
  Y_all <- rbind(Y[[1]],Y[[2]]); X_all <- rbind(X[[1]],X[[2]])
  # Multivariate linear regression
  lm_fit <- lm(Y_all~X_all-1)
  # dim = (G+1) * G
  L0 <- unname(lm_fit$coefficients)
  V0 <- diag(100,nrow=nrow(L0))
  # Prior for Sigma_j
  Phi0 <- var(lm_fit$residuals)/5^(2/G)
  # Prior mean is Phi0/(omega0-G-1)
  omega0 <- G+2

  # Initial L_j and Sigma_j drawn from full conditional
  unique_output <- unique_params_update(Y, X, J, Z_initial, L0, V0, Phi0, omega0)
  L_1_J_initial <- unique_output$L_new
  Sigma_1_J_initial <- unique_output$Sigma_new


  #--------------------------------------------------------------------------------

  # Acceptance probability
  acceptance_count_avg <- data.frame(P_accept = rep(0,niter), alpha_accept = rep(0,niter),
                                     alpha_0_accept = rep(0,niter), mu_accept=rep(0,niter),
                                     lambda_accept=rep(0,niter), h_accept = rep(0, niter), m_2_accept = rep(0, niter))

  #----------------------- Step 3: Set the initial_values as new values ----------------------------
  L_1_J_new <- L_1_J_initial
  Sigma_1_J_new <- Sigma_1_J_initial

  mu_J_D_new <- mu_J_D_initial
  sigma_2_J_D_new <- sigma_2_J_D_initial
  lambda_J_D_new <- lambda_J_D_initial

  r_J_new <- r_J_initial
  s_2_new <- s_2_initial
  h_J_new <- h_J_initial
  m_2_new <- m_2_initial

  Q_J_D_new <- Q_J_D_initial
  P_C_J_D_new <- P_C_J_D_initial

  Xi_C_D_new <- Xi_C_D_update_pkernel(Q_J_D = Q_J_D_new, C_d = C_d, t = t, mu_J_D = mu_J_D_new, lambda_J_D=lambda_J_D_new,
                              sigma_2_J_D = sigma_2_J_D_new)

  U_C_J_D_new <- U_C_J_D_update_pkernel(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C_d = C_d, t = t,
                                mu_J_D = mu_J_D_new, lambda_J_D=lambda_J_D_new, sigma_2_J_D = sigma_2_J_D_new,
                                pkernel = pkernel)

  Z_new <- Z_initial
  P_new <- P_initial
  alpha_new <- alpha_initial
  alpha_0_new <- alpha_0_initial

  # Total count of acceptance
  P_count <- 0; alpha_count <- 0; alpha_0_count <- 0; mu_count <- 0; lambda_count <- 0
  h_count <- 0; m_2_count <- 0

  #----------------------- Step 4: Prepare for the adaptive covariance update -----------------------------

  # 0) For h_J and m_2
  sd_h_new <- rep(0.01,J)
  sd_m_2_new <- 0.01

  # 1) For mu_J_D, matrix
  sd_mu_new <- matrix(0.01, nrow=J, ncol=D)

  # 2) For lambda_J_D, matrix
  sd_lambda_new <- matrix(0.01, nrow=J, ncol=D)

  # 3) For Component probabilities
  sd_P_new <- 0.001
  mean_X_component_new <- ln(matrix(P_new[1:(J-1)]/P_new[J], nrow = 1))
  tilde_s_component_new <- t(mean_X_component_new)%*%mean_X_component_new
  covariance_component_new <- matrix(0, nrow = J-1, ncol = J-1)

  # 4) For alpha
  mean_X_alpha_new <- ln(alpha_new)
  M_2_alpha_new <- 0
  variance_alpha_new <- 0

  # 5) For alpha_0
  mean_X_alpha_0_new <- ln(alpha_0_new)
  M_2_alpha_0_new <- 0
  variance_alpha_0_new <- 0

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
    # In the case of thinning, if we require to save some samples that are not divisible by thinning (for consesnsus clustering) or before burn-in
    if(is.null(save_ind)) {
      criterion1 <- FALSE
    }else{
      criterion1 <- any((iter-1)==save_ind)
    }

    criterion2 <- (iter-1 > burn_in & (iter-1-burn_in)%%thinning == 0)

    if(criterion1 | criterion2){
      output_index <- output_index + 1
      update <- TRUE
    }else{
      update <- FALSE
    }


    # 0) Update latent/auxiliary variables Xi_C_D
    Xi_C_D_new <- Xi_C_D_update_pkernel(Q_J_D = Q_J_D_new, C_d = C_d, t = t, mu_J_D = mu_J_D_new, lambda_J_D=lambda_J_D_new,
                                sigma_2_J_D = sigma_2_J_D_new)

    # 1) Update dataset-specific vector q_j_d
    Q_J_D_new <- Q_J_D_update_pkernel(Z = Z_new, alpha = alpha_new, P = P_new, Xi_C_D = Xi_C_D_new,
                              t = t, mu_J_D = mu_J_D_new, lambda_J_D=lambda_J_D_new,
                              sigma_2_J_D = sigma_2_J_D_new)

    # 2) Update the allocation variable, if Z_fix is not provided
    if(is.null(Z_fix)){
      Z_new <- allocation_variables_update_pkernel(Y = Y, X = X, t = t, L_1_J = L_1_J_new, Sigma_1_J = Sigma_1_J_new,
                                           Q_J_D = Q_J_D_new, mu_J_D = mu_J_D_new,lambda_J_D=lambda_J_D_new,
                                           sigma_2_J_D = sigma_2_J_D_new)
    }

    # 3) Update cluster-specific parameters
    unique_params_output_sim <- unique_params_update(Y = Y, X = X, J = J, Z = Z_new, L0 = L0, V0 = V0,
                                                     Phi0 = Phi0, omega0 = omega0)

    L_1_J_new <- unique_params_output_sim$L_new
    Sigma_1_J_new <- unique_params_output_sim$Sigma_new


    # 4) Update kernel parameters
    # 4-1) mu_J_D
    mu_J_D_output_sim <- mu_J_D_update(mu_J_D = mu_J_D_new, Z = Z_new, t = t, lambda_J_D=lambda_J_D_new,
                                       sigma_2_J_D = sigma_2_J_D_new, Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new,
                                       pkernel = pkernel, sd_mu = sd_mu_new, iter_num = iter, target_accept)

    mu_J_D_new <- mu_J_D_output_sim$mu_J_D_new
    mu_count <- mu_count + mu_J_D_output_sim$accept_count
    acceptance_count_avg$mu_accept[iter-1] <- mu_count/((iter-1)*J*D)
    sd_mu_new <- mu_J_D_output_sim$sd_mu_new

    sd_mu_output[[iter-1]] <- sd_mu_new

    # 4-2) lambda_J_D
    lambda_J_D_output_sim <- lambda_J_D_update(lambda_J_D = lambda_J_D_new, r_J = r_J_new, s_2 = s_2_new, Z = Z_new, t = t,
                                               mu_J_D=mu_J_D_new, sigma_2_J_D = sigma_2_J_D_new, Xi_C_D = Xi_C_D_new,
                                               Q_J_D = Q_J_D_new, pkernel = pkernel, sd_lambda = sd_lambda_new,
                                               iter_num = iter, target_accept)

    lambda_J_D_new <- lambda_J_D_output_sim$lambda_J_D_new
    lambda_count <- lambda_count + lambda_J_D_output_sim$accept_count
    acceptance_count_avg$lambda_accept[iter-1] <- lambda_count/((iter-1)*J*D)
    sd_lambda_new <- lambda_J_D_output_sim$sd_lambda_new

    sd_lambda_output[[iter-1]] <- sd_lambda_new

    # 4-3) U_C_J_D
    U_C_J_D_new <- U_C_J_D_update_pkernel(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C_d = C_d, t = t,
                                  mu_J_D = mu_J_D_new, lambda_J_D=lambda_J_D_new, sigma_2_J_D = sigma_2_J_D_new,
                                  pkernel = pkernel)

    # 4-4) sigma_2_J_D
    sigma_2_J_D_new <- sigma_2_J_D_update(h_J = h_J_new, m_2 = m_2_new, Z = Z_new,
                                          t = t, mu_J_D = mu_J_D_new, lambda_J_D = lambda_J_D_new, U_C_J_D = U_C_J_D_new,
                                          Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, iter_num = iter)


    # Update P_C_J_D
    P_C_J_D_new <- NULL
    P_C_J_D_new[[1]] <- t(sapply(t[[1]],function(time) {
      Q_J_D_new[,1]*pkernel(time,mu_J_D_new[,1],lambda_J_D_new[,1],sigma_2_J_D_new[,1])/(
        sum(Q_J_D_new[,1]*pkernel(time,mu_J_D_new[,1],lambda_J_D_new[,1],sigma_2_J_D_new[,1])))
    }))
    P_C_J_D_new[[2]] <- t(sapply(t[[2]],function(time) {
      Q_J_D_new[,2]*pkernel(time,mu_J_D_new[,2],lambda_J_D_new[,2],sigma_2_J_D_new[,2])/(
        sum(Q_J_D_new[,2]*pkernel(time,mu_J_D_new[,2],lambda_J_D_new[,2],sigma_2_J_D_new[,2])))
    }))

    # 5) Update hyper parameters in priors for lambda, sigma_2
    # 5-1) r_J
    r_J_new <- r_J_update_pkernel(lambda_J_D = lambda_J_D_new, mu_r = mu_r, sigma_r = sigma_r, s_2 = s_2_new)

    # 5-2) s_2
    s_2_new <- s_2_update_pkernel(lambda_J_D = lambda_J_D_new, eta_1 = eta_1, eta_2 = eta_2, r_J = r_J_new)

    # 5-3) h_J
    h_J_output_sim <- h_J_update_pkernel(h_J = h_J_new, sigma_2_J_D = sigma_2_J_D_new, mu_h = mu_h, sigma_h = sigma_h, m_2 = m_2_new,
                                 sd_h = sd_h_new, iter_num = iter, target_accept)

    h_J_new <- h_J_output_sim$h_J_new
    h_count <- h_count + h_J_output_sim$accept_count
    acceptance_count_avg$h_accept[iter-1] <- h_count/((iter-1)*J)
    sd_h_new <- h_J_output_sim$sd_h_new

    sd_h_output[[iter-1]] <- sd_h_new

    # 5-4) m_2
    m_2_output_sim <- m_2_update_pkernel(m_2 = m_2_new, sigma_2_J_D = sigma_2_J_D_new, kappa_1 = kappa_1, kappa_2 = kappa_2, h_J = h_J_new,
                                 sd_m_2 = sd_m_2_new, iter_num = iter, target_accept)

    m_2_new <- m_2_output_sim$m_2_new
    m_2_count <- m_2_count + m_2_output_sim$accept
    acceptance_count_avg$m_2_accept[iter-1] <- m_2_count/(iter-1)
    sd_m_2_new <- m_2_output_sim$sd_m_2_new

    sd_m_2_output[iter-1] <- sd_m_2_new


    # 6) Update the component probabilities P
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

    # 7) Update alpha
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

    # 8) update alpha_0
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

    #-------------------------- Step 6: Update simulated values ------------------------

    if(update == TRUE){

      if(is.null(Z_fix)){
        Z_output[[output_index]] <- Z_new
      }

      P_C_J_D_output[[1]][[output_index]] <- P_C_J_D_new[[1]]
      P_C_J_D_output[[2]][[output_index]] <- P_C_J_D_new[[2]]

      P_output[[output_index]] <- P_new
      alpha_output[output_index] <- alpha_new
      alpha_0_output[output_index] <- alpha_0_new

      L_1_J_output[[output_index]] <- L_1_J_new
      Sigma_1_J_output[[output_index]] <- Sigma_1_J_new

      Q_J_D_output[[output_index]] <- Q_J_D_new
      Xi_C_D_output[[output_index]] <- Xi_C_D_new

      mu_J_D_output[[output_index]] <- mu_J_D_new
      sigma_2_J_D_output[[output_index]] <- sigma_2_J_D_new
      lambda_J_D_output[[output_index]] <- lambda_J_D_new

      U_C_J_D_output[[1]][[output_index]] <- U_C_J_D_new[[1]]
      U_C_J_D_output[[2]][[output_index]] <- U_C_J_D_new[[2]]

      r_J_output[[output_index]] <- as.vector(unname(r_J_new))
      s_2_output[output_index] <- s_2_new

      h_J_output[[output_index]] <- as.vector(unname(h_J_new))
      m_2_output[output_index] <- m_2_new
    }

    if((iter-1) %% save_frequency == 0 && auto.save == TRUE){
      my_list <- list('Z_output' = Z_output,'P_C_J_D_output' = P_C_J_D_output,  'P_output' = P_output,
                      'alpha_output' = alpha_output,'alpha_0_output' = alpha_0_output, 'L_1_J_output' = L_1_J_output,
                      'Sigma_1_J_output' = Sigma_1_J_output,
                      'Q_J_D_output' = Q_J_D_output, 'Xi_C_D_output' = Xi_C_D_output, 'mu_J_D_output' = mu_J_D_output,
                      'lambda_J_D_output' = lambda_J_D_output, 'sigma_2_J_D_output' = sigma_2_J_D_output,
                      'U_C_J_D_output' = U_C_J_D_output,
                      'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output,
                      'acceptance_count_avg' = acceptance_count_avg,
                      'output_index' = output_index,
                      'L0' = L0, 'V0' = V0, 'Phi0' = Phi0, 'omega0' = omega0)
      save(my_list, file=partial.save.name)
    }

  } # End for loop

  ## Return the list
  my_list <- list('Z_output' = Z_output,'P_C_J_D_output' = P_C_J_D_output,  'P_output' = P_output,
                  'alpha_output' = alpha_output,'alpha_0_output' = alpha_0_output, 'L_1_J_output' = L_1_J_output,
                  'Sigma_1_J_output' = Sigma_1_J_output,
                  'Q_J_D_output' = Q_J_D_output, 'Xi_C_D_output' = Xi_C_D_output, 'mu_J_D_output' = mu_J_D_output,
                  'lambda_J_D_output' = lambda_J_D_output, 'sigma_2_J_D_output' = sigma_2_J_D_output,
                  'U_C_J_D_output' = U_C_J_D_output,
                  'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output,
                  'acceptance_count_avg' = acceptance_count_avg,
                  'output_index' = output_index,
                  'L0' = L0, 'V0' = V0, 'Phi0' = Phi0, 'omega0' = omega0)

  close(pb)

  cat('\n')
  end_time <- Sys.time()
  print(paste('End:',end_time))
  cat('\n')
  diff_time <- difftime(end_time,start_time_mcmc)
  print(paste('MCMC running time:', round(diff_time, digits = 3),units(diff_time)))
  return(my_list)

}

