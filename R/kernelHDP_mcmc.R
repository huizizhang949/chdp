#' Perform MCMC for covariate-dependent hierarchical Dirichlet process on single-cell data
#'
#' @import stats
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom SciViews ln
#' @param Y a
#' @param t a
#' @param niter a
#' @param J a
#' @param burn_in a
#' @param thinning a
#' @param empirical a
#' @param empirical_z a
#' @param Z_fix a
#' @param b_initial a
#' @param alpha_initial a
#' @param alpha_0_initial a
#' @param quadratic a
#' @param MH.variance a
#' @param mu_r a
#' @param sigma_r a
#' @param eta_1 a
#' @param eta_2 a
#' @param mu_h a
#' @param sigma_h a
#' @param kappa_1 a
#' @param kappa_2 a
#' @param beta.mean a
#' @param alpha_mu_2 a
#' @param partial.save.name a
#' @param save_frequency a
#' @param auto.save a
#' @param partial_pca a
#' @param save_ind a
#'
#' @return MCMC output for all parameters
#' @export
#'
kernelHDP_mcmc <- function(Y, t, niter, J, burn_in = 1000, thinning = 5,
                           empirical = TRUE, empirical_z = NULL, Z_fix = NULL,
                           b_initial = NULL, alpha_initial = 1, alpha_0_initial = 1,
                           quadratic=FALSE, MH.variance = 0.01,
                           mu_r = 0.5, sigma_r = 0.5, eta_1 = 5, eta_2 = 1, mu_h = -5,
                           sigma_h = 0.5, kappa_1 = 5, kappa_2 = 1,
                           beta.mean = 0.06, alpha_mu_2 = NULL,
                           partial.save.name = NULL, save_frequency = 100, auto.save = FALSE,
                           partial_pca = FALSE, save_ind = NULL){

  # Show time
  start_time <- Sys.time()
  print(paste('Start Initialization:', start_time))
  cat('\n')

  # Add cell names if inputted dataset has no cell names
  if(is.null(colnames(Y[[1]]))){
    colnames(Y[[1]]) <- paste0('d1-',1:ncol(Y[[1]]))
  }

  if(is.null(colnames(Y[[2]]))){
    colnames(Y[[2]]) <- paste0('d2-',1:ncol(Y[[2]]))
  }

  # Sample size after burn-in and thinning
  # Apply burn-in first, then apply thinning to the rest samples
  # output_size <- floor((niter-burn_in)/thinning)

  # Define Gaussian kernel (Radial basis function)
  rbf <- function(t_c_d, t_j_d, sigma_j_d_2) {exp(-(t_c_d-t_j_d)^2/(2*sigma_j_d_2))}

  D <- length(Y); G <- nrow(Y[[1]])

  # Data size
  C_d <- rep(0,2); C_d[1] <- ncol(Y[[1]]); C_d[2] <- ncol(Y[[2]])

  # --------------- Step 0: Functions for each step of the Gibb sampling ----------

  ##--------------- Simulation of latent Xi ----------------

  Xi_C_D_update <- function(Q_J_D, C_d, t, t_star_J_D, sigma_star_2_J_D){
    J <- nrow(Q_J_D)
    D <- ncol(Q_J_D)

    # Set up the list to save updated values
    Xi <- NULL
    Xi[[1]] <- rep(NA, C_d[1])
    Xi[[2]] <- rep(NA, C_d[2])

    # Draw from Gamma (full conditional distribution)
    for (d in 1:D) {

      rates <- vapply(1:C_d[d], function(c) {

        log_rbf_value <- -(t[[d]][c]-t_star_J_D[,d])^2/2/sigma_star_2_J_D[,d]
        log_temp <- ln(Q_J_D[,d])+log_rbf_value
        log_K <- max(log_temp)

        return(exp(log_K)*sum(exp(log_temp-log_K)))
      }, FUN.VALUE = numeric(1))

      Xi[[d]] <- rgamma(C_d[d], shape = 1, rate = rates)
    }

    return(Xi)
  }

  ##--------------- Simulation of dataset-specific vector q_j_d ----------------

  Q_J_D_update <- function(Z, alpha, P, Xi_C_D, t, t_star_J_D, sigma_star_2_J_D){
    J <- nrow(t_star_J_D)
    D <- ncol(t_star_J_D)

    # Set up the matrix to save updated values
    Q <- matrix(NA, nrow = J, ncol = D)

    for (d in 1:D) {
      shape_rate_param <- vapply(1:J, function(j) {
        shape <- sum(Z[[d]]==j)+alpha*P[j]

        log_rbf_value <- -(t[[d]]-t_star_J_D[j,d])^2/2/sigma_star_2_J_D[j,d]
        log_temp <- ln(Xi_C_D[[d]])+log_rbf_value
        log_K <- max(log_temp)
        rate <- 1+exp(log_K)*sum(exp(log_temp-log_K))

        return(c(shape, rate))
      }, FUN.VALUE = numeric(2))

      Q[,d] <- rgamma(J, shape = shape_rate_param[1,], rate = shape_rate_param[2,])

      if(any(Q[,d]==0)) {
        Q[Q[,d]==0,d] <- shape_rate_param[1,Q[,d]==0]/shape_rate_param[2,Q[,d]==0]
      }
    }

    return(Q)
  }

  ##--------------- Simulation of allocations ----------------

  allocation_variables_update <- function(Y, t, mu_star_1_J, phi_star_1_J, Beta, Q_J_D,
                                          t_star_J_D, sigma_star_2_J_D){

    D <- length(Y)
    C_d <- unlist(lapply(Y, ncol))
    J <- nrow(mu_star_1_J)

    # Set up the list to save updated values
    Z <- NULL
    for (d in 1:D){
      Z[[d]] <- rep(0,C_d[d])
    }

    for (d in 1:D){

      loop.result <- vapply(1:C_d[d], function(c) {

        ## Set log probability
        LP <- vapply(1:J, function(j) {

          sum(dnbinom(Y[[d]][,c],mu=mu_star_1_J[j,]*Beta[[d]][c], size = phi_star_1_J[j,], log = TRUE)) +
            ln(Q_J_D[j,d]) - (t[[d]][c]-t_star_J_D[j,d])^2/2/sigma_star_2_J_D[j,d]

        }, FUN.VALUE = numeric(1))

        ## Compute the normalizing constant
        nc <- -max(LP)
        P <- exp(LP+nc)/sum(exp(LP+nc)) # LP is a vector of length J
        Z <- sample(1:J, 1, prob=P)

        return(Z)

      }, FUN.VALUE = numeric(1))

      Z[[d]] <- loop.result
    }

    return(Z)
  }

  ##----------- Simulation of latent variables U_C_J_D for updating kernel parameters----------

  U_C_J_D_update <- function(Xi_C_D, Q_J_D, C_d, t, t_star_J_D, sigma_star_2_J_D, rbf){
    J <- nrow(Q_J_D)
    D <- ncol(Q_J_D)

    # Set up the list to save updated values
    U_C_J_D <- NULL
    U_C_J_D[[1]] <- matrix(NA, nrow = C_d[1], ncol = J)
    U_C_J_D[[2]] <- matrix(NA, nrow = C_d[2], ncol = J)

    for (d in 1:D) {
      for (j in 1:J){
        K_C_J_D <- exp(-Xi_C_D[[d]]*Q_J_D[j,d]*rbf(t[[d]],t_star_J_D[j,d],sigma_star_2_J_D[j,d]))

        U_C_J_D[[d]][,j] <- runif(C_d[d], min = 0, max = K_C_J_D)
      }
    }

    return(U_C_J_D)

  }

  ##----------- Simulation of kernel parameters t_star_J_D-------------

  t_star_J_D_update <- function(r_J, s_2, Z, t, C_d, sigma_star_2_J_D, U_C_J_D, Xi_C_D, Q_J_D){
    J <- nrow(Q_J_D)
    D <- ncol(Q_J_D)

    # Set up the matrix to save updated values
    t_star <- matrix(NA, nrow = J, ncol = D)

    # Find the truncation region

    for (d in 1:D) {
      for (j in 1:J){
        # Find which cells to truncate regions
        ind <- c(1:C_d[d])[-ln(U_C_J_D[[d]][,j])<Xi_C_D[[d]]*Q_J_D[j,d]]
        if(length(ind)==0) {
          truncate <- FALSE
        }else{
          # If truncate, find truncated regions
          truncate <- TRUE
          i_all_complement <- lapply(ind, function(c) {
            temp <- -2*sigma_star_2_J_D[j,d]*(ln(-ln(U_C_J_D[[d]][c,j]))-ln(Xi_C_D[[d]][c])-ln(Q_J_D[j,d]))

            i_complement <- sets::interval(t[[d]][c]-sqrt(temp),t[[d]][c]+sqrt(temp),'[]')
            return(i_complement)
          })

          i_all <- sets::interval_complement(sets::interval_union(i_all_complement))
        }

        # Parameters in posterior distribution (Normal)
        ind1 <- (Z[[d]]==j)
        sum_t_j_d <- sum(t[[d]][ind1])
        N_j_d <- sum(ind1)
        r_j_hat <- (r_J[j]*sigma_star_2_J_D[j,d]+s_2*sum_t_j_d)/(sigma_star_2_J_D[j,d]+N_j_d*s_2)
        s_2_hat <- s_2*sigma_star_2_J_D[j,d]/(sigma_star_2_J_D[j,d]+N_j_d*s_2)

        if(!truncate){
          # No truncation, draw from the Normal posterior directly
          t_star[j,d] <- rnorm(1,mean = r_j_hat, sd = sqrt(s_2_hat))
        }else{
          N_interval <- length(i_all)
          log_p <- c()

          # Compute the probability of lying in each interval
          for (l in 1:N_interval) {
            i <- i_all[[l]]
            lower <- as.numeric(unlist(i)[1])
            upper <- as.numeric(unlist(i)[2])

            if(is.infinite(lower)){
              log_p <- c(log_p,pnorm(upper,mean = r_j_hat, sd = sqrt(s_2_hat),log.p = TRUE))
            }else if(is.infinite(upper)){
              log_p <- c(log_p,pnorm(lower,mean = r_j_hat, sd = sqrt(s_2_hat),log.p = TRUE, lower.tail = FALSE))
            }else{
              log_p <- c(log_p,ln(pnorm(upper,mean = r_j_hat, sd = sqrt(s_2_hat))-
                                    pnorm(lower,mean = r_j_hat, sd = sqrt(s_2_hat))))
            }
          }

          log_K <- -max(log_p)

          # Select one truncated region (one interval)
          ind2 <- sample(1:N_interval, size = 1, prob = exp(log_p+log_K)) #breakpoint
          i_chosen <- i_all[[ind2]]
          t_star[j,d] <- truncnorm::rtruncnorm(1,a=as.numeric(unlist(i_chosen)[1]),
                                    b=as.numeric(unlist(i_chosen)[2]),
                                    mean = r_j_hat, sd = sqrt(s_2_hat))
        }
      }
    }

    return(t_star)
  }

  ##----------- Simulation of kernel parameters sigma_star_2_J_D-------------

  # log(full conditional distribution)
  sigma_star_2_log_prob <- function(sigma_star_2_j_d, m_2, h_j, t_star_j_d, t_c_d_sub){
    lprod <- -ln(sigma_star_2_j_d)-(ln(sigma_star_2_j_d)-h_j)^2/2/m_2-sum((t_c_d_sub-t_star_j_d)^2)/2/sigma_star_2_j_d
    return(lprod)
  }

  sigma_star_2_J_D_update <- function(sigma_star_2_J_D_old,h_J, m_2, Z, t, C_d, t_star_J_D, U_C_J_D, Xi_C_D, Q_J_D,
                                      X_mean, M_2, variance, iter_num, MH.variance){

    J <- nrow(Q_J_D)
    D <- ncol(Q_J_D)

    # All are JxD matrices
    X_mean_old <- X_mean
    M_2_old <- M_2
    variance_old <- variance

    # Save updated covariance matrices
    X_new <- matrix(NA, nrow=J, ncol = D)
    X_mean_new <- matrix(NA, nrow=J, ncol = D)
    M_2_new <- matrix(NA, nrow=J, ncol = D)
    variance_new <- matrix(NA, nrow=J, ncol = D)

    # Save updated sigma_star_2_J_D
    sigma_star_2_J_D_new <- matrix(NA, nrow=J, ncol = D)

    # Iteration index
    n <- iter_num

    # Total count of accepted proposals
    accept_count <- 0

    for (d in 1:D) {
      for (j in 1:J) {
        # Compute truncation region if needed
        # Find which cells to truncate regions
        ind <- c(1:C_d[d])[-ln(U_C_J_D[[d]][,j])<Xi_C_D[[d]]*Q_J_D[j,d]]
        if(length(ind)==0) {
          truncate <- FALSE
        }else{
          truncate <- TRUE
          uppers <- -(t[[d]]-t_star_J_D[j,d])^2/2/(ln(-ln(U_C_J_D[[d]][,j]))-ln(Xi_C_D[[d]])-ln(Q_J_D[j,d]))
          upper <- min(uppers[ind])
          if(upper<=sigma_star_2_J_D_old[j,d]){
            print(upper)
          }
        }

        # For empty clusters in dataset d, sample from the log-normal prior (may be truncated), always accept
        if(sum(Z[[d]]==j)==0) {
          if(!truncate){
            # No truncation
            sigma_star_2_J_D_new[j,d] <- rlnorm(1, meanlog = h_J[j], sdlog = sqrt(m_2))
            X_new[j,d] <- ln(sigma_star_2_J_D_new[j,d])

            accept_count <- accept_count+1
          }else{
            # Truncate
            # X = -log(1/sigma_star - 1/upper)
            p.upper <- plnorm(upper, meanlog = h_J[j], sdlog = sqrt(m_2))
            u <- runif(1, 0, p.upper)
            sigma_star_2_J_D_new[j,d] <- qlnorm(u, meanlog = h_J[j], sdlog = sqrt(m_2))

            X_new[j,d] <- -ln(1/sigma_star_2_J_D_new[j,d] - 1/upper)

            accept_count <- accept_count+1
          }
        }else{
          # Non-empty clusters, use adaptive MH
          if(!truncate){
            # No truncation, perform adaptive MH directly, propose sigma_star_j_d from log-normal, X=log(sigma_star_2)
            if(n <= 100){
              X_new[j,d] <- rnorm(n = 1, mean = ln(sigma_star_2_J_D_old[j,d]), sd = 0.1)
            }else{
              X_new[j,d] <- rnorm(n = 1, mean = ln(sigma_star_2_J_D_old[j,d]), sd = sqrt(2.4^2*variance_old[j,d]+2.4^2*MH.variance))
            }

            # Transform the new value of X back to new value of sigma_star_2_j_d
            sigma_star_2_J_D_new[j,d] <- exp(X_new[j,d])
            if(is.infinite(sigma_star_2_J_D_new[j,d])) {
              print(paste('No truncate:','cluster',j,'dataset',d,": sigma==Inf"))
              sigma_star_2_J_D_new[j,d] <- .Machine$double.xmax
            }
            if(sigma_star_2_J_D_new[j,d]==0) {
              print(paste('No truncate:','cluster',j,'dataset',d,": sigma==0"))
              sigma_star_2_J_D_new[j,d] <- .Machine$double.eps
            }

            # Subset of t such that Z[[d]][c]=j
            t_c_d_sub <- t[[d]][Z[[d]]==j]

            # Compute log acceptance probability
            log_acceptance <- sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_new[j,d], m_2 = m_2,
                                                    h_j = h_J[j], t_star_j_d = t_star_J_D[j,d],
                                                    t_c_d_sub = t_c_d_sub) -
              sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_old[j,d], m_2,
                                    h_J[j], t_star_J_D[j,d], t_c_d_sub) -
              ln(sigma_star_2_J_D_old[j,d]) + ln(sigma_star_2_J_D_new[j,d])
            acceptance_sigma <- exp(log_acceptance)
            acceptance_sigma <- min(1,acceptance_sigma)

            # Update sigma_star_2_J_D_new, X_new, X_mean_new, variance_new
            outcome <- rbinom(n = 1, size = 1, prob = acceptance_sigma)
            if(is.na(outcome) == TRUE | outcome == 0){
              X_new[j,d] <- ln(sigma_star_2_J_D_old[j,d])
              sigma_star_2_J_D_new[j,d] <- sigma_star_2_J_D_old[j,d]
            }else{
              accept_count <- accept_count+1
            }

          }else{
            # Truncate
            # X = -log(1/sigma_star - 1/upper)
            if(n <= 100){
              X_new[j,d] <- rnorm(n = 1, mean = -ln(1/sigma_star_2_J_D_old[j,d] - 1/upper), sd = 0.1)
            }else{
              X_new[j,d] <- rnorm(n = 1, mean = -ln(1/sigma_star_2_J_D_old[j,d] - 1/upper), sd = sqrt(2.4^2*variance_old[j,d]+2.4^2*MH.variance))
            }

            # Transform the new value of X back to new value of sigma_star_2_j_d
            sigma_star_2_J_D_new[j,d] <- 1/(exp(-X_new[j,d])+1/upper)
            if(sigma_star_2_J_D_new[j,d]==0) {
              print(paste('Truncate:','cluster',j,'dataset',d,": sigma==0"))
              sigma_star_2_J_D_new[j,d] <- .Machine$double.eps
            }

            # Subset of t such that Z[[d]][c]=j
            t_c_d_sub <- t[[d]][Z[[d]]==j]

            # Compute log acceptance probability
            log_acceptance <- sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_new[j,d], m_2 = m_2,
                                                    h_j = h_J[j], t_star_j_d = t_star_J_D[j,d],
                                                    t_c_d_sub = t_c_d_sub) -
              sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_old[j,d], m_2,
                                    h_J[j], t_star_J_D[j,d], t_c_d_sub) -
              ln(sigma_star_2_J_D_old[j,d]) - ln(upper-sigma_star_2_J_D_old[j,d]) +
              ln(sigma_star_2_J_D_new[j,d]) + ln(upper-sigma_star_2_J_D_new[j,d])

            acceptance_sigma <- exp(log_acceptance)
            acceptance_sigma <- min(1,acceptance_sigma)

            # Update sigma_star_2_J_D_new, X_new, X_mean_new, variance_new
            outcome <- rbinom(n = 1, size = 1, prob = acceptance_sigma)
            if(is.na(outcome) == TRUE | outcome == 0){
              X_new[j,d] <- -ln(1/sigma_star_2_J_D_old[j,d] - 1/upper)
              sigma_star_2_J_D_new[j,d] <- sigma_star_2_J_D_old[j,d]
            }else{
              accept_count <- accept_count+1
            }

          } # End sampling in truncation case (non-empty clusters)

        } # End non-empty clusters

        # Update variance in adaptive MH
        X_mean_new[j,d] <- (1-1/n)*X_mean_old[j,d] + 1/n*X_new[j,d]
        M_2_new[j,d] <- M_2_old[j,d] + (X_new[j,d]-X_mean_old[j,d])*(X_new[j,d]-X_mean_new[j,d])
        variance_new[j,d] <- 1/(n-1)*M_2_new[j,d]

      } # End for j in 1:J

    }

    return(list(sigma_star_2_J_D_new=sigma_star_2_J_D_new, X_mean_new=X_mean_new, M_2_new=M_2_new,
                variance_new=variance_new, accept=accept_count))
  }

  ##----------- Simulation of hyper parameters r_J (in prior for t_star) ----------------
  r_J_update <- function(t_star_J_D, mu_r, sigma_r, s_2){
    J <- nrow(t_star_J_D)
    D <- ncol(t_star_J_D)

    # Update
    r_J <- sapply(1:J, function(j) {
      mu_r_hat <- (mu_r*s_2+sigma_r^2*sum(t_star_J_D[j,]))/(s_2+D*sigma_r^2)
      sigma_r_2_hat <- sigma_r^2*s_2/(s_2+D*sigma_r^2)
      return(rnorm(1, mean = mu_r_hat, sd = sqrt(sigma_r_2_hat)))
    })
  }

  ##----------- Simulation of hyper parameters s_2 (in prior for t_star) ----------------
  s_2_update <- function(t_star_J_D, eta_1, eta_2, r_J){
    J <- nrow(t_star_J_D)
    D <- ncol(t_star_J_D)

    # Update
    shape <- J*D/2+eta_1

    temp <- sapply(1:D, function(d) {
      return(t_star_J_D[,d]-r_J)
    })

    rate <- eta_2+sum(temp^2)/2

    s_2 <- extraDistr::rinvgamma(1,alpha = shape,beta = rate)

    return(s_2)
  }

  ##----------- Simulation of hyper parameters h_J (in prior for sigma_star_2) ----------------
  h_J_update <- function(sigma_star_2_J_D, mu_h, sigma_h, m_2){
    J <- nrow(sigma_star_2_J_D)
    D <- ncol(sigma_star_2_J_D)

    # Update
    h_J <- sapply(1:J, function(j) {
      mu_h_hat <- (mu_h*m_2+sigma_h^2*sum(ln(sigma_star_2_J_D[j,])))/(m_2+D*sigma_h^2)
      sigma_h_2_hat <- sigma_h^2*m_2/(m_2+D*sigma_h^2)
      return(rnorm(1, mean = mu_h_hat, sd = sqrt(sigma_h_2_hat)))
    })
  }

  ##----------- Simulation of hyper parameters m_2 (in prior for sigma_star_2) ----------------
  m_2_update <- function(sigma_star_2_J_D, kappa_1, kappa_2, h_J){
    J <- nrow(sigma_star_2_J_D)
    D <- ncol(sigma_star_2_J_D)

    # Update
    shape <- J*D/2+kappa_1

    temp <- sapply(1:D, function(d) {
      return(ln(sigma_star_2_J_D[,d])-h_J)
    })

    rate <- kappa_2+sum(temp^2)/2

    m_2 <- extraDistr::rinvgamma(1,alpha = shape,beta = rate)

    return(m_2)
  }

  ##----------- Simulation of Component probabilities P----------------

  component_log_prob <- function(P, Q_J_D, alpha_0, alpha){

    J <- nrow(Q_J_D)
    D <- ncol(Q_J_D)

    lprod <- sum((alpha_0/J-1)*ln(P))
    for(d in 1:D){
      lprod <- lprod + sum(alpha*P*ln(Q_J_D[,d])-lgamma(alpha*P))
    }

    return(lprod)
  }

  component_probabilities <- function(P, Q_J_D, alpha_0, alpha, covariance,
                                      mean_x, tilde_s, iter_num, s_d_P, MH.variance){

    J <- nrow(Q_J_D)
    D <- ncol(Q_J_D)


    # Define the inputs
    P_old <- P
    covariance_old <- covariance
    mean_x_old <- mean_x
    tilde_s_old <- tilde_s
    s_d_P_old <- s_d_P

    X_old <- ln(P_old[1:(J-1)]/P_old[J]) # Length = J-1

    # Iteration index
    n <- iter_num

    # Adaptive step
    if(n <= 100){
      X_new <- mvtnorm::rmvnorm(n = 1, mean = X_old, sigma = 0.01*diag(x=1,nrow = J-1, ncol = J-1))
    }else{
      X_new <- mvtnorm::rmvnorm(n = 1, mean = X_old, sigma = s_d_P_old*(covariance_old + MH.variance*diag(1, nrow = J-1, ncol = J-1)))
    }

    # Compute P_new (Length = J) from X_new
    P_new <- c(exp(X_new)/(1+sum(exp(X_new))),1/(1+sum(exp(X_new))))

    if(any(P_new==0)){
      print(P_new)
    }

    # Compute acceptance probability
    log_acceptance <- component_log_prob(P_new, Q_J_D, alpha_0, alpha) -
      component_log_prob(P_old, Q_J_D, alpha_0, alpha) +
      sum(ln(P_new)-ln(P_old))
    acceptance_P <- exp(log_acceptance)
    acceptance_P <- min(1,acceptance_P)
    if(is.na(acceptance_P)) {
      print(acceptance_P)
    }

    # Update adaptive scaling parameter
    s_d_P_new <- exp(log(s_d_P_old)+n^(-0.7)*(acceptance_P-0.234))
    if(s_d_P_new>exp(50)) {s_d_P_new <- exp(50)}
    if(s_d_P_new<exp(-50)) {s_d_P_new <- exp(-50)}

    outcome <- rbinom(n = 1, size = 1, prob=acceptance_P)
    if(is.na(outcome) == TRUE | outcome == 0){
      X_new <- X_old
      P_new <- P_old
      accept <- 0
    }else{
      accept <- 1
    }

    # Update covariance, mean_x and tilde_s
    tilde_s_new <- tilde_s_old + matrix(X_new, ncol = 1)%*%matrix(X_new, nrow = 1)
    mean_x_new <- mean_x_old*(1-1/n) + 1/n*matrix(X_new, nrow = 1)
    covariance_new <- 1/(n-1)*tilde_s_new - n/(n-1)*t(mean_x_new)%*%mean_x_new

    return(list(P_new=P_new, tilde_s_new=tilde_s_new, mean_x_new=mean_x_new,
                covariance_new=covariance_new, accept=accept, accept_prob=acceptance_P, s_d_P_new=s_d_P_new))
  }

  ## ----------------------- Simulation of Concentration parameter alpha -----------------------
  alpha_log_prob <- function(Q_J_D, P, alpha){

    D <- ncol(Q_J_D)
    # Construct the log-probability
    lprod <- -alpha + sum(vapply(1:D, function(d) {
      sum(alpha*P*ln(Q_J_D[,d])-lgamma(alpha*P))
    }, FUN.VALUE = numeric(1)))

    return(lprod)
  }

  alpha_sim <- function(Q_J_D, P, alpha, X_mean, M_2, variance, iter_num, MH.variance){

    D <- ncol(Q_J_D)

    # Define the inputs
    alpha_old <- alpha; X_old <- ln(alpha_old)
    X_mean_old <- X_mean
    M_2_old <- M_2
    variance_old <- variance

    # Iteration index
    n <- iter_num

    # Apply AMH based on the iterative number of the current iteration
    # to simulated new value of X
    if(n <= 100){
      X_new <- rnorm(n = 1, mean = X_old, sd = 0.1)
    }else{
      X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(2.4^2*variance_old + 2.4^2*MH.variance))
    }

    # Transform the new value of X back to new value of alpha, namely alpha_new
    alpha_new <- exp(X_new)

    # Compute log acceptance probability
    log_acceptance <- alpha_log_prob(Q_J_D, P, alpha = alpha_new) -
      alpha_log_prob(Q_J_D, P, alpha = alpha_old) +
      ln(alpha_new) - ln(alpha_old)
    acceptance_alpha <- exp(log_acceptance)
    acceptance_alpha <- min(1,acceptance_alpha)

    # Update X_alpha
    outcome <- rbinom(n = 1, size = 1, prob = acceptance_alpha)
    if(is.na(outcome) == TRUE | outcome == 0){
      X_new <- X_old
      alpha_new <- alpha_old
      accept <- 0
    }else{
      accept <- 1
    }

    X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
    M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
    variance_new <- 1/(n-1)*M_2_new

    return(list(alpha_new=alpha_new, X_mean_new=X_mean_new, M_2_new=M_2_new,
                variance_new=variance_new, accept=accept, accept_prob=acceptance_alpha))
  }

  ##---------------------- Simulation of concentration parameter alpha_0 -----------------
  alpha_0_log_prob <- function(P,alpha_0){

    J <- length(P)

    lprob <- -alpha_0 + lgamma(alpha_0) - J*lgamma(alpha_0/J) + sum(alpha_0/J*ln(P))

    # Return the log-probability
    return(lprob)
  }

  alpha_0_sim <- function(P, alpha_0, X_mean, M_2, variance, iter_num, MH.variance){

    J <- length(P)

    # Define the inputs
    alpha_0_old <- alpha_0; X_old <- ln(alpha_0_old)
    variance_old <- variance
    M_2_old <- M_2
    X_mean_old <- X_mean

    n <- iter_num
    # Apply AMH
    if(n <= 100){
      X_new <- rnorm(n = 1, mean = X_old, sd = 0.1)
    }else{
      X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(2.4^2*variance_old + 2.4^2*MH.variance))
    }

    # Obtain the new simulated value for alpha_0
    alpha_0_new <- exp(X_new)

    # Compute acceptance probability
    log_acceptance <- alpha_0_log_prob(P,alpha_0 = alpha_0_new) -
      alpha_0_log_prob(P, alpha_0 = alpha_0_old) +
      ln(alpha_0_new) - ln(alpha_0_old)
    acceptance_alpha0 <- exp(log_acceptance)
    acceptance_alpha0 <- min(1,acceptance_alpha0)

    # Update X_alpha_0 and output alpha_0_new
    outcome <- rbinom(n = 1, size = 1, prob = acceptance_alpha0)
    if(is.na(outcome) == TRUE | outcome == 0){
      X_new <- X_old
      alpha_0_new <- alpha_0_old
      accept <- 0
    }else{
      accept <- 1
    }


    X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
    M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
    variance_new <- 1/(n-1)*M_2_new

    return(list(alpha_0_new=alpha_0_new, X_mean_new=X_mean_new, M_2_new=M_2_new,
                variance_new=variance_new, accept=accept, accept_prob=acceptance_alpha0))
  }

  ##--------------- Simulation of regression parameters -------------
  mean_dispersion <- function(mu_star_1_J, phi_star_1_J,m_b,v_1,v_2,quadratic=FALSE){
    # mu_star_1_J: matrix JxG: current values of mu_star_j_g
    # phi_star_1_J: matrix JxG: current values of phi_star_j_g

    J <- nrow(mu_star_1_J)
    G <- ncol(mu_star_1_J)

    # Let m_b be a vector
    num_unknown <- length(m_b)
    # Set initial values
    parameter0 <- diag(x=1, nrow=num_unknown, ncol=num_unknown)
    parameter1 <- matrix(m_b,ncol=1)
    parameter2 <- 0

    # For each j, update the parameters
    if(quadratic==FALSE){
      for(j in 1:J) {
        mu_tilde_j <- matrix(c(rep(1,G),log(mu_star_1_J[j,])), ncol=2, nrow=G)
        # For posterior variance of b
        parameter0 <- parameter0 + t(mu_tilde_j)%*%(mu_tilde_j)
        # For posterior mean of b
        parameter1 <- parameter1 + t(mu_tilde_j)%*%(matrix(ln(phi_star_1_J[j,]), nrow=G, ncol=1))
        # For scale parameter of IG prior for alpha_phi^2
        parameter2 <- parameter2 + (matrix(ln(phi_star_1_J[j,]), nrow=1, ncol=G))%*%
          (matrix(ln(phi_star_1_J[j,]), nrow=G, ncol=1))
      }
    }else{
      for(j in 1:J){
        mu_tilde_j <- matrix(c(rep(1,G),log(mu_star_1_J[j,]),(log(mu_star_1_J[j,]))^2), ncol=3, nrow=G)
        # For posterior variance of b
        parameter0 <- parameter0 + t(mu_tilde_j)%*%(mu_tilde_j)
        # For posterior mean of b
        parameter1 <- parameter1 + t(mu_tilde_j)%*%(matrix(ln(phi_star_1_J[j,]), nrow=G, ncol=1))
        # For scale parameter of IG prior for alpha_phi^2
        parameter2 <- parameter2 + (matrix(ln(phi_star_1_J[j,]), nrow=1, ncol=G))%*%
          (matrix(ln(phi_star_1_J[j,]), nrow=G, ncol=1))
      }
    }


    # Setting parameters
    V_tilde_b <- Matrix::solve(parameter0)
    m_tilde_b <- V_tilde_b%*%parameter1
    v_tilde_1 <- v_1 + J*G/2
    v_tilde_2 <- v_2 + 1/2*(parameter2-t(m_tilde_b)%*%parameter0%*%(m_tilde_b) + sum(m_b^2))

    # Simulate alpha_phi_2
    alpha_phi_2 <- extraDistr::rinvgamma(n = 1, alpha = v_tilde_1, beta = v_tilde_2)

    # Simulate b
    b <- mvtnorm::rmvnorm(n=1, mean = m_tilde_b, sigma = alpha_phi_2*(V_tilde_b))

    # Return alpha_phi_2 and b
    return(list(alpha_phi_2=alpha_phi_2,b=b))
  }

  ##--------------- Simulation of unique parameters --------------
  unique_parameters_log_prob <- function(mu_star, phi_star, Z, b, alpha_phi_2, Y, Beta, j, g, alpha_mu_2, quadratic=FALSE){

    # Define B_c_d_j to be the vector of Beta_c_d, such that Z_c_d = j
    # The algorithm select the cells separately for the 1st and 2nd dataset
    # Beta: a list 2, each element is a vector of length=cell size
    B_c_d_j <- c(Beta[[1]][which(Z[[1]]==j)],Beta[[2]][which(Z[[2]]==j)])

    # Define Y_c_g_d_j to be the vector of Y_c_g_d, such that Z_c_d = j and g given
    # The algorithm select the cells separately for the 1st and 2nd dataset
    Y_c_g_d_j <- c(Y[[1]][g,which(Z[[1]]==j)], Y[[2]][g,which(Z[[2]]==j)])

    if(quadratic==FALSE){
      lprod1 <- -ln(mu_star*phi_star) - 1/(2*alpha_mu_2)*(ln(mu_star))^2 -
        (ln(phi_star)-(b[1]+b[2]*ln(mu_star)))^2/(2*alpha_phi_2)
    }else{
      lprod1 <- -ln(mu_star*phi_star) - 1/(2*alpha_mu_2)*(ln(mu_star))^2 -
        (ln(phi_star)-(b[1]+b[2]*ln(mu_star)+b[3]*(ln(mu_star))^2))^2/(2*alpha_phi_2)
    }

    if(length(Y_c_g_d_j) != 0){
      # Compute the log-probability

      lprod2 <- sum(phi_star*ln(phi_star/(mu_star*B_c_d_j+phi_star)) +
                      lgamma(Y_c_g_d_j + phi_star) - lgamma(phi_star) - lgamma(Y_c_g_d_j+1) +
                      Y_c_g_d_j*ln(mu_star/(mu_star*B_c_d_j+phi_star)))
      lprod <- lprod1 + lprod2

    }else{
      lprod <- lprod1
    }

    # Return the log-probability
    return(list(lprod1,lprod2))
  }

  unique_parameters_sim <- function(mu_star_1_J, phi_star_1_J, mean_X_mu_phi, tilde_s_mu_phi,
                                    J, G, Z, b, alpha_phi_2, Beta, alpha_mu_2, covariance, iter_num, quadratic=FALSE, Y,
                                    MH.variance){

    # Create a matrix for mu_star_1_J to store simulated values for output
    # Create a matrix for phi_star_1_J to store simulated values for output
    mu_star_1_J_new <- matrix(0, nrow = J, ncol = G); phi_star_1_J_new <- matrix(0, nrow = J, ncol = G)

    mean_X_mu_phi_old <- mean_X_mu_phi; covariance_old <- covariance; tilde_s_mu_phi_old <- tilde_s_mu_phi

    # Iteration index
    n <- iter_num

    # Prepare for outputs
    covariance_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
    tilde_s_mu_phi_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
    mean_X_mu_phi_new <- rep(list(rep(list(matrix(0,nrow=1,ncol=2)),G)),J)

    accept_count <- 0

    # Loop over each j and each g
    for(j in 1:J){
      loop.result <- lapply(1:G, function(g) {

        # Transform into X
        X_mu_phi_star_old <- c(ln(mu_star_1_J[j,g]),ln(phi_star_1_J[j,g]))
        mu_star_old <- mu_star_1_J[j,g]; phi_star_old <- phi_star_1_J[j,g]

        # If cluster j is empty with respect to all datasets, sample mu and phi from their priors
        if(sum(unlist(sapply(Z,function(l) l==j)))==0){
          mu_star_new <- rlnorm(1, meanlog = 0, sdlog = sqrt(alpha_mu_2))
          phi_star_new <- rlnorm(1, meanlog = b[1]+b[2]*ln(mu_star_new), sdlog = sqrt(alpha_phi_2))
          X_mu_phi_star_new <- c(ln(mu_star_new),ln(phi_star_new))
          if(is.na(phi_star_new)){
            print(phi_star_new)
          }

          accept <- 1

        }else{
          # Adaptive MH
          # Simulate new value of X, based on the previous value and covariance structure
          if(n <= 100){
            X_mu_phi_star_new <- mvtnorm::rmvnorm(n = 1, mean = X_mu_phi_star_old, sigma = 0.01*diag(1,nrow = 2, ncol = 2))
          }else{
            X_mu_phi_star_new <- mvtnorm::rmvnorm(n = 1, mean = X_mu_phi_star_old,
                                         sigma = 1*(covariance_old[[j]][[g]] + MH.variance*diag(1,nrow = 2, ncol = 2)))
          }

          # Convert back to mu_new and phi_new and based on these values
          # to compute the acceptance probability
          mu_star_new <- exp(X_mu_phi_star_new[1]); phi_star_new <- exp(X_mu_phi_star_new[2])

          if(is.finite(phi_star_new)){
            # If the new phi is not infinite, then compute for probability
            acceptance_prob_log <- sum(unlist(unique_parameters_log_prob(mu_star = mu_star_new, phi_star = phi_star_new,
                                                                         Z, b, alpha_phi_2, Y, Beta, j, g, alpha_mu_2, quadratic = quadratic))) -
              sum(unlist(unique_parameters_log_prob(mu_star = mu_star_old, phi_star = phi_star_old,
                                                    Z,b,alpha_phi_2,Y,Beta,j,g,alpha_mu_2, quadratic = quadratic))) - ln(mu_star_old) -ln(phi_star_old) +
              ln(mu_star_new) + ln(phi_star_new)

            acceptance_unique <- min(1, exp(acceptance_prob_log))

          }else{
            acceptance_unique <- 0
          }

          outcome <- rbinom(n = 1, size = 1, prob = acceptance_unique)

          if(is.na(outcome) == TRUE | outcome == 0){
            X_mu_phi_star_new <- X_mu_phi_star_old
            mu_star_new <- mu_star_old
            phi_star_new <- phi_star_old
            accept <- 0
          }else{
            accept <- 1
          }

        } # End doing adaptive MH


        # Return outputs
        return(list(mu_star_new, phi_star_new, X_mu_phi_star_new, accept))
      }) # End lapply

      ## Update AMH information

      for(g in 1:G){

        mu_star_1_J_new[j,g] <- loop.result[[g]][[1]]
        phi_star_1_J_new[j,g] <- loop.result[[g]][[2]]

        X_mu_phi_star_new <- loop.result[[g]][[3]]

        accept_count <- accept_count + loop.result[[g]][[4]]

        tilde_s_mu_phi_new_j_g <- tilde_s_mu_phi_old[[j]][[g]] + matrix(X_mu_phi_star_new,ncol=1) %*%
          matrix(X_mu_phi_star_new,nrow=1)
        mean_X_mu_phi_new_j_g <- mean_X_mu_phi_old[[j]][[g]]*(1-1/n) + 1/n*matrix(X_mu_phi_star_new,nrow=1)
        covariance_new_j_g <- 1/(n-1)*tilde_s_mu_phi_new_j_g - n/(n-1)*t(mean_X_mu_phi_new_j_g)%*%mean_X_mu_phi_new_j_g

        covariance_new[[j]][[g]] <- covariance_new_j_g
        tilde_s_mu_phi_new[[j]][[g]] <- tilde_s_mu_phi_new_j_g
        mean_X_mu_phi_new[[j]][[g]] <- mean_X_mu_phi_new_j_g
      }


    } # End for j in 1:J


    # Returning outputs
    return(list(mu_star_1_J_new=mu_star_1_J_new, phi_star_1_J_new=phi_star_1_J_new,
                accept_count=accept_count, tilde_s_mu_phi_new=tilde_s_mu_phi_new,
                mean_X_mu_phi_new=mean_X_mu_phi_new,covariance_new=covariance_new))
  }

  ##------------- Simulation of capture efficiencies --------------------

  capture_efficiencies_log_prob <- function(Beta_single, Y, Z, mu_star_1_J, phi_star_1_J, c, d, a_d_beta, b_d_beta){

    j <- Z[[d]][c]
    # Construct the log-probability
    lprod1 <- (a_d_beta[d]-1)*ln(Beta_single) + (b_d_beta[d]-1)*ln(1-Beta_single)
    lprod2 <- sum((phi_star_1_J[j,]+Y[[d]][,c])*ln(phi_star_1_J[j,]+mu_star_1_J[j,]*Beta_single) - Y[[d]][,c]*ln(Beta_single))
    lprod <- lprod1 - lprod2

    # Return the log-probability
    return(lprod)
  }


  capture_efficiencies_sim <- function(Beta, Y, Z, mu_star_1_J, phi_star_1_J, a_d_beta, b_d_beta,
                                       iter_num, M_2, mean_X, variance,
                                       MH.variance){

    D <- length(Z)
    C_d <- unlist(lapply(Z, length))

    M_2_old <- M_2; mean_X_old <- mean_X; variance_old <- variance

    # To store outputs
    Beta_new <- NULL; Beta_new[[1]] <- rep(0,C_d[1]); Beta_new[[2]] <- rep(0,C_d[2])

    variance_new <- NULL; variance_new[[1]] <- rep(0,C_d[1]); variance_new[[2]] <- rep(0,C_d[2])
    mean_X_new <- NULL; mean_X_new[[1]] <- rep(0,C_d[1]); mean_X_new[[2]] <- rep(0,C_d[1])
    M_2_new <- NULL; M_2_new[[1]] <- rep(0,C_d[1]); M_2_new[[2]] <- rep(0,C_d[2])

    n <- iter_num

    accept_count <- 0

    for(d in 1:D){
      loop.result <- lapply(1:C_d[d], function(c) {

        Beta_single_old <- Beta[[d]][c]
        X_old <- ln(Beta_single_old/(1-Beta_single_old))

        if(n <= 100){
          X_new <- rnorm(n = 1, mean = X_old, sd = 0.1)
        }else{
          X_new <- rnorm(n = 1, mean = X_old, sd = sqrt((2.4^2)*variance_old[[d]][c] + 2.4^2*MH.variance))
        }

        # Transform back to beta
        Beta_single_new <- 1/(1+exp(-X_new))
        if(Beta_single_new==1) {
          Beta_single_new <- 1-.Machine$double.eps
        }

        # Compute the acceptance probability
        acceptance_prob_log <- capture_efficiencies_log_prob(Beta_single = Beta_single_new, Y, Z, mu_star_1_J,
                                                             phi_star_1_J,c,d,a_d_beta,b_d_beta) -
          capture_efficiencies_log_prob(Beta_single = Beta_single_old, Y, Z, mu_star_1_J, phi_star_1_J,
                                        c,d,a_d_beta,b_d_beta) + ln(Beta_single_new) + ln(1-Beta_single_new) -
          ln(Beta_single_old) - ln(1-Beta_single_old)

        acceptance_beta <- min(1,exp(acceptance_prob_log))

        outcome <- rbinom(n = 1, size = 1, prob = acceptance_beta)
        if(is.na(outcome) == TRUE | outcome == 0){
          X_new <- X_old
          Beta_single_new <- Beta_single_old
          accept <- 0
        }else{
          accept <- 1
        }

        ## Output result
        list(Beta_single_new, X_new, accept)
      })

      for(c in 1:C_d[d]){
        Beta_single_new <- loop.result[[c]][[1]]
        X_new <- loop.result[[c]][[2]]
        accept_count <- accept_count+loop.result[[c]][[3]]

        Beta_new[[d]][c] <- Beta_single_new
        mean_X_new[[d]][c] <- (1-1/n)*mean_X_old[[d]][c]+(1/n)*X_new
        M_2_new[[d]][c] <- M_2_old[[d]][c] + (X_new-mean_X_old[[d]][c])*(X_new-mean_X_new[[d]][c])
        variance_new[[d]][c] <- 1/(n-1)*M_2_new[[d]][c]
      }
    }

    return(list(Beta_new=Beta_new, accept_count=accept_count, mean_X_new=mean_X_new,
                M_2_new=M_2_new, variance_new=variance_new))
  }

  #------------------------ Step 1: Prepare for outputs -----------------
  b_output <- NULL
  alpha_phi_2_output <- c()
  Z_output <- NULL

  P_C_J_D_output <- list(list(NULL),list(NULL))

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

  U_C_J_D_output <- list(list(NULL),list(NULL))

  r_J_output <- NULL
  s_2_output <- c()
  h_J_output <- NULL
  m_2_output <- c()


  #----------------------- Step 2: Initial values in MCMC-----------

  # Use K-means for initialization of Z, t_star_J_D, sigma_star_2_J_D, r_J, s^2, h_j, m^2
  # after dimension reduction (t-sne)
  if(is.null(Z_fix)){
    if(empirical_z==TRUE){
      Y_all <- t(cbind(Y[[1]],Y[[2]]))
      tsne_results <- Rtsne::Rtsne(Y_all, perplexity=30, check_duplicates = FALSE, partial_pca = partial_pca)
      km_cluster <- kmeans(tsne_results$Y,centers=J)$cluster

      # Allocation variables Z
      Z_initial <- NULL
      Z_initial[[1]] <- km_cluster[1:C_d[1]]
      Z_initial[[2]] <- km_cluster[(C_d[1]+1):(C_d[2]+C_d[1])]

    }else{
      Z_initial <- NULL
      Z_initial[[1]] <- sample(1:J, size = C_d[[1]],replace = TRUE,prob = rep(1/J,J))
      Z_initial[[2]] <- sample(1:J, size = C_d[[2]],replace = TRUE,prob = rep(1/J,J))
    }

  }else{
    Z_initial <- Z_fix
  }


  # Compute the mean of t within each cluster in each dataset to initialize t_star_J_D
  # matrix of [J, D]
  t_star_J_D_initial <- matrix(NA, nrow = J, ncol = D)

  # Mean of t within each cluster of each dataset
  t_star_J_1 <- tapply(t[[1]], Z_initial[[1]], mean)
  t_star_J_2 <- tapply(t[[2]], Z_initial[[2]], mean)

  # There may be empty clusters
  # Fill the non-empty values into the intialization matrix, impute NA values with dataset-specific mean
  t_star_J_D_initial[as.numeric(names(t_star_J_1)),1] <- t_star_J_1
  t_star_J_D_initial[is.na(t_star_J_D_initial[,1]),1] <- mean(t[[1]])

  t_star_J_D_initial[as.numeric(names(t_star_J_2)),2] <- t_star_J_2
  t_star_J_D_initial[is.na(t_star_J_D_initial[,2]),2] <- mean(t[[2]])

  # Initials for r_J as an empirical mean from t_star_J_D_initial
  r_J_initial <- apply(t_star_J_D_initial, 1, mean)

  # Initials for s^2
  s_2_initial <- mean(apply(t_star_J_D_initial, 1, var))

  # Initials for sigma_star_2_J_D
  sigma_star_2_J_D_initial <- matrix(NA, nrow = J, ncol = D)

  # sd of t within each cluster of each dataset
  sigma_star_2_J_1 <- tapply(t[[1]], Z_initial[[1]], var)
  sigma_star_2_J_2 <- tapply(t[[2]], Z_initial[[2]], var)

  # There may be singleton clusters (NA sd) and empty clusters
  # Fill the non-empty values into the intialization matrix, impute NA values with dataset-specific mean
  sigma_star_2_J_D_initial[as.numeric(names(sigma_star_2_J_1)),1] <- sigma_star_2_J_1
  sigma_star_2_J_D_initial[is.na(sigma_star_2_J_D_initial[,1]),1] <- var(t[[1]])

  sigma_star_2_J_D_initial[as.numeric(names(sigma_star_2_J_2)),2] <- sigma_star_2_J_2
  sigma_star_2_J_D_initial[is.na(sigma_star_2_J_D_initial[,2]),2] <- var(t[[2]])

  # Initials for h_J from log(sigma_star_2_J_D)
  h_J_initial <- apply(ln(sigma_star_2_J_D_initial), 1, mean)

  # Initials for m^2
  m_2_initial <- mean(apply(ln(sigma_star_2_J_D_initial), 1, var))


  # Component probabilities p_j, add 1 to avoid zero p if cluster is empty
  P_initial <- sapply(1:J, function(j) mean(unlist(Z_initial)==j))
  if(any(P_initial==0)){
    # Add 1 to avoid zero p if cluster is empty
    P_initial <- sapply(1:J, function(j) sum(unlist(Z_initial)==j)+1)/(sum(C_d)+J)
  }

  # Dataset-specific vector q_j_d, and compute covariate-dependent probabilities p_j_d
  Q_J_D_initial <- matrix(NA,nrow = J, ncol = D)
  Q_J_D_initial[,1] <- sapply(1:J,function(j) sum(Z_initial[[1]]==j))+1
  Q_J_D_initial[,2] <- sapply(1:J,function(j) sum(Z_initial[[2]]==j))+1

  P_C_J_D_initial <- NULL
  P_C_J_D_initial[[1]] <- t(sapply(t[[1]],function(time) {
    Q_J_D_initial[,1]*rbf(time,t_star_J_D_initial[,1],sigma_star_2_J_D_initial[,1])/(
      sum(Q_J_D_initial[,1]*rbf(time,t_star_J_D_initial[,1],sigma_star_2_J_D_initial[,1])))
  }))
  P_C_J_D_initial[[2]] <- t(sapply(t[[2]],function(time) {
    Q_J_D_initial[,2]*rbf(time,t_star_J_D_initial[,2],sigma_star_2_J_D_initial[,2])/(
      sum(Q_J_D_initial[,2]*rbf(time,t_star_J_D_initial[,2],sigma_star_2_J_D_initial[,2])))
  }))


  # Initialize mu, phi, alpha_phi_2 and beta based on bayNorm supplied with allocations
  # then we random choose k clusters with initial values
  # given by bayNorm without supplying the allocations (global initial).
  # This is to ensure the initial mu, phi have some differences and also similarities across clusters

  # bayNorm without Z
  baynorm_HET <- bayNorm::bayNorm(Data = Y[[1]], BETA_vec = NULL, mode_version = TRUE, mean_version = FALSE,
                         BB_SIZE = FALSE, verbose = FALSE)
  baynorm_HOM <- bayNorm::bayNorm(Data = Y[[2]], BETA_vec = NULL, mode_version = TRUE, mean_version = FALSE,
                         BB_SIZE = FALSE, verbose = FALSE)
  baynorm_tot <- bayNorm::bayNorm(Data = cbind(Y[[1]], Y[[2]]), BETA_vec = NULL, mode_version = TRUE,
                         mean_version = FALSE, BB_SIZE = FALSE, verbose = FALSE)

  mu.estimate <- baynorm_tot$PRIORS$MME_prior$MME_MU
  # For mu=0, it is replaced by the smallest non-zero value
  mu.estimate <- ifelse(mu.estimate == 0, min(baynorm_tot$PRIORS$MME_prior$MME_MU[baynorm_tot$PRIORS$MME_prior$MME_MU!=0]),
                        mu.estimate)
  # Adjust mu according to the given beta.mean because beta estimate from baynorm has a mean of 0.06
  mu.estimate <- mu.estimate*mean(baynorm_tot$PRIORS$BETA_vec)/beta.mean


  phi.estimate <- baynorm_tot$PRIORS$MME_prior$MME_SIZE
  # For phi=Inf, it is replaced by the maximum finite value
  phi.estimate <- ifelse(is.infinite(phi.estimate), max(baynorm_tot$PRIORS$MME_prior$MME_SIZE[is.finite(baynorm_tot$PRIORS$MME_prior$MME_SIZE)]),
                         phi.estimate)

  # bayNorm with Z
  baynorm_HET_by_z <- bayNorm::bayNorm(Data = Y[[1]], BETA_vec = NULL, mode_version = TRUE, mean_version = FALSE,
                              BB_SIZE = FALSE, verbose = FALSE, Conditions = Z_initial[[1]], Prior_type = 'LL')
  baynorm_HOM_by_z <- bayNorm::bayNorm(Data = Y[[2]], BETA_vec = NULL, mode_version = TRUE, mean_version = FALSE,
                              BB_SIZE = FALSE, verbose = FALSE, Conditions = Z_initial[[2]], Prior_type = 'LL')

  baynorm_tot_by_z <- bayNorm::bayNorm(Data = cbind(Y[[1]], Y[[2]]), BETA_vec = NULL, mode_version = TRUE,
                              mean_version = FALSE, BB_SIZE = FALSE, verbose = FALSE, Conditions = c(Z_initial[[1]],
                                                                                                     Z_initial[[2]]),
                              Prior_type = 'LL')

  # Cluster-specific parameters in NB likelihood
  mu_star_1_J_initial <- matrix(0, nrow = J, ncol = G)
  mu_all <- unlist(lapply(baynorm_tot_by_z$PRIORS_LIST, function(l) l$MME_prior$MME_MU))
  for (j in 1:J){
    mu_star_1_J_initial[j,] <- ifelse(baynorm_tot_by_z$PRIORS_LIST[[paste('Group',j)]]$MME_prior$MME_MU == 0,
                                      min(mu_all[mu_all!=0]),
                                      baynorm_tot_by_z$PRIORS_LIST[[paste('Group',j)]]$MME_prior$MME_MU)
  }

  # Adjust mu according to the given beta.mean because beta estimate from baynorm has a mean of 0.06
  mu_star_1_J_initial <- mu_star_1_J_initial*mean(unlist(lapply(baynorm_tot_by_z$PRIORS_LIST, function(l) l$BETA_vec)))/beta.mean

  phi_star_1_J_initial <- matrix(0, nrow = J, ncol = G)
  phi_all <- unlist(lapply(baynorm_tot_by_z$PRIORS_LIST, function(l) l$MME_prior$MME_SIZE))
  for(j in 1:J){
    phi_star_1_J_initial[j,] <- ifelse(is.infinite(baynorm_tot_by_z$PRIORS_LIST[[paste('Group',j)]]$MME_prior$MME_SIZE),
                                       max(phi_all[is.finite(phi_all)]),
                                       baynorm_tot_by_z$PRIORS_LIST[[paste('Group',j)]]$MME_prior$MME_SIZE)
  }

  # Save the adjusted cluster-specific initial mu and phi for later use in initializing alpha_phi_2
  mu.estimate.z <- c(mu_star_1_J_initial)
  phi.estimate.z <- c(phi_star_1_J_initial)

  # Select the number of clusters with global initial mu and phi
  k <- sample(0:J,size=1)
  # If decide to use global initials, random select some clusters labels jj to assign global initials
  if(k!=0){
    jj <- sample(1:J, size = k)
    for (j in jj) {
      mu_star_1_J_initial[j,] <- mu.estimate
      phi_star_1_J_initial[j,] <- phi.estimate
    }
  }
  print(paste(k,'clusters have global mu/phi initialization'))

  # Cell-specific capture efficiency beta_c_d from bayNorm estimates
  Beta_initial <- NULL
  Beta_initial[[1]] <- unlist(lapply(baynorm_HET_by_z$PRIORS_LIST, function(l) l$BETA_vec))
  # Order Beta_initial according to cell orders rather than by group
  # Currently the beta's from bayNorm are named as 'Group 1.cellname'
  # extract everything after the first dot
  group1_name <- sub("^[^\\.]*\\.","",names(Beta_initial[[1]]))
  Beta_initial[[1]] <- unname(Beta_initial[[1]][match(colnames(Y[[1]]),group1_name)])
  Beta_initial[[1]] <- ifelse(Beta_initial[[1]]/mean(Beta_initial[[1]])*beta.mean >= 1,
                              0.99, Beta_initial[[1]]/mean(Beta_initial[[1]])*beta.mean)

  Beta_initial[[2]] <- unlist(lapply(baynorm_HOM_by_z$PRIORS_LIST, function(l) l$BETA_vec))
  # Order Beta_initial according to cell orders rather than by group
  group2_name <- sub("^[^\\.]*\\.","",names(Beta_initial[[2]]))
  Beta_initial[[2]] <- unname(Beta_initial[[2]][match(colnames(Y[[2]]),group2_name)])
  Beta_initial[[2]] <- ifelse(Beta_initial[[2]]/mean(Beta_initial[[2]])*beta.mean >= 1,
                              0.99, Beta_initial[[2]]/mean(Beta_initial[[2]])*beta.mean)

  # Initialize alpha_phi_2 from mu and phi given by bayNorm with Z
  x.1 <- log(mu.estimate.z)
  y.1 <- log(phi.estimate.z)

  if(quadratic==FALSE){
    lm.1 <- lm(y.1 ~ x.1)
  }else{
    lm.1 <- lm(y.1 ~ x.1+I(x.1^2))
  }

  # Estimated sigma^2 (of the linear model): variance of log-phi, used as the initial value of alpha_phi^2
  alpha_phi_2_initial <- deviance(lm.1)/df.residual(lm.1)
  print(paste("alpha_phi_2_initial:",alpha_phi_2_initial))

  # Initialize b
  if(is.null(b_initial)) {
    b_initial <- as.numeric(coef(lm.1))
  }
  print(paste("b_initial:",b_initial))

  # ---------Prior settings-------------
  # Only rely on mu, phi and beta from bayNorm without Z, so the priors for different runs are the same as long as the dataset is fixed
  # even for different J

  # alpha_mu_2
  if(is.null(alpha_mu_2)){
    # prior log(mu) ~ N(0, alpha_mu^2), so alpha_mu^2 = mean(log(mu)^2)
    alpha_mu_2 <- mean(c(log(mu.estimate)^2))
  }
  print(paste("alpha_mu_2:",alpha_mu_2))

  # m_b, v1 and v2
  x.2 <- log(mu.estimate)
  y.2 <- log(phi.estimate)

  if(quadratic==FALSE){
    lm.2 <- lm(y.2 ~ x.2)
  }else{
    lm.2 <- lm(y.2 ~ x.2+I(x.2^2))
  }

  rse.lm.2.squared <- deviance(lm.2)/df.residual(lm.2)
  # Variance of alpha_phi^2 is 1, mean is rse.lm.2.squared in the inverse gamma prior for alpha_phi^2
  variance <- 1
  v_1_empirical <- rse.lm.2.squared^2/variance + 2
  v_2_empirical <- (v_1_empirical - 1)*rse.lm.2.squared

  if(empirical == TRUE){
    v_1 <- v_1_empirical
    v_2 <- v_2_empirical
  }else{
    v_1 <- 2; v_2 <- 1
  }
  print(paste("v1,v2:",v_1,v_2))

  if(empirical == TRUE){
    m_b <- as.numeric(coef(lm.2))
  }else{
    if(quadratic==TRUE){
      m_b <- c(-1,2,0)
    }else{
      m_b <- c(-1,2)
    }
  }
  print(paste("m_b:",m_b))

  # Function to set empirical values for hyper parameters in prior of beta ~ Beta(a_d,b_d)
  baynorm_estimate_beta_param <- function(baynorm_output, additional_variance, beta.mean){

    ## Normlized beta
    beta_vec <- baynorm_output$PRIORS$BETA_vec
    beta.normalized <- beta_vec/mean(beta_vec)*beta.mean
    beta.normalized <- ifelse(beta.normalized >= 1, 0.99, beta.normalized)

    baynorm_mean_capeff <- mean(beta.normalized) # Mean of the beta prior for beta
    baynorm_var_capeff <- var(beta.normalized) + additional_variance # Variance of the beta prior
    a_beta <- ((1-baynorm_mean_capeff)/baynorm_var_capeff - 1/baynorm_mean_capeff)*baynorm_mean_capeff^2
    b_beta <- a_beta*(1/baynorm_mean_capeff - 1)

    # to ensure unimodal prior
    while(a_beta<1 | b_beta<1) {
      baynorm_var_capeff <- baynorm_var_capeff/2
      a_beta <- ((1-baynorm_mean_capeff)/baynorm_var_capeff - 1/baynorm_mean_capeff)*baynorm_mean_capeff^2
      b_beta <- a_beta*(1/baynorm_mean_capeff - 1)

    }

    return(c(a_beta,b_beta))

  }

  # Emprical values for for a_d_beta and b_d_beta (parameters in prior for beta ~ Beta(a_d, b_d))
  bay_summary_HET <- baynorm_estimate_beta_param(baynorm_output = baynorm_HET, additional_variance = 0.01, beta.mean = beta.mean)
  bay_summary_HOM <- baynorm_estimate_beta_param(baynorm_output = baynorm_HOM, additional_variance = 0.01, beta.mean = beta.mean)

  a_d_beta <- rep(0,2); b_d_beta <- rep(0,2)
  a_d_beta[1] <- bay_summary_HET[1]; b_d_beta[1] <- bay_summary_HET[2]
  a_d_beta[2] <- bay_summary_HOM[1]; b_d_beta[2] <- bay_summary_HOM[2]

  if(any(c(a_d_beta,b_d_beta)<0)) {
    stop('Prior setting for beta inappropriate!!')
  }
  print(paste('betas:',a_d_beta,b_d_beta))


  #--------------------------------------------------------------------------------

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
  P_C_J_D_new <- P_C_J_D_initial

  Xi_C_D_new <- Xi_C_D_update(Q_J_D = Q_J_D_new, C_d = C_d, t = t, t_star_J_D = t_star_J_D_new,
                              sigma_star_2_J_D = sigma_star_2_J_D_new)

  U_C_J_D_new <- U_C_J_D_update(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C_d = C_d, t = t,
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
  mean_X_sigma_star_2_new <- ln(sigma_star_2_J_D_new)
  M_2_sigma_star_2_new <- matrix(0, nrow = J, ncol = D)
  variance_sigma_star_2_new <- matrix(0, nrow = J, ncol = D)

  for (d in 1:D) {
    for (j in 1:J) {
      ind <- c(1:C_d[d])[-ln(U_C_J_D_new[[d]][,j])<Xi_C_D_new[[d]]*Q_J_D_new[j,d]]
      if(length(ind)!=0) {
        uppers <- -(t[[d]]-t_star_J_D_new[j,d])^2/2/(ln(-ln(U_C_J_D_new[[d]][,j]))-ln(Xi_C_D_new[[d]])-ln(Q_J_D_new[j,d]))
        upper <- min(uppers[ind])
        mean_X_sigma_star_2_new[j,d] <- -ln(1/sigma_star_2_J_D_new[j,d] - 1/upper)
      }
    }
  }

  # 1) For Component probabilities
  s_d_P_new <- 0.001
  #X_n (n: indicate iterations)
  mean_X_component_new <- ln(matrix(P_initial[1:(J-1)]/P_initial[J], nrow = 1)) # 1x(J-1)
  tilde_s_component_new <- t(mean_X_component_new)%*%mean_X_component_new
  #At 1st iteration, the covariance based on the intial values are 0
  covariance_component_new <- matrix(0, nrow = J-1, ncol = J-1)

  # 2) For alpha
  mean_X_alpha_new <- ln(alpha_new)
  M_2_alpha_new <- 0
  variance_alpha_new <- 0

  # 3) For alpha_0
  mean_X_alpha_0_new <- ln(alpha_0_new)
  M_2_alpha_0_new <- 0
  variance_alpha_0_new <- 0

  # 4) Unique parameters
  # JxG ge 2x2 matrices
  covariance_unique_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
  tilde_s_unique_new <- rep(list(rep(list(matrix(0,nrow=2,ncol=2)),G)),J)
  mean_X_unique_new <- rep(list(rep(list(matrix(0,nrow=1,ncol=2)),G)),J)
  for(j in 1:J){
    for(g in 1:G){
      mean_X_unique_new[[j]][[g]] <- matrix(c(ln(mu_star_1_J_new[j,g]),ln(phi_star_1_J_new[j,g])),nrow=1)
      tilde_s_unique_new[[j]][[g]] <- t(mean_X_unique_new[[j]][[g]])%*%mean_X_unique_new[[j]][[g]]
    }
  }

  # 5) Capture efficiency
  mean_X_capture_new <- NULL
  mean_X_capture_new[[1]] <- ln(Beta_new[[1]]/(1-Beta_new[[1]]))
  mean_X_capture_new[[2]] <- ln(Beta_new[[2]]/(1-Beta_new[[2]]))

  M_2_capture_new <- NULL
  M_2_capture_new[[1]] <- rep(0,C_d[1])
  M_2_capture_new[[2]] <- rep(0,C_d[2])

  variance_capture_new <- NULL
  variance_capture_new[[1]] <- rep(0,C_d[1])
  variance_capture_new[[2]] <- rep(0,C_d[2])

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
    # In the case of thinning, if we require to save some samples that are not divisible by thinning (for consensus clustering) or before burn-in
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
    Xi_C_D_new <- Xi_C_D_update(Q_J_D = Q_J_D_new, C_d = C_d, t = t, t_star_J_D = t_star_J_D_new,
                                sigma_star_2_J_D = sigma_star_2_J_D_new)

    # 1) Update dataset-specific vector q_j_d
    Q_J_D_new <- Q_J_D_update(Z = Z_new, alpha = alpha_new, P = P_new, Xi_C_D = Xi_C_D_new,
                              t = t, t_star_J_D = t_star_J_D_new,
                              sigma_star_2_J_D = sigma_star_2_J_D_new)

    # 2) Update the allocation variable, if Z_fix is not provided
    if(is.null(Z_fix)){
      Z_new <- allocation_variables_update(Y = Y, t = t, mu_star_1_J = mu_star_1_J_new, phi_star_1_J = phi_star_1_J_new,
                                           Beta = Beta_new, Q_J_D = Q_J_D_new, t_star_J_D = t_star_J_D_new,
                                           sigma_star_2_J_D = sigma_star_2_J_D_new)
    }


    # 3) Update kernel parameters
    # 3-1) latent variable U_C_J_D
    U_C_J_D_new <- U_C_J_D_update(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C_d = C_d, t = t,
                                  t_star_J_D = t_star_J_D_new, sigma_star_2_J_D = sigma_star_2_J_D_new,
                                  rbf = rbf)

    # 3_2) t_star_J_D
    t_star_J_D_new <- t_star_J_D_update(r_J = r_J_new, s_2 = s_2_new, Z = Z_new, t = t, C_d = C_d,
                                        sigma_star_2_J_D = sigma_star_2_J_D_new, U_C_J_D = U_C_J_D_new,
                                        Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new)

    # Update U again
    U_C_J_D_new <- U_C_J_D_update(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C_d = C_d, t = t,
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


    # Update P_C_J_D
    P_C_J_D_new <- NULL
    P_C_J_D_new[[1]] <- t(sapply(t[[1]],function(time) {
      Q_J_D_new[,1]*rbf(time,t_star_J_D_new[,1],sigma_star_2_J_D_new[,1])/(
        sum(Q_J_D_new[,1]*rbf(time,t_star_J_D_new[,1],sigma_star_2_J_D_new[,1])))
    }))
    P_C_J_D_new[[2]] <- t(sapply(t[[2]],function(time) {
      Q_J_D_new[,2]*rbf(time,t_star_J_D_new[,2],sigma_star_2_J_D_new[,2])/(
        sum(Q_J_D_new[,2]*rbf(time,t_star_J_D_new[,2],sigma_star_2_J_D_new[,2])))
    }))


    # 4) Update hyper parameters in priors for t_star, sigma_star_2
    # 4-1) r_J
    r_J_new <- r_J_update(t_star_J_D = t_star_J_D_new, mu_r = mu_r, sigma_r = sigma_r, s_2 = s_2_new)

    # 4-2) s_2
    s_2_new <- s_2_update(t_star_J_D = t_star_J_D_new, eta_1 = eta_1, eta_2 = eta_2, r_J = r_J_new)

    # 4-3) h_J
    h_J_new <- h_J_update(sigma_star_2_J_D = sigma_star_2_J_D_new, mu_h = mu_h, sigma_h = sigma_h, m_2 = m_2_new)

    # 4_4) m_2
    m_2_new <- m_2_update(sigma_star_2_J_D = sigma_star_2_J_D_new, kappa_1 = kappa_1, kappa_2 = kappa_2, h_J = h_J_new)


    # 4) Update the component probabilities P
    component_output <- component_probabilities(P = P_new, Q_J_D = Q_J_D_new, alpha_0 = alpha_0_new,
                                                alpha = alpha_new, covariance = covariance_component_new,
                                                mean_x = mean_X_component_new,
                                                tilde_s = tilde_s_component_new,
                                                iter_num = iter, s_d_P = s_d_P_new[iter-1],
                                                MH.variance = MH.variance)
    P_new <- component_output$P_new
    tilde_s_component_new <- component_output$tilde_s_new
    mean_X_component_new <- component_output$mean_x_new
    covariance_component_new <- component_output$covariance_new
    s_d_P_new[iter] <- component_output$s_d_P_new
    P_count <- P_count + component_output$accept
    acceptance_count_avg$P_accept[iter-1] <- P_count/(iter-1)

    # 5) Update alpha
    alpha_output_sim <- alpha_sim(Q_J_D = Q_J_D_new, P = P_new, alpha = alpha_new,
                                  X_mean = mean_X_alpha_new, M_2 = M_2_alpha_new,
                                  variance = variance_alpha_new, iter_num = iter,
                                  MH.variance = MH.variance)

    alpha_new <- alpha_output_sim$alpha_new
    mean_X_alpha_new <- alpha_output_sim$X_mean_new
    M_2_alpha_new <- alpha_output_sim$M_2_new
    variance_alpha_new <- alpha_output_sim$variance_new
    alpha_count <- alpha_count + alpha_output_sim$accept
    acceptance_count_avg$alpha_accept[iter-1] <- alpha_count/(iter-1)

    # 6) Update alpha_0
    alpha_0_output_sim <- alpha_0_sim(P = P_new, alpha_0 = alpha_0_new,
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

    # 7) Update mean_dispersion
    mean_dispersion_output <- mean_dispersion(mu_star_1_J = mu_star_1_J_new,
                                              phi_star_1_J = phi_star_1_J_new,
                                              v_1 = v_1, v_2 = v_2, m_b = m_b,
                                              quadratic = quadratic)

    alpha_phi_2_new <- mean_dispersion_output$alpha_phi_2
    b_new <- mean_dispersion_output$b

    # 8) Update unique parameters
    unique_output_sim <- unique_parameters_sim(mu_star_1_J = mu_star_1_J_new,
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

    # 9) Update capture efficiency
    capture_output_sim <- capture_efficiencies_sim(Beta = Beta_new, Y = Y, Z = Z_new,
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

    if((iter-1) %% save_frequency == 0 && auto.save == TRUE){
      my_list <- list('b_output' = b_output, 'alpha_phi2_output' = alpha_phi_2_output, 'Z_output' = Z_output,
                      'P_C_J_D_output' = P_C_J_D_output,  'P_output' = P_output, 'alpha_output' = alpha_output,
                      'alpha_0_output' = alpha_0_output, 'mu_star_1_J_output' = mu_star_1_J_output,
                      'phi_star_1_J_output' = phi_star_1_J_output, 'Beta_output' = Beta_output,
                      'Q_J_D_output' = Q_J_D_output, 'Xi_C_D_output' = Xi_C_D_output, 't_star_J_D_output' = t_star_J_D_output,
                      'sigma_star_2_J_D_output' = sigma_star_2_J_D_output, 'U_C_J_D_output' = U_C_J_D_output,
                      'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output,
                      'acceptance_count_avg' = acceptance_count_avg,
                      'output_index' = output_index,
                      's_d_P_new' = s_d_P_new,
                      'alpha_mu_2' = alpha_mu_2, 'v_1' = v_1, 'v_2' = v_2, 'm_b' = m_b,
                      'a_d_beta' = a_d_beta, 'b_d_beta' = b_d_beta)
      save(my_list, file=partial.save.name)
    }

  }

  ## Return the list
  my_list <- list('b_output' = b_output, 'alpha_phi2_output' = alpha_phi_2_output, 'Z_output' = Z_output,
                  'P_C_J_D_output' = P_C_J_D_output,  'P_output' = P_output, 'alpha_output' = alpha_output,
                  'alpha_0_output' = alpha_0_output, 'mu_star_1_J_output' = mu_star_1_J_output,
                  'phi_star_1_J_output' = phi_star_1_J_output, 'Beta_output' = Beta_output,
                  'Q_J_D_output' = Q_J_D_output, 'Xi_C_D_output' = Xi_C_D_output, 't_star_J_D_output' = t_star_J_D_output,
                  'sigma_star_2_J_D_output' = sigma_star_2_J_D_output, 'U_C_J_D_output' = U_C_J_D_output,
                  'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output,
                  'acceptance_count_avg' = acceptance_count_avg,
                  'output_index' = output_index,
                  's_d_P_new' = s_d_P_new,
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

