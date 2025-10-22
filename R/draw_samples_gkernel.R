# --------------- Step 0: Functions for each step of the Gibb sampling ----------
# ---------------- For single-cell clustering ---------------------------------

##--------------- Simulation of latent Xi ----------------
Xi_C_D_update_gkernel <- function(Q_J_D, C_d, t, t_star_J_D, sigma_star_2_J_D){

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
      log_temp <- log(Q_J_D[,d])+log_rbf_value
      log_K <- max(log_temp)

      return(exp(log_K)*sum(exp(log_temp-log_K)))
    }, FUN.VALUE = numeric(1))

    Xi[[d]] <- rgamma(C_d[d], shape = 1, rate = rates)
  }

  return(Xi)
}

##--------------- Simulation of dataset-specific vector q_j_d ----------------
Q_J_D_update_gkernel <- function(Z, alpha, P, Xi_C_D, t, t_star_J_D, sigma_star_2_J_D){

  J <- nrow(t_star_J_D)
  D <- ncol(t_star_J_D)

  # Set up the matrix to save updated values
  Q <- matrix(NA, nrow = J, ncol = D)

  for (d in 1:D) {
    # dim = 2*J, first row for shape, second row for rate
    shape_rate_param <- vapply(1:J, function(j) {
      shape <- sum(Z[[d]]==j)+alpha*P[j]

      log_rbf_value <- -(t[[d]]-t_star_J_D[j,d])^2/2/sigma_star_2_J_D[j,d]
      log_temp <- log(Xi_C_D[[d]])+log_rbf_value
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
allocation_variables_update_gkernel <- function(Y, t, mu_star_1_J, phi_star_1_J, Beta, Q_J_D,
                                                t_star_J_D, sigma_star_2_J_D){

  D <- length(Y)
  C_d <- unlist(lapply(Y, ncol))
  J <- nrow(mu_star_1_J)

  # Set up the list to save updated values
  Z <- NULL
  for(d in 1:D){
    Z[[d]] <- rep(0,C_d[d])
  }

  for(d in 1:D){

    P_d <- compute_probability_cpp(Y_d = Y[[d]], Beta = Beta[[d]], mu = mu_star_1_J, phi = phi_star_1_J, Q_J = Q_J_D[,d],
                                   t_d = t[[d]], t_star_J = t_star_J_D[,d], sigma_star_2_J = sigma_star_2_J_D[,d])

    Z_d <- apply(P_d, 1, function(x) extraDistr::rcat(n = 1, prob = x))
    Z[[d]] <- Z_d
  }
  return(Z)
}

##----------- Simulation of latent variables U_C_J_D for updating kernel parameters----------

U_C_J_D_update_gkernel <- function(Xi_C_D, Q_J_D, C_d, t, t_star_J_D, sigma_star_2_J_D, rbf){

  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)

  # Set up the list to save updated values
  U_C_J_D <- NULL
  U_C_J_D[[1]] <- matrix(NA, nrow = C_d[1], ncol = J)
  U_C_J_D[[2]] <- matrix(NA, nrow = C_d[2], ncol = J)

  for (d in 1:D) {
    for (j in 1:J){
      K_C_J_D <- exp(-Xi_C_D[[d]]*Q_J_D[j,d]*rbf(t[[d]], t_star_J_D[j,d], sigma_star_2_J_D[j,d]))
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
      ind <- c(1:C_d[d])[-log(U_C_J_D[[d]][,j])<Xi_C_D[[d]]*Q_J_D[j,d]]
      if(length(ind)==0) {
        truncate <- FALSE
      }else{
        # If truncate, find truncated regions
        truncate <- TRUE
        i_all_complement <- lapply(ind, function(c) {
          temp <- -2*sigma_star_2_J_D[j,d]*(log(-log(U_C_J_D[[d]][c,j]))-log(Xi_C_D[[d]][c])-log(Q_J_D[j,d]))

          i_complement <- matrix(c(t[[d]][c]-sqrt(temp),t[[d]][c]+sqrt(temp)),nrow=1)
          return(i_complement)
        })

        i_all <- intervals::interval_complement(intervals::interval_union(intervals::Intervals(do.call(rbind,i_all_complement))))@.Data
      }

      # Parameters in posterior distribution (Normal)
      ind1 <- (Z[[d]]==j)
      sum_t_j_d <- sum(t[[d]][ind1])
      N_j_d <- sum(ind1)
      r_j_hat <- (r_J[j]*sigma_star_2_J_D[j,d]+s_2*sum_t_j_d)/(sigma_star_2_J_D[j,d]+N_j_d*s_2)
      s_2_hat <- s_2*sigma_star_2_J_D[j,d]/(sigma_star_2_J_D[j,d]+N_j_d*s_2)

      if(!truncate){
        # No truncation, draw from the Normal posterior directly
        t_star[j,d] <- rnorm(1, mean = r_j_hat, sd = sqrt(s_2_hat))
      }else{
        N_interval <- nrow(i_all)
        # Compute the probability of lying in each interval
        log_p <- unlist(lapply(1:N_interval,function(l) {
          i <- i_all[l,]
          lower <- i[1]
          upper <- i[2]
          if(is.infinite(lower)){
            val <- pnorm(upper,mean = r_j_hat, sd = sqrt(s_2_hat),log.p = TRUE)
          }else if(is.infinite(upper)){
            val <- pnorm(lower,mean = r_j_hat, sd = sqrt(s_2_hat),log.p = TRUE, lower.tail = FALSE)
          }else{
            lp1 <- pnorm(upper,mean = r_j_hat, sd = sqrt(s_2_hat), log.p = TRUE)
            lp2 <- pnorm(lower,mean = r_j_hat, sd = sqrt(s_2_hat), log.p = TRUE)
            val <- lp1+log(1-exp(lp2-lp1))
          }
          return(val)
        }))
        log_K <- -max(log_p)

        # Select one truncated region (one interval)
        ind2 <- sample(1:N_interval, size = 1, prob = exp(log_p+log_K)) #breakpoint
        i_chosen <- i_all[ind2,]
        t_star[j,d] <- truncnorm::rtruncnorm(1, a = i_chosen[1],
                                             b = i_chosen[2],
                                             mean = r_j_hat, sd = sqrt(s_2_hat))
      }
    }
  }

  return(t_star)
}

##----------- Simulation of kernel parameters sigma_star_2_J_D-------------

# log(full conditional distribution)
sigma_star_2_log_prob <- function(sigma_star_2_j_d, m_2, h_j, t_star_j_d, t_c_d_sub){
  lprod <- -log(sigma_star_2_j_d)-(log(sigma_star_2_j_d)-h_j)^2/2/m_2-sum((t_c_d_sub-t_star_j_d)^2)/2/sigma_star_2_j_d
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
      ind <- c(1:C_d[d])[-log(U_C_J_D[[d]][,j])<Xi_C_D[[d]]*Q_J_D[j,d]]
      if(length(ind)==0) {
        truncate <- FALSE
      }else{
        truncate <- TRUE
        uppers <- -(t[[d]]-t_star_J_D[j,d])^2/2/(log(-log(U_C_J_D[[d]][,j]))-log(Xi_C_D[[d]])-log(Q_J_D[j,d]))
        upper <- min(uppers[ind])
      }

      # For empty clusters in dataset d, sample from the log-normal prior (may be truncated), always accept
      if(sum(Z[[d]]==j)==0) {
        if(!truncate){
          # No truncation
          sigma_star_2_J_D_new[j,d] <- rlnorm(1, meanlog = h_J[j], sdlog = sqrt(m_2))
          X_new[j,d] <- log(sigma_star_2_J_D_new[j,d])

          accept_count <- accept_count+1
        }else{
          # Truncate
          # X = -log(1/sigma_star - 1/upper)
          p.upper <- plnorm(upper, meanlog = h_J[j], sdlog = sqrt(m_2))
          u <- runif(1, 0, p.upper)
          sigma_star_2_J_D_new[j,d] <- qlnorm(u, meanlog = h_J[j], sdlog = sqrt(m_2))

          X_new[j,d] <- -log(1/sigma_star_2_J_D_new[j,d] - 1/upper)

          accept_count <- accept_count+1
        }
      }else{
        # Non-empty clusters, use adaptive MH
        if(!truncate){
          # No truncation, perform adaptive MH directly, propose sigma_star_j_d from log-normal, X=log(sigma_star_2)
          if(n <= 100){
            X_new[j,d] <- rnorm(n = 1, mean = log(sigma_star_2_J_D_old[j,d]), sd = 0.1)
          }else{
            X_new[j,d] <- rnorm(n = 1, mean = log(sigma_star_2_J_D_old[j,d]), sd = sqrt(2.4^2*variance_old[j,d]+2.4^2*MH.variance))
          }

          # Transform the new value of X back to new value of sigma_star_2_j_d
          sigma_star_2_J_D_new[j,d] <- exp(X_new[j,d])

          # Subset of t such that Z[[d]][c]=j
          t_c_d_sub <- t[[d]][Z[[d]]==j]

          # Compute log acceptance probability
          log_acceptance <- sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_new[j,d], m_2 = m_2,
                                                  h_j = h_J[j], t_star_j_d = t_star_J_D[j,d],
                                                  t_c_d_sub = t_c_d_sub) -
            sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_old[j,d], m_2,
                                  h_J[j], t_star_J_D[j,d], t_c_d_sub) -
            log(sigma_star_2_J_D_old[j,d]) + log(sigma_star_2_J_D_new[j,d])
          acceptance_sigma <- exp(log_acceptance)
          acceptance_sigma <- min(1, acceptance_sigma)

          # Update sigma_star_2_J_D_new, X_new, X_mean_new, variance_new
          outcome <- rbinom(n = 1, size = 1, prob = acceptance_sigma)
          if(is.na(outcome) == TRUE | outcome == 0){
            X_new[j,d] <- log(sigma_star_2_J_D_old[j,d])
            sigma_star_2_J_D_new[j,d] <- sigma_star_2_J_D_old[j,d]
          }else{
            accept_count <- accept_count+1
          }

        }else{
          # Truncate
          # X = -log(1/sigma_star - 1/upper)
          if(n <= 100){
            X_new[j,d] <- rnorm(n = 1, mean = -log(1/sigma_star_2_J_D_old[j,d] - 1/upper), sd = 0.1)
          }else{
            X_new[j,d] <- rnorm(n = 1, mean = -log(1/sigma_star_2_J_D_old[j,d] - 1/upper), sd = sqrt(2.4^2*variance_old[j,d]+2.4^2*MH.variance))
          }

          # Transform the new value of X back to new value of sigma_star_2_j_d
          sigma_star_2_J_D_new[j,d] <- 1/(exp(-X_new[j,d])+1/upper)
          if(sigma_star_2_J_D_new[j,d]==0) {
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
            log(sigma_star_2_J_D_old[j,d]) - log(upper-sigma_star_2_J_D_old[j,d]) +
            log(sigma_star_2_J_D_new[j,d]) + log(upper-sigma_star_2_J_D_new[j,d])

          acceptance_sigma <- exp(log_acceptance)
          acceptance_sigma <- min(1, acceptance_sigma)

          # Update sigma_star_2_J_D_new, X_new, X_mean_new, variance_new
          outcome <- rbinom(n = 1, size = 1, prob = acceptance_sigma)
          if(is.na(outcome) == TRUE | outcome == 0){
            X_new[j,d] <- -log(1/sigma_star_2_J_D_old[j,d] - 1/upper)
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

  return(list(sigma_star_2_J_D_new=sigma_star_2_J_D_new,
              X_mean_new=X_mean_new, M_2_new=M_2_new,
              variance_new=variance_new, accept=accept_count))
}

##----------- Simulation of hyper parameters r_J (in prior for t_star) ----------------
r_J_update_gkernel <- function(t_star_J_D, mu_r, sigma_r, s_2){
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
s_2_update_gkernel <- function(t_star_J_D, eta_1, eta_2, r_J){
  J <- nrow(t_star_J_D)
  D <- ncol(t_star_J_D)

  # Update
  shape <- J*D/2+eta_1

  temp <- sapply(1:D, function(d) {
    return(t_star_J_D[,d]-r_J)
  })

  rate <- eta_2+sum(temp^2)/2

  s_2 <- extraDistr::rinvgamma(1, alpha = shape, beta = rate)

  return(s_2)
}

##----------- Simulation of hyper parameters h_J (in prior for sigma_star_2) ----------------
h_J_update_gkernel <- function(sigma_star_2_J_D, mu_h, sigma_h, m_2){
  J <- nrow(sigma_star_2_J_D)
  D <- ncol(sigma_star_2_J_D)

  # Update
  h_J <- sapply(1:J, function(j) {
    mu_h_hat <- (mu_h*m_2+sigma_h^2*sum(log(sigma_star_2_J_D[j,])))/(m_2+D*sigma_h^2)
    sigma_h_2_hat <- sigma_h^2*m_2/(m_2+D*sigma_h^2)
    return(rnorm(1, mean = mu_h_hat, sd = sqrt(sigma_h_2_hat)))
  })
}

##----------- Simulation of hyper parameters m_2 (in prior for sigma_star_2) ----------------
m_2_update_gkernel <- function(sigma_star_2_J_D, kappa_1, kappa_2, h_J){
  J <- nrow(sigma_star_2_J_D)
  D <- ncol(sigma_star_2_J_D)

  # Update
  shape <- J*D/2+kappa_1

  temp <- sapply(1:D, function(d) {
    return(log(sigma_star_2_J_D[,d])-h_J)
  })

  rate <- kappa_2+sum(temp^2)/2

  m_2 <- extraDistr::rinvgamma(1, alpha = shape, beta = rate)

  return(m_2)
}

##----------- Simulation of Component probabilities P----------------

component_log_prob <- function(P, Q_J_D, alpha_0, alpha){

  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)

  lprod <- sum((alpha_0/J-1)*log(P))
  for(d in 1:D){
    lprod <- lprod + sum(alpha*P*log(Q_J_D[,d])-lgamma(alpha*P))
  }

  return(lprod)
}

component_probabilities_update <- function(P, Q_J_D, alpha_0, alpha, covariance,
                                           mean_x, tilde_s, iter_num, sd_P, MH.variance){

  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)


  # Define the inputs
  P_old <- P
  covariance_old <- covariance
  mean_x_old <- mean_x
  tilde_s_old <- tilde_s
  sd_P_old <- sd_P

  X_old <- log(P_old[1:(J-1)]/P_old[J]) # Length = J-1

  # Iteration index
  n <- iter_num

  # Adaptive step
  if(n <= 100){
    X_new <- mvtnorm::rmvnorm(n = 1, mean = X_old, sigma = 0.01*diag(x=1,nrow = J-1, ncol = J-1))
  }else{
    X_new <- mvtnorm::rmvnorm(n = 1, mean = X_old, sigma = sd_P_old*(covariance_old + MH.variance*diag(1, nrow = J-1, ncol = J-1)))
  }

  # Compute P_new (Length = J) from X_new
  P_new <- c(exp(X_new)/(1+sum(exp(X_new))),1/(1+sum(exp(X_new))))

  # Compute acceptance probability
  log_acceptance <- component_log_prob(P_new, Q_J_D, alpha_0, alpha) -
    component_log_prob(P_old, Q_J_D, alpha_0, alpha) +
    sum(log(P_new)-log(P_old))
  acceptance_P <- exp(log_acceptance)
  acceptance_P <- min(1, acceptance_P)

  # Update adaptive scaling parameter
  sd_P_new <- exp(log(sd_P_old)+n^(-0.7)*(acceptance_P-0.234))
  if(sd_P_new>exp(50)) {sd_P_new <- exp(50)}
  if(sd_P_new<exp(-50)) {sd_P_new <- exp(-50)}

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
              covariance_new=covariance_new, accept=accept,
              accept_prob=acceptance_P, sd_P_new=sd_P_new))
}

## ----------------------- Simulation of Concentration parameter alpha -----------------------
alpha_log_prob <- function(Q_J_D, P, alpha){

  D <- ncol(Q_J_D)
  # Construct the log-probability
  lprod <- -alpha + sum(vapply(1:D, function(d) {
    sum(alpha*P*log(Q_J_D[,d])-lgamma(alpha*P))
  }, FUN.VALUE = numeric(1)))

  return(lprod)
}

alpha_update <- function(Q_J_D, P, alpha, X_mean, M_2, variance, iter_num, MH.variance){

  D <- ncol(Q_J_D)

  # Define the inputs
  alpha_old <- alpha; X_old <- log(alpha_old)
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
    log(alpha_new) - log(alpha_old)
  acceptance_alpha <- exp(log_acceptance)
  acceptance_alpha <- min(1, acceptance_alpha)

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

  lprob <- -alpha_0 + lgamma(alpha_0) - J*lgamma(alpha_0/J) + sum(alpha_0/J*log(P))

  # Return the log-probability
  return(lprob)
}

alpha_0_update <- function(P, alpha_0, X_mean, M_2, variance, iter_num, MH.variance){

  J <- length(P)

  # Define the inputs
  alpha_0_old <- alpha_0; X_old <- log(alpha_0_old)
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
    log(alpha_0_new) - log(alpha_0_old)
  acceptance_alpha0 <- exp(log_acceptance)
  acceptance_alpha0 <- min(1, acceptance_alpha0)

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
  if(!quadratic){
    for(j in 1:J) {
      mu_tilde_j <- matrix(c(rep(1,G),log(mu_star_1_J[j,])), ncol=2, nrow=G)
      # For posterior variance of b
      parameter0 <- parameter0 + t(mu_tilde_j)%*%(mu_tilde_j)
      # For posterior mean of b
      parameter1 <- parameter1 + t(mu_tilde_j)%*%(matrix(log(phi_star_1_J[j,]), nrow=G, ncol=1))
      # For scale parameter of IG prior for alpha_phi^2
      parameter2 <- parameter2 + (matrix(log(phi_star_1_J[j,]), nrow=1, ncol=G))%*%
        (matrix(log(phi_star_1_J[j,]), nrow=G, ncol=1))
    }
  }else{
    for(j in 1:J){
      mu_tilde_j <- matrix(c(rep(1,G),log(mu_star_1_J[j,]),(log(mu_star_1_J[j,]))^2), ncol=3, nrow=G)
      # For posterior variance of b
      parameter0 <- parameter0 + t(mu_tilde_j)%*%(mu_tilde_j)
      # For posterior mean of b
      parameter1 <- parameter1 + t(mu_tilde_j)%*%(matrix(log(phi_star_1_J[j,]), nrow=G, ncol=1))
      # For scale parameter of IG prior for alpha_phi^2
      parameter2 <- parameter2 + (matrix(log(phi_star_1_J[j,]), nrow=1, ncol=G))%*%
        (matrix(log(phi_star_1_J[j,]), nrow=G, ncol=1))
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
  b <- mvtnorm::rmvnorm(n = 1, mean = m_tilde_b, sigma = alpha_phi_2*(V_tilde_b))

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

  if(!quadratic){
    lprod1 <- -log(mu_star*phi_star) - 1/(2*alpha_mu_2)*(log(mu_star))^2 -
      (log(phi_star)-(b[1]+b[2]*log(mu_star)))^2/(2*alpha_phi_2)
  }else{
    lprod1 <- -log(mu_star*phi_star) - 1/(2*alpha_mu_2)*(log(mu_star))^2 -
      (log(phi_star)-(b[1]+b[2]*log(mu_star)+b[3]*(log(mu_star))^2))^2/(2*alpha_phi_2)
  }

  if(length(Y_c_g_d_j) != 0){
    # Compute the log-probability

    lprod2 <- sum(phi_star*log(phi_star/(mu_star*B_c_d_j+phi_star)) +
                    lgamma(Y_c_g_d_j + phi_star) - lgamma(phi_star) - lgamma(Y_c_g_d_j+1) +
                    Y_c_g_d_j*log(mu_star/(mu_star*B_c_d_j+phi_star)))
    lprod <- lprod1 + lprod2

  }else{
    lprod <- lprod1
  }

  # Return the log-probability
  return(list(lprod1,lprod2))
}

unique_parameters_update <- function(mu_star_1_J, phi_star_1_J, mean_X_mu_phi, tilde_s_mu_phi,
                                     J, G, Z, b, alpha_phi_2, Beta, alpha_mu_2, covariance, iter_num, quadratic, Y,
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
      X_mu_phi_star_old <- c(log(mu_star_1_J[j,g]),log(phi_star_1_J[j,g]))
      mu_star_old <- mu_star_1_J[j,g]; phi_star_old <- phi_star_1_J[j,g]

      # If cluster j is empty with respect to all datasets, sample mu and phi from their priors
      if(!any(unlist(Z)==j)){
        mu_star_new <- rlnorm(1, meanlog = 0, sdlog = sqrt(alpha_mu_2))
        if(quadratic){
          phi_star_new <- rlnorm(1, meanlog = b[1]+b[2]*log(mu_star_new)+b[3]*log(mu_star_new)^2, sdlog = sqrt(alpha_phi_2))
        }else{
          phi_star_new <- rlnorm(1, meanlog = b[1]+b[2]*log(mu_star_new), sdlog = sqrt(alpha_phi_2))
        }
        X_mu_phi_star_new <- c(log(mu_star_new), log(phi_star_new))
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
                                                  Z, b, alpha_phi_2, Y, Beta, j, g, alpha_mu_2,  quadratic = quadratic))) -
            log(mu_star_old) -log(phi_star_old) +
            log(mu_star_new) + log(phi_star_new)

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

capture_efficiencies_log_prob <- function(Beta_d, Y, Z, mu_star_1_J, phi_star_1_J, d, a_d_beta, b_d_beta){

  j <- Z[[d]]
  # Construct the log-probability
  lprod1 <- (a_d_beta[d]-1)*log(Beta_d) + (b_d_beta[d]-1)*log(1-Beta_d)

  # mu x beta: C x G
  mu_star_beta <- apply(mu_star_1_J[j,], 2, function(x) x*Beta_d)

  # Y x log(beta): C x G
  Y_log_beta <- apply(t(Y[[d]]), 2, function(x) x*log(Beta_d))

  lprod2 <- rowSums((phi_star_1_J[j,]+t(Y[[d]]))*log(phi_star_1_J[j,] + mu_star_beta) - Y_log_beta)

  lprod <- lprod1 - lprod2

  # Return the log-probability
  return(lprod)
}


capture_efficiencies_update <- function(Beta, Y, Z, mu_star_1_J, phi_star_1_J, a_d_beta, b_d_beta,
                                        iter_num, M_2, mean_X, variance,
                                        MH.variance){

  D <- length(Z)
  C_d <- unlist(lapply(Z, length))

  M_2_old <- M_2; mean_X_old <- mean_X; variance_old <- variance

  # To store outputs
  Beta_new <- NULL
  variance_new <- NULL
  mean_X_new <- NULL
  M_2_new <- NULL
  X_new <- NULL

  n <- iter_num

  accept_count <- 0

  for(d in 1:D){
    Beta_d_old <- Beta[[d]]
    X_d_old <- log(Beta_d_old/(1-Beta_d_old))

    if(n <= 100){
      X_d_new <- rnorm(n = length(X_d_old), mean = X_d_old, sd = 0.1)
    }else{
      X_d_new <- rnorm(n = length(X_d_old), mean = X_d_old, sd = sqrt((2.4^2)*variance_old[[d]] + 2.4^2*MH.variance))
    }

    Beta_d_new <- 1/(1+exp(-X_d_new))
    if(any(Beta_d_new==1)) {
      Beta_d_new[Beta_d_new==1] <- 1-.Machine$double.eps
    }
    # Compute the acceptance probability
    acceptance_prob_log <- capture_efficiencies_log_prob(Beta_d = Beta_d_new, Y, Z, mu_star_1_J,
                                                         phi_star_1_J,d,a_d_beta,b_d_beta) -
      capture_efficiencies_log_prob(Beta_d = Beta_d_old, Y, Z, mu_star_1_J, phi_star_1_J,
                                    d,a_d_beta,b_d_beta) + log(Beta_d_new) + log(1-Beta_d_new) -
      log(Beta_d_old) - log(1-Beta_d_old)


    outcome <- sapply(1:C_d[d], function(c) {
      rbinom(n = 1,size = 1,prob = min(1, exp(acceptance_prob_log[c])))
    })

    accept.num <- 0
    for(c in 1:C_d[d]){

      if(is.na(outcome[c]) == TRUE | outcome[c] == 0){
        X_d_new[c] <- X_d_old[c]
        Beta_d_new[c] <- Beta_d_old[c]
      }else{
        accept.num <- accept.num + 1
      }
    }

    Beta_new[[d]] <- Beta_d_new
    X_new[[d]] <- X_d_new
    accept_count <- accept_count+accept.num

    mean_X_new[[d]] <- (1-1/n)*mean_X_old[[d]] + (1/n)*X_new[[d]]
    M_2_new[[d]] <- M_2_old[[d]] + (X_new[[d]]-mean_X_old[[d]])*(X_new[[d]]-mean_X_new[[d]])
    variance_new[[d]] <- 1/(n-1)*M_2_new[[d]]
  }

  return(list(Beta_new=Beta_new, accept_count=accept_count, mean_X_new=mean_X_new,
              M_2_new=M_2_new, variance_new=variance_new))
}
