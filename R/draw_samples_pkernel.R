# ---------------- For clustering calcium imaging data ---------------------------------
# update alpha, alpha_0, component probabilities P is exactly the same as the steps in clustering single-cell
##--------------- Simulation of latent Xi ----------------

Xi_C_D_update_pkernel <- function(Q_J_D, C_d, t, mu_J_D, lambda_J_D, sigma_2_J_D){
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)

  # Set up the list to save updated values
  Xi <- NULL
  Xi[[1]] <- rep(NA, C_d[1])
  Xi[[2]] <- rep(NA, C_d[2])

  # Draw from Gamma (full conditional distribution)
  for (d in 1:D) {

    rates <- vapply(1:C_d[d], function(cc) {
      # vector of length = J
      val1 <- -2/sigma_2_J_D[,d]
      val2 <- sin((t[[d]][cc]-mu_J_D[,d])/lambda_J_D[,d])
      log_pkernel_value <- val1*val2^2

      # vector of length = J
      log_temp <- log(Q_J_D[,d])+log_pkernel_value
      log_K <- max(log_temp)

      return(exp(log_K)*sum(exp(log_temp-log_K)))
    }, FUN.VALUE = numeric(1))

    Xi[[d]] <- rgamma(C_d[d], shape = 1, rate = rates)
  }

  return(Xi)
}

##--------------- Simulation of dataset-specific vector q_j_d ----------------

Q_J_D_update_pkernel <- function(Z, alpha, P, Xi_C_D, t, mu_J_D, lambda_J_D, sigma_2_J_D){
  J <- nrow(mu_J_D)
  D <- ncol(mu_J_D)

  # Set up the matrix to save updated values
  Q <- matrix(NA, nrow = J, ncol = D)

  for (d in 1:D) {
    # dim = 2*J, first row for shape, second row for rate
    shape_rate_param <- vapply(1:J, function(j) {
      shape <- sum(Z[[d]]==j)+alpha*P[j]

      # vector of length = C_d[[d]]
      val1 <- -2/sigma_2_J_D[j,d]
      val2 <- sin((t[[d]]-mu_J_D[j,d])/lambda_J_D[j,d])
      log_pkernel_value <- val1*val2^2

      # vector of length = C_d[[d]]
      log_temp <- log(Xi_C_D[[d]])+log_pkernel_value
      log_K <- max(log_temp)

      rate <- 1+exp(log_K)*sum(exp(log_temp-log_K))

      return(c(shape, rate))
    }, FUN.VALUE = numeric(2))

    Q[,d] <- rgamma(J, shape = shape_rate_param[1,], rate = shape_rate_param[2,])

    # replace zero values with mean
    if(any(Q[,d]==0)) {
      Q[Q[,d]==0,d] <- shape_rate_param[1,Q[,d]==0]/shape_rate_param[2,Q[,d]==0]
    }
  }

  return(Q)
}

##--------------- Simulation of allocations ----------------
allocation_variables_update_pkernel <- function(Y, X, t, L_1_J, Sigma_1_J, Q_J_D,
                                        mu_J_D, lambda_J_D, sigma_2_J_D){

  D <- length(Y)
  C_d <- unlist(lapply(Y, nrow))
  J <- nrow(Q_J_D)

  # Set up the list to save updated values
  Z <- NULL
  for(d in 1:D){
    Z[[d]] <- rep(0,C_d[d])
  }

  means_by_j <- NULL
  for(d in 1:D){
    means_by_j[[d]] <- lapply(1:J,function(j) {
      X[[d]] %*% L_1_J[[j]]
    })
  }

  for(d in 1:D){

    # n *J
    Q_mat <- matrix(rep(Q_J_D[,d],C_d[d]),nrow=C_d[d],byrow = TRUE)
    mu_mat <- matrix(rep(mu_J_D[,d],C_d[d]),nrow=C_d[d],byrow = TRUE)
    lambda_mat <- matrix(rep(lambda_J_D[,d],C_d[d]),nrow=C_d[d],byrow = TRUE)
    sigma_2_mat <- matrix(rep(sigma_2_J_D[,d],C_d[d]),nrow=C_d[d],byrow = TRUE)
    t_mat <- matrix(rep(t[[d]],J),nrow=C_d[d],byrow = FALSE)

    # n * J
    LP1 <- -2/sigma_2_mat*(sin((t_mat-mu_mat)/lambda_mat))^2+log(Q_mat)
    LP2 <- lapply(1:J,function(j) {
      temp <- mniw::dmNorm(as.matrix(Y[[d]]), mu=as.matrix(means_by_j[[d]][[j]]), Sigma = Sigma_1_J[[j]], log = TRUE)
      return(temp)
    })
    LP2 <- do.call(cbind,LP2)
    LP <- LP1+LP2

    nc <- -apply(LP,1,max) # length=n
    # n * J
    LP_plus_nc <- LP+matrix(rep(nc,J),ncol=J,byrow = FALSE)


    P <- t(apply(LP_plus_nc, 1, function(x) {
      exp(x)/sum(exp(x))
    }))

    Z[[d]] <- apply(P, 1, function(x) extraDistr::rcat(1, prob=x))
  }


  return(Z)
}


##--------------- Simulation of unique parameters --------------
unique_params_update <- function(Y,X,J,Z,L0,V0,Phi0,omega0,iter_num){

  D <- length(Y)

  V0_inv <- solve(V0)
  temp <- Matrix::crossprod(L0,V0_inv)%*%L0

  loop.result <- lapply(1:J,function(j) {

    N_j <- sum(unlist(Z)==j)

    if(N_j!=0){

      Y_j <- do.call(rbind,lapply(1:D,function(d) Y[[d]][Z[[d]]==j,]))
      X_j <- do.call(rbind,lapply(1:D,function(d) X[[d]][Z[[d]]==j,]))

      Vn <- Matrix::solve(Matrix::crossprod(X_j,X_j)+V0_inv)
      Ln <- Vn%*%(Matrix::crossprod(X_j,Y_j)+V0_inv%*%L0)

      Phi_n <- Phi0+Matrix::crossprod(Y_j,Y_j)+temp-Matrix::crossprod(Ln,solve(Vn))%*%Ln
      omega_n <- N_j+omega0

      Sigma_j_new <- MCMCpack::riwish(omega_n, Phi_n)
      L_j_new <- mniw::rMNorm(n=1,Lambda=Ln,SigmaR=Vn,SigmaC=Sigma_j_new)

    }else{
      # for empty clusters, draw from the prior
      Sigma_j_new <- MCMCpack::riwish(omega0, Phi0)
      L_j_new <- mniw::rMNorm(n=1,L0,V0,Sigma_j_new)

    }

    return(list(L_j_new=L_j_new,Sigma_j_new=Sigma_j_new))

  })


  L_new <- lapply(loop.result, function(l) l$L_j_new)
  Sigma_new <- lapply(loop.result, function(l) l$Sigma_j_new)

  return(list(L_new=L_new, Sigma_new=Sigma_new))
}



##----------- Simulation of latent variables U_C_J_D for updating kernel parameters----------

U_C_J_D_update_pkernel <- function(Xi_C_D, Q_J_D, C_d, t, mu_J_D, lambda_J_D, sigma_2_J_D, pkernel){
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)

  # Set up the list to save updated values
  U_C_J_D <- NULL
  U_C_J_D[[1]] <- matrix(NA, nrow = C_d[1], ncol = J)
  U_C_J_D[[2]] <- matrix(NA, nrow = C_d[2], ncol = J)

  for (d in 1:D) {
    for (j in 1:J){
      K_C_J_D <- exp(-Xi_C_D[[d]]*Q_J_D[j,d]*pkernel(t[[d]], mu_J_D[j,d], lambda_J_D[j,d], sigma_2_J_D[j,d]))
      U_C_J_D[[d]][,j] <- runif(C_d[d], min = 0, max = K_C_J_D)
    }
  }

  return(U_C_J_D)

}


##----------- Simulation of kernel parameters mu_J_D-------------
mu_log_prob <- function(mu_j_d, lambda_j_d, sigma_2_j_d, t_c_d_sub, t_d, pkernel, xi_d, q_j_d){

  # From likelihood
  if(any(is.na(t_c_d_sub))){
    lp1 <- 0
  }else{
    val1 <- -2/sigma_2_j_d
    val2 <- sin((t_c_d_sub-mu_j_d)/lambda_j_d)
    log_pkernel_value <- val1*val2^2

    lp1 <- sum(log_pkernel_value)
  }

  lp2 <- sum(-xi_d*q_j_d*pkernel(t_d,mu_j_d,lambda_j_d,sigma_2_j_d))

  lprod <- lp1+lp2

  return(lprod)
}

mu_J_D_update <- function(mu_J_D, Z, t, lambda_J_D, sigma_2_J_D, Xi_C_D, Q_J_D, pkernel,
                          sd_mu, iter_num, target_accept){
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)

  # All are JxD matrices
  mu_J_D_old <- mu_J_D
  mu_lower <- -pi * lambda_J_D/2; mu_upper <- pi * lambda_J_D/2
  X_old <- log(mu_J_D_old - mu_lower) - log(mu_upper - mu_J_D_old)
  sd_mu_old <- sd_mu

  # Save updated mu_J_D
  X_new <- matrix(NA, nrow=J, ncol=D)
  mu_J_D_new <- matrix(NA, nrow=J, ncol=D)
  sd_mu_new <- matrix(NA, nrow=J, ncol=D)

  # Iteration index
  n <- iter_num

  accept_count <- 0

  for (d in 1:D) {
    for (j in 1:J) {
      # Apdative MH (algorithm 5)
      X_new[j,d] <- rnorm(n=1, mean=X_old[j,d], sd=sqrt(sd_mu_old[j,d]))

      # Transform the new value of X back to new value of mu_j_d
      mu_J_D_new[j,d] <- mu_upper[j,d] + (mu_lower[j,d] - mu_upper[j,d]) / (1 + exp(X_new[j,d]))

      if(mu_J_D_new[j,d] <= mu_lower[j,d]) {
        mu_J_D_new[j,d] <- mu_lower[j,d] + .Machine$double.eps
      }

      if(mu_J_D_new[j,d] >= mu_upper[j,d]) {
        mu_J_D_new[j,d] <- mu_upper[j,d] - .Machine$double.eps
      }

      # Subset of t such that Z[[d]][c]=j
      if(any(Z[[d]] == j)) {
        t_c_d_sub <- t[[d]][Z[[d]] == j]
      } else{
        t_c_d_sub <- NA
      }

      # Compute log acceptance probability
      log_acceptance <- mu_log_prob(mu_j_d=mu_J_D_new[j,d],lambda_j_d=lambda_J_D[j,d],
                                    sigma_2_j_d=sigma_2_J_D[j,d],t_c_d_sub=t_c_d_sub,t_d=t[[d]],pkernel=pkernel,
                                    xi_d=Xi_C_D[[d]],q_j_d=Q_J_D[j,d]) -
        mu_log_prob(mu_j_d=mu_J_D_old[j,d], lambda_J_D[j,d], sigma_2_J_D[j,d], t_c_d_sub, t[[d]], pkernel, Xi_C_D[[d]], Q_J_D[j,d]) +
        log(mu_J_D_new[j,d] - mu_lower[j,d]) + log(mu_upper[j,d] - mu_J_D_new[j,d]) -
        log(mu_J_D_old[j,d] - mu_lower[j,d]) - log(mu_upper[j,d] - mu_J_D_old[j,d])

      acceptance_mu <- exp(log_acceptance)
      acceptance_mu <- min(1, acceptance_mu)

      # Decision
      outcome <- rbinom(n=1, size=1, prob=acceptance_mu)
      if(is.na(outcome) == TRUE | outcome == 0) {
        X_new[j,d] <- X_old[j,d]
        mu_J_D_new[j,d] <- mu_J_D_old[j,d]
      } else{
        accept_count <- accept_count + 1
      }

      # Update scale parameter in the proposal distribution
      sd_mu_new[j,d] <- exp(log(sd_mu_old[j,d]) + n^(-0.7)*(acceptance_mu - target_accept))
      if(sd_mu_new[j,d]>exp(50)) {sd_mu_new[j,d] <- exp(50)}
      if(sd_mu_new[j,d]<exp(-50)) {sd_mu_new[j,d] <- exp(-50)}

    } # End for j in 1:J
  } # End for d in 1:D

  return(list(mu_J_D_new=mu_J_D_new, accept_count=accept_count,sd_mu_new=sd_mu_new))

}

##----------- Simulation of kernel parameters lambda_J_D-------------
lambda_log_prob <- function(lambda_j_d, r_j, s_2, mu_j_d, sigma_2_j_d, t_c_d_sub, t_d, pkernel, xi_d, q_j_d){

  # From likelihood
  if(any(is.na(t_c_d_sub))){
    lp1 <- 0
  }else{
    val1 <- -2/sigma_2_j_d
    val2 <- sin((t_c_d_sub-mu_j_d)/lambda_j_d)
    log_pkernel_value <- val1*val2^2

    lp1 <- sum(log_pkernel_value)
  }

  lp2 <- sum(-xi_d*q_j_d*pkernel(t_d,mu_j_d,lambda_j_d,sigma_2_j_d))

  lprod <- lp1-2*log(lambda_j_d)-(log(lambda_j_d)-r_j)^2/2/s_2+lp2

  return(lprod)
}

lambda_J_D_update <- function(lambda_J_D, r_J, s_2, Z, t, mu_J_D, sigma_2_J_D, Xi_C_D, Q_J_D, pkernel,
                              sd_lambda, iter_num, target_accept){
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)

  # All are JxD matrices
  lambda_J_D_old <-lambda_J_D
  lambda_lower <- 2*abs(mu_J_D)/pi
  X_old <- log(lambda_J_D_old-lambda_lower)
  sd_lambda_old <- sd_lambda

  # Save updated lambda_J_D
  X_new <- matrix(NA, nrow=J, ncol=D)
  lambda_J_D_new <- matrix(NA, nrow=J, ncol=D)
  sd_lambda_new <- matrix(NA, nrow=J, ncol=D)

  # Iteration index
  n <- iter_num

  accept_count <- 0

  for (d in 1:D) {
    for (j in 1:J) {
      # Apdative MH (algorithm 5)
      X_new[j,d] <- rnorm(n=1, mean=X_old[j,d], sd=sqrt(sd_lambda_old[j,d]))

      # Transform the new value of X back to new value of lambda_j_d
      lambda_J_D_new[j,d] <- lambda_lower[j,d] + exp(X_new[j,d])

      if(lambda_J_D_new[j,d] <= lambda_lower[j,d]) {
        lambda_J_D_new[j,d] <- lambda_lower[j,d] + .Machine$double.eps
      }

      # Subset of t such that Z[[d]][c]=j
      if(any(Z[[d]] == j)) {
        t_c_d_sub <- t[[d]][Z[[d]] == j]
      } else{
        t_c_d_sub <- NA
      }

      # Compute log acceptance probability
      log_acceptance <- lambda_log_prob(lambda_j_d=lambda_J_D_new[j,d], r_j=r_J[j], s_2=s_2,
                                        mu_j_d=mu_J_D[j,d],sigma_2_j_d=sigma_2_J_D[j,d],t_c_d_sub=t_c_d_sub,
                                        t_d=t[[d]],pkernel=pkernel,xi_d=Xi_C_D[[d]],q_j_d=Q_J_D[j,d]) -
        lambda_log_prob(lambda_j_d=lambda_J_D_old[j,d], r_J[j], s_2, mu_J_D[j,d], sigma_2_J_D[j,d], t_c_d_sub,
                        t[[d]], pkernel, Xi_C_D[[d]], Q_J_D[j,d]) +
        log(lambda_J_D_new[j,d] - lambda_lower[j,d]) - log(lambda_J_D_old[j,d] - lambda_lower[j,d])

      acceptance_lambda <- exp(log_acceptance)
      acceptance_lambda <- min(1, acceptance_lambda)

      # Decision
      outcome <- rbinom(n=1, size=1, prob=acceptance_lambda)
      if(is.na(outcome) == TRUE | outcome == 0) {
        X_new[j,d] <- X_old[j,d]
        lambda_J_D_new[j,d] <- lambda_J_D_old[j,d]
      } else{
        accept_count <- accept_count + 1
      }

      # Update scale parameter in the proposal distribution
      sd_lambda_new[j,d] <- exp(log(sd_lambda_old[j,d]) + n^(-0.7)*(acceptance_lambda - target_accept))
      if(sd_lambda_new[j,d]>exp(50)) {sd_lambda_new[j,d] <- exp(50)}
      if(sd_lambda_new[j,d]<exp(-50)) {sd_lambda_new[j,d] <- exp(-50)}


    } # End for j in 1:J
  } # End for d in 1:D


  return(list(lambda_J_D_new=lambda_J_D_new, accept_count=accept_count,sd_lambda_new=sd_lambda_new))

}

##----------- Simulation of kernel parameters sigma_2_J_D-------------
sigma_2_J_D_update <- function(h_J, m_2, Z, t, mu_J_D, lambda_J_D, U_C_J_D, Xi_C_D, Q_J_D, iter_num){

  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)
  C_d <- unlist(lapply(Z,length))

  # Save updated sigma_2_J_D
  sigma_2_J_D_new <- matrix(NA, nrow=J, ncol = D)

  for (d in 1:D) {
    for (j in 1:J) {
      # Compute truncation region if needed
      # Find which cells to truncate regions
      ind <- c(1:C_d[d])[-log(U_C_J_D[[d]][,j])<Xi_C_D[[d]]*Q_J_D[j,d]]
      if(length(ind)==0) {
        truncate <- FALSE
      }else{
        truncate <- TRUE

        val <- sin((t[[d]]-mu_J_D[j,d])/lambda_J_D[j,d])

        uppers <- -2*val^2/(log(-log(U_C_J_D[[d]][,j]))-log(Xi_C_D[[d]])-log(Q_J_D[j,d]))
        upper <- min(uppers[ind])
      }

      # Compute the updated parameter in IG
      if(any(Z[[d]] == j)) {
        t_c_d_sub <- t[[d]][Z[[d]] == j]
        val <- sin((t_c_d_sub-mu_J_D[j,d])/lambda_J_D[j,d])
        temp <- 2*sum(val^2)

      }else{
        temp <- 0
      }

      a_j <- 2+h_J[j]^2/m_2; b_j <- h_J[j]+h_J[j]^3/m_2
      if(truncate){
        # Draw from truncated IG, instead of traditional inverse CDF, use a improved inverse CDF to generate
        # generate X from truncated gamma (1/upper, +Inf), then get 1/X

        # Inverse CDF does not work when truncated region is of low probability, which may give p.upper==0 exactly

        a=a_j;b=b_j+temp
        p <- pgamma(1/upper,shape=a,rate=b,lower.tail = F,log.p = T)
        y <- rexp(1,1); v <- y-p
        x <- 1/qgamma(p=-v,shape=a,rate=b,lower.tail = F,log.p = T)
        sigma_2_J_D_new[j,d] <- x

        # Replace zero values with a very small value
        if(sigma_2_J_D_new[j,d]==0){
          sigma_2_J_D_new[j,d] <- .Machine$double.eps
        }

        if(sigma_2_J_D_new[j,d]>=upper){
          sigma_2_J_D_new[j,d] <- upper-.Machine$double.eps
        }
      }else{
        # Draw from IG directly
        sigma_2_J_D_new[j,d] <- extraDistr::rinvgamma(n=1, alpha = a_j, beta = b_j+temp)

        # Replace zero values with a very small value
        if(sigma_2_J_D_new[j,d]==0){
          sigma_2_J_D_new[j,d] <- .Machine$double.eps
        }
      }

    } # End for j in 1:J

  } # End for d in 1:D

  return(sigma_2_J_D_new)
}


##----------- Simulation of hyper parameters r_J (in prior for lambda) ----------------
r_J_update_pkernel <- function(lambda_J_D, mu_r, sigma_r, s_2){
  J <- nrow(lambda_J_D)
  D <- ncol(lambda_J_D)

  # Update
  r_J <- sapply(1:J, function(j) {
    mu_r_hat <- (mu_r*s_2+sigma_r^2*sum(log(lambda_J_D[j,])))/(s_2+D*sigma_r^2)
    sigma_r_2_hat <- sigma_r^2*s_2/(s_2+D*sigma_r^2)
    return(rnorm(1, mean = mu_r_hat, sd = sqrt(sigma_r_2_hat)))
  })

  return(r_J)
}

##----------- Simulation of hyper parameters s_2 (in prior for lambda) ----------------
s_2_update_pkernel <- function(lambda_J_D, eta_1, eta_2, r_J){
  J <- nrow(lambda_J_D)
  D <- ncol(lambda_J_D)

  # Update
  shape <- J*D/2+eta_1

  temp <- sapply(1:D, function(d) {
    return(log(lambda_J_D[,d])-r_J)
  })

  rate <- eta_2+sum(temp^2)/2

  s_2 <- extraDistr::rinvgamma(1, alpha = shape, beta = rate)

  return(s_2)
}

##----------- Simulation of hyper parameters h_J (in prior for sigma_2) ----------------
h_log_prob <- function(h_j,m_2,sigma_2_d, mu_h, sigma_h){

  a_j <- 2+h_j^2/m_2; b_j <- h_j+h_j^3/m_2
  D <- length(sigma_2_d)
  lp <- a_j*D*log(b_j)-D*lgamma(a_j)-(a_j+1)*sum(log(sigma_2_d))-b_j*sum(1/sigma_2_d)-log(h_j)-(log(h_j)-mu_h)^2/2/sigma_h^2

  return(lp)
}

h_J_update_pkernel <- function(h_J, sigma_2_J_D, mu_h, sigma_h, m_2, sd_h, iter_num, target_accept){

  J <- nrow(sigma_2_J_D)
  D <- ncol(sigma_2_J_D)

  h_J_old <- h_J
  X_old <- log(h_J_old)
  sd_h_old <- sd_h

  # Save updated h_J
  X_new <- rep(NA,J)
  h_J_new <- rep(NA,J)
  sd_h_new <- rep(NA,J)

  # Iteration index
  n <- iter_num

  accept_count <- 0

  for (j in 1:J) {
    # Apdative MH (algorithm 5)
    X_new[j] <- rnorm(n=1, mean=X_old[j], sd=sqrt(sd_h_old[j]))

    # Transform the new value of X back to new value of h_J
    h_J_new[j] <- exp(X_new[j])

    # Compute log acceptance probability
    log_acceptance <- h_log_prob(h_j=h_J_new[j], m_2=m_2, sigma_2_d=sigma_2_J_D[j,], mu_h=mu_h, sigma_h=sigma_h) -
      h_log_prob(h_j=h_J_old[j],m_2,sigma_2_J_D[j,],mu_h,sigma_h) +
      log(h_J_new[j]) - log(h_J_old[j])


    acceptance_h <- exp(log_acceptance)
    acceptance_h <- min(1, acceptance_h)

    # Decision
    outcome <- rbinom(n=1, size=1, prob=acceptance_h)
    if(is.na(outcome) == TRUE | outcome == 0) {
      X_new[j] <- X_old[j]
      h_J_new[j] <- h_J_old[j]
    } else{
      accept_count <- accept_count + 1
    }

    # Update scale parameter in the proposal distribution
    sd_h_new[j] <- exp(log(sd_h_old[j]) + n^(-0.7)*(acceptance_h - target_accept))
    if(sd_h_new[j]>exp(50)) {sd_h_new[j] <- exp(50)}
    if(sd_h_new[j]<exp(-50)) {sd_h_new[j] <- exp(-50)}


  } # End for j in 1:J


  return(list(h_J_new=h_J_new, accept_count=accept_count,sd_h_new=sd_h_new))
}

##----------- Simulation of hyper parameters m_2 (in prior for sigma_2) ----------------
m_2_log_prob <- function(m_2,h_J,sigma_2_J_D,kappa_1,kappa_2){

  J <- nrow(sigma_2_J_D)
  D <- ncol(sigma_2_J_D)

  temp <- sapply(1:J,function(j) {
    h_j <- h_J[j]; sigma_2_d <- sigma_2_J_D[j,]
    a_j <- 2+h_j^2/m_2; b_j <- h_j+h_j^3/m_2

    a_j*D*log(b_j)-D*lgamma(a_j)-(a_j+1)*sum(log(sigma_2_d))-b_j*sum(1/sigma_2_d)
  })

  lp <- sum(temp)-(kappa_1+1)*log(m_2)-kappa_2/m_2

}

m_2_update_pkernel <- function(m_2, sigma_2_J_D, kappa_1, kappa_2, h_J, sd_m_2, iter_num, target_accept){


  m_2_old <- m_2
  X_old <- log(m_2_old)
  sd_m_2_old <- sd_m_2

  # Iteration index
  n <- iter_num

  # Apdative MH (algorithm 5)
  X_new <- rnorm(n=1, mean=X_old, sd=sqrt(sd_m_2_old))

  # Transform the new value of X back to new value of m_2
  m_2_new <- exp(X_new)

  # Compute log acceptance probability
  log_acceptance <- m_2_log_prob(m_2=m_2_new,h_J=h_J,sigma_2_J_D=sigma_2_J_D,kappa_1=kappa_1,kappa_2=kappa_2) -
    m_2_log_prob(m_2=m_2_old,h_J,sigma_2_J_D,kappa_1,kappa_2) +
    log(m_2_new) - log(m_2_old)

  acceptance_m_2 <- exp(log_acceptance)
  acceptance_m_2 <- min(1, acceptance_m_2)

  # Decision
  outcome <- rbinom(n=1, size=1, prob=acceptance_m_2)
  if(is.na(outcome) == TRUE | outcome == 0) {
    X_new <- X_old
    m_2_new <- m_2_old
    accept <- 0
  } else{
    accept <- 1
  }

  # Update scale parameter in the proposal distribution
  sd_m_2_new <- exp(log(sd_m_2_old) + n^(-0.7)*(acceptance_m_2 - target_accept))
  if(sd_m_2_new>exp(50)) {sd_m_2_new <- exp(50)}
  if(sd_m_2_new<exp(-50)) {sd_m_2_new <- exp(-50)}

  return(list(m_2_new=m_2_new, accept=accept,sd_m_2_new=sd_m_2_new))

}





