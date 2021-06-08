find_theta_y <- function(y, alpha = 1.001, beta = 1.001){
  return(list(interval = pscl::betaHPD(alpha + sum(y), beta + length(y) - sum(y)),
              mean = (sum(y) + alpha) / (alpha + length(y) + beta ) ))
}

find_theta_z <- function(z_mat, num_mcmc, step_size, alpha = 1.001, beta = 1.001){
  ## initialize MCMC
  num_pooled <- dim(z_mat)[2]
  num_samples <- dim(z_mat)[1]
  burn <- num_mcmc / 10
  theta_samples <- rep(.5, num_mcmc)
  acceptance <- rep(0, num_mcmc)
  theta_star <- -1
  prob_star <- 1
  pooled_success <- sum(rowSums(z_mat) > 0)
  
  ## run MCMC
  for (iter in 2:num_mcmc){
    # propose
    while(theta_star < 0 | theta_star > 1 ){
      theta_star <- theta_samples[iter - 1] + rnorm(1,0,step_size)
      prob_star <- 1 - (1 - theta_star) ^ num_pooled
    }
    
    if (log(1 - prob_star) == -Inf & (num_samples - pooled_success) > 0){
      log_pi_star <- -Inf
    } else if (log(1 - prob_star) == -Inf & (num_samples - pooled_success) == 0){
      log_pi_star <- (alpha - 1) * log(theta_star)  + (beta - 1) * log(1 - theta_star)
    }  else{
      log_pi_star <- pooled_success * log(prob_star) + (num_samples - pooled_success) * 
        log(1 - prob_star) + (alpha - 1) * log(theta_star)  + (beta - 1) * log(1 - theta_star)
    }
    
    prob_old <- 1 - (1 - theta_samples[iter - 1]) ^ num_pooled
    theta_old <- theta_samples[iter - 1]
    if (log(1 - prob_old) == -Inf & (num_samples - pooled_success) > 0){
      log_pi_old <- -Inf
    } else if (log(1 - prob_old) == -Inf & (num_samples - pooled_success) == 0){
      log_pi_old <- (alpha - 1) * log(theta_old)  + (beta - 1) * log(1 - theta_old)
    }  else{
      log_pi_old <- pooled_success * log(prob_old) + (num_samples - pooled_success) * 
        log(1 - prob_old) + (alpha - 1) * log(theta_old)  + (beta - 1) * log(1 - theta_old)
    }
    
    # acceptance
    if ((log_pi_star - log_pi_old) > log(runif(1))){
      acceptance[iter] <- 1
      theta_samples[iter] <- theta_star
      theta_star <- -1
    } else {
      theta_samples[iter] <- theta_samples[iter - 1]
      theta_star <- -1
    }
  }
  return(list(interval = hdi(theta_samples[burn:num_mcmc]),theta = theta_samples[burn:num_mcmc], acceptance_rate = mean(acceptance)))
}

find_theta_z2 <- function(z_pooled, n_pooled, num_mcmc, step_size, alpha = 1.001, beta = 1.001){
  ## initialize MCMC
  burn <- num_mcmc / 10
  theta_samples <- rep(.5, num_mcmc)
  acceptance <- rep(0, num_mcmc)
  theta_star <- -1
  prob_star <- 1
  
  ## run MCMC
  for (iter in 2:num_mcmc){
    # propose
    while(theta_star < 0 | theta_star > 1 ){
      theta_star <- theta_samples[iter - 1] + rnorm(1,0,step_size)
    }
    prob_star <-  1 - (1 - theta_star) ^ n_pooled
    
    log_pi_star <- sum(z_pooled * log(prob_star) + (1 - z_pooled) * 
                         log(1 - prob_star) + (alpha - 1) * log(theta_star)  + (beta - 1) * log(1 - theta_star))
    if (is.na(log_pi_star)){
      log_pi_star <- (alpha - 1) * log(theta_star)  + (beta - 1) * log(1 - theta_star)
      for(indiv_loop in 1:length(z_pooled)){
        like_inc <- z_pooled[indiv_loop] * log(prob_star[indiv_loop]) + (1 - z_pooled[indiv_loop]) * 
          log(1 - prob_star[indiv_loop])
        log_pi_star <- log_pi_star + ifelse(is.na(like_inc),
                                            ifelse(z_pooled[indiv_loop] == 1, 0, -Inf),
                                            like_inc)
      }
    }
    
    prob_old <- 1 - (1 - theta_samples[iter - 1]) ^ n_pooled
    theta_old <- theta_samples[iter - 1]
    
    log_pi_old <- sum(z_pooled * log(prob_old) + (1 - z_pooled) * 
                        log(1 - prob_old) + (alpha - 1) * log(theta_old)  + (beta - 1) * log(1 - theta_old))
    
    if (is.na(log_pi_old)){
      log_pi_old <- (alpha - 1) * log(theta_old)  + (beta - 1) * log(1 - theta_old)
      for(indiv_loop in 1:length(z_pooled)){
        like_inc <- z_pooled[indiv_loop] * log(prob_old[indiv_loop]) + (1 - z_pooled[indiv_loop]) * 
          log(1 - prob_old[indiv_loop])
        log_pi_old <- log_pi_old + ifelse(is.na(like_inc),
                                          ifelse(z_pooled[indiv_loop] == 1, 0, -Inf),
                                          like_inc)
      }
    }
    
    # acceptance
    if ((log_pi_star - log_pi_old) > log(runif(1))){
      acceptance[iter] <- 1
      theta_samples[iter] <- theta_star
      theta_star <- -1
    } else {
      theta_samples[iter] <- theta_samples[iter - 1]
      theta_star <- -1
    }
  }
  return(list(interval = hdi(theta_samples[burn:num_mcmc]),theta = theta_samples[burn:num_mcmc], acceptance_rate = mean(acceptance)))
}

find_theta_yz <- function(y, z_mat, num_mcmc, step_size, alpha = 1.001, beta = 1.001){
  ## initialize MCMC
  num_pooled <- dim(z_mat)[2]
  num_indiv <- length(y)
  num_samples <- dim(z_mat)[1]
  burn <- num_mcmc / 10
  theta_samples <- rep(.5, num_mcmc)
  acceptance <- rep(0, num_mcmc)
  theta_star <- -1
  prob_star <- 1
  pooled_success <- sum(rowSums(z_mat) > 0)
  indiv_success <- sum(y) 
  
  ## run MCMC
  for (iter in 2:num_mcmc){
    # propose
    while(theta_star < 0 | theta_star > 1){
      theta_star <- theta_samples[iter - 1] + rnorm(1,0,step_size)
      prob_star <- 1 - (1 - theta_star) ^ num_pooled
    }
    
    if (log(1 - prob_star) == -Inf & (num_samples - pooled_success) > 0){ # prob one with zeros
      log_pi_star <- -Inf
    } else if (log(1 - prob_star) == -Inf & (num_samples - pooled_success) == 0){ # prob one with all ones
      log_pi_star <- indiv_success * log(theta_star) + (num_indiv - indiv_success) * log(1 - theta_star) + 
        (alpha - 1) * log(theta_star)  + (beta - 1) * log(1 - theta_star)
    }  else{
      log_pi_star <- indiv_success * log(theta_star) + (num_indiv - indiv_success) * log(1 - theta_star) + 
        pooled_success * log(prob_star) + (num_samples - pooled_success) * 
        log(1 - prob_star) + (alpha - 1) * log(theta_star)  + (beta - 1) * log(1 - theta_star)
    }
    
    prob_old <- 1 - (1 - theta_samples[iter - 1]) ^ num_pooled
    theta_old <- theta_samples[iter - 1]
    
    log_pi_old <- indiv_success * log(theta_old) + (num_indiv - indiv_success) * log(1 - theta_old) +
      pooled_success * log(prob_old) + (num_samples - pooled_success) * 
      log(1 - prob_old) + (alpha - 1) * log(theta_old)  + (beta - 1) * log(1 - theta_old)
    
    if (log(1 - prob_old) == -Inf & (num_samples - pooled_success) > 0){
      log_pi_old <- -Inf
    } else if (log(1 - prob_old) == -Inf & (num_samples - pooled_success) == 0){
      log_pi_old <- indiv_success * log(theta_old) + (num_indiv - indiv_success) * log(1 - theta_old) + 
        (alpha - 1) * log(theta_old)  + (beta - 1) * log(1 - theta_old)
    }  else{
      log_pi_old <- indiv_success * log(theta_old) + (num_indiv - indiv_success) * log(1 - theta_old) +
        pooled_success * log(prob_old) + (num_samples - pooled_success) * 
        log(1 - prob_old) + (alpha - 1) * log(theta_old)  + (beta - 1) * log(1 - theta_old)
    }
    
    # acceptance
    if ((log_pi_star - log_pi_old) > log(runif(1))){
      acceptance[iter] <- 1
      theta_samples[iter] <- theta_star
      theta_star <- -1
    } else {
      theta_samples[iter] <- theta_samples[iter - 1]
      theta_star <- -1
    }
  }
  return(list(interval = hdi(theta_samples[burn:num_mcmc]),theta = theta_samples[burn:num_mcmc], acceptance_rate = mean(acceptance)))
}
