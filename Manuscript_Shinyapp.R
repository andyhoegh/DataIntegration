library(shiny)
library(tibble)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(pscl)
library(HDInterval)

find_theta_y <- function(y, alpha = 1.1, beta = 1.1){
  return(list(interval = pscl::betaHPD(alpha + sum(y), beta + length(y) - sum(y)),
              mean = (sum(y) + alpha) / (alpha + length(y) + beta ) ))
}

find_theta_z <- function(z_mat, num_mcmc, step_size, alpha = 1, beta = 1){
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

find_theta_z2 <- function(z_pooled, n_pooled, num_mcmc, step_size, alpha = 1, beta = 1){
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

find_theta_yz <- function(y, z_mat, num_mcmc, step_size, alpha = 1, beta = 1){
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

ui <- navbarPage("Estimating Viral Prevalence with Data Integration for Adaptive Two-Phase Pooled Sampling",
                 tabPanel("Overview",
                          titlePanel("Shiny App Overview"),
                              br(),
                              mainPanel(h5("Efficacy of Data Integration"),
                              p("This tab shows how data integration can improve either pooled or individual tests."),
                              h5("Pooling Strategy"),
                              p("This tab shows how credible interval width and posterior mean differ by the number of pooled samples."),
                              h5("Adaptive Strategy"),
                              p("This tab compares  adaptive techniques for the second phase of a sampling framework.")
                              )),
########################################################################
                         tabPanel("Data Integration",
                          titlePanel("Data Integration"),
                          
                          # Sidebar with a slider input for number of bins
                          sidebarLayout(
                              sidebarPanel(
                                  sliderInput("num_pooled_samples",
                                              "Number of Pooled Samples:",
                                              min = 10,
                                              max = 200,
                                              value = 50,
                                              step = 10),
                                  sliderInput("num_indiv_samples",
                                              "Number of Individual Samples:",
                                              min = 10,
                                              max = 200,
                                              value = 50,
                                              step = 10),
                                  sliderInput("pool_size",
                                              "Number of Individual Samples in a Pool:",
                                              min = 1,
                                              max = 10,
                                              value = 3,
                                              step = 1),
                                  sliderInput("theta_resolution",
                                              "Prevalence Resolution,  \n Note: fine resolution will take longer to render",
                                              min = .01,
                                              max = .2,
                                              value = .1,
                                              step = .01),
                                  sliderInput("theta_max",
                                              "Maximum Prevalence: ",
                                              min = 0,
                                              max = 1,
                                              value = .9,
                                              step = .01),
                                  sliderInput("theta_min",
                                              "Minimum Prevalence: ",
                                              min = 0,
                                              max = 1,
                                              value = .1,
                                              step = .01),
                                  numericInput("num_replicates",
                                              "Number of Simulation Replicates: ",
                                              min = 1,
                                              max = 100,
                                              value = 2,
                                              step = 1)
                              ),
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                  h5("Data Integration"),
                                  p("These calculations show that data integration, using both the pooled and individual samples, results in better estimates than either alone."),
                                  plotOutput("di_plot")))),
########################################################################
tabPanel("Pooling Strategy",  titlePanel("Pooling Strategy"),
         
         # Sidebar with a slider input for number of bins
         sidebarLayout(
           sidebarPanel(
             selectizeInput(
               "pool_size2", "Pool Size",
               choices = 1:10,
               selected = c(1,5,10),
               multiple = TRUE
             ),
             sliderInput("n_pooled_samples",
                         "Number of Pooled Samples:",
                         min = 10,
                         max = 200,
                         value = 50,
                         step = 10),
             numericInput("pool_replicates",
                          "Number of Simulation Replicates: ",
                          min = 1,
                          max = 100,
                          value = 4,
                          step = 1),
             sliderInput("theta_resolution_pool",
                         "Prevalence Resolution,  \n Note: fine resolution will take longer to render",
                         min = .01,
                         max = .2,
                         value = .2,
                         step = .01),
             sliderInput("theta_max_pool",
                         "Maximum Prevalence: ",
                         min = 0,
                         max = 1,
                         value = .9,
                         step = .01),
             sliderInput("theta_min_pool",
                         "Minimum Prevalence: ",
                         min = 0,
                         max = 1,
                         value = .1,
                         step = .01)
           ),
           
           # Show a plot of the generated distribution
           mainPanel(
             h5("Pooling Strategy"),
             p("This tab shows how credible interval width and posterior mean differ by the number of pooled samples."),
             plotOutput("pool_plot")))),
########################################################################
tabPanel("Adaptive Strategy",  titlePanel("Adaptive Strategy"),
         # Sidebar with a slider input for number of bins
         sidebarLayout(
           sidebarPanel(
             numericInput("phase1_samples",
                         "Number of Pooled Phase 1 Samples",
                         min = 1,
                         max = 100,
                         value = 20,
                         step = 1),
             numericInput("phase1_pool_size",
                         "Phase 1 Pool Size",
                         min = 1,
                         max = 10,
                         value = 5,
                         step = 1),
             numericInput("phase2_samples",
                         "Number of Phase 2 Samples",
                         min = 1,
                         max = 200,
                         value = 80,
                         step = 1),
             numericInput("adaptive_replicates",
                          "Number of Simulation Replicates:",
                          min = 1,
                          max = 50,
                          value = 2,
                          step = 1)),
           
           # Show a plot of the generated distribution
           mainPanel(
             h5("Adaptive Strategy"),
             p("Allows users to evaluate a two-phase strategy. With the default values, the plot takes about 40 seconds to render."),
             plotOutput("adaptive_plot")))))
# Define server logic required to draw a histogram
server <- function(input, output) {
  ########################################################################
  output$di_plot <- renderPlot({
      
        num_pooled_samples <- input$num_pooled_samples
        num_indiv_samples <- input$num_indiv_samples
        pool_size <- input$pool_size
        theta_resolution <- input$theta_resolution
        theta_max <- input$theta_max
        theta_min <- input$theta_min
        num_replicates <- input$num_replicates
        
        
        theta_seq <- seq(theta_min, theta_max, by = theta_resolution)
        num_theta <- length(theta_seq)
        # Data observations: individual bats
        y_replicates <- array(0, dim = c(num_theta,num_replicates, 2))
        z_replicates <- array(0, dim = c(num_theta,num_replicates, 2))
        yz_replicates <- array(0, dim = c(num_theta,num_replicates, 2))
        y_mean <- y_width <- matrix(0, num_theta, num_replicates)
        z_mean <- z_width <- matrix(0, num_theta, num_replicates)
        yz_mean <- yz_width <- matrix(0, num_theta, num_replicates)
        
        for (j in 1:num_theta){
          theta <- theta_seq[j]
          for (i in 1:num_replicates){
            y <- rbinom(num_indiv_samples, 1, theta)
            z_mat <- matrix(rbinom(num_pooled_samples * pool_size, 1, theta), nrow = num_pooled_samples, ncol = pool_size)
            z <- as.numeric(rowSums(z_mat) > 0)
            
            y_samples <- find_theta_y(y)
            y_replicates[j,i,] <- y_samples$interval
            y_mean[j,i] <- y_samples$mean
            y_width[j,i] <- diff(y_samples$interval)
            
            z_samples <- find_theta_z(z_mat, 10000, .1)
            z_replicates[j,i,] <- z_samples$interval
            z_mean[j,i] <- mean(z_samples$theta)
            z_width[j,i] <- diff(z_samples$interval)
            
            yz_samples <- find_theta_yz(y,z_mat, 1000, .1)
            yz_replicates[j,i,] <- yz_samples$interval
            yz_mean[j,i] <- mean(yz_samples$theta)
            yz_width[j,i] <- diff(yz_samples$interval)
          }
        }
        
        f1a <- tibble(ci_width = c(y_mean, z_mean, yz_mean),
                      prevalence = rep(theta_seq, 3 * num_replicates),
                      Method = rep(c('individual samples','pooled only','data integration'), each = num_theta * num_replicates)) %>% 
          ggplot(aes(y=ci_width, x = prevalence, color = Method)) + geom_point(alpha = .2) + theme_bw() + 
          geom_smooth(method = "loess", formula = y ~ x) + ggtitle('Posterior Mean by Prevalence')  + ylab('Width') + xlab(expression(paste("prevalence (",theta,')'))) + facet_wrap(.~Method) + theme(legend.position='none')
        
        f1b <- tibble(ci_width = c(y_width, z_width, yz_width),
                      prevalence = rep(theta_seq, 3 * num_replicates),
                      Method = rep(c('individual samples','pooled only','data integration'), each = num_theta * num_replicates)) %>% 
          ggplot(aes(y=ci_width, x = prevalence, color = Method)) + geom_point(alpha = .2) + theme_bw() + 
          geom_smooth(method = "loess", formula = y ~ x) + ggtitle('Credible Interval Width by Prevalence') + ylab('Width') + xlab(expression(paste("prevalence (",theta,')'))) + theme(legend.position='bottom')
        
        grid.arrange(f1a,f1b)
     
    })
  ########################################################################
    output$pool_plot <- renderPlot({
      pool_replicates <- input$pool_replicates
      num_pooled_samples <- input$n_pooled_samples
      pool_size <- as.numeric(input$pool_size2)
      pool_length <- length(pool_size)
      max_samples <- max(pool_size) * num_pooled_samples
      theta_resolution <- input$theta_resolution_pool
      theta_max <- input$theta_max_pool
      theta_min <- input$theta_min_pool
      theta_seq <- seq(theta_min, theta_max, by = theta_resolution)
      num_theta <- length(theta_seq)
      
      ## Simulate and Fit Data
      z_mean <- z_width <- array(0, dim = c(pool_length,num_theta, pool_replicates))
      
      for (j in 1:num_theta){
        theta <- theta_seq[j]
        for (i in 1:pool_replicates){
          y <- rbinom(max_samples, 1, theta)
          for (k in 1:pool_length){
            z_mat <- matrix(y[1:(pool_size[k]*num_pooled_samples)], nrow = num_pooled_samples, ncol = pool_size[k])
            z_samples <- find_theta_z(z_mat, 10000, .1)
            z_mean[k,j,i] <- mean(z_samples$theta)
            z_width[k,j,i] <- diff(z_samples$interval)
          }
        }
      }
       
      f2a <- tibble(ci_mean = c(z_mean),
                     replicate = rep(1:pool_replicates, each = pool_length * num_theta),
                     num_pooled = factor(rep(pool_size, pool_replicates * num_theta)),
                     prevalence = rep(rep(theta_seq, each = pool_length), pool_replicates)) %>%
        ggplot(aes(y=ci_mean, x = prevalence, color = num_pooled)) +
        geom_point(alpha = .2) + geom_smooth(method = "loess", formula = y ~ x) +  theme_bw() + ggtitle('Posterior Mean by Prevalence')  +
        ylab(expression(paste("Posterior mean for prevalence (",theta,')'))) + xlab(expression(paste("prevalence (",theta,')'))) + theme(legend.position='none')
      
      f2b <- tibble(ci_mean = c(z_width),
                   replicate = rep(1:pool_replicates, each = pool_length * num_theta),
                   `pool size` = factor(rep(pool_size, pool_replicates * num_theta)),
                   prevalence = rep(rep(theta_seq, each = pool_length), pool_replicates)) %>%
        ggplot(aes(y=ci_mean, x = prevalence, color = `pool size`)) +
        geom_point(alpha = .2) + geom_smooth(method = "loess", formula = y ~ x) +  theme_bw() + ggtitle('Credible Interval Width by Prevalence')  +
        ylab("Credible Interval Width") + xlab(expression(paste("prevalence (",theta,')'))) + theme(legend.position='bottom')
      
      grid.arrange(f2a,f2b)
    })
    ########################################################################
    output$adaptive_plot <- renderPlot({
      num_replicates <- input$adaptive_replicates
      theta_seq <- seq(0, 1, by = .05)
      theta_length <- length(theta_seq)
      num_theta <- length(theta_seq)
      phase1_samples <- input$phase1_samples
      pool_size <- input$phase1_pool_size
      total_samples <- phase1_samples * pool_size
      z <- array(0, dim = c(num_replicates, theta_length, total_samples / pool_size, pool_size))
      # Take Phase 1 Sample
      z_mean <- z_width <- array(0, dim = c(num_replicates, theta_length))
      for (i in 1:num_replicates){
        for (j in 1:theta_length){
          z[i,j,,] <- matrix(rbinom(total_samples,1,theta_seq[j]), nrow = total_samples / pool_size, ncol = pool_size)
          z_samples <- find_theta_z(z[i,j,,], 10000, .1)
          z_mean[i,j] <- mean(z_samples$theta)
          z_width[i,j] <- diff(z_samples$interval)
        }
      }
      ################################################
      ### Phase 2. Take an additional 80 samples.
      ################################################
      # Strategies
      ## a. individual samples
      ## b. pooled size 2
      ## c. pooled size 3
      ## d. pooled size 4
      ## e. pooled size 5
      
      # a. Add individual samples
      phase2_samples <- input$phase2_samples
      
      y_width <- y_mean <- array(0, dim = c(num_replicates, theta_length))
      
      for (i in 1:num_replicates){
        for (j in 1:theta_length){
          y <- rbinom(phase2_samples,1,theta_seq[j])
          y_samples <- find_theta_yz(y,z[i,j,,], 10000, .1)
          y_mean[i,j] <- mean(y_samples$theta)
          y_width[i,j] <- diff(y_samples$interval)
        }
      }
 
      # # b. Add pooled samples of size 2
      # pool_size <- 2
      # z2_width <- z2_mean <- array(0, dim = c(num_replicates, theta_length))
      # 
      # for (i in 1:num_replicates){
      #   for (j in 1:theta_length){
      #     z2 <- matrix(rbinom(phase2_samples * pool_size,1,theta_seq[j]), nrow = phase2_samples, ncol = pool_size)
      #     z_pooled <- as.numeric(c(rowMeans(z[i,j,,]) > 0, rowMeans(z2) > 0))
      #     n_pooled <- c(rep(5, nrow(z[i,j,,])), rep(pool_size, nrow(z2)))
      #     z2_samples <- find_theta_z2(z_pooled, n_pooled, 10000, .1)
      #     z2_mean[i,j] <- mean(z2_samples$theta)
      #     z2_width[i,j] <- diff(z2_samples$interval)
      #   }
      # }
 
      ## c. Add pooled samples of size 3
      # pool_size <- 3
      # z3_width <- z3_mean <- array(0, dim = c(num_replicates, theta_length))
      # 
      # for (i in 1:num_replicates){
      #   for (j in 1:theta_length){
      #     z3 <- matrix(rbinom(phase2_samples * pool_size,1,theta_seq[j]), nrow = phase2_samples, ncol = pool_size)
      #     z_pooled <- as.numeric(c(rowMeans(z[i,j,,]) > 0, rowMeans(z3) > 0))
      #     n_pooled <- c(rep(5, nrow(z[i,j,,])), rep(pool_size, nrow(z3)))
      #     z3_samples <- find_theta_z2(z_pooled, n_pooled, 10000, .1)
      #     z3_mean[i,j] <- mean(z3_samples$theta)
      #     z3_width[i,j] <- diff(z3_samples$interval)
      #   }
      # }
      # 
      # # d. Add pooled samples of size 4
      # pool_size <- 4
      # z4_width <- z4_mean <- array(0, dim = c(num_replicates, theta_length))
      # 
      # for (i in 1:num_replicates){
      #   for (j in 1:theta_length){
      #     z4 <- matrix(rbinom(phase2_samples * pool_size,1,theta_seq[j]), nrow = phase2_samples, ncol = pool_size)
      #     z_pooled <- as.numeric(c(rowMeans(z[i,j,,]) > 0, rowMeans(z4) > 0))
      #     n_pooled <- c(rep(5, nrow(z[i,j,,])), rep(pool_size, nrow(z4)))
      #     z4_samples <- find_theta_z2(z_pooled, n_pooled, 10000, .1)
      #     z4_mean[i,j] <- mean(z4_samples$theta)
      #     z4_width[i,j] <- diff(z4_samples$interval)
      #   }
      # }

      # e. Add pooled samples of size 5
      pool_size <- 5
      z5_width <- z5_mean <- array(0, dim = c(num_replicates, theta_length))

      for (i in 1:num_replicates){
        for (j in 1:theta_length){
          z5 <- matrix(rbinom(phase2_samples * pool_size,1,theta_seq[j]), nrow = phase2_samples, ncol = pool_size)
          z_pooled <- as.numeric(c(rowMeans(z[i,j,,]) > 0, rowMeans(z5) > 0))
          n_pooled <- c(rep(5, nrow(z[i,j,,])), rep(pool_size, nrow(z5)))
          z5_samples <- find_theta_z2(z_pooled, n_pooled, 10000, .1)
          z5_mean[i,j] <- mean(z5_samples$theta)
          z5_width[i,j] <- diff(z5_samples$interval)
        }
      }
      
      f3a <- tibble( vals = c(c(z_mean),c(y_mean), c(z5_mean)),
                     theta = rep(rep(theta_seq, each = num_replicates),3),
                     method = rep(c('initial pool','phase 2 with pools of size 1', 'phase 2 with pools of size 5'), each = num_replicates * theta_length)) %>%
        ggplot(aes(y = vals, x = theta, color = method)) + geom_point(alpha = .1) + geom_smooth(method = 'gam') +
        theme_bw() + ylim(0,1) + ggtitle('Posterior Mean by Strategy') + theme(legend.position = "none")
      
      f3b <- tibble( vals = c(c(z_width),c(y_width), c(z5_width)),
                     theta = rep(rep(theta_seq, each = num_replicates),3),
                     method = rep(c('initial pool','phase 2 with pools of size 1', 'phase 2 with pools of size 5'), each = num_replicates * theta_length)) %>%
        ggplot(aes(y = vals, x = theta, color = method)) + geom_point(alpha = .1) + geom_smooth(method = 'gam') +
        theme_bw() + ylim(0,1) + ggtitle('Credible Interval Width by Strategy') + theme(legend.position='bottom')
      
      #   f3a <- tibble( vals = c(c(z_mean),c(y_mean), c(z2_mean), c(z3_mean), c(z4_mean), c(z5_mean)),
      #                  theta = rep(rep(theta_seq, each = num_replicates),6),
      #                  method = rep(c('initial pool','phase 2 with 1', 'phase2 with 2','phase2 with 3', 'phase2 with 4', 'phase2 with 5'), each = num_replicates * theta_length)) %>%
      #   ggplot(aes(y = vals, x = theta, color = method)) + geom_point(alpha = .1) + geom_smooth(method = 'gam') +
      #   theme_bw() + ylim(0,1) + ggtitle('Posterior Mean by Strategy') + theme(legend.position = "none")
      # 
      # f3b <- tibble( vals = c(c(z_width),c(y_width), c(z2_width), c(z3_width), c(z4_width), c(z5_width)),
      #                theta = rep(rep(theta_seq, each = num_replicates),6),
      #                method = rep(c('initial pool','phase 2 with 1', 'phase2 with 2','phase2 with 3', 'phase2 with 4', 'phase2 with 5'), each = num_replicates * theta_length)) %>%
      #   ggplot(aes(y = vals, x = theta, color = method)) + geom_point(alpha = .1) + geom_smooth(method = 'gam') +
      #   theme_bw() + ylim(0,1) + ggtitle('Credible Interval Width by Strategy') + theme(legend.position='bottom')
      # 
      grid.arrange(f3a, f3b)
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
