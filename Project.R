# N
# days 
# y = recovery rate
# c = potentially infectious contacts picked randomly
# a = transmission probability
rm(list=ls())
set.seed(565)

sir_model <- function(N, days, gam, c, alph, I0){
  
  # (1) Initialize population matrix
  # 1 = susceptible | 2 = infected | 3 = recovered
  
  # Starting off, everyone is susceptible
  population_matrix <- matrix(1L, ncol = 1, nrow = N)
  colnames(population_matrix) <- "state"
  
  # Infected individuals are sampled
  infected_idx <- sample.int(N, I0)
  population_matrix[infected_idx, "state"] <- 2L
  
  # Store Day 0
  sumStats <- NULL
  sumStats <- rbind(sumStats, c(day=0, 
                                S = sum(population_matrix[ , "state"] == 1L), 
                                I = sum(population_matrix[ , "state"] == 2L), 
                                R = sum(population_matrix[ , "state"] == 3L)))
  
  colnames(sumStats) <- c("Day", "S", "I", "R")
  
  # Daily Events
  
  for(t in 1:days){
    
  
  #----------------------------------------------
  
  # (2) Recovery Phase
    infected_idx <- which(population_matrix[, "state"] == 2L)
    
    if(length(infected_idx) > 0){
      recovery_events <- rbinom(n = length(infected_idx), size = 1, prob = gam)
      
      recovered_positions <- which(recovery_events == 1)
      
      recovered_idx <- infected_idx[recovered_positions]
      
      population_matrix[recovered_idx, "state"] <- 3L
    }
  
  
  #----------------------------------------------
  
  # (3) Infection Phase
    infected_idx <- which(population_matrix[, "state"] == 2L)
    
    if(length(infected_idx) > 0){
      number_of_contacts <- length(infected_idx)*c
      
      sampled_contacts <- sample.int(N, number_of_contacts, TRUE)
      
      sampled_contacts_states <- population_matrix[sampled_contacts, "state"]
      
      susceptible_contacts_idx <- sampled_contacts[sampled_contacts_states == 1L]
      
      
      if(length(susceptible_contacts_idx) > 0){
        infection_events <- rbinom(length(susceptible_contacts_idx), 1, alph)
        
        infected_positions <- which(infection_events == 1)
        
        newly_infected_idx <- susceptible_contacts_idx[infected_positions]
        
        population_matrix[newly_infected_idx, "state"] <- 2L
      }
      
    }
  
  #----------------------------------------------
    
  # (4) Store Summary Data
    sumStats <- rbind(sumStats,
                      c(day = t,
                        S = sum(population_matrix[, "state"] == 1L),
                        I = sum(population_matrix[, "state"] == 2L),
                        R = sum(population_matrix[, "state"] == 3L)))
  }

  return(as.data.frame(sumStats))
  
}


summary <- sir_model(100, 100, 0.3, 4, 0.2, 10)


#c) Time until Peak

#d) Duration of pandemic

# Run x replicates, store results in DF, compute diagnostics (a-d), plot mean trajectory
explore_sir_model <- function(N, days, gam, c, alph, I0, repl){
  susceptible_summary <- matrix(NA, nrow=days+1, ncol=repl)
  infected_summary <- matrix(NA, nrow=days+1, ncol=repl)
  recovered_summary <- matrix(NA, nrow=days+1, ncol=repl)
  
  day_vec <- 0:days
  rownames(susceptible_summary) <- day_vec
  rownames(infected_summary) <- day_vec
  rownames(recovered_summary) <- day_vec
  
  colnames(susceptible_summary) <- paste0("rep_", 1:repl)
  colnames(infected_summary) <- paste0("rep_", 1:repl)
  colnames(recovered_summary) <- paste0("rep_", 1:repl)
  
  for(r in 1:repl){
    summary <- sir_model(N = N, day = days, gam = gam, c = c, alph = alph, I0 = I0)
    #store summary data in 
    susceptible_summary[, r] <- summary$S
    infected_summary[, r] <- summary$I
    recovered_summary[, r] <- summary$R
  }
  # mean infected proportion across replicates
  mean_final_susceptible <- mean(susceptible_summary[nrow(susceptible_summary), ])
  mean_infected_proportion <- (N - mean_final_susceptible) / N
  
  # mean peak infected across replicates
  peak_infected_per_rep <- apply(infected_summary, 2, max)
  mean_infected_peak <- mean(peak_infected_per_rep)

  old_par <- par(mar = c(6, 4, 4, 2) + 0.1)
  
  plot(day_vec, infected_summary[, 1],
       type = "n",
       xlab = "Day",
       ylab = "Number of individuals",
       main = "SIR trajectories (all replicates)",
       ylim = range(c(susceptible_summary,
                      infected_summary,
                      recovered_summary),
                    na.rm = TRUE)
  )
  
  col_S <- rgb(0,   0,   1,   0.3)  
  col_I <- rgb(1,   0,   0,   0.3)  
  col_R <- rgb(0, 0.7,   0,   0.3)  
  
  for (r in 1:repl) {
    lines(day_vec, susceptible_summary[, r], col = col_S)
    lines(day_vec, infected_summary[, r],    col = col_I)
    lines(day_vec, recovered_summary[, r],   col = col_R)
  }
  
  legend("topright",
         legend = c("Susceptible", "Infected", "Recovered"),
         col    = c("blue", "red", "darkgreen"),
         lwd    = 2,
         cex    = 0.7,
         bty    = "n")
  
  txt1 <- sprintf("Mean infected proportion: %.3f", mean_infected_proportion)
  txt2 <- sprintf("Mean peak infected number: %.1f", mean_infected_peak)
  mtext(txt1, side = 1, line = 4, adj = 0)
  mtext(txt2, side = 1, line = 5, adj = 0)
  
  par(old_par) 
  
  susceptible_df <- data.frame(day = day_vec, susceptible_summary, row.names = NULL)
  infected_df    <- data.frame(day = day_vec, infected_summary, row.names = NULL)
  recovered_df   <- data.frame(day = day_vec, recovered_summary, row.names = NULL)
  
  metrics <- list(
    mean_infected_proportion = mean_infected_proportion,
    mean_infected_peak = mean_infected_peak
  )
  
  return(list(
    susceptible = susceptible_df,
    infected = infected_df,
    recovered = recovered_df,
    metrics = metrics
  ))
}

explore_sir_model(100, 100, 0.3, 4, 0.2, 10, 10)

  


