


# ----- Genetic Optimizer -----



# A genetic optimizer for the sensitivity analysis described in Daniel McCaffrey's 
# Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies


create_population <- function(pop_size, test_data, G){
  population <- matrix(runif(nrow(test_data), 1/G, G), 
                       nrow = pop_size, 
                       ncol = nrow(test_data))
  return(population)
}

crossover <- function(X, Y, crossover_rate, G){
  p <- matrix(runif(length(X), 0, 1), nrow = nrow(X), ncol = ncol(X))
  doCross <- matrix((p < crossover_rate), nrow = nrow(X), ncol = ncol(X))
  gamma <- matrix(runif(length(X), 0, 1), nrow = nrow(X), ncol = ncol(X))
  crosses <- gamma * (X - Y) + X
  Z <- X
  Z <- ifelse(doCross, crosses, Z)
  Z <- ifelse(Z > G, G, Z)
  Z <- ifelse(Z < (1/G), (1/G), Z)
  return(Z)
}

mutate <- function(offspring, mutation_rate = 0.01, G){
  p <- matrix(runif(length(offspring), 0, 1), nrow = nrow(offspring),
              ncol = ncol(offspring))
  
  m <- matrix(runif(length(offspring), 1/G, G), nrow = nrow(offspring),
              ncol = ncol(offspring))
  
  offspring <- ifelse(p < mutation_rate, m, offspring)
  return(offspring)
}

evaluate_chromosomes <- function(population, v, y, prob_vaccine){
  v_0 <- matrix(ifelse(v == 0, 1, 0), dim(population)[1], ncol = length(v), byrow = TRUE)
  y <- matrix(y, dim(population)[1], ncol = length(v), byrow = TRUE)
  prob_vaccine <- matrix(prob_vaccine, dim(population)[1], ncol = length(v), byrow = TRUE)
  p_star <- (prob_vaccine * population) / ((prob_vaccine * (population - 1)) + 1) 
  w <- 1 / p_star
  num <- rowSums((v_0 * y * w))
  denom <- rowSums(v_0)
  return(num/denom)
}

# Run the optimizer
genetic_optimizer <- function(data, prob_vaccine, G) {

  # Optimization parameters
  pop_size <- 1000
  crossover_rate <- 0.5
  mutation_rate <- 0.01
  iters <- 500
  epochs <- 50
  fitness_record <- rep(NA, epochs)
  population <- create_population(
    pop_size, data, G)
  
  # Initialize the record
  fitness_record <- rep(NA, epochs)
  population <- create_population(pop_size, data, G)
  
  # Optimize
  for (i in 1:epochs) {
    
    # Evaluate fitness
    pop_fitness <- evaluate_chromosomes(
      population, 
      v = data$vaccine, 
      y = data$outcome,
      prob_vaccine = prob_vaccine
    )
    pop_fitness <- unlist(pop_fitness)
    
    # Selection
    fitness_median <- median(pop_fitness)
    population <- population[pop_fitness >= fitness_median,]
    
    # Crossover
    num_offspring <- pop_size-dim(population)[1]
    if(num_offspring == 1){
      num_offspring_eff <- 2
    } else {
      num_offspring_eff <- num_offspring
    }
    x <- population[sample(1:dim(population)[1], num_offspring_eff, replace = TRUE),]
    y <- population[sample(1:dim(population)[1], num_offspring_eff, replace = TRUE),]
    offspring <- crossover(x, y, crossover_rate, G)
    
    # Mutate
    offspring <- mutate(offspring, mutation_rate, G)
    
    # Update population
    if(num_offspring == 1){
      population <- rbind(population, offspring[1,])
    } else {
      population <- rbind(population, offspring)
    }
    
    # Update fitness
    fitness_record[i] <- mean(pop_fitness)
    fittest_ind <- which.max(pop_fitness)
    
  }
  
  # Return
  return(population[fittest_ind,])
  
}
