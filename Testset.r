


# ----- Test Functions for the Sensitivity Analysis -----



# Set the file paths
here::i_am("sensitivity_analysis/Testset.R")

# Load the functions for the causal testbed and generate some data
# Generating a prime number just so that we can check against the recycle rule
source(here::here("sensitivity_analysis", "CausalTestbed.R"))

# Load the genetic operators
source(here::here("sensitivity_analysis", "GeneticOperators.R"))



# Test the population function
# We'll do this inside an enclosing environment so that tests don't
# interfere with each other
test_create_population <- function(){
  pop_size <- 37
  test_data <- obs_dgp(n = 1229)
  G = 4
  test_pop <- create_population(pop_size = pop_size, test_data = test_data, G = G)
  num_rows <- dim(test_pop)[1]
  num_cols <- dim(test_pop)[2]
  if(!(num_rows == pop_size)){
    stop("Number of rows in population is not correct.")
  }
  if(!(num_cols == nrow(test_data))){
    stop("Number of columns in population is not correct.")
  }
  invisible()
}
test_create_population()


# Test evaluate chromosomes
# population, v, y, prob_vaccine
test_evaluate_chromosome <- function(){
  test_pop <- create_population(pop_size = pop_size, test_data = test_data, G = G)
  test_data <- obs_dgp(n = 1229)
  test_v = test_data$vaccine
  test_y = test_data$outcome
  test_prob_vaccine <- ifelse(v == 1, runif(n = 1, min = 0, max = 1), runif(n = 1, min = 0, max = 1))
  test_chromosomes <- evaluate_chromosomes(test_pop, test_data, test_v, test_y)
}


# Test the crossover function
test_crossover <- function(){
  
  # Parameters
  num_offspring <- 37
  num_cols <- 51
  cross_rate <- 0.19
  test_G <- 4
  
  # Dummy crossover
  parent_x <- matrix(
    rep(seq(1 / test_G, test_G, length.out = num_offspring), each = num_cols),
    num_offspring,
    num_cols,
    byrow = TRUE
  )
  parent_y <- matrix(
    rep(-seq(1 / test_G, test_G, length.out = num_offspring), each = num_cols),
    num_offspring,
    num_cols,
    byrow = TRUE
  )
  test_offspring <- crossover(parent_x, parent_y, cross_rate, G = test_G)
  
  # Check dimensions
  num_rows <- num_offspring
  num_cols <- num_cols
  if(!(num_rows == dim(test_offspring)[1])){
    stop("Number of rows in offspring is not correct.")
  }
  if(!(num_cols == dim(test_offspring)[2])){
    stop("Number of columns in offspring is not correct.")
  }
  
  # Check if values are in range (between G & 1/G)
  if(any(test_offspring < (1 / test_G))){
    stop("There is at least 1 a-value below 1/G in the offspring.")
  }
  if(any(test_offspring > test_G)){
    stop("There is at least 1 a-value above G in the offspring.")
  }
  
  # Check if offspring is a duplicate of X
  est_prop_cross <- sapply(1:dim(test_offspring)[1], function(i){
    return(1 - sum(test_offspring[i, ] == parent_x[i, ]) / dim(parent_x)[2])
  })
  if(abs(mean(est_prop_cross) - cross_rate) > 0.1){
    stop("Offspring did not crossover at the correct rate.")
  }
  
  invisible()
}
test_crossover()


# Test the mutation function
test_mutation <- function(){
  
  # Parameters
  num_offspring <- 37
  num_cols <- 51
  mut_rate <- 0.07
  test_G <- 3
  
  # Dummy crossover
  offspring_x <- matrix(
    rep(seq(1 / test_G, test_G, length.out = num_offspring), each = num_cols),
    num_offspring,
    num_cols,
    byrow = TRUE
  )
  test_offspring <- mutate(offspring = offspring_x, mut_rate, G = test_G)
  
  # Check dimensions
  num_rows <- num_offspring
  num_cols <- dim(offspring_x)[2]
  if(!(num_rows == dim(test_offspring)[1])){
    stop("Number of rows in offspring is not correct.")
  }
  if(!(num_cols == dim(test_offspring)[2])){
    stop("Number of columns in offspring is not correct.")
  }
  
  # Check if values are in range (between G & 1/G)
  if(any(test_offspring < (1 / test_G))){
    stop("There is at least 1 a-value below 1/G in the offspring.")
  }
  if(any(test_offspring > test_G)){
    stop("There is at least 1 a-value above G in the offspring.")
  }
  
  # Check if offspring is a duplicate of X
  est_prop_mut <- sapply(1:dim(test_offspring)[1], function(i){
    return(1 - sum(test_offspring[i, ] == offspring_x[i, ]) / dim(offspring_x)[2])
  })
  if(abs(mean(est_prop_mut) - mut_rate) > 0.1){
    stop("Offspring did not mutate at the correct rate.")
  }
  
  invisible()
  
}
test_mutation()
