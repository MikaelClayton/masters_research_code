# libraries
file_to_source <- 'C:/Users/Mikae/OneDrive/Documents/Bronwyn_research_code/masters_research_code/Simulation_code/data_generation.R'
normalized_path <- normalizePath(file_to_source, mustWork = TRUE)
source(normalized_path)
file_to_source <- 'C:/Users/Mikae/OneDrive/Documents/Bronwyn_research_code/masters_research_code/Simulation_code/all_models_use.R'
normalized_path <- normalizePath(file_to_source, mustWork = TRUE)
source(normalized_path)

## parameters
n <- 1000 # data points 
m <- 4 # number of nodes
mu1=c(0,0,0,0)
sigma1=matrix(c(3,0,0,0,0,3,0,0,0,0,3,0,0,0,0,3),ncol=4,nrow=4, byrow=TRUE)
mu2=c(7,7,7,7)
sigma2=sigma1
mu3=c(3,3,3,3)
sigma3=sigma1
probs <- c(0.3, 0.3, 0.4)
means <- list(mu1, mu2, mu3)
covs <- list(sigma1, sigma2, sigma3)
params <- list('means' = means,
               'covs' = covs,
               'mixing_probs' = probs)

## Function to calculate global statistics for parallel approach
calculate_global_statistics_multivariate <- function(node_local_statistics) {
  k <- length(node_local_statistics[[1]]$suff_stat_two)  # Number of components
  num_nodes <- length(node_local_statistics)  # Number of nodes
  
  global_suff_stat_one <- NULL
  global_suff_stat_two <- list()
  global_suff_stat_three <- list()
  
  # Initialize sums for suff_stat_two and suff_stat_three
  for (i in 1:k) {
    global_suff_stat_two[[i]] <- 0
    global_suff_stat_three[[i]] <- 0
  }
  
  # Sum local statistics across all nodes for each component
  for (node in 1:num_nodes) {
    local_stats <- node_local_statistics[[node]]
    
    # First sufficient statistic (just summing across nodes)
    if (is.null(global_suff_stat_one)) {
      global_suff_stat_one <- local_stats$suff_stat_one
    } else {
      global_suff_stat_one <- global_suff_stat_one + local_stats$suff_stat_one
    }
    
    # Second and third sufficient statistics (sum across nodes for each component)
    for (i in 1:k) {
      global_suff_stat_two[[i]] <- global_suff_stat_two[[i]] + local_stats$suff_stat_two[[i]]
      global_suff_stat_three[[i]] <- global_suff_stat_three[[i]] + local_stats$suff_stat_three[[i]]
    }
  }
  
  return(list(
    "global_suff_stat_one" = global_suff_stat_one,
    "global_suff_stat_two" = global_suff_stat_two,
    "global_suff_stat_three" = global_suff_stat_three
  ))
}

## Generate data for nodes
data_chunks <- lapply(1:m, function(x) generate_multivariate_mixture_gaussian_data(means, 
                                                           covs,
                                                           probs,
                                                           n/m))
## Check code for parameter estimation
local_node_statistics <- lapply(1:m, function(i) calculate_node_local_sufficient_statistics_multivariate(
  data_chunks[[i]], params
))
global_sufficient_statistics <- calculate_global_statistics_multivariate(local_node_statistics)
params <- estimate_parameters_with_sufficient_statistics_multivariate(global_sufficient_statistics)
params
