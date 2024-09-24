#####################File info#####################
# Code that stores logic for Parallel expectation maximisation model.

# The functions stored is:
#     1. calculate_global_statistics 
#     2. iem_distributed_model

# libraries
file_to_source <- 'C:/Users/Mikae/OneDrive/Documents/Bronwyn_research_code/masters_research_code/Simulation_code/data_generation.R'
normalized_path <- normalizePath(file_to_source, mustWork = TRUE)
source(normalized_path)
file_to_source <- 'C:/Users/Mikae/OneDrive/Documents/Bronwyn_research_code/masters_research_code/Simulation_code/all_models_use.R'
normalized_path <- normalizePath(file_to_source, mustWork = TRUE)
source(normalized_path)
library(parallel)

calculate_global_statistics <- function(...){
    local_statistics <- list(...)
    
    global_suff_stat_one <- NULL
    global_suff_stat_two <- NULL
    global_suff_stat_three <- NULL
    
    for(local_stat in local_statistics){
        if(is.null(global_suff_stat_one)){
            global_suff_stat_one <- local_stat$suff_stat_one
            global_suff_stat_two <- local_stat$suff_stat_two
            global_suff_stat_three <- local_stat$suff_stat_three
        } else{
            global_suff_stat_one <- global_suff_stat_one + local_stat$suff_stat_one
            global_suff_stat_two <- global_suff_stat_two + local_stat$suff_stat_two
            global_suff_stat_three <- global_suff_stat_three + local_stat$suff_stat_three
        }
    }
    
    if(sum(global_suff_stat_one) <= 0){
        stop("The sum of the first global sufficient statistics must be positi`ve")
    }
    
    if(any(global_suff_stat_three < 0)){
        stop("Third global sufficient statistics is negative")
    }
    return(list("global_suff_stat_one" = global_suff_stat_one,
                "global_suff_stat_two" = global_suff_stat_two,
                "global_suff_stat_three" = global_suff_stat_three))
}

GMM_using_parallel <- function(data, params_init,
                               tol = 0.05,
                               max_iter = 1000,
                               worker_nodes = length(data)){
    params_old <- params_init
    
    # Determine the number of components (k)
    k <- length(params_old$means)
    
    # Create column names dynamically
    column_names <- c("iteration")
    column_names <- c(column_names,
                      paste0("component_", 1:k, "_mean"),
                      paste0("component_", 1:k, "_stds"),
                      paste0("component_", 1:k, "_mixing_probs"))
    
    iteration <- 0
    
    # Initialize results matrix with the correct number of columns
    results <- matrix(ncol = length(column_names), nrow = 0)
    colnames(results) <- column_names
    
    # Set up cluster
    cl <- makeCluster(worker_nodes, type = "PSOCK")
    clusterExport(cl, list("calculate_node_local_sufficient_statistics", 
                           "calculate_data_belongings"
                           ))
    param_diff <- 1000
    while(param_diff > tol){
        iteration <- iteration + 1
        if (iteration > max_iter){
            break
        }
        
        # Export current parameters to the workers
        clusterExport(cl, varlist = "params_old", envir = environment())
        
        # E-Step (done on workers)
        local_sufficient_statistics <- parLapply(cl, data, function(d){
            calculate_node_local_sufficient_statistics(data = d, 
                                                       params = params_old)
        })
        
        # Get global sufficient statistics (manager node)
        global_sufficient_statistics <- do.call(calculate_global_statistics, 
                                                local_sufficient_statistics)
        
        # Update the parameters
        params_new <- estimate_parameters_with_sufficient_statistics(
            global_sufficient_statistics = global_sufficient_statistics)
        # log_likelihood_new <- calculate_log_likelihood_with_full_data(
        #     data = full_data, params = params_new)

        param_diff <- max(abs(unlist(params_new) - unlist(params_old)))
        new_row <- c(iteration, 
                     params_new$means, params_new$stds, params_new$mixing_probs)
        
        # Append the new row to results
        results <- rbind(results, new_row)
        params_old <- params_new
    }
    stopCluster(cl)
    
    return(list("params" = params_old,
                "performance" = results))
}

######### MULTIVARIATE ########
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

GMM_using_parallel_multivariate <- function(data, params_init,
                               tol = 0.05,
                               max_iter = 1000,
                               worker_nodes = length(data)){
  params_old <- params_init
  
  # Determine the number of components (k)
  k <- length(params_old$means)
  results <- list()
  iteration <- 0
  
  # Set up cluster
  cl <- makeCluster(worker_nodes, type = "PSOCK")
  clusterExport(cl, list("calculate_node_local_sufficient_statistics_multivariate", 
                         "calculate_data_belongings_multivariate",
                         "dmvnorm"
  ))
  param_diff <- 1000
  while(param_diff > tol){
    iteration <- iteration + 1
    if (iteration > max_iter){
      break
    }
    
    # Export current parameters to the workers
    clusterExport(cl, varlist = "params_old", envir = environment())
    
    # E-Step (done on workers)
    local_sufficient_statistics <- parLapply(cl, data, function(d){
      calculate_node_local_sufficient_statistics_multivariate(data = d, 
                                                 params = params_old)
    })
    
    # Get global sufficient statistics (manager node)
    global_sufficient_statistics <- calculate_global_statistics_multivariate(local_sufficient_statistics)
    
    # Update the parameters
    params_new <- estimate_parameters_with_sufficient_statistics_multivariate(
      global_sufficient_statistics = global_sufficient_statistics)
    
    param_diff <- max(abs(unlist(params_new) - unlist(params_old)))
    
    # Append the new row to results
    results[[iteration]] <- params_new
    params_old <- params_new
  }
  stopCluster(cl)
  
  return(list("params" = params_new,
              "performance" = results))
}
params_new <- GMM_using_parallel_multivariate(data_chunks, 
                                              params)

