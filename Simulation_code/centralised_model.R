file_to_source <- 'C:/Users/Mikae/OneDrive/Documents/Bronwyn_research_code/masters_research_code/Simulation_code/data_generation.R'
normalized_path <- normalizePath(file_to_source, mustWork = TRUE)
source(normalized_path)
file_to_source <- 'C:/Users/Mikae/OneDrive/Documents/Bronwyn_research_code/masters_research_code/Simulation_code/all_models_use.R'
normalized_path <- normalizePath(file_to_source, mustWork = TRUE)
source(normalized_path)

estimate_parameters_with_full_data <- function(data, data_belongings){
    n_k <- colSums(data_belongings)
    n <- length(data)
    k <- ncol(data_belongings)
    
    probs <- c()
    means <- c()
    stds <- c()
    for(i in 1:k){
        probs <- cbind(probs, n_k[i]/n)
        means <- cbind(means, sum(data_belongings[,i]*data)/n_k[i])
        stds <- cbind(stds, 
                      sqrt((t(data_belongings[,i]*(data - means[i]))%*%(data - means[i]))/n_k[i]))
    }
    return(list('mixing_probs' = probs,
                 'means' = means,
                 'stds' = stds))
}

centralised_model <- function(data,
                              initial_params,
                              tol = 0.05,
                              max_iter = 1000){
    params_old <- initial_params
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
    params_diff <- 10000
    while(params_diff > tol){
        iteration <- iteration + 1
        if (iteration > max_iter){
            break
        }
        gamma <- calculate_data_belongings(params_old, data)
        params_new <- estimate_parameters_with_full_data(data,
                                                         gamma)
        # log_likelihood_new <- calculate_log_likelihood_with_full_data(
        #     data = data, params = params_new)
        
        params_diff <- max(abs(unlist(params_new) - unlist(params_old)))
        new_row <- c(iteration, 
                     params_new$means, params_new$stds, params_new$mixing_probs)
        
        # Append the new row to results
        results <- rbind(results, new_row)
        params_old <- params_new
    }
    return(list("params" = params_old,
                "performance" = results))
}


######## MULTIVARIATE ##########
estimate_parameters_with_full_data_multivariate <- function(data, data_belongings){
  n_k <- colSums(data_belongings)
  n <- nrow(data)
  k <- ncol(data_belongings)
  
  probs <- c()
  means <- list()
  covs <- list()
  for(i in 1:k){
    probs <- cbind(probs, n_k[i]/n)
    means[[i]] <- colSums(data_belongings[, i] * data) / n_k[i]
    covs[[i]] <- (t(data_belongings[,i]*(data - means[[i]]))%*%(data - means[[i]]))/n_k[i]
  }
  return(list('mixing_probs' = probs,
              'means' = means,
              'covs' = covs))
}

centralised_model_multivariate <- function(data,
                              initial_params,
                              tol = 0.05,
                              max_iter = 1000){
  params_old <- initial_params
  # Determine the number of components (k)
  k <- length(params_old$means)
  results <- list()
  iteration <- 0
  
  params_diff <- 10000
  while(params_diff > tol){
    iteration <- iteration + 1
    if (iteration > max_iter){
      break
    }
    gamma <- calculate_data_belongings_multivariate(params_old, data)
    params_new <- estimate_parameters_with_full_data_multivariate(data,
                                                     gamma)
    
    params_diff <- max(abs(unlist(params_new) - unlist(params_old)))
    
    # Append the new row to results
    results[[iteration]] <- params_new
    params_old <- params_new
  }
  return(list("params" = params_old,
              "performance" = results))
}
params_centralised <- centralised_model_multivariate(data,
                                                     params)



