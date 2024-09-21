source('C:/Users/mccal/OneDrive/Desktop/Masters/Research/masters_research/data_generation.R')
source('C:/Users/mccal/OneDrive/Desktop/Masters/Research/masters_research/all_models_use.R')

calculate_data_belongings <- function(params, data){
    #' Calculate component responsibilities for Univariate Gaussian Mixture Model
    #' 
    #' @description
    #' The function calculates the probability that a data point belongs to
    #' to specified component.
    #' 
    #' This is based on the data and the current parameter estimates.
    #'
    #'@param params list. A list of the different components params. Mean, sd and pi
    #'@param data matrix or vector. Matrix if multivariate and vector if univariate
    #'
    #'@return data_belongings matrix. A column for each component and row for data points.
    
    ### Extracting parameters and checking if they are they same.
    means <- params$means
    stds <- params$stds
    probs <- params$mixing_probs
    if (length(means) != length(stds) || length(means) != length(probs)){
        stop("Number of means, standard deviations, and weights must be equal.")
    }
    
    ### Get the data sample and number of components
    n_components <- length(means) 
    n <- length(data)
    
    ### Calculating the belongings 
    gamma_num <- matrix(nrow = n, ncol = n_components)
    for(i in 1:n_components){
        gamma_num[,i] <- probs[i] * dnorm(data, mean = means[i], sd = stds[i])
    }
    gamma_den <- rowSums(gamma_num)
    gamma <- gamma_num / gamma_den
    
    ### Check if data belongings sum to 1 (have grace for rounding errors).
    if (any(abs(rowSums(gamma) - 1) > .Machine$double.eps * n)){
        stop("Calculated belongings for each data point do not sum to 1.")
    }
    
    ### Check if any data belonging is negative
    if (any(gamma < 0)){
        stop("Calculated belongings contain negative values.")
    }
    return(gamma)
}


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




