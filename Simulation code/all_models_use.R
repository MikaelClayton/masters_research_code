#####################File info#####################
# Code that stores all functions used by all models.
# The functions stored is:
#     1. calculate_data_belongings
#     2. calculate_node_local_sufficient_statistics
#     3. estimate_parameters_with_sufficient_statistics
#     4. TODO: estimate_parameters_with_full_data
#     5. calculate_log_likelihood_with_full_data
#     6. TODO: calculate_log_likelihood_with_sufficient_statistics

# Libraries
library(docstring)

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

calculate_node_local_sufficient_statistics <- function(data, params){
    #' Calculate sufficient statistics for Univariate Gaussian Mixture Modelling
    #' 
    #' @description
    #' The function calculates the local sufficient statistics for a univariate
    #' Gaussian mixture model. This is:
    #' 1. Sum of component belonginings
    #' 2. Sum of the product of component belongings and data points.
    #' 3. Sum of the product of component belongins and data points squared.
    #'
    #' This is based on the data and the current parameter estimates.
    #'
    #'@param params list. A list of the different components params. Mean, sd and pi
    #'@param data matrix or vector. Matrix if multivariate and vector if univariate
    #'
    #'@return sufficient_statistics list. Each element is a set of one of the the ss.
    
    ### Calculating the data belongings
    belongings <- calculate_data_belongings(params = params, data = data)
    
    ### Calculating the local sufficient statistics for the node
    #### Number 1: sum of data belongings
    sufficient_stat_one <- colSums(belongings)
    
    #### Number 2: sum of product (data and data belongings)
    sufficient_stat_two <- colSums(belongings * matrix(data, nrow = length(data), ncol = ncol(belongings), byrow = FALSE))
    
    #### Number 3: sum of product (data squared and data belongings)
    sufficient_stat_three <- colSums(belongings * matrix(data^2, nrow = length(data), ncol = ncol(belongings), byrow = FALSE))
    
    if(sum(sufficient_stat_one) != length(data)){
        stop("The sum of the first sufficient statistic does not equal node sample size")
    }
    
    if(any(sufficient_stat_three < 0)){
        stop("Third sufficient statistic is negative")
    }
    return(list("suff_stat_one" = sufficient_stat_one, 
                "suff_stat_two" = sufficient_stat_two,
                "suff_stat_three" = sufficient_stat_three)
    )
}

estimate_parameters_with_sufficient_statistics <- function(global_sufficient_statistics){
    #' Reestimate parameters for univariate Gaussian Mixture Model.
    #' 
    #' @description
    #' The function re-estimates univariate GMM parameters. This is:
    #' 1. The components' mixing probabilities.
    #' 2. The components' means.
    #' 3. The components' standard deviations.
    #'
    #' This is based on the sufficient statistics for the full dataset.
    #'
    #'@param params list. A list of the different components params. Mean, sd and pi
    #'@param data matrix or vector. Matrix if multivariate and vector if univariate
    #'
    #'@return new_parameter_estimates list. List with mixing_probs, means, sds.
    
    # Initialize vectors for means, standard deviations, and probabilities
    means <- numeric(length(global_sufficient_statistics$global_suff_stat_one))
    stds <- numeric(length(global_sufficient_statistics$global_suff_stat_one))
    probs <- numeric(length(global_sufficient_statistics$global_suff_stat_one))
    
    for (k in 1:length(global_sufficient_statistics$global_suff_stat_one)) {
        mean <- global_sufficient_statistics$global_suff_stat_two[k] / global_sufficient_statistics$global_suff_stat_one[k]
        var <- (global_sufficient_statistics$global_suff_stat_three[k] / global_sufficient_statistics$global_suff_stat_one[k]) - (mean^2)
        std <- sqrt(var)
        prob <- global_sufficient_statistics$global_suff_stat_one[k] / sum(global_sufficient_statistics$global_suff_stat_one)
        
        # Assign calculated values to the respective vectors
        means[k] <- mean
        stds[k] <- std
        probs[k] <- prob
    }
    
    return(list("mixing_probs" = probs,
                "means" = means,
                "stds" = stds))
}


calculate_log_likelihood_with_full_data <- function(data, params){
    #' Calculate log-likelihood for univariate Gaussian Mixture Model.
    #' 
    #' @description
    #' The function calculates the log-likelihood for a univariate 
    #' Gaussian mixture model.
    #'
    #' This is based on the data and the current parameter estimates.
    #'
    #'@param params list. A list of the different components params. Mean, sd and pi
    #'@param data matrix or vector. Matrix if multivariate and vector if univariate
    #'
    #'@return log-likelihood float.
    
    # Initialize vectors for means, standard deviations, and probabilities
    mu <- params$means
    sigma <- params$stds
    weight <- params$mixing_probs
    k <- length(mu)
    n <- length(data)
    likelihood <- numeric(n)
    for (i in 1:n) {
        point_likelihood <- 0
        for (j in 1:k) {
            component_density <- weight[j] * dnorm(data[i], mean = mu[j], sd = sigma[j])
            point_likelihood <- point_likelihood + component_density
        }
        likelihood[i] <- point_likelihood
    }
    return(sum(log(likelihood)))
}

