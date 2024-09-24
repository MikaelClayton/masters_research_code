# Centralised multivariate gaussian mixture modeling
## parameters
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
               'probs' = probs)

## libaries
library(ggplot2)
library(MASS)
library(reshape2)
library(mvtnorm)

generate_multivariate_mixture_gaussian_data <- function(means, covariances, probs, n) {
  # Ensure lengths match
  if (length(means) != length(covariances) || length(means) != length(probs)) {
    stop("Number of means, covariance matrices, and component probabilities must be equal.")
  }
  
  data <- matrix(nrow = n, ncol = length(means[[1]]))
  components <- numeric(n)  # Initialize components vector to store component indices
  
  for (i in 1:n) {
    component <- which.max(rmultinom(1, 1, probs))  # Get the component index
    components[i] <- component                      # Store component number
    data[i, ] <- rmvnorm(1, mean = means[[component]], sigma = covariances[[component]])
  }
  
  return(data)  # Return both data and components
}

plot_mixture_histograms_with_density <- function(data, means, covariances, probs) {
  num_dims <- ncol(data)
  num_components <- length(means)
  n <- nrow(data)
  
  # Calculate component assignment based on probabilities
  components <- numeric(n)
  
  for (i in 1:n) {
    components[i] <- which.max(rmultinom(1, 1, probs))  # Assign component based on probabilities
  }
  
  components <- as.factor(components)  # Convert components to factor for discrete coloring
  
  # Create a data frame for plotting
  df <- as.data.frame(data)
  df$components <- components  # Add the calculated components
  
  # Create a long-form data for easier plotting
  df_long <- melt(df, id.vars = "components", variable.name = "dimension")
  
  # Initialize plot
  plot <- ggplot(df_long, aes(x = value, fill = components)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6, position = "identity") +
    facet_wrap(~dimension, scales = "free", ncol = 2) +  # Create a grid for each dimension
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    theme_minimal() +
    labs(title = "Histogram of Mixture Components with Density Overlays", x = "Value", y = "Density")
  
  # Overlay density curves for each component
  for (comp in 1:num_components) {
    for (dim in 1:num_dims) {
      x_vals <- seq(min(data[, dim]), max(data[, dim]), length.out = 100)
      density_vals <- dnorm(x_vals, mean = means[[comp]][dim], sd = sqrt(covariances[[comp]][dim, dim]))
      df_density <- data.frame(x = x_vals, density = density_vals, dimension = paste0("V", dim))
      
      plot <- plot + 
        geom_line(data = df_density, aes(x = x, y = density), color = c("#F8766D", "#00BA38", "#619CFF")[comp], linewidth = 1)
    }
  }
  
  # Display the plot
  print(plot)
}

calculate_node_local_sufficient_statistics_multivariate <- function(data, params){
  ### Initiate needed values
  sufficient_stat_two <- list()
  sufficient_stat_three <- list()
  means <- params$means
  
  ### Calculating the data belongings
  belongings <- calculate_data_belongings(params = params, data = data)
  k <- ncol(belongings)
  
  ### Calculating the local sufficient statistics for the node
  #### Number 1: sum of data belongings
  sufficient_stat_one <- colSums(belongings)
  
  for(i in 1:k){
    #### Number 2: sum of product (data and data belongings)
    sufficient_stat_two[[i]] <- colSums(belongings[, i] * data)
    #### Number 3: sum of product (data squared and data belongings)
    sufficient_stat_three[[i]] <- t(belongings[,i]*(data - means[[i]]))%*%(data - means[[i]])
  }
  
  if(sum(sufficient_stat_one) != nrow(data)){
    stop("The sum of the first sufficient statistic does not equal node sample size")
  }
  
  #if(any(diag(sufficient_stat_three) < 0)){
  #  stop("Third sufficient statistic is negative")
  #}
  
  return(list("suff_stat_one" = sufficient_stat_one, 
              "suff_stat_two" = sufficient_stat_two,
              "suff_stat_three" = sufficient_stat_three)
  )
}

calculate_global_statistics_multivariate <- function(node_local_statistics, gamma) {
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

estimate_parameters_with_sufficient_statistics_multivariate <- function(global_sufficient_statistics){
  
  # Initialize vectors for means, standard deviations, and probabilities
  means <- list()
  covs <- list()
  probs <- c()
  k <- length(global_sufficient_statistics$global_suff_stat_one)
  
  for (i in 1:k) {
    means[[i]] <- global_sufficient_statistics$global_suff_stat_two[[i]] / global_sufficient_statistics$global_suff_stat_one[i]
    covs[[i]] <- global_sufficient_statistics$global_suff_stat_three[[i]] / global_sufficient_statistics$global_suff_stat_one[i]
    probs <- cbind(probs,
                   global_sufficient_statistics$global_suff_stat_one[i] / sum(global_sufficient_statistics$global_suff_stat_one))
  }
  
  return(list("mixing_probs" = probs,
              "means" = means,
              "covs" = covs))
}

# Example call
data_chunks <- lapply(1:4, function(x) generate_multivariate_mixture_gaussian_data(
  means, covs, probs, 100/4
))
node_local_statistics <- lapply(1:4, function(x) calculate_node_local_sufficient_statistics_multivariate(
  data_chunks[[x]], params
))
global_sufficient_statistics <- calculate_global_statistics_multivariate(node_local_statistics, gamma)
params_new_dist <- estimate_parameters_with_sufficient_statistics_multivariate(global_sufficient_statistics)
params_new_dist
calculate_log_likelihood_with_full_data_multivariate(data, params_new_dist)
calculate_log_likelihood_with_full_data_multivariate(data, centralised_estimates)

data <- NULL
for(i in 1:length(data_chunks)){
  data <- rbind(data, data_chunks[[i]])
}
gamma <- calculate_data_belongings(params, data)
centralised_estimates <- estimate_parameters_with_full_data(data, gamma)

