#####################File info#####################
# Code that stores the logic for data generation 

# The functions stored is:
#     1. generate_univariate_mixture_gaussian_data
#     2. generate_multivariate_mixture_gaussian_data
#     3. plot_univariate_gaussian_data
#     4. TODO: plot_bivariate_gaussian_data

#libraries
library(plotly)
library(reshape2)
library(ggplot2)
library(mvtnorm)
generate_univariate_mixture_gaussian_data <- function(means, stds, probs, n){
    
    ## Ensure different parameter lengths match
    if(length(means) != length(stds) || length(means) != length(probs)){
        stop("Number of means, standard deviations and components probabilities 
             must be equal.")
    }
    data <- matrix(nrow = n)
    for(i in 1:n){
        component <- rmultinom(1, 1, probs) + 1
        data[i] <- rnorm(1, mean = means[component], sd = stds[component])
    }
    return(data)
}

plot_univariate_gaussian_data <- function(data, n_components, params, 
                                          overlay = T, freq = F){
    
    ## Extracting parameters
    means <- params$means
    stds <- params$stds
    probs <- params$mixing_probs
    ## Histogram
    hist(data, breaks = 50, col = "lightblue", xlab = "Value", 
         main = "Histogram of two-component univariate Gaussian mixture model",
         freq = freq)
    ## Overlayed line plot
    if (overlay) {
        colors <- rainbow(n_components)
        for (i in 1:n_components){
            curve(dnorm(x, mean = means[i], sd = stds[i])*probs[i], add = TRUE, 
                  col = colors[i], lwd = 2)
        }
    }
}

######## Multivariate #########
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
    data[i, ] <- round(rmvnorm(1, mean = means[[component]], sigma = covs[[component]]), 6)
  }
  
  return(data)  # Return both data and components
}

plot_mixture_histograms_with_density <- function(data, means, covariances, probs) {
  num_dims <- ncol(data)  # Number of dimensions (features)
  num_components <- length(means)  # Number of mixture components
  n <- nrow(data)  # Number of data points
  
  # Assign components based on probabilities
  components <- numeric(n)
  for (i in 1:n) {
    components[i] <- which.max(rmultinom(1, 1, probs))  # Assign component based on probabilities
  }
  
  components <- as.factor(components)  # Convert components to a factor for coloring
  
  # Create a data frame for plotting
  df <- as.data.frame(data)
  df$components <- components  # Add components to the data frame
  
  # Reshape the data for easier plotting
  df_long <- melt(df, id.vars = "components", variable.name = "dimension")
  
  # Initialize the plot
  plot <- ggplot(df_long, aes(x = value, fill = components)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6, position = "identity") +
    facet_wrap(~dimension, scales = "free", ncol = 2) +  # Facet by dimension
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    theme_minimal() +
    labs(title = "Histogram of Mixture Components with Density Overlays", x = "Value", y = "Density")
  
  # Overlay density curves for each component in each dimension
  for (comp in 1:num_components) {
    for (dim in 1:num_dims) {
      # Generate x values for density calculation
      x_vals <- seq(min(data[, dim]), max(data[, dim]), length.out = 100)
      
      # Compute multivariate normal density for each component and dimension
      density_vals <- dnorm(x_vals, mean = means[[comp]][dim], sd = sqrt(covariances[[comp]][dim, dim]))
      
      # Create a data frame for the density curve
      df_density <- data.frame(x = x_vals, density = density_vals, dimension = paste0("V", dim))
      
      # Add density curves to the plot
      plot <- plot + 
        geom_line(data = df_density, aes(x = x, y = density), color = c("#F8766D", "#00BA38", "#619CFF")[comp], linewidth = 1)
    }
  }
  
  # Display the plot
  print(plot)
}

                               


