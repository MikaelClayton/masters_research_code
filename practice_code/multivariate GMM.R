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


#################Code for EM#################
calculate_data_belongings <- function(params, data) {
  
  # Extracting parameters
  means <- params$means
  covs <- params$covs
  probs <- params$probs
  
  # Ensure the number of means, covariances, and probabilities are the same
  if (length(means) != length(covs) || length(means) != length(probs)) {
    stop("Number of means, covariance matrices, and component probabilities must be equal.")
  }
  
  # Number of components and data points
  n_components <- length(means)
  n <- nrow(data)
  
  # Initialize a matrix to store the numerator for gamma
  gamma_num <- matrix(nrow = n, ncol = n_components)
  
  # Calculate the numerator for each component
  for (i in 1:n_components) {
    gamma_num[, i] <- probs[i] * dmvnorm(data, mean = means[[i]], sigma = covs[[i]])  
  }
  
  # Calculate the denominator (sum across all components for each data point)
  gamma_den <- rowSums(gamma_num)
  
  # Calculate gamma (posterior probabilities for each component)
  gamma <- gamma_num / gamma_den
  
  # Check if gamma rows sum to 1 (allowing small rounding errors)
  if (any(abs(rowSums(gamma) - 1) > .Machine$double.eps * n)) {
    stop("Calculated belongings for each data point do not sum to 1.")
  }
  
  # Check for negative values in gamma
  if (any(gamma < 0)) {
    stop("Calculated belongings contain negative values.")
  }
  
  return(gamma)
}

estimate_parameters_with_full_data <- function(data, data_belongings){
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
  return(list('probs' = probs,
              'means' = means,
              'covs' = covs))
}

# Generate data and count components
data <- generate_multivariate_mixture_gaussian_data(means, covs, probs, 100)
gamma <- calculate_data_belongings(params, data)
estimate_parameters_with_full_data(data, gamma)
plot_mixture_histograms_with_density(data, means, covs, probs)
params
