# libraries
file_to_source <- 'C:/Users/Mikae/OneDrive/Documents/Bronwyn_research_code/masters_research_code/Simulation_code/data_generation.R'
normalized_path <- normalizePath(file_to_source, mustWork = TRUE)
source(normalized_path)

n <- 100 # data points 
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

## data_belongings
calculate_data_belongings_multivariate <- function(params, data) {
  
  # Extracting parameters
  means <- params$means
  covs <- params$covs
  probs <- params$mixing_probs
  
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
    # The covariance is using a newer mean value. 
    covs[[i]] <- (t(data_belongings[,i]*(data - means[[i]]))%*%(data - means[[i]]))/n_k[i]
  }
  print(paste('ss_1: ', n_k))
  print(colSums(data_belongings[, i] * data))
  print(t(data_belongings[,i]*(data - means[[i]]))%*%(data - means[[i]]))
  return(list('mixing_probs' = probs,
              'means' = means,
              'covs' = covs))
}

calculate_node_local_sufficient_statistics_multivariate <- function(data, params){
  ### Initiate needed values
  sufficient_stat_two <- list()
  sufficient_stat_three <- list()
  means <- params$means
  
  ### Calculating the data belongings
  belongings <- calculate_data_belongings_multivariate(params = params, data = data)
  k <- ncol(belongings)
  
  ### Calculating the local sufficient statistics for the node
  #### Number 1: sum of data belongings
  sufficient_stat_one <- colSums(belongings)
  
  for(i in 1:k){
    #### Number 2: sum of product (data and data belongings)
    sufficient_stat_two[[i]] <- colSums(belongings[, i] * data)
    #### Number 3: sum of product (data squared and data belongings)
    sufficient_stat_three[[i]] <-t(belongings[,i]*(data - means[[i]]))%*%(data - means[[i]])
    #sufficient_stat_three[[i]] <-t(belongings[,i]*(data - means[[i]]))%*%(data - means[[i]])
  }
  
  
  if(sum(sufficient_stat_one) != nrow(data)){
    stop("The sum of the first sufficient statistic does not equal node sample size")
  }
  
  return(list("suff_stat_one" = sufficient_stat_one, 
              "suff_stat_two" = sufficient_stat_two,
              "suff_stat_three" = sufficient_stat_three)
  )
}

calculate_global_statistics_multivariate <- function(node_local_statistics) {
  k <- length(node_local_statistics[[1]]$suff_stat_two)  # Number of components
  
  # Summing the first sufficient statistic across nodes
  global_suff_stat_one <- Reduce(`+`, lapply(node_local_statistics, `[[`, "suff_stat_one"))
  
  # Summing the second and third sufficient statistics across nodes
  global_suff_stat_two <- lapply(1:k, function(i) {
    Reduce(`+`, lapply(node_local_statistics, function(node) node$suff_stat_two[[i]]))
  })
  
  global_suff_stat_three <- lapply(1:k, function(i) {
    Reduce(`+`, lapply(node_local_statistics, function(node) node$suff_stat_three[[i]]))
  })
  
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
    print(global_sufficient_statistics$global_suff_stat_two[[i]])
    print(global_sufficient_statistics$global_suff_stat_one[i])
    covs[[i]] <- global_sufficient_statistics$global_suff_stat_three[[i]] / global_sufficient_statistics$global_suff_stat_one[i]
    print(global_sufficient_statistics$global_suff_stat_three[[i]])
    probs <- cbind(probs,
                   global_sufficient_statistics$global_suff_stat_one[i] / sum(global_sufficient_statistics$global_suff_stat_one))
    print(sum(global_sufficient_statistics$global_suff_stat_one))
  }
  
  return(list("mixing_probs" = probs,
              "means" = means,
              "covs" = covs))
}

### Check
data_chunks <- lapply(1:m, function(x) generate_multivariate_mixture_gaussian_data(
  means, covs, probs, n/m
))
data <- NULL
for(i in 1:m){
  data <- rbind(data,data_chunks[[i]])
}
data_chunks
data

data_belongings_chunks <- lapply(1:m, function(x) calculate_data_belongings_multivariate(
  params, data_chunks[[x]]
))

data_belongings_chunks

data_belongings <- calculate_data_belongings_multivariate(params, data)
data_belongings

estimate_parameters_with_full_data_multivariate(data, data_belongings)

local_node_statistics <- lapply(1:m, function(x) calculate_node_local_sufficient_statistics_multivariate(
  data_chunks[[x]], params
))
local_node_statistics[[1]]


calculate_global_statistics_multivariate(local_node_statistics)


estimate_parameters_with_full_data_multivariate(data_chunks[[1]], data_belongings_chunks[[1]])
calculate_node_local_sufficient_statistics_multivariate(data_chunks[[1]], params)
