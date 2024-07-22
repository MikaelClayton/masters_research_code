library(promises)
library(future)
library(later)

plan(multisession) # Use multiple R sessions for parallel processing

# Initialize an empty results data frame
results <- data.frame(
    Task = integer(),
    simulated_sleep = integer(),
    start_time = as.POSIXct(character()),
    end_time = as.POSIXct(character()),
    duration = numeric(),
    random_number = numeric(),
    stringsAsFactors = FALSE
)

generate_random_number <- function(task_name = Sys.getpid(), wait_time = 10) {
    future({
        start_time <- Sys.time()
        Sys.sleep(wait_time)
        random_number <- rnorm(1, mean = 0, sd = 1)
        end_time <- Sys.time()
        duration <- round(difftime(end_time, start_time, units = "secs"))
        result <- data.frame(
            Task = task_name,
            simulated_sleep = wait_time,
            start_time = start_time,
            end_time = end_time,
            duration = duration,
            random_number = random_number,
            stringsAsFactors = FALSE
        )
        
        return(result)
    }, seed = TRUE)
}

diff_sleep_times <- c(20, 30, 40)

promises <- lapply(diff_sleep_times, function(wait_time) generate_random_number(wait_time = wait_time))

update_results <- function(promise){
    promise %...>% (function(result){
        # Safely update results data frame
        results <<- rbind(results, result)
    }) %...!% (function(e) {
        cat("Error:", e$message, "\n")
    })
}

# Apply update_results to each promise
lapply(promises, function(promise) update_results(promise))

# Function to check and print the table periodically
check_table <- function(interval = 5) {
    # Print the current state of the results table
    print(results)
    
    # Check if all promises are resolved
    if (nrow(results) < length(diff_sleep_times)) {
        # Schedule the next check
        later(function() check_table(interval), interval)
    } else {
        print("Complete")
    }
}

# Start the periodic checking
check_table(interval = 5)
