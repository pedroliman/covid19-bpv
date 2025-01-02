tFinal <- c(50)
R0step = .05
R0 <- seq(3,3,by=R0step)
i0 <- c(0.0001)
sFinalStep <- 0.01
sFinal <- seq(0.01,1,by=sFinalStep)

exper.design.mu <- expand.grid(tFinal,R0,i0,sFinal)

param.list <- c("tFinal","R0","i0","sFinal")

names(exper.design.mu) <- param.list

exper.design.mu$case.id.mu <- 1:dim(exper.design.mu)[1]

# Create a copy of the experimental design for neighbor calculations
data.steps <- exper.design.mu
data.steps$R0 <- data.steps$R0 / R0step
data.steps$sFinal <- data.steps$sFinal / sFinalStep

# Function to find neighbors one step away along each parameter dimension
find_one_step_neighbors <- function(y) {
  current_row <- y[1:4]  # Exclude the 'case.id.mu' column for the current row
  
  # Find neighbors one step away along each parameter dimension
  neighbors <- apply(data.steps[, param.list], 1, function(x) {
    diff = sum(abs(current_row - x))
    is_one_step_away <- abs(diff - 1) < 0.0000001
    return(is_one_step_away)
  })
  
  # Find the indices of the neighbors
  neighbor_indices <- which(neighbors)
  
  # Get case IDs of neighbors and pad with NA if less than 4 neighbors found
  out <- data.steps[neighbor_indices, "case.id.mu"]
  if (length(out) < 4) {
    pad_length <- 4 - length(out)
    out <- c(out, rep(NA, pad_length))
  }
  return(out)
}

# Apply the function to each row of data.steps
data.steps.with.neighbors <- apply(data.steps, 1, find_one_step_neighbors)

# Create column names for the neighbor columns
exper.design.mu$solved <- FALSE
neighbor_col_names <- paste0("neighbor_", 1:4)



# Add the neighbor columns to the original data frame
exper.design.mu[, neighbor_col_names] <- t(data.steps.with.neighbors)



# Save the data frame with the added neighbor information
write.csv(exper.design.mu, "INTERMEDIATE/exper.design.mu.csv")