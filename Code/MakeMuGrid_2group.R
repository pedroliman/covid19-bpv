require(tidyverse)
tFinal <- c(50)
R0step = .05
R0 <- seq(1.5,2,by=R0step)
R0 <- 3
i0 <- c(0.0001)
sFinalStep <- 0.01
sFinal <- c(seq(.0001,.001,by=.0001),seq(.001,.01,by=.001),seq(0.01,1,by=sFinalStep))
#delta <- seq(0,.1,by=sFinalStep)
delta <- c(seq(0,.001,by=.0001),seq(.001,.01,by=.001),seq(0.01,1,by=sFinalStep))
exper.design.mu <- expand.grid(tFinal,R0,i0,sFinal,delta)

param.list <- c("tFinal","R0","i0","sFinal","delta")

names(exper.design.mu) <- param.list

exper.design.mu$case.id.mu <- 1:dim(exper.design.mu)[1]
exper.design.mu <- subset(exper.design.mu,sFinal-delta>0)

# Create a copy of the experimental design for neighbor calculations
data.steps <- exper.design.mu
data.steps$R0 <- data.steps$R0 / R0step
data.steps$sFinal <- data.steps$sFinal / sFinalStep
data.steps$delta <- data.steps$delta / sFinalStep

# Function to find neighbors one step away along each parameter dimension
find_one_step_neighbors <- function(y) {
  current_row <- y[1:5]  # Exclude the 'case.id.mu' column for the current row
  
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
  if (length(out) < 6) {
    pad_length <- 6 - length(out)
    out <- c(out, rep(NA, pad_length))
  }
  return(out)
}


# Apply the function to each row of data.steps
data.steps.with.neighbors <- apply(data.steps, 1, find_one_step_neighbors)

# Create column names for the neighbor columns
exper.design.mu$solved <- FALSE
neighbor_col_names <- paste0("neighbor_", 1:6)



# Add the neighbor columns to the original data frame
exper.design.mu[, neighbor_col_names] <- do.call('rbind',data.steps.with.neighbors)[,1:6]



# Save the data frame with the added neighbor information
write.csv(exper.design.mu, "INTERMEDIATE/exper.design.mu0.2groups.csv")

sols <- read.csv("INTERMEDIATE/mu.solutions.csv")

aux <- subset(exper.design.mu,delta==0)
aux$solved <- NULL
aux2 <- inner_join(aux,sols)
aux2$s1 <- aux2$s/2
aux2$s2 <- aux2$s/2
aux2$i1 <- aux2$i/2
aux2$i2 <- aux2$i/2
aux2$mu1 <- aux2$mu_s
aux2$mu2 <- aux2$mu_s
aux2$nu1 <- aux2$mu_i
aux2$nu2 <- aux2$mu_i
aux2$rho11 <- aux2$R_tau
aux2$rho12 <- aux2$R_tau
aux2$rho21 <- aux2$R_tau
aux2$rho22 <- aux2$R_tau

aux2$case.id <- NULL
aux2$s <- NULL
aux2$i <- NULL
aux2$mu_s <- NULL
aux2$mu_i <- NULL
aux2$R_tau <- NULL
aux2$Rt <- NULL

neighbors <- exper.design.mu %>%
  pivot_longer(cols = starts_with("neighbor_"),
               names_to = "name",
               values_to = "neighbor_num") 
neighbors$neighbor_solved <- FALSE

write.csv(aux2,"INTERMEDIATE/mu.solutions.2groups.csv")

write.csv(neighbors,"INTERMEDIATE/mu.neighbors0.2groups.csv")
