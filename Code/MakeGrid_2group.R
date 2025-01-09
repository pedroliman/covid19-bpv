require(tidyverse)
tFinal <- c(50)
R0step = .05
R0 <- seq(1.5,2,by=R0step)
R0 <- 3
i0 <- c(0.0001)
cStep <- 0.01
c1.vec <- seq(0,1,by=cStep)
c2.vec <- seq(0,1,by=cStep)
num.neighbors <- 10
exper.design <- expand.grid(tFinal,R0,i0,c1.vec,c2.vec)

param.list <- c("tFinal","R0","i0","c1","c2")

names(exper.design) <- param.list

exper.design$case.id <- 1:dim(exper.design)[1]
exper.design <- subset(exper.design,c1>=c2)

# Create a copy of the experimental design for neighbor calculations
data.steps <- exper.design
data.steps$R0 <- data.steps$R0 / R0step
data.steps$c1 <- data.steps$c1 / cStep
data.steps$c2 <- data.steps$c2 / cStep

# Function to find neighbors one step away along each parameter dimension
find_one_step_neighbors <- function(y) {
  current_row <- y[1:5]  # Exclude the 'case.id.mu' column for the current row
  
  # Find neighbors one step away along each parameter dimension
  diffs <- apply(data.steps[, param.list], 1, function(x) {
    diff = sum(abs(current_row - x))

    return(diff)
  })
  smallest_diffs <- sort(diffs)[1:num.neighbors]
  # Find the indices of the neighbors
  neighbor_indices <- which(diffs %in% smallest_diffs)[1:num.neighbors]
  
  # Get case IDs of neighbors and pad with NA if less than 4 neighbors found
  out <- data.steps[neighbor_indices, "case.id.mu"]
  if (length(out) < num.neighbors) {
    pad_length <- num.neighbors - length(out)
    out <- c(out, rep(NA, pad_length))
  }
  return(out)
}


# Apply the function to each row of data.steps
data.steps.with.neighbors <- apply(data.steps, 1, find_one_step_neighbors)

# Create column names for the neighbor columns
exper.design$solved <- FALSE
neighbor_col_names <- paste0("neighbor_", 1:num.neighbors)



# Add the neighbor columns to the original data frame
#exper.design.mu[, neighbor_col_names] <- do.call('rbind',data.steps.with.neighbors)[,1:6]

exper.design[, neighbor_col_names] <- t(data.steps.with.neighbors)


# Save the data frame with the added neighbor information
write.csv(exper.design, "INTERMEDIATE/exper.design_2groups_lambda.csv")

sols <- read.csv("INTERMEDIATE/mu.solutions.2groups.final.csv")

exper.mu <- subset(sols,tau==tFinal)
exper.mu <- exper.mu[!duplicated(exper.mu$case.id.mu),]
exper.mu$c1 <- round(1/exper.mu$mu1,2)
exper.mu$c2 <- round(1/exper.mu$mu2,2)
exper.mu$b2 <- exper.mu$c1/exper.mu$c2
exper.mu$Cost.full <- exper.mu$Cost+(1-exper.mu$s1)+exper.mu$b2*(1-exper.mu$s2)
exper.mu <- subset(exper.mu,c1 %in% c1.vec & c2 %in% c2.vec)
exper.mu$c1.c2 <- paste(exper.mu$c1,exper.mu$c2,sep="_")
exper.mu.split <- split(exper.mu,exper.mu$c1.c2)
find.min <- function(x){
  j <- which(x$Cost.full==min(x$Cost.full))
  return(x[j,])
}
exper.mu.split2 <- lapply(exper.mu.split,find.min)
exper.mu2 <- do.call('rbind',exper.mu.split2)


# Find matches and add a column to df_A
exper.design$param.id <- apply(exper.design[, c("tFinal", "R0", "i0", "c1", "c2")], 1, paste, collapse = "-")
exper.mu2$param.id <- apply(exper.mu2[, c("tFinal", "R0", "i0", "c1", "c2")], 1, paste, collapse = "-")
row.names(exper.mu2) <- exper.mu2$param.id
row.names(exper.design) <- exper.design$param.id
exper.design[,"case.id.mu"] <- exper.mu2[exper.design$param.id,"case.id.mu"]
exper.design[!(exper.design$param.id %in% exper.mu2$param.id),"case.id.mu"] <- NA
aux <- subset(sols,case.id.mu %in% exper.design$case.id.mu)
sub <- subset(exper.design,!is.na(exper.design$case.id.mu))
row.names(sub) <- sub$case.id.mu
aux$c <- sub[aux$case.id.mu,"c1"]
aux$case.id <- sub[aux$case.id.mu,"case.id"]
aux$b2 <- sub[aux$case.id.mu,"c1"]/sub[aux$case.id.mu,"c2"]
aux$lambda1 <- aux$c*aux$mu1-1
aux$lambda2 <- aux$c*aux$mu2-aux$b2
aux$eta1 <- aux$c*aux$nu1
aux$eta2 <- aux$c*aux$nu2


# Function to find neighbors one step away along each parameter dimension
find_one_step_neighbors <- function(y) {
  current_row <- y[1:5]  # Exclude the 'case.id.mu' column for the current row
  
  # Find neighbors one step away along each parameter dimension
  diffs <- apply(data.steps[, param.list], 1, function(x) {
    diff = sum(abs(current_row - x))
    
    return(diff)
  })
  smallest_diffs <- sort(diffs)[1:num.neighbors]
  # Find the indices of the neighbors
  neighbor_indices <- which(diffs %in% smallest_diffs)[1:num.neighbors]
  
  # Get case IDs of neighbors and pad with NA if less than 4 neighbors found
  out <- data.steps[neighbor_indices, "case.id"]
  if (length(out) < num.neighbors) {
    pad_length <- num.neighbors - length(out)
    out <- c(out, rep(NA, pad_length))
  }
  return(out)
}


# Apply the function to each row of data.steps
data.steps.with.neighbors <- apply(data.steps, 1, find_one_step_neighbors)
exper.design[, neighbor_col_names] <- t(data.steps.with.neighbors)

neighbors <- exper.design %>%
  pivot_longer(cols = starts_with("neighbor_"),
               names_to = "name",
               values_to = "neighbor_num") 
neighbors$neighbor_solved <- FALSE

aux$c1 <- aux$c
aux$c2 <- aux$c*aux$b2

write.csv(aux,"INTERMEDIATE/SeedSolutionsFromMuForLambda_2group.csv")

write.csv(neighbors,"INTERMEDIATE/neighbors02groups_lambda.csv")
iter.num <- 0
save(iter.num,file="iter.num.2groups_lambda.Rdata")

