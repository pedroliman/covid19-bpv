FindUniqueSols <- function(mult,costs){
  

  costs$cost <- round(as.numeric(costs$cost),digits=4)
  costs$cost2 <- round(as.numeric(costs$cost2),digits=4)
  
dup <- duplicated(costs[,c("cost","cost2","R0","c","i0","tFinal")])

costs2 <- costs[!dup,]

costs2$split.var <- as.factor(paste(costs2$i0,costs2$R0,costs2$c,costs2$tFinal,sep="_"))
costs2$sol.id <- paste(costs2$case.id,costs2$case.id.c,sep="_")
costs.split <- split(costs2,costs2$split.var)

findMin <- function(x){
  if(dim(x)[1]==1){
    return(x)
  } else {
    min <- which(x$cost2==min(x$cost2))
   return(x[min,])
  }
}

min.costs <- lapply(costs.split,findMin)
min.costs <- do.call(rbind,min.costs)

not.min <- subset(costs2,!(sol.id %in% min.costs$sol.id))

mult$sol.id <- paste(mult$case.id,mult$case.id.c,sep="_")
mult <- subset(mult,sol.id %in% costs2$sol.id)
mult$is.min <- FALSE
mult[mult$sol.id %in% min.costs$sol.id,"is.min"] <- TRUE
mult <- merge(mult,costs2,all.x = TRUE)
return(mult)
}