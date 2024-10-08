
# Needed packages to perform the analysis
solved.count <- 1
while(solved.count !=0){
  rm(list=ls())
  gc()
load <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} 

packages <- c("knitr","dplyr",
              "ggplot2","bvpSolve","deSolve","class","pracma","tidyr")
load(packages)

knitr::opts_chunk$set(warning = FALSE, message = FALSE) 

### get the arguments
args = commandArgs(TRUE)
pos <- as.numeric(args[1])+as.numeric(args[2])
R0.cases <- args[3]
if(is.na(pos)){
  pos <- 1
  R0.cases <- 6
}

oneComparmentSIRmu <- function(x, y, parms) {
  
  with(as.list(c(y, parms)),{
    s <- y[1] 
    i <- y[2] 
    mu_s <- y[3]
    mu_i <- y[4]
    
    R_tau = R_0/(1+R_0*s*i*(mu_s-mu_i))
    
    ds <- -R_tau *s *i
    
    di <- i * (R_tau*s-1)
    
    dmu_s <- R_tau*i*(mu_s-mu_i)
    
    dmu_i <- R_tau*s*(1+mu_s-mu_i)+mu_i
    
    dCost <- (-log(R_tau/R_0)+R_tau/R_0-1)
    
    return(list(c(ds,di,dmu_s,dmu_i,dCost)))
  })
}

iter.num <- get(base::load("iter.num.Rdata"))
exper.design <- read.csv(paste0("INTERMEDIATE/exper.design.mu",iter.num,".csv"),row.names = NULL)
neighbors <- read.csv(paste0("INTERMEDIATE/mu.neighbors",iter.num,".csv"),row.names = NULL)
sol.explore <- read.csv(paste0("INTERMEDIATE/mu.solutions",iter.num,".csv"),row.names = NULL)

sol.explore$X <- NULL
    solved.count <- 0
    
aux <- subset(neighbors,solved==TRUE&neighbor_solved==FALSE&!is.na(neighbor_num))
aux <- aux[-1,]
for(i in 1:dim(aux)[1]){
print(Sys.time())
  pos <- which(exper.design$case.id.mu==aux[i,"neighbor_num"])
      
  case.id <- exper.design$case.id.mu[pos]
  aux[i,"case.id.mu"]
  print(case.id)
   tFinal <- exper.design$tFinal[pos]
   donor.case <- aux[i,"case.id.mu"]
   # c needs to become s (select s final from exper.design)
   sFinal <- exper.design$sFinal[pos]
   R_0 <- exper.design$R0[pos]
   i0 <- exper.design$i0[pos]
    t <- seq(0,tFinal,1/7)
    t.temp <- t
    
    #s <- y[1] 
    #i <- y[2] 
    #mu_s <- y[3]
    #mu_i <- y[4]
    yini = c(1-i0,i0,NA,NA,0)
    yend = c(sFinal,NA,NA,0,NA)
    solved <- FALSE


  unique.costs <- c()
 
    solved <- FALSE
guess <- subset(sol.explore,case.id == donor.case)

    res <- try(bvptwp(yini = yini,yend=yend,x=t,parms = list(R_0=R_0),func = oneComparmentSIRmu,
              xguess = guess[,1],yguess = t(guess[,2:6]),nmax=length(t)*5))
    
      solved <- !inherits(res, "try-error")
      print(solved)
     # print(gc())
    
      
      if(solved==TRUE){
        print(solved.count)
        solved.count <- solved.count +1
        sol <- as.data.frame(res)
        sol <- sol[,1:6]
        sol.cost <- round(tail(sol[,6],n=1),6)
        names(sol) <- c("x",as.character(1:5))
        
        sol <- as.data.frame(sol) %>% 
          rename(tau=x,s='1',i='2',mu_s='3',mu_i='4',Cost='5') %>% 
          mutate(R_tau = R_0/(1+R_0*s*i*(mu_s-mu_i)),
                 Rt= R_tau*s,
                 tFinal=tFinal,
                 case.id = case.id)
        sol$case.id <- case.id
        sol$R0 <- R_0
        sol$tFinal <- tFinal
        sol$sFinal <- sFinal
        sol$i0 <- i0
        sol$solved <- solved
        exper.design[pos,"solved"] <- TRUE
        neighbors[which(neighbors$neighbor_num==case.id),"neighbor_solved"] <- TRUE
        
        neighbors[which(neighbors$case.id.mu==case.id),"solved"] <- TRUE

      
          sol.explore  <- bind_rows(sol.explore,sol)
        
      } else {
        row <- intersect(which(neighbors$neighbor_num==case.id),which(neighbors$case.id.mu==donor.case))
        neighbors[row,"neighbor_solved"] <- "TF" 
      }
 
    
          

    print(solved)
print(Sys.time())

}
print(solved.count)
iter.num <- iter.num+1
write.csv(sol.explore,paste0("INTERMEDIATE/exper.design.mu",iter.num,".csv"),row.names=FALSE)
write.csv(exper.design,paste0("INTERMEDIATE/mu.neighbors",iter.num,".csv"),row.names=FALSE)
write.csv(neighbors,paste0("INTERMEDIATE/mu.solutions",iter.num,".csv"),row.names=FALSE)

save(iter.num,file="iter.num.Rdata")

}

