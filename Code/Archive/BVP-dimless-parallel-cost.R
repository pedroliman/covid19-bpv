
# Needed packages to perform the analysis
load <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} 

packages <- c("knitr","dplyr",
              "ggplot2","bvpSolve","deSolve","class","pracma")
load(packages)

### get the arguments
args = commandArgs(TRUE)
pos <- as.numeric(args[1])+as.numeric(args[2])
R0.cases <- args[3]
if(is.na(pos)){
  pos <- 1
  R0.cases <- 6
}
oneComparmentSIR <- function(x, y, parms) {

    with(as.list(c(y, parms)),{
        s <- y[1] 
        i <- y[2] 
        lambda_s <- y[3]
        lambda_i <- y[4]
        
        R_tau = R_0*c/(c+R_0*s*i*(1+lambda_s-lambda_i))
        
        ds <- -R_tau *s *i
        
        di <- i * (R_tau*s-1)
        
        dlambda_s <- R_tau*i*(1+lambda_s-lambda_i)
        
        dlambda_i <- R_tau*s*(1+lambda_s-lambda_i)+lambda_i
        
        dCost <- R_tau*s*i+c*(-log(R_tau/R_0)+R_tau/R_0-1)

    return(list(c(ds,di,dlambda_s,dlambda_i,dCost)))
    })
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


sol.explore <- read.csv("INTERMEDIATE/sol.explore13.csv")
sol.explore$case.id <- paste0(sol.explore$case.id,sol.explore$case.id.c)
sol.explore <- sol.explore[,c("tau","s","i","lambda_s","lambda_i","R_tau","Rt","tFinal","case.id")]
sol.explore.p <- read.csv("INTERMEDIATE/sol.explore.all.csv")
sol.explore.p <- subset(sol.explore.p,solved == TRUE)[,c("tau","s","i","lambda_s","lambda_i","R_tau","Rt","tFinal","case.id")]
sol.explore.p2 <- read.csv("INTERMEDIATE/parallel.output3.csv")
sol.explore.p2 <- subset(sol.explore.p2,solved == TRUE)[,c("tau","s","i","lambda_s","lambda_i","R_tau","Rt","tFinal","case.id")]
sol.explore.p3 <- read.csv("INTERMEDIATE/output-R03.csv")
sol.explore.p3 <- subset(sol.explore.p3,solved == TRUE)[,c("tau","s","i","lambda_s","lambda_i","R_tau","Rt","tFinal","case.id")]
sol.explore <- rbind(sol.explore.p,sol.explore,sol.explore.p3)
case.id.p <- unique(sol.explore$case.id)
exper.design <- read.csv("INTERMEDIATE/exper.design-round13.csv")
exper.design$bvp.success[which(exper.design$case.id %in% case.id.p)] <- TRUE
exper.design.s <- subset(exper.design,bvp.success==TRUE)
exper.design.u <- subset(exper.design,R0 == R0.cases & case.id>417)
sol.explore <- subset(sol.explore,!is.na(tau)&!grepl("NA",case.id))


if(!("Cost" %in% names(sol.explore))){
  sol.explore$Cost <- 0
}
print(Sys.time())
        c.tried <- c()
        
  case.id <- exper.design.u$case.id[pos]
   tFinal <- exper.design.u$tFinal[pos]
   c <- exper.design.u$c[pos]
   R_0 <- exper.design.u$R0[pos]
   i0 <- exper.design.u$i0[pos]
    t <- seq(0,tFinal,1/7)
    t.temp <- t
    yini = c(1-i0,i0,NA,NA,0)
    yend = c(NA,NA,0,0,NA)
    solved <- FALSE

# exper.design.s$R0.dist <- (exper.design.s$R0-R_0)^2/R_0^2
# exper.design.s$i0.dist <- (exper.design.s$i0-i0)^2/i0^2
# exper.design.s$c.dist <- (exper.design.s$c-c)^2/c^2
# exper.design.s$tFinal.dist <- (exper.design.s$tFinal-tFinal)^2/tFinal^2
# exper.design.s$dist <- exper.design.s$R0.dist+exper.design.s$i0.dist+exper.design.s$c.dist+exper.design.s$tFinal.dist
# same.c <- exper.design.s[order(exper.design.s$dist),]

  unique.costs <- c()
 
    solved <- FALSE
    solved.count <- 0
    
### Start the loop
    for(case.id.c in case.id.p){
      print(paste("Case id c is",case.id.c))
      guess <- subset(sol.explore,case.id == case.id.c)
      
     if((solved == FALSE | solved.count < 20|length(unique.costs<3))&dim(guess)[1]>0){
      worked <- TRUE
if(tFinal>=max(guess$tau)){
    # Here we try the two-point method. Note we silenced the error so if it doesn't converge
    # we use the shooting method solution. 
  tFinal.start <- max(guess$tau)
  tFinal.end <- tFinal
  tFinal.vec <- seq(tFinal.start,tFinal.end,1)
} else {
  tFinal.start <- tFinal
  tFinal.end <- tFinal
  tFinal.vec <- seq(tFinal.start,tFinal.end,1)
}
  #print(paste("tFinal vec is",tFinal.vec))
        for(tFinal in tFinal.vec){
          if(worked){
                      print(paste("tFinal is",tFinal))

      t <- seq(0,tFinal,1/7)
guess <- as.data.frame(guess)
     if(length(t) < dim(guess)[1]){
        guess <- guess[1:length(t),]
      } else if(length(t) > dim(guess)[1]){
        guess[(dim(guess)[1]+1):length(t),]<- guess[dim(guess)[1],]
      }
    res <- try(guess <- bvptwp(yini = yini,yend=yend,x=t,parms = list(R_0 = R_0,c=c),
                         func = oneComparmentSIR,xguess = t,yguess = t(guess[,c("s","i","lambda_s","lambda_i","Cost")]),nmax = (7*tFinal+3)),
               silent=FALSE)
        t.temp <- guess[,1]
    print(length(t.temp))
    worked <- !inherits(res,"try-error")
    print(paste("did it work?",worked))
          if(worked){
              c.tried <- c(c.tried,tFinal.start)
              guess <- as.data.frame(guess)
              names(guess) <- c("tau","s","i","lambda_s","lambda_i","Cost")
              sol <- guess
                    }
          }
        }
      solved <- !inherits(res, "try-error")
      if(solved==TRUE){
        print(solved.count)
        solved.count <- solved.count +1
        sol <- as.data.frame(sol)
        sol <- sol[,1:6]
        sol.cost <- round(tail(sol[,6],n=1),6)
        names(sol) <- c("x",as.character(1:5))
        
        sol <- as.data.frame(sol) %>% 
          rename(tau=x,s='1',i='2',lambda_s='3',lambda_i='4',Cost='5') %>% 
          mutate(R_tau = R_0*c/(c+R_0*s*i*(1+lambda_s-lambda_i)),
                 Rt= R_tau*s,
                 tFinal=tFinal,
                 case.id = case.id)
        sol$case.id <- case.id
        sol$R0 <- R_0
        sol$tFinal <- tFinal
        sol$c <- c
        sol$i0 <- i0
        sol$solved <- solved
        sol$case.id.c <- case.id.c
        if(solved.count == 1){
          sol.explore2 <- sol
        } else {
          sol.explore2  <- bind_rows(sol.explore2,sol)
        } 
        if(!(sol.cost %in% unique.costs)){
          write.csv(sol,paste0("INTERMEDIATE/sol.exploreR0",R0.cases,"seed",case.id.c,".case-",case.id,".csv"),row.names = FALSE)
          unique.costs <- unique(c(unique.costs,sol.cost))
        }
      }
      
  }

   

 
}
    
          

    print(solved)
print(Sys.time())



