
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

    return(list(c(ds,di,dlambda_s,dlambda_i)))
    })
}


sol.explore <- read.csv("INTERMEDIATE/sol.explore-step12.csv")
sol.explore <- sol.explore[,c("tau","s","i","lambda_s","lambda_i","R_tau","Rt","tFinal","case.id")]
sol.explore.p <- read.csv("INTERMEDIATE/sol.explore.all.csv")
sol.explore.p <- subset(sol.explore.p,solved == TRUE)[,c("tau","s","i","lambda_s","lambda_i","R_tau","Rt","tFinal","case.id")]
sol.explore <- rbind(sol.explore.p,sol.explore)
case.id.p <- unique(sol.explore.p$case.id)
exper.design <- read.csv("INTERMEDIATE/exper.design-round13.csv")
exper.design$bvp.success[which(exper.design$case.id %in% case.id.p)] <- TRUE
exper.design.s <- subset(exper.design,bvp.success==TRUE)
exper.design.u <- subset(exper.design,(is.na(bvp.success)|bvp.success==FALSE) & case.id>417)
print(Sys.time())
        c.tried <- c()
        
  case.id <- exper.design.u$case.id[pos]
   tFinal <- exper.design$tFinal[pos]
   c <- exper.design$c[pos]
   R_0 <- exper.design$R0[pos]
   i0 <- exper.design$i0[pos]
    t <- seq(0,tFinal,1/7)
    t.temp <- t
    yini = c(1-i0,i0,NA,NA)
    yend = c(NA,NA,0,0)
    solved <- FALSE
    
 
  same.c <- exper.design.s[exper.design.s$tFinal <= tFinal&exper.design.s$R0==R_0,]
 # dists <- dist.matrix[same.c$case.id,case.run]
#  cases <- names(sort(dists)[1:5])
#  same.c <- subset(exper.design.s,case.id %in% cases)
    solved <- FALSE
    for(case.id.c in same.c$case.id){
     if(solved == FALSE){
      guess <- subset(sol.explore,case.id == case.id.c)
      worked <- TRUE
if(tFinal>=same.c[which(same.c$case.id==case.id.c),"tFinal"]&!(same.c[which(same.c$case.id==case.id.c),"tFinal"] %in% c.tried)){
    # Here we try the two-point method. Note we silenced the error so if it doesn't converge
    # we use the shooting method solution. 
  tFinal.start <- same.c[which(same.c$case.id==case.id.c),"tFinal"]
  tFinal.end <- tFinal
  tFinal.vec <- seq(tFinal.start,tFinal.end,1)
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
                         func = oneComparmentSIR,xguess = t,yguess = t(guess[,c("s","i","lambda_s","lambda_i")]),nmax = 20*tFinal),
               silent=FALSE)
        t.temp <- guess[,1]
    print(length(t.temp))
    worked <- !inherits(res,"try-error")
    print(paste("did it work?",worked))
    if(worked){
      c.tried <- c(c.tried,tFinal.start)
    names(guess) <- c("tau","s","i","lambda_s","lambda_i")
    sol <- guess
    }
          }
    }
}
     }
    

    if(!inherits(res, "try-error")){
   solved <- TRUE
    }
          
    if(solved==TRUE){
      sol.explore <- subset(sol.explore,case.id!=case.run)
    sol <- as.data.frame(sol)
    sol <- sol[,1:5]
    names(sol) <- c("x",as.character(1:4))
    
    sol <- as.data.frame(sol) %>% 
    rename(tau=x,s='1',i='2',lambda_s='3',lambda_i='4') %>% 
    mutate(R_tau = R_0*c/(c+R_0*s*i*(1+lambda_s-lambda_i)),
           Rt= R_tau*s,
           tFinal=tFinal,
           case.id = case.id)
    
    sol.explore  <- bind_rows(sol.explore,sol)
    } 
     }
    print(solved)
print(Sys.time())
if(solved == TRUE){
sol$case.id <- case.id
sol$R0 <- R_0
sol$tFinal <- tFinal
sol$c <- c
sol$i0 <- i0
sol$solved <- solved
write.csv(sol,paste0("INTERMEDIATE/sol.explore.case-",case.id,".csv"),row.names = FALSE)
write.csv(exper.design,paste0("INTERMEDIATE/exper.design.case-",case.id,".csv"),row.names = FALSE)
}

