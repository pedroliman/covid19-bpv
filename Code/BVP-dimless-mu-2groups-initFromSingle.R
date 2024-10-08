
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
              "ggplot2","bvpSolve","class","pracma","tidyr")
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

twoComparmentSIRmu <- function(x, y, parms) {
  
with(as.list(c(y, parms)),{
  s1 <- y[1]
  s2 <- y[2]
  i1 <- y[3]
  i2 <- y[4]
  mu1 <- y[5]
  mu2 <- y[6]
  nu1 <- y[7]
  nu2 <- y[8]
  
  
  rho11 <- (r11*a11)/(a11+r11*s1*i1*(mu1-nu1))
  rho12 <- (r12*a12)/(a12+r12*s1*i2*(mu1-nu1))
  rho21 <- (r21*a21)/(a21+r21*s2*i1*(mu2-nu2))
  rho22 <- (r22*a22)/(a22+r22*s2*i2*(mu2-nu2))
  
  ds1 <- -s1*rho11*i1-s1*rho12*i2
  ds2 <- -s2*rho21*i1-s2*rho22*i2
  di1 <- s1*rho11*i1+s1*rho12*i2-i1
  di2 <- s2*rho21*i1+s2*rho22*i2-i2
  
  dmu1 <- (rho11*i1+rho12*i2)*(mu1-nu1)
  dmu2 <- (rho21*i1+rho22*i2)*(mu2-nu2)
  dnu1 <- s1*rho11*(mu1-nu1)+s2*rho21*(mu2-nu2)+nu1
  dnu2 <- s1*rho12*(mu1-nu1)+s2*rho22*(mu2-nu2)+nu2
  
  g <- function(x){
    y <- -log(x)+x-1
    return(y)
  }
  dCost <- a11*rho11*g(rho11/r11)+a12*rho12*g(rho12/r12)+a21*rho21*g(rho21/r21)+
    a22*rho22*g(rho22/r22)
  
  if(is.nan(dCost)){
    dCost <- 0
  }
  
  return(list(c(ds1,ds2,di1,di2,dmu1,dmu2,dnu1,dnu2,dCost)))
})
}

iter.num <- get(base::load("iter.num.2groups.Rdata"))
#exper.design <- read.csv(paste0("INTERMEDIATE/exper.design.mu",iter.num,".2groups.csv"),row.names = NULL)
exper.design <- read.csv(paste0("INTERMEDIATE/exper.design.mu0.2groups.csv"),row.names = NULL)

neighbors <- read.csv(paste0("INTERMEDIATE/mu.neighbors",iter.num,".2groups.csv"),row.names = NULL)
sol.explore <- read.csv(paste0("INTERMEDIATE/mu.solutions.2groups.csv"),row.names = NULL)

sol.explore$X <- NULL
    solved.count <- 0
    
#aux <- subset(neighbors,solved==TRUE&neighbor_solved==FALSE&!is.na(neighbor_num))
aux <- subset(neighbors,case.id.mu %in% sol.explore$case.id.mu)
aux <- aux[!duplicated(aux$case.id.mu),]
#aux <- aux[-1,]
for(i in 1:dim(aux)[1]){
print(Sys.time())
  pos <- which(exper.design$case.id.mu==aux[i,"neighbor_num"])
  pos <- which(exper.design$case.id.mu==aux[i,"case.id.mu"])
  
  case.id <- exper.design$case.id.mu[pos]
  aux[i,"case.id.mu"]
  print(case.id)
   tFinal <- exper.design$tFinal[pos]
   donor.case <- aux[i,"case.id.mu"]
   # c needs to become s (select s final from exper.design)
   sFinal <- exper.design$sFinal[pos]
   R_0 <- exper.design$R0[pos]
   i0 <- exper.design$i0[pos]
   delta <- exper.design$delta[pos]
    t <- seq(0,tFinal,.2)
    t.temp <- t
    
    #s <- y[1] 
    #i <- y[2] 
    #mu_s <- y[3]
    #mu_i <- y[4]
    s.init <- (1-i0)/2
    i.init <- i0/2
    yini = c(s.init,s.init,i.init,i.init,NA,NA,NA,NA,0)
    yend = c(sFinal/2,(sFinal-delta)/2,NA,NA,NA,NA,0,0,NA)
    solved <- FALSE


  unique.costs <- c()
 
    solved <- FALSE
     guess <- subset(sol.explore,case.id.mu == donor.case)
     guess <- guess[!duplicated(guess),]
     xguess <- guess[,"tau"]
     yguess <- t(guess[,c("s1","s2","i1","i2","mu1","mu2","nu1","nu2","Cost")])
     res <- try(bvptwp(yini = yini,yend=yend,x=t,parms = list(r11=R_0,r12=R_0,r21=R_0,r22=R_0,a11=.25,a12=.25,a21=.25,a22=.25),func = twoComparmentSIRmu,
                       xguess = xguess,yguess = yguess,nmax=length(t)*5))


      solved <- !inherits(res, "try-error")
      #gc()
      print(solved)
     # print(gc())
   
      if(solved==TRUE){
        print(solved.count)
        solved.count <- solved.count +1
        sol <- as.data.frame(res)
        sol <- sol[,1:10]
        sol <- sol[!duplicated(sol),]
        sol.cost <- round(tail(sol[,10],n=1),10)
        names(sol) <- c("x",as.character(1:9))
        
        sol <- as.data.frame(sol) %>% 
          rename(tau=x,s1='1',s2='2',i1='3',i2='4',mu1='5',mu2='6',nu1='7',nu2='8',Cost='9') %>% 
          mutate(  rho11 = (.25*R_0)/(.25+R_0*s1*i1*(mu1-nu1)),
                   rho12 = (.25*R_0)/(.25+R_0*s1*i2*(mu1-nu1)),
                   rho21 = (.25*R_0)/(.25+R_0*s2*i1*(mu2-nu2)),
                   rho22 = (.25*R_0)/(.25+R_0*s2*i2*(mu2-nu2)),
                 tFinal=tFinal,
                 case.id = case.id)
        sol <- subset(sol,tau %in% t)
        sol$case.id <- case.id
        sol$R0 <- R_0
        sol$tFinal <- tFinal
        sol$sFinal <- sFinal
        sol$delta <- delta
        sol$i0 <- i0
        sol$solved <- solved
        exper.design[pos,"solved"] <- TRUE
        neighbors[which(neighbors$neighbor_num==case.id),"neighbor_solved"] <- TRUE
        
        neighbors[which(neighbors$case.id.mu==case.id),"solved"] <- TRUE
        write.table(
          sol,
          file = paste0("INTERMEDIATE/mu.solutions.2groups.final.csv"),
          sep = ",",         # Set the separator according to your CSV format
          col.names = FALSE, # Don't write column names if the file already exists
          row.names = FALSE, # Don't write row names
          append = TRUE      # Append to the existing file
        )
        
      } else {
        row <- intersect(which(neighbors$neighbor_num==case.id),which(neighbors$case.id.mu==donor.case))
        neighbors[row,"neighbor_solved"] <- "TF" 
      }
 
    
}          

    print(solved)
print(Sys.time())


print(solved.count)
iter.num <- iter.num+1
write.csv(exper.design,paste0("INTERMEDIATE/exper.design.mu",iter.num,".csv"),row.names=FALSE)
write.csv(neighbors,paste0("INTERMEDIATE/mu.neighbors",iter.num,"2groups.csv"),row.names=FALSE)

save(iter.num,file="iter.num.Rdata")

}


