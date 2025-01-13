
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


twoComparmentSIR <- function(x, y, parms) {
  
with(as.list(c(y, parms)),{
  s1 <- y[1]
  s2 <- y[2]
  i1 <- y[3]
  i2 <- y[4]
  lambda1 <- y[5]
  lambda2 <- y[6]
  eta1 <- y[7]
  eta2 <- y[8]
  c <- c1
  b <- c2/c1
  
  rho11 <- (c*r11*a11)/(c*a11+r11*s1*i1*(1+lambda1-eta1))
  rho12 <- (c*r12*a12)/(c*a12+r12*s1*i2*(1+lambda1-eta1))
  rho21 <- (c*r21*a21)/(c*a21+r21*s2*i1*(c2/c1+lambda2-eta2))
  rho22 <- (c*r22*a22)/(c*a22+r22*s2*i2*(c2/c1+lambda2-eta2))
  
  ds1 <- -s1*rho11*i1-s1*rho12*i2
  ds2 <- -s2*rho21*i1-s2*rho22*i2
  di1 <- s1*rho11*i1+s1*rho12*i2-i1
  di2 <- s2*rho21*i1+s2*rho22*i2-i2
  
  dlambda1 <- (1+rho11*i1+rho12*i2)*(lambda1-eta1)
  dlambda2 <- (c2/c1+rho21*i1+rho22*i2)*(lambda2-eta2)
  deta1 <- s1*rho11*(1+lambda1-eta1)+s2*rho21*(b+lambda2-eta2)+eta1
  deta2 <- s1*rho12*(1+lambda1-eta1)+s2*rho22*(b+lambda2-eta2)+eta2
  
  g <- function(x){
    y <- -log(x)+x-1
    return(y)
  }
  dCost <- a11*rho11*g(rho11/r11)+a12*rho12*g(rho12/r12)+a21*rho21*g(rho21/r21)+
    a22*rho22*g(rho22/r22)+s1*rho11*i1+s1*rho12*i2+c2/c1*(s2*rho21*i1+s2*rho22*i2)
  
  if(is.nan(dCost)){
    dCost <- 0
  }
  
  return(list(c(ds1,ds2,di1,di2,dlambda1,dlambda2,deta1,deta2,dCost)))
})
}

iter.num <- get(base::load("iter.num.2groups_lambda.Rdata"))
#exper.design <- read.csv(paste0("INTERMEDIATE/exper.design.mu",iter.num,".2groups.csv"),row.names = NULL)
exper.design <- read.csv(paste0("INTERMEDIATE/exper.design_2groups_lambda.csv"),row.names = NULL)

neighbors <- read.csv(paste0("INTERMEDIATE/neighbors",iter.num,"2groups_lambda.csv"),row.names = NULL)
sol.explore.try <- read.csv(paste0("INTERMEDIATE/SeedSolutionsFromMuForLambda_2group.csv"),row.names = NULL)
sol.explore<- sol.explore.try[0,]
sol.explore$X <- NULL
    solved.count <- 0
    
aux <- subset(sol.explore.try,tau==tFinal&!is.na(case.id))
#aux <- subset(neighbors,case.id.mu %in% sol.explore$case.id.mu)
#aux <- aux[!duplicated(aux$case.id.mu),]
#aux <- aux[-1,]
print(dim(aux))
solved.cases <- c()
for(i in 1:dim(aux)[1]){
#print(Sys.time())
  pos <- which(exper.design$case.id==aux[i,"case.id"])

  case.id <- exper.design$case.id[pos]
if(!(case.id %in% solved.cases)){
   tFinal <- exper.design$tFinal[pos]
   donor.case.mu <- aux[i,"case.id.mu"]
   # c needs to become s (select s final from exper.design)
   R_0 <- exper.design$R0[pos]
   i0 <- exper.design$i0[pos]
   c1 <- exper.design$c1[pos]
   c2 <- exper.design$c2[pos]
   
    t <- seq(0,tFinal,.2)
    t.temp <- t
    
    #s <- y[1] 
    #i <- y[2] 
    #mu_s <- y[3]
    #mu_i <- y[4]
    s.init <- (1-i0)/2
    i.init <- i0/2
    yini = c(s.init,s.init,i.init,i.init,NA,NA,NA,NA,0)
    yend = c(NA,NA,NA,NA,0,0,0,0,NA)


  unique.costs <- c()
 
    solved <- FALSE
   
     guess <- subset(sol.explore.try,case.id.mu %in% donor.case.mu)
     guess <- guess[!duplicated(guess$tau),]
     xguess <- guess[,"tau"]
     yguess <- t(guess[,c("s1","s2","i1","i2","lambda1","lambda2","eta1","eta2","Cost")])
     res <- try(bvptwp(yini = yini,yend=yend,x=t,parms = list(r11=R_0,r12=R_0,r21=R_0,r22=R_0,a11=.25,a12=.25,a21=.25,a22=.25,c1=c1,c2=c2),func = twoComparmentSIR,
                       xguess = xguess,yguess = yguess,nmax=length(t)*5))


      solved <- !inherits(res, "try-error")
      
      if(inherits(res,"try-error")){
      error_message <- conditionMessage(attr(res, "condition"))
      
      # Check if the error is 'subscript out of bounds'
      if (grepl("subscript out of bounds", error_message)) {
        # Print the specific error message and stop execution
        stop("Execution stopped due to 'subscript out of bounds' error: ", error_message)
      } else {
        # Print other errors and allow continuation
        print("Non-critical error detected:")
        print(error_message)
      }
}


      #gc()
     # print(solved)
     # print(gc())
   
      if(solved==TRUE){
        solved.cases <- c(solved.cases,case.id)
        #print(solved.count)
        solved.count <- solved.count +1
        sol <- as.data.frame(res)
        sol <- sol[,1:10]
        sol <- sol[!duplicated(sol),]
        sol.cost <- round(tail(sol[,10],n=1),10)
        names(sol) <- c("x",as.character(1:9))
        
        sol <- as.data.frame(sol) %>% 
          rename(tau=x,s1='1',s2='2',i1='3',i2='4',lambda1='5',lambda2='6',eta1='7',eta2='8',Cost='9') %>% 
          mutate(  rho11 = (c1*.25*R_0)/(c1*.25+.25*s1*i1*(1+lambda1-eta1)),
                   rho12 = (c1*.25*R_0)/(c1*.25+.25*s1*i2*(1+lambda1-eta1)),
                   rho21 = (c1*.25*R_0)/(c1*.25+.25*s2*i1*(c2/c1+lambda2-eta2)),
                   rho22 = (c1*.25*R_0)/(c1*.25+.25*s2*i2*(c2/c1+lambda2-eta2)),
                 tFinal=tFinal,
                 case.id.mu = case.id)
        sol <- subset(sol,tau %in% t)
        sol$R0 <- R_0
        sol$tFinal <- tFinal
        sol$c1 <- c1
        sol$c2 <- c2
        sol$i0 <- i0
        sol$solved <- solved
        exper.design[pos,"solved"] <- TRUE
        neighbors[which(neighbors$neighbor_num==case.id),"neighbor_solved"] <- TRUE
        neighbors[which(neighbors$case.id.mu==case.id),"solved"] <- TRUE
        write.table(
          sol,
          file = paste0("INTERMEDIATE/lambda.solutions.2groups.final.csv"),
          sep = ",",         # Set the separator according to your CSV format
          col.names = FALSE, # Don't write column names if the file already exists
          row.names = FALSE, # Don't write row names
          append = TRUE      # Append to the existing file
        )
        
      } else {
        row <- intersect(which(neighbors$neighbor_num==case.id),which(neighbors$case.id.mu==donor.case))
        neighbors[row,"neighbor_solved"] <- "TF" 
        print("tried, failed")
      }
 
    
} else {
  print("already solved")
}
}

 #   print(solved)
#print(Sys.time())


#print(solved.count)
iter.num <- iter.num+1
#write.csv(exper.design,paste0("INTERMEDIATE/exper.design.mu",iter.num,".csv"),row.names=FALSE)
write.csv(neighbors,paste0("INTERMEDIATE/neighbors",iter.num,"2groups_lambda.csv"),row.names=FALSE)

save(iter.num,file="iter.num.2groups_lambda.Rdata")
stop()
}


