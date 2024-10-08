CalculateCost <- function(all.sol){
  
  if("case.id.c" %in% names(all.sol)){
    all.sol$split.var <- as.factor(paste(all.sol$i0,all.sol$R0,all.sol$c,all.sol$tFinal,all.sol$case.id.c,sep="_"))
  } else {
    all.sol$split.var <- as.factor(paste(all.sol$i0,all.sol$R0,all.sol$c,all.sol$tFinal,sep="_"))
    
}
  all.sol.split <- split(all.sol,as.factor(all.sol$split.var))

  CostSummary <- function(x){
   # print("Case id is:")
   # print(x$case.id[1])
   # print("Case id c is:")
   # print(x$case.id.c[1])
    
    x <- x[order(as.numeric(x$tau)),]
    x <- x[!duplicated(x),]
    x$c <- as.numeric(x$c)
    x$s <- as.numeric(x$s)
    x$R_tau <- as.numeric(x$R_tau)
    x$i <- as.numeric(x$i)
    x$dcost <- -x$c*log(as.numeric(x$R_tau)/as.numeric(x$R0))+x$c*(as.numeric(x$R_tau)/as.numeric(x$R0))+x$R_tau*x$s*x$i-as.numeric(x$c)
    x$dcost2 <- -x$c*log(as.numeric(x$R_tau)/as.numeric(x$R0))+x$c*(as.numeric(x$R_tau)/as.numeric(x$R0))-as.numeric(x$c)
    x$total.epi <- 1-x$s[length(x$s)]
    cost <- AUC(as.numeric(x$tau),x$dcost)
    cost2 <- AUC(as.numeric(x$tau),x$dcost2)+1-(tail(as.numeric(x$s),1)+tail(as.numeric(x$i),1))
    if("case.id.c" %in% names(x)){
    return(as.data.frame(cbind(cost=cost,cost2=cost2,R0 = x$R0[1],c=x$c[1],i0 = x$i0[1],tFinal=x$tFinal[1],case.id=x$case.id[1],case.id.c=x$case.id.c[1],
                               split.var=x$split.var[1])))
    } else {
      return(as.data.frame(cbind(cost=cost,cost2=cost2,R0 = x$R0[1],c=x$c[1],i0 = x$i0[1],tFinal=x$tFinal[1],case.id=x$case.id[1],split.var=x$split.var[1],
                                 total.epi=x$total.epi[1])))
      
    }
  }

sol.costs <- lapply(all.sol.split,CostSummary)
sol.costs <- as.data.frame(do.call(rbind,sol.costs))
row.names(sol.costs) <- sol.costs$split.var
return(sol.costs)
}