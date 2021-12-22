tFinal <- seq(8,80,by=4)
R0 <- seq(1.5,6,by=1.5)
i0 <- c(0.0001,0.005)
sFinal <- seq(0.01,1,by=0.01)

exper.design.mu <- expand.grid(tFinal,R0,i0,sFinal)

names(exper.design.mu) <- c("tFinal","R0","i0","sFinal")

exper.design.mu$case.id.mu <- 1:dim(exper.design.mu)[1]

write.csv(exper.design.mu,"INTERMEDIATE/exper.design.mu.csv")
