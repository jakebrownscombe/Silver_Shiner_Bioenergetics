#functions for silver shiner bioenergetics model

grow <- function(data, meta, ndays,  p){
  
  for (i in 1:ndays){
    
    #consumption
    data$V[i] <- (meta$CTM - data$temp[i]) / (meta$CTM - meta$CTO)
    data$ft[i] <- ifelse(data$temp[i] < meta$CTM, data$V[i]^CX * exp(CX * (1 - data$V[i])), 0)
    data$Cmax[i] <- meta$CA * data$weight[i] ^ meta$CB
    data$Cons.p[i] <- p
    data$C[i] <- data$Cmax[i] * data$Cons.p[i] * data$ft[i]
    data$Cons.g[i] <- data$C[i]*data$weight[i]
    data$Cons.J[i] <- data$Cons.g[i]*meta$EDP 
    
    #wastes
    data$Eg[i] = meta$FA*data$Cons.J[i] 
    data$Ex[i] = meta$UA*(data$Cons.J[i]-data$Eg[i])
    
    #SDA
    data$SDA[i] <- meta$SDA *(data$Cons.J[i]-data$Eg[i])  
    
    #metabolism
    data$R.ft[i] <- exp(meta$RQ*data$temp[i])
    data$Rmax[i] <- meta$RA * data$weight[i] ^ meta$RB
    data$Met.J[i] <- data$Rmax[i] * data$R.ft[i] * data$weight[i] * meta$oxycal
    data$Met.act[i] <- data$Met.J[i]+(data$Met.J[i]*(meta$BACT*data$temp[i]))
    
    #growth
    data$Growth.J[i] <- data$Cons.J[i] - (data$Met.act[i]+data$SDA[i]+data$Eg[i]+data$Ex[i])
    data$Growth.g[i] <- data$Growth.J[i]/meta$ED 
    
    #update weight
    data$E[i+1] <- data$E[i]+data$Growth.J[i]
    data$weight[i+1] <- data$E[i+1]/meta$ED 
  }
  
  #outputs
  data <- data %>% filter(day<(nrow(data)-1))
  FW <- data$weight[nrow(data)] 
  outputs <- list(data, FW)
  return(outputs)
}


fit.p <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow(data=SS.grow, meta=shiner.meta, ndays=nrow(SS.grow)-1,  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow(data=SS.grow, meta=shiner.meta, ndays=nrow(SS.grow)-1,  p)
    W.p <- g[[2]]
  }
  return(p)
}  






