LH.Sampling <- function(n.samp, parameters, constraints = NULL){ 
  if(!is.null(constraints)){
    if(any(constraints %in% "Mode")){
      r <- data.frame(NULL)
      r.add <- Mode.Sample(parameters)
      for(i in 1:n.samp){       
      r <- rbind(r,r.add) }
      colnames(r) <- names(r.add)  
      rownames(r) <- c(1:n.samp)  
      return(r)   
    } 
    constraints.string <- constraints[1]
    if(length(constraints)>1)for(i in 2: length(constraints)){
    constraints.string <- paste(constraints.string, " & ",constraints[i], sep="")}
  } 
  
   
  lhs.tab <- t(apply(parameters,1,UniTri.Inputs))
    
  r <- lhs(n.samp, lhs.tab[,c("Lower","Upper")], lhs.tab[,"Shape"], lhs.tab[,"Mode"]) 
  colnames(r) <- rownames(parameters)
  ## To Do: Trim sample removing constraints
  r <- as.data.frame(r)
  if(!is.null(constraints)){ 
  r <- subset(r,eval(parse(text= constraints.string)))   
  ## Repopulate the lost samples.  
  while(nrow(r)<n.samp){
    r.add <- lhs(n.samp-nrow(r), lhs.tab[,c("Lower","Upper")], lhs.tab[,"Shape"], lhs.tab[,"Mode"])
    colnames(r.add) <- rownames(parameters)
    r <- rbind(r,as.data.frame(r.add)) 
    r <- subset(r,eval(parse(text= constraints.string)))
  }
  }
  rownames(r) <- c(1:n.samp)  
  return(r)                           
}