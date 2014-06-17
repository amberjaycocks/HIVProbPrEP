Mode.Sample <- function(param){
 r <- as.array(param$Mode)
 names(r) <- rownames(param) 
 return(r)
}