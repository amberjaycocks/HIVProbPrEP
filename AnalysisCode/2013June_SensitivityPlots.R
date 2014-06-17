
library(tgp) 
library(RODBC)
library(MASS)
library(rpart) # for CART routines: rpart = recursive partitioning and regression trees
library(randomForest)
library(gam) # for gam...

if(.Platform$OS.type == "windows"){
  OSdir = "C:/Documents and Settings/rvardava/My Documents/Projects_2011/HIV_Wagner/"
  OSdir = "C:/Documents and Settings/jaycocks/My Documents/20110830_HIV_Wagner/"
  setwd(OSdir)
  
  library.dir   <- "Library/"
  library <- file.path(library.dir, "library.R")
  source(library)
  
  data.dir   <- "Data/"
  figure.dir   <- "Figures/"
  
}else{
  OSdir = "/Users/jaycocks/Documents/20110830_HIV_Wagner/R_model"
  setwd(OSdir)
  
  library.dir   <- "Library"
  library <- file.path(library.dir, "library.R")
  source(library)
  
  data.dir   <- "Data/"
  figure.dir   <- "Figures/"
}



########################################
###                                  ###
###       GET MODEL PARAMETES        ###
###                                  ###
########################################  

pregnancy <- read.csv(paste(data.dir,"pregnancy_table",".csv",sep=""))
colnames(pregnancy) = c("age","Prob.of.conception","Prob.of.delivery")
rownames(pregnancy) <- pregnancy$age

parameters <- read.csv(paste(data.dir,"model.parameters",".csv",sep=""))
parameters <- parameters[1:9,1:5]
rownames(parameters) <- parameters$parameter 

#parameters["h.PrEP",]        <- parameters["eff.PrEP",]
#parameters["h.PrEP","mode"]  <- 1-parameters["eff.PrEP","mode"]
#parameters["h.PrEP","lower"] <- 1-parameters["eff.PrEP","upper"]
#parameters["h.PrEP","upper"] <- 1-parameters["eff.PrEP","lower"]
#parameters <- parameters[!(rownames(parameters) %in% "eff.PrEP") , !(colnames(parameters) %in% "parameter")]
parameters = parameters[,-1]

age = 40 #change this for plots at different ages
sex.acts.N = 6
N = 50
alpha = parameters["alpha","mode"]
h.PrEP = parameters["h.PrEP","mode"]
h.other.STIs = parameters["h.other.STIs","mode"]
h.tx = parameters["h.tx","mode"]
h.late = parameters["h.late","mode"]
h.tx.MTCT = parameters["h.tx.MTCT","mode"]
p.MTCT = parameters["p.MTCT","mode"]

p.conception.act <- pregnancy[pregnancy$age==age,"Prob.of.conception"] / sex.acts.N 
p.delivery <- pregnancy[pregnancy$age==age,"Prob.of.delivery"]

#############
# Treatment #
#############
xa = c(parameters["h.tx",2],parameters["h.tx",3])
xd = (xa[2]-xa[1])/10
tx.use = xa[1]
for(i in 2:10) tx.use = c(tx.use,tx.use[i-1]+xd) 
tx.use = c(tx.use,xa[2],h.tx)

outcome = NULL
for(i in 1:length(tx.use)){
  h.tx = tx.use[i]
  alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late                      

  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)

  outcome.temp = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  
  alpha.overall <- alpha*h.tx    
  h.tx.MTCT = 1
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp2 = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  outcome.temp3 = c(outcome.temp,outcome.temp2)
  outcome = rbind(outcome,outcome.temp3)
}
tx.outcome = cbind(tx.use,outcome)
colnames(tx.outcome) = c("tx.val","neg_succ.all","pos_succ_neg.all","pos_succ_pos.all","neg_unsucc.all","pos_unsucc.all","neg_succ.tx","pos_succ_neg.tx","pos_succ_pos.tx","neg_unsucc.tx","pos_unsucc.tx")
tx.outcome = tx.outcome[order(tx.outcome[,1]),]

########
# PrEP #
########
h.tx = parameters["h.tx","mode"]

xa = c(parameters["h.PrEP",2],parameters["h.PrEP",3])
xd = (xa[2]-xa[1])/10
tx.use = xa[1]
for(i in 2:10) tx.use = c(tx.use,tx.use[i-1]+xd) 
tx.use = c(tx.use,xa[2],h.PrEP)

outcome = NULL
for(i in 1:length(tx.use)){
  h.PrEP = tx.use[i]
  alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late                      
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  
  alpha.overall <- alpha*h.PrEP  
  h.tx.MTCT = 1
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp2 = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  outcome.temp3 = c(outcome.temp,outcome.temp2)
  outcome = rbind(outcome,outcome.temp3)
}
p.outcome = cbind(tx.use,outcome)
colnames(p.outcome) = c("prep.val","neg_succ.all","pos_succ_neg.all","pos_succ_pos.all","neg_unsucc.all","pos_unsucc.all","neg_succ.prep","pos_succ_neg.prep","pos_succ_pos.prep","neg_unsucc.prep","pos_unsucc.prep")
p.outcome = p.outcome[order(p.outcome[,1]),]

########
# STIs #
########
h.PrEP = parameters["h.PrEP","mode"]

xa = c(parameters["h.other.STIs",2],parameters["h.other.STIs",3])
xd = (xa[2]-xa[1])/10
tx.use = xa[1]
for(i in 2:10) tx.use = c(tx.use,tx.use[i-1]+xd) 
tx.use = c(tx.use,xa[2],h.other.STIs)

outcome = NULL
for(i in 1:length(tx.use)){
  h.other.STIs = tx.use[i]
  alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late                      
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  
  alpha.overall <- alpha*h.other.STIs   
  h.tx.MTCT = 1
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp2 = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  outcome.temp3 = c(outcome.temp,outcome.temp2)
  outcome = rbind(outcome,outcome.temp3)
}
sti.outcome = cbind(tx.use,outcome)
colnames(sti.outcome) = c("sti.val","neg_succ.all","pos_succ_neg.all","pos_succ_pos.all","neg_unsucc.all","pos_unsucc.all","neg_succ.sti","pos_succ_neg.sti","pos_succ_pos.sti","neg_unsucc.sti","pos_unsucc.sti")
sti.outcome = sti.outcome[order(sti.outcome[,1]),]

########
# Late #
########
h.other.STIs = parameters["h.other.STIs","mode"]

xa = c(parameters["h.late",2],parameters["h.late",3])
xd = (xa[2]-xa[1])/10
tx.use = xa[1]
for(i in 2:10) tx.use = c(tx.use,tx.use[i-1]+xd) 
tx.use = c(tx.use,xa[2],h.late)

outcome = NULL
for(i in 1:length(tx.use)){
  h.late = tx.use[i]
  alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late                      
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  
  alpha.overall <- alpha*h.late
  h.tx.MTCT = 1
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp2 = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  outcome.temp3 = c(outcome.temp,outcome.temp2)
  outcome = rbind(outcome,outcome.temp3)
}
late.outcome = cbind(tx.use,outcome)
colnames(late.outcome) = c("late.val","neg_succ.all","pos_succ_neg.all","pos_succ_pos.all","neg_unsucc.all","pos_unsucc.all","neg_succ.late","pos_succ_neg.late","pos_succ_pos.late","neg_unsucc.late","pos_unsucc.late")
late.outcome = late.outcome[order(late.outcome[,1]),]

#########
# alpha #
#########
h.late = parameters["h.late","mode"]

xa = c(parameters["alpha",2],parameters["alpha",3])
xd = (xa[2]-xa[1])/10
tx.use = xa[1]
for(i in 2:10) tx.use = c(tx.use,tx.use[i-1]+xd) 
tx.use = c(tx.use,xa[2],alpha)

outcome = NULL
for(i in 1:length(tx.use)){
  alpha = tx.use[i]
  alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late                      
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  
  alpha.overall <- alpha     
  h.tx.MTCT = 1
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp2 = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  outcome.temp3 = c(outcome.temp,outcome.temp2)
  outcome = rbind(outcome,outcome.temp3)
}
alpha.outcome = cbind(tx.use,outcome)
colnames(alpha.outcome) = c("late.val","neg_succ.all","pos_succ_neg.all","pos_succ_pos.all","neg_unsucc.all","pos_unsucc.all","neg_succ.alpha","pos_succ_neg.alpha","pos_succ_pos.alpha","neg_unsucc.alpha","pos_unsucc.alpha")
alpha.outcome = alpha.outcome[order(alpha.outcome[,1]),]

##################
# Haart, tx.MTCT #
##################
alpha = parameters["alpha","mode"]

xa = c(parameters["h.tx.MTCT",2],parameters["h.tx.MTCT",3])
xd = (xa[2]-xa[1])/10
tx.use = xa[1]
for(i in 2:10) tx.use = c(tx.use,tx.use[i-1]+xd) 
tx.use = c(tx.use,xa[2],parameters["h.tx.MTCT",1])

outcome = NULL
for(i in 1:length(tx.use)){
  h.tx.MTCT = tx.use[i]
  alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late                      
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  
  alpha.overall <- alpha                
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp2 = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  outcome.temp3 = c(outcome.temp,outcome.temp2)
  outcome = rbind(outcome,outcome.temp3)
}
haart.outcome = cbind(tx.use,outcome)
colnames(haart.outcome) = c("late.val","neg_succ.all","pos_succ_neg.all","pos_succ_pos.all","neg_unsucc.all","pos_unsucc.all","neg_succ.haart","pos_succ_neg.haart","pos_succ_pos.haart","neg_unsucc.haart","pos_unsucc.haart")
haart.outcome = haart.outcome[order(haart.outcome[,1]),]

########
# age #
######

tx.use = pregnancy$age
outcome = NULL
for(i in 1:length(tx.use)){
  age.use = tx.use[i]
  p.conception.act <- pregnancy[pregnancy$age==age.use,"Prob.of.conception"] / sex.acts.N 
  p.delivery <- pregnancy[pregnancy$age==age.use,"Prob.of.delivery"]
  
  
  alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late                      
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  
  alpha.overall <- alpha
  h.tx.MTCT = 1
  
  # scenario 1: HIV- wife is succesful in delivering baby
  neg_success <- p.delivery*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in delivering baby
  neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-alpha.overall)^N
  # scenario 4: HIV- wife is not succesful in conception
  neg_notsuccess_conception <- ((1-p.conception.act)^N)*(1-alpha.overall)^N 
  # scenario 2: HIV+ wife is succesful in delivering HIV- baby
  pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 3: HIV+ wife is succesful in delivering HIV+ baby
  pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)  
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^N)*(1-(1-alpha.overall)^N)
  # scenario 5: HIV+ wife is not succesful
  pos_notscuccess_conception <- (1-p.conception.act)^N*(1-(1-alpha.overall)^N)
  
  outcome.temp2 = c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception)
  outcome.temp3 = c(outcome.temp,outcome.temp2)
  outcome = rbind(outcome,outcome.temp3)
}
age.outcome = cbind(tx.use,outcome)
colnames(age.outcome) = c("late.val","neg_succ.all","pos_succ_neg.all","pos_succ_pos.all","neg_unsucc.all","pos_unsucc.all","neg_succ.haart","pos_succ_neg.haart","pos_succ_pos.haart","neg_unsucc.haart","pos_unsucc.haart")
age.outcome = age.outcome[order(age.outcome[,1]),]

pdf("AgeVary2_07Sept2012.pdf",height=6,width=6)
par (fig=c(0,1,0,1), omi=c(.1,0.1,0,0), mai=c(0.5,0.5,0.2,0.1),cex.lab=1.2)
layout(matrix(c(1,2), 2, 1, byrow = TRUE),heights=c(1,.2))

matplot(age.outcome[,1],age.outcome[,2:6],xlab="",ylab="",cex.lab=1.5,type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "age", line = 2.15, cex=1)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=1)

# setup for no margins on the legend
par(mai=c(0, .2, 0, .05))
# c(bottom, left, top, right)
plot.new()
legend("center",c("Female HIV-,Child HIV-","Female HIV+,Child HIV-","Female HIV+,Child HIV+","Female HIV-,No Child","Female HIV+,No Child"),lty=c(1,2,3,3,1),lwd=c(2.5,2.5,2.5,1.5,1.5), col=c("dark grey","grey","dark grey","black","black"),ncol=3,cex=0.75,pt.cex=1,seg.len=3)
dev.off()

#########
# Plots #
#########
pdf("ModeVaryParam_07Sept2012.pdf",height=6,width=6)
par (fig=c(0,1,0,1), omi=c(.1,0.1,0,0), mai=c(0.5,0.5,0.2,0.1),cex.lab=1.2)
layout(matrix(c(1,2,3,4,5,6,7,7,7), 3, 3, byrow = TRUE),heights=c(1,1,.2))

matplot(tx.outcome[,1],tx.outcome[,2:6],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.Tx", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(p.outcome[,1],p.outcome[,2:6],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.PrEP", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(sti.outcome[,1],sti.outcome[,2:6],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.other.STIs", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(late.outcome[,1],late.outcome[,2:6],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.late", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(haart.outcome[,1],haart.outcome[,2:6],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.tx.MTCT", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(alpha.outcome[,1],alpha.outcome[,2:6],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "alpha", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

# setup for no margins on the legend
par(mai=c(0, .2, 0, .05))
# c(bottom, left, top, right)
plot.new()
legend("center",c("Female HIV-,Child HIV-","Female HIV+,Child HIV-","Female HIV+,Child HIV+","Female HIV-,No Child","Female HIV+,No Child"),lty=c(1,2,3,3,1),lwd=c(2.5,2.5,2.5,1.5,1.5), col=c("dark grey","grey","dark grey","black","black"),ncol=3,cex=1,pt.cex=1.5,seg.len=3.5)
dev.off()

#only varying var, others off
pdf("OnlyVaryParam_Age40_06Sept2012.pdf",height=6,width=6)
par (fig=c(0,1,0,1), omi=c(.1,0.1,0,0), mai=c(0.5,0.5,0.2,0.1),cex.lab=1.2)
layout(matrix(c(1,2,3,4,5,6,7,7,7), 3, 3, byrow = TRUE),heights=c(1,1,.2))

matplot(tx.outcome[,1],tx.outcome[,7:11],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.Tx", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(p.outcome[,1],p.outcome[,7:11],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.PrEP", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(sti.outcome[,1],sti.outcome[,7:11],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.other.STIs", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(late.outcome[,1],late.outcome[,7:11],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.late", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(haart.outcome[,1],haart.outcome[,7:11],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "h.tx.MTCT", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

matplot(alpha.outcome[,1],alpha.outcome[,7:11],xlab="",ylab="",type=rep("l",5),col=c("dark grey","grey","dark grey","black","black"),lwd=c(2.5,2.5,2.5,1.5,1.5),lty=c(1,2,3,3,1))
mtext(side = 1, text = "alpha", line = 2.15, cex=0.7)
mtext(side = 2,text = 'Outcome Probability', line=2.15,cex=0.7)

# setup for no margins on the legend
par(mai=c(0, .2, 0, .05))
# c(bottom, left, top, right)
plot.new()
legend("center",c("Female HIV-,Child HIV-","Female HIV+,Child HIV-","Female HIV+,Child HIV+","Female HIV-,No Child","Female HIV+,No Child"),lty=c(1,2,3,3,1),lwd=c(2.5,2.5,2.5,1.5,1.5), col=c("dark grey","grey","dark grey","black","black"),ncol=3,cex=1,pt.cex=1.5,seg.len=3.5)
dev.off()