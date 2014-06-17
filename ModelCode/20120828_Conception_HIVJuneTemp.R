########################################
###                                  ###
###           MODEL SETUP            ###
###                                  ###
########################################  


rm(list=ls())
if(.Platform$OS.type == "windows") memory.limit(4000) 
options(width=135)

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

### Number of LH samples
n.samp <-  10^4
outcome.case <- "suboptimal" #"suboptimal" #"optimal" #"less"
Tx.Overwhelming.Late.Assumption <- FALSE

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
if (outcome.case == "optimal") {
  parameters["N","mode"] = 6
  parameters["N","lower"] = 3
  parameters["N", "upper"] = 12
} else { #less optimal or suboptimal
  parameters["N","mode"] = 1
  parameters["N","lower"] = 60
  parameters["N", "upper"] = 30 
}

lh.samp <- LH.Sampling(n.samp,parameters,constraints=NULL)
lh.samp$N <- as.integer(lh.samp$N)
lh.samp$age <- as.integer(lh.samp$age)

#binary options
statuses <- integer.base.b(0:(2^5-1),2)                 
colnames(statuses) <- c("Late","Tx","o.STDs","PrEP","Tx.Preg")




full.lh.samp <- data.frame(NULL)
for(n in 91183:nrow(statuses)){
dummy <- lh.samp
if(statuses[n,"Late"]==0)    dummy[,"h.late"] <- 1
if(statuses[n,"Tx"]==0)      dummy[,"h.tx"] <- 1     
if(Tx.Overwhelming.Late.Assumption & statuses[n,"Tx"]==1) dummy[,"h.late"] <- 1
if(statuses[n,"o.STDs"]==0)  dummy[,"h.other.STIs"] <- 1
if(statuses[n,"PrEP"]==0)    dummy[,"h.PrEP"] <- 1
if(statuses[n,"Tx.Preg"]==0) dummy[,"h.tx.MTCT"] <- 1

dummy.state <- t(matrix(statuses[n,],nrow=5 ,ncol= nrow(lh.samp)))
colnames(dummy.state) <- colnames(statuses)

dummy <- cbind(dummy,dummy.state)  

full.lh.samp <- rbind(full.lh.samp,dummy)

}
full.lh.samp$sex.acts.N = full.lh.samp$N
if (outcome.case == "suboptimal") full.lh.samp$sex.acts.N = round((full.lh.samp$N/30)*3)
if (outcome.case == "less") full.lh.samp$sex.acts.N = round((full.lh.samp$N/30)*6)

#number of sex act options
#lh.temp = cbind(full.lh.samp,rep(sex.acts[1],dim(full.lh.samp)[1]))
#colnames(lh.temp)[15] <- "sex.acts.N"
#if (length(sex.acts) > 1) {
#  for (i in sex.acts[2:length(sex.acts)]) {
#    sex.acts.N = rep(i,dim(full.lh.samp)[1])
#    lh.temp2 = cbind(full.lh.samp,sex.acts.N)
#    lh.temp = rbind(lh.temp,lh.temp2) 
#  }
#}

#lh.samp <- lh.temp
lh.samp <- full.lh.samp

########################################
###                                  ###
###   INTEGRATE MODEL OVER ALL       ###
##     LATIN HYPERCUBE SAMPLES       ###
###                                  ###
########################################  


outcome <- data.frame(NULL)
for(n in 160437:nrow(lh.samp)){

alpha     <- lh.samp[n,"alpha"]
h.tx      <- lh.samp[n,"h.tx"]
h.late    <- lh.samp[n,"h.late"]
h.PrEP    <- lh.samp[n,"h.PrEP"]     
h.other.STIs <- lh.samp[n,"h.other.STIs"]   
p.MTCT    <- lh.samp[n,"p.MTCT"]
h.tx.MTCT <- lh.samp[n,"h.tx.MTCT"]

age <- lh.samp[n,"age"]

N <- lh.samp[n,"N"]
n.sex <- lh.samp[n,"sex.acts.N"]

p.conception.act <- pregnancy[pregnancy$age==age,"Prob.of.conception"] / n.sex
p.delivery <- pregnancy[pregnancy$age==age,"Prob.of.delivery"]

alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late                      

# scenario 1: HIV- wife is succesful in delivering baby
if(outcome.case=="optimal"){
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
  
} else {
neg_success <- p.delivery*(1-(1-p.conception.act)^n.sex)*(1-alpha.overall)^N
# scenario 4: HIV- wife is not succesful in delivering baby
neg_notsuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^n.sex)*(1-alpha.overall)^N
# scenario 4: HIV- wife is not succesful in conception
neg_notsuccess_conception <- ((1-p.conception.act)^n.sex)*(1-alpha.overall)^N 
# scenario 2: HIV+ wife is succesful in delivering HIV- baby
pos_success_neg <- p.delivery*(1-h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^n.sex)*(1-(1-alpha.overall)^N)
# scenario 3: HIV+ wife is succesful in delivering HIV+ baby
pos_success_pos <- p.delivery*(h.tx.MTCT*p.MTCT)*(1-(1-p.conception.act)^n.sex)*(1-(1-alpha.overall)^N)  
# scenario 5: HIV+ wife is not succesful
pos_notscuccess_delivery <- (1-p.delivery)*(1-(1-p.conception.act)^n.sex)*(1-(1-alpha.overall)^N)
# scenario 5: HIV+ wife is not succesful
pos_notscuccess_conception <- (1-p.conception.act)^n.sex*(1-(1-alpha.overall)^N)
}
#outcome <- rbind(outcome,c(neg_success,pos_success_neg,pos_success_pos,
#neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception,alpha.overall,p.conception.act,p.delivery))
outcome.temp <- t(as.data.frame(c(neg_success,pos_success_neg,pos_success_pos,neg_notsuccess_delivery+neg_notsuccess_conception, pos_notscuccess_delivery+pos_notscuccess_conception,alpha.overall,p.conception.act,p.delivery)))
colnames(outcome.temp) <- c("Neg_Succ","Pos_Succ_Neg","Pos_Succ_Pos","Neg_NonSucc","Pos_NonSucc","alpha.overall","p.conception","p.delivery")
outcome.temp <- cbind(lh.samp[n,],outcome.temp)
outcome <- rbind(outcome,outcome.temp)
#if(n==1) write.csv(outcome,"Data/20120823_Outcomes.csv",append=FALSE,row.names=FALSE,col.names=TRUE)
#if(n!=1) write.csv(outcome,"Data/20120823_Outcomes.csv",append=TRUE,row.names=FALSE,col.names=FALSE)
}

write.csv(outcome,"Data/20130108_Outcomes_260927to320000.csv")
write.csv(outcome,"Data/20130108_Outcomes_213960to260927.csv")
write.csv(outcome,"Data/20130108_Outcomes_170324to213960.csv")
write.csv(outcome,"Data/20130108_Outcomes_129554to170324.csv")
write.csv(outcome,"Data/20130108_Outcomes_80903to129554.csv")
write.csv(outcome,"Data/20120823_Outcomes_1to245924.csv",append=FALSE,row.names=FALSE,col.names=TRUE)


write.csv(lh.samp,file=paste(data.dir,"lh.samp.outcomes",".csv",sep=""))


########################################
###                                  ###
###  DO A CART SENSITIVITY ANALYSIS  ###
###                                  ###
########################################
load("~/Documents/20110830_HIV_Wagner/R_Model/20120321_n406553_unfinished.RData")
rm(lh.samp)
lh.samp.all = read.csv("/Users/jaycocks/Documents/20110830_HIV_Wagner/R_model/Data/lh.samp.outcomes.csv")
setwd("/Users/jaycocks/Documents/20110830_HIV_Wagner/R_model/")
diff.outcomes = c("Neg_Succ","Pos_Succ_Neg","Pos_Succ_Pos","Neg_NonSucc","Pos_NonSucc")
sub.variables = c("PrEP","o.STDs")
sub.sexacts = unique(lh.samp.all$sex.acts.N)
age.groups = c(49,42,34,26,18)

pdf(paste(figure.dir,"Sensitivity_Analysis_Output2.pdf",sep=""))
for (j in sub.sexacts) {
  lh.samp.temp = lh.samp.all[which(lh.samp.all$sex.acts.N==j),]
for (a in 2:length(age.groups)) {
  lh.samp = lh.samp.temp[which(lh.samp.temp$age <= age.groups[a-1] & lh.samp.temp$age > age.groups[a]),]
for (i in diff.outcomes) {
 #use bootstrapping to generate forest of trees
  par(mfrow = c(2, 1))
  fit.temp = randomForest(lh.samp[,i] ~ alpha + h.late + h.tx + h.other.STIs + p.MTCT + h.tx.MTCT + N + age + h.PrEP + Late + Tx + o.STDs + PrEP + Tx.Preg,data = lh.samp,importance = TRUE,ntree=50)
  bp = barplot(sort(fit.temp$importance[,1], dec = TRUE), width = 1, space = NULL,
        names.arg = NULL, legend.text = NULL, beside = FALSE,
        horiz = FALSE, density = NULL, angle = 45,
        col = NULL, border = par("fg"),ylab = "importance metric",
        axes = TRUE, axisnames = TRUE,xaxt="n",cex.axis=0.8)
        title(paste("mean decrease in accuracy:",i,j,a),font.main=1)
  text(bp,.1*max(fit.temp$importance[,1])+sort(fit.temp$importance[,1], dec = TRUE),names(sort(fit.temp$importance[,1], dec = TRUE)),xpd=TRUE,cex=0.6)
  bp2 = barplot(sort(fit.temp$importance[,2], dec = TRUE), width = 1, space = NULL,
        names.arg = NULL, legend.text = NULL, beside = FALSE,
        horiz = FALSE, density = NULL, angle = 45,
        col = NULL, border = par("fg"),ylab = "importance metric",
        axes = TRUE, axisnames = TRUE,xaxt="n",cex.axis=0.8)
  title(paste("mean decrease in MSE:",i,j,a),font.main=1)
  text(bp2,.05*max(fit.temp$importance[,2])+sort(fit.temp$importance[,2], dec = TRUE),names(sort(fit.temp$importance[,1], dec = TRUE)),xpd=TRUE,cex=0.6)
  
  #use bootstrapping to generate forest of trees on reduced parameter set
  fit.temp2 = randomForest(lh.samp[,i] ~ alpha + h.late + h.tx + h.other.STIs + p.MTCT + h.tx.MTCT + h.PrEP,data = lh.samp,importance = TRUE,ntree=50)
  bp3 = barplot(sort(fit.temp2$importance[,1], dec = TRUE), width = 1, space = NULL,
        names.arg = NULL, legend.text = NULL, beside = FALSE,
        horiz = FALSE, density = NULL, angle = 45,
        col = NULL, border = par("fg"),ylab = "importance metric",
        axes = TRUE, axisnames = TRUE,xaxt="n",cex.axis=0.8)
        title(paste("mean decrease in accuracy:",i,j,a),font.main=1)
  text(bp3,.1*max(fit.temp2$importance[,1])+sort(fit.temp2$importance[,1], dec = TRUE),names(sort(fit.temp2$importance[,1], dec = TRUE)),xpd=TRUE,cex=0.6)
  bp4 = barplot(sort(fit.temp2$importance[,2], dec = TRUE), width = 1, space = NULL,
        names.arg = NULL, legend.text = NULL, beside = FALSE,
        horiz = FALSE, density = NULL, angle = 45,
        col = NULL, border = par("fg"),ylab = "importance metric",
        axes = TRUE, axisnames = TRUE,xaxt="n",cex.axis=0.8)
  title(paste("mean decrease in MSE:",i,j,a),font.main=1)
  text(bp4,.1*max(fit.temp2$importance[,2])+sort(fit.temp2$importance[,2], dec = TRUE),names(sort(fit.temp2$importance[,1], dec = TRUE)),xpd=TRUE,cex=0.6)
  
  #plot the tree using one iteration of rpart
  #par(mfrow = c(1, 1))
  #fit = rpart(lh.samp[,i] ~ alpha + h.late + h.tx + h.other.STIs + p.MTCT + h.tx.MTCT + N + age + h.PrEP + Late + Tx + o.STDs + PrEP + Tx.Preg,data = lh.samp,,method="anova")
  #pfit = prune(fit,fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
  #plot(pfit,uniform=TRUE,main=paste("Classification Tree (pruned):",i))
  #text(pfit,use.n=TRUE,all=TRUE,cex=0.8)
  
  #plot the, reduced parameter, tree using one iteration of rpart
  #par(mfrow = c(1, 1))
  #fit2 = rpart(lh.samp[,i] ~ alpha + h.late + h.tx + h.other.STIs + p.MTCT + h.tx.MTCT + h.PrEP,data = lh.samp,method="anova")
  #pfit2 = prune(fit2,fit2$cptable[which.min(fit2$cptable[,"xerror"]),"CP"])
  #plot(pfit2,uniform=TRUE,main=paste("Redcued Parameter Classification Tree (pruned):",i))
  #text(pfit2,use.n=TRUE,all=TRUE,cex=0.7)   
} }}#End of sensitivity loop
dev.off()

#To DO: create additional category
#Use classification based upon the category

#To DO: try derived variables
#alpha.overall <- alpha*h.PrEP*h.other.STIs*h.tx*h.late 
