
#optimal case
lh.samp.o = read.csv("/Users/jaycocks/Documents/20110830_HIV_Wagner/R_Model/Data/20130109_Optimal_Outcomes_1to320000.csv")
lh.samp.s = read.csv("/Users/jaycocks/Documents/20110830_HIV_Wagner/R_Model/Data/20130113_Suboptimal_Outcome_1to320000.csv")

lh.samp = lh.samp.o
alpha.new = lh.samp$alpha*lh.samp$h.late*lh.samp$h.tx*lh.samp$h.other.STIs*lh.samp$h.PrEP
pMTCT.new = lh.samp$p.MTCT*lh.samp$h.tx.MTCT
ds = cbind(lh.samp$h.late*lh.samp$Late,lh.samp$h.tx*lh.samp$Tx,lh.samp$h.other.STIs*lh.samp$o.STDs,lh.samp$h.PrEP*lh.samp$PrEP,lh.samp$Tx.Preg*lh.samp$h.tx.MTCT)
colnames(ds) = c("LATE.N", "TX.N", "STIs.N","PrEP.N","HAART.N")
ds.out = cbind(ds,lh.samp[,-1])
colnames(ds.out)[15:19] = c("LATE.I","TX.I","STIs.I","PrEP.I","HAART.I")
#The treatment overwhelming late was on so need to account for this
#ds.out[which(ds.out$LATE==1),"LATE.I"] = 0
ds.out$alpha.new = alpha.new
ds.out$pMTCT.new = pMTCT.new

ds.out.o = ds.out
ds.out.o$Succ = ds.out.o$Neg_Succ + ds.out.o$Pos_Succ_Neg + ds.out.o$Pos_Succ_Pos
ds.out.o$NonSucc = ds.out.o$Neg_NonSucc + ds.out.o$Pos_NonSucc
ds.out.o$Pos = ds.out.o$Pos_Succ_Neg + ds.out.o$Pos_Succ_Pos + ds.out.o$Pos_NonSucc
ds.out.o$Neg = ds.out.o$Neg_Succ + ds.out.o$Neg_NonSucc
##############
# Suboptimal #
##############
lh.samp = lh.samp.s
alpha.new = lh.samp$alpha*lh.samp$h.late*lh.samp$h.tx*lh.samp$h.other.STIs*lh.samp$h.PrEP
pMTCT.new = lh.samp$p.MTCT*lh.samp$h.tx.MTCT
ds = cbind(lh.samp$h.late*lh.samp$Late,lh.samp$h.tx*lh.samp$Tx,lh.samp$h.other.STIs*lh.samp$o.STDs,lh.samp$h.PrEP*lh.samp$PrEP,lh.samp$Tx.Preg*lh.samp$h.tx.MTCT)
colnames(ds) = c("LATE.N", "TX.N", "STIs.N","PrEP.N","HAART.N")
ds.out = cbind(ds,lh.samp[,-1])
colnames(ds.out)[15:19] = c("LATE.I","TX.I","STIs.I","PrEP.I","HAART.I")
#The treatment overwhelming late was on so need to account for this
#ds.out[which(ds.out$LATE==1),"LATE.I"] = 0
ds.out$alpha.new = alpha.new
ds.out$pMTCT.new = pMTCT.new

ds.out$Succ = ds.out$Neg_Succ + ds.out$Pos_Succ_Neg + ds.out$Pos_Succ_Pos
ds.out$NonSucc = ds.out$Neg_NonSucc + ds.out$Pos_NonSucc
ds.out$Pos = ds.out$Pos_Succ_Neg + ds.out$Pos_Succ_Pos + ds.out$Pos_NonSucc
ds.out$Neg = ds.out$Neg_Succ + ds.out$Neg_NonSucc

####################################
# Get Average Annual Probabilities #
####################################
ns.o = ds.out.o[,"Succ"]
ave.neg_succ1 = cbind(ns.o,(1-ns.o)*ns.o,(1-ns.o)^2*ns.o,(1-ns.o)^3*ns.o,(1-ns.o)^4*ns.o,(1-ns.o)^5*ns.o,(1-ns.o)^6*ns.o,(1-ns.o)^7*ns.o,(1-ns.o)^8*ns.o,(1-ns.o)^9*ns.o,(1-ns.o)^10*ns.o,(1-ns.o)^11*ns.o)
ave.neg_succ2 = cbind(ave.neg_succ1[,1],rowSums(ave.neg_succ1),rowSums(ave.neg_succ1[,1:11]),rowSums(ave.neg_succ1[,1:10]),rowSums(ave.neg_succ1[,1:9]),rowSums(ave.neg_succ1[,1:8]),rowSums(ave.neg_succ1[,1:7]),rowSums(ave.neg_succ1[,1:6]),rowSums(ave.neg_succ1[,1:5]),rowSums(ave.neg_succ1[,1:4]),rowSums(ave.neg_succ1[,1:3]),rowSums(ave.neg_succ1[,1:2]))
ave.neg_succ3 = cbind(ave.neg_succ1[,1],rowMeans(ave.neg_succ2))
ds.out.o$AveSucc = ave.neg_succ3[,2]
ds.out.o$AveNonSucc = 1 - ave.neg_succ3[,2]

ns.o = ds.out.o[,"Pos"]
ave.neg_succ1 = cbind(ns.o,(1-ns.o)*ns.o,(1-ns.o)^2*ns.o,(1-ns.o)^3*ns.o,(1-ns.o)^4*ns.o,(1-ns.o)^5*ns.o,(1-ns.o)^6*ns.o,(1-ns.o)^7*ns.o,(1-ns.o)^8*ns.o,(1-ns.o)^9*ns.o,(1-ns.o)^10*ns.o,(1-ns.o)^11*ns.o)
ave.neg_succ2 = cbind(ave.neg_succ1[,1],rowSums(ave.neg_succ1),rowSums(ave.neg_succ1[,1:11]),rowSums(ave.neg_succ1[,1:10]),rowSums(ave.neg_succ1[,1:9]),rowSums(ave.neg_succ1[,1:8]),rowSums(ave.neg_succ1[,1:7]),rowSums(ave.neg_succ1[,1:6]),rowSums(ave.neg_succ1[,1:5]),rowSums(ave.neg_succ1[,1:4]),rowSums(ave.neg_succ1[,1:3]),rowSums(ave.neg_succ1[,1:2]))
ave.neg_nonsucc3 = cbind(ave.neg_succ1[,1],rowMeans(ave.neg_succ2))
ds.out.o$AvePos = ave.neg_nonsucc3[,2]
ds.out.o$AveNeg = 1-ave.neg_nonsucc3[,2]

ds.out.o$SuccNeg = ds.out.o$AveSucc*ds.out.o$AveNeg
ds.out.o$NonSuccNeg = ds.out.o$AveNonSucc*ds.out.o$AveNeg
ds.out.o$NonSuccPos = ds.out.o$AveNonSucc*ds.out.o$AvePos

ns.o = ds.out[,"Succ"]
ave.neg_succ1 = cbind(ns.o,(1-ns.o)*ns.o,(1-ns.o)^2*ns.o,(1-ns.o)^3*ns.o,(1-ns.o)^4*ns.o,(1-ns.o)^5*ns.o,(1-ns.o)^6*ns.o,(1-ns.o)^7*ns.o,(1-ns.o)^8*ns.o,(1-ns.o)^9*ns.o,(1-ns.o)^10*ns.o,(1-ns.o)^11*ns.o)
ave.neg_succ2 = cbind(ave.neg_succ1[,1],rowSums(ave.neg_succ1),rowSums(ave.neg_succ1[,1:11]),rowSums(ave.neg_succ1[,1:10]),rowSums(ave.neg_succ1[,1:9]),rowSums(ave.neg_succ1[,1:8]),rowSums(ave.neg_succ1[,1:7]),rowSums(ave.neg_succ1[,1:6]),rowSums(ave.neg_succ1[,1:5]),rowSums(ave.neg_succ1[,1:4]),rowSums(ave.neg_succ1[,1:3]),rowSums(ave.neg_succ1[,1:2]))
ave.neg_succ3 = cbind(ave.neg_succ1[,1],rowMeans(ave.neg_succ2))
ds.out$AveSucc = ave.neg_succ3[,2]
ds.out$AveNonSucc = 1 - ave.neg_succ3[,2]

ns.o = ds.out[,"Pos"]
ave.neg_succ1 = cbind(ns.o,(1-ns.o)*ns.o,(1-ns.o)^2*ns.o,(1-ns.o)^3*ns.o,(1-ns.o)^4*ns.o,(1-ns.o)^5*ns.o,(1-ns.o)^6*ns.o,(1-ns.o)^7*ns.o,(1-ns.o)^8*ns.o,(1-ns.o)^9*ns.o,(1-ns.o)^10*ns.o,(1-ns.o)^11*ns.o)
ave.neg_succ2 = cbind(ave.neg_succ1[,1],rowSums(ave.neg_succ1),rowSums(ave.neg_succ1[,1:11]),rowSums(ave.neg_succ1[,1:10]),rowSums(ave.neg_succ1[,1:9]),rowSums(ave.neg_succ1[,1:8]),rowSums(ave.neg_succ1[,1:7]),rowSums(ave.neg_succ1[,1:6]),rowSums(ave.neg_succ1[,1:5]),rowSums(ave.neg_succ1[,1:4]),rowSums(ave.neg_succ1[,1:3]),rowSums(ave.neg_succ1[,1:2]))
ave.neg_nonsucc3 = cbind(ave.neg_succ1[,1],rowMeans(ave.neg_succ2))
ds.out$AvePos = ave.neg_nonsucc3[,2]
ds.out$AveNeg = 1-ave.neg_nonsucc3[,2]

ds.out$SuccNeg = ds.out$AveSucc*ds.out$AveNeg
ds.out$NonSuccNeg = ds.out$AveNonSucc*ds.out$AveNeg
ds.out$NonSuccPos = ds.out$AveNonSucc*ds.out$AvePos

##################
# random forests #
##################
rfit1.o = rpart(SuccNeg~alpha +N +TX.I+PrEP.I+LATE.I+STIs.I+age,data=ds.out.o,method="anova")
pdf("/Users/jaycocks/Documents/20110830_HIV_Wagner/Analysis/Tree_Jan2013.pdf")
pfit1.o=prune(rfit1.o,rfit1.o$cptable[8])
plot(pfit1.o,uniform=TRUE,main=paste("Classification Tree (pruned):Female HIV-,Child HIV-"))
text(pfit1.o,use.n=F,all=F,cex=0.8)
plot(as.party(pfit1.o),tp_args=list(id=FALSE))

rfit1.o2 = rpart(SuccNeg~alpha +N +TX.I+PrEP.I+age,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0),],method="anova")
pfit1.o2=prune(rfit1.o2,rfit1.o2$cptable[6])
plot(pfit1.o2,uniform=TRUE,main=paste("Classification Tree (pruned):Female HIV-,Child HIV- (no STIs or LATE"))
text(pfit1.o2,use.n=F,all=F,cex=0.8)
plot(as.party(pfit1.o2),tp_args=list(id=FALSE))

rfit1 = rpart(SuccNeg~alpha +N +TX.I+PrEP.I+LATE.I+STIs.I+age,data=ds.out,method="anova")
plotcp(rfit1)
pfit1=prune(rfit1,rfit1$cptable[8])
plot(pfit1,uniform=TRUE,main=paste("Classification Tree (pruned):Female HIV-,Child HIV- (sub)" ))
text(pfit1,use.n=F,all=F,cex=0.8)
plot(as.party(pfit1),tp_args=list(id=FALSE))

rfit2.o = rpart(NonSuccNeg~alpha +N +TX.I+PrEP.I+LATE.I+STIs.I+age,data=ds.out.o,method="anova")
plotcp(rfit2.o)
pfit2.o=prune(rfit2.o,rfit2.o$cptable[8])
plot(pfit2.o,uniform=TRUE,main=paste("Classification Tree (pruned):Female HIV-,No Child" ))
text(pfit2.o,use.n=F,all=F,cex=0.8)
plot(as.party(pfit2.o),tp_args=list(id=FALSE))

rfit2 = rpart(NonSuccNeg~alpha +N +TX.I+PrEP.I+LATE.I+STIs.I+age,data=ds.out,method="anova")
plotcp(rfit2)
pfit2=prune(rfit2,rfit2$cptable[7])
plot(pfit2,uniform=TRUE,main=paste("Classification Tree (pruned):Female HIV-,Child HIV- (sub)" ))
text(pfit2,use.n=F,all=F,cex=0.8)
plot(as.party(pfit2),tp_args=list(id=FALSE))

rfit3.o = rpart(NonSuccPos~alpha +N +TX.I+PrEP.I+LATE.I+STIs.I+age,data=ds.out.o,method="anova")
plotcp(rfit3.o)
pfit3.o=prune(rfit3.o,rfit3.o$cptable[6])
plot(pfit3.o,uniform=TRUE,main=paste("Classification Tree (pruned):Female HIV+,No Child" ))
text(pfit3.o,use.n=F,all=F,cex=0.8)
plot(as.party(pfit3.o),tp_args=list(id=FALSE))


rfit3 = rpart(NonSuccPos~alpha +N +TX.I+PrEP.I+LATE.I+STIs.I+age,data=ds.out,method="anova")
plotcp(rfit3)
pfit3=prune(rfit3,rfit3$cptable[7])
plot(pfit3,uniform=TRUE,main=paste("Classification Tree (pruned):Female HIV+,No Child (sub)" ))
text(pfit3,use.n=F,all=F,cex=0.8)
plot(as.party(pfit3),tp_args=list(id=FALSE))
dev.off()


################################################
# regressions for monthly probability: Optimal #
################################################
b1.o = lm(Neg_Succ ~ alpha + N + age + LATE.I + TX.I + STIs.I + PrEP.I,data=ds.out, family=)
summary(b1.o)
b1.o = glm(Neg_Succ ~ alpha + N + age + LATE.I + TX.I + STIs.I + PrEP.I,data=ds.out)
summary(b1.o)
b2.o = lm(Neg_NonSucc ~ alpha + N + age + LATE.I + TX.I + STIs.I + PrEP.I ,data=ds.out)
summary(b2.0)
b3.o = lm(Pos_NonSucc ~ alpha + N + age + LATE.I + TX.I + STIs.I + PrEP.I, data=ds.out)
summary(b3.o)

###################################################
# regressions for monthly probability: Suboptimal #
###################################################
b1.s = lm(Neg_Succ ~ alpha + N + age + LATE.I + TX.I + STIs.I + PrEP.I,data=ds.out)
summary(b1.s)
b1.s2 = lm(Neg_Succ ~ alpha.new + N + age + LATE.I + TX.I + STIs.I + PrEP.I,data=ds.out)
summary(b1.s2)
b2.s = lm(Neg_NonSucc ~ alpha + N + age + LATE.I + TX.I + STIs.I + PrEP.I ,data=ds.out)
summary(b2.s)
b2.s2 = lm(Neg_NonSucc ~ alpha.new + N + age + LATE.I + TX.I + STIs.I + PrEP.I ,data=ds.out)
summary(b2.s2)
b3.s = lm(Pos_NonSucc ~ alpha + N + age + LATE.I + TX.I + STIs.I + PrEP.I, data=ds.out)
summary(b3.s)

#######################################
# regressions for annual probability #
######################################
bns.o1 = lm(SuccNeg~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I,data=ds.out.o)
bns1 = lm(SuccNeg~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I,data=ds.out)
bnns.o1 = lm(NonSuccNeg~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I,data=ds.out.o)
bnns1 = lm(NonSuccNeg~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I,data=ds.out)
bpns.o1 = lm(NonSuccPos~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I,data=ds.out.o)
bpns1 = lm(NonSuccPos~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I,data=ds.out)

bns.o = lm(SuccNeg~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I + TX.I*PrEP.I,data=ds.out.o)
bns = lm(SuccNeg~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I + TX.I*PrEP.I,data=ds.out)
bnns.o = lm(NonSuccNeg~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I + TX.I*PrEP.I,data=ds.out.o)
bnns = lm(NonSuccNeg~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I + TX.I*PrEP.I,data=ds.out)
bpns.o = lm(NonSuccPos~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I + TX.I*PrEP.I,data=ds.out.o)
bpns = lm(NonSuccPos~alpha + age +N + TX.I + STIs.I + PrEP.I + LATE.I + TX.I*PrEP.I,data=ds.out)

bns.o2 = lm(SuccNeg~alpha + age +N + TX.I + PrEP.I +TX.I*PrEP.I,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0),])
bns2 = lm(SuccNeg~alpha + age +N + TX.I +  PrEP.I+TX.I*PrEP.I,data=ds.out[which(ds.out$STIs.I==0 & ds.out$LATE.I==0),])
bpns.o2 = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I + TX.I*PrEP.I,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0),])
bpns2 = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I + TX.I*PrEP.I,data=ds.out[which(ds.out$STIs.I==0 & ds.out$LATE.I==0),])

ds.out.o$age30 = 0
ds.out.o[which(ds.out.o$age<30),"age30"]=1
bns.o2b = lm(SuccNeg~alpha + age +N + TX.I + PrEP.I + PrEP.I*TX.I+age30,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0),])
bns2b = lm(SuccNeg~alpha + age +N + TX.I +  PrEP.I + TX.I*PrEP.I+age30,data=ds.out[which(ds.out$STIs.I==0 & ds.out$LATE.I==0),])
bpns.o2b = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I +TX.I*PrEP.I+age30,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0),])
bpns2b = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I + TX.I*PrEP.I+age30,data=ds.out[which(ds.out$STIs.I==0 & ds.out$LATE.I==0),])

#30 and younger
bns.o2c = lm(SuccNeg~alpha + age +N + TX.I + PrEP.I + PrEP.I*TX.I,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0 & ds.out.o$age<31),])
bns2c = lm(SuccNeg~alpha + age +N + TX.I +  PrEP.I + TX.I*PrEP.I,data=ds.out[which(ds.out$STIs.I==0 & ds.out$LATE.I==0 & ds.out.o$age<31),])
bpns.o2c = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I +TX.I*PrEP.I,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0 & ds.out.o$age<31),])
bpns2c = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I + TX.I*PrEP.I,data=ds.out[which(ds.out$STIs.I==0 & ds.out$LATE.I==0 & ds.out.o$age<31),])
#older than 30
bns.o2d = lm(SuccNeg~alpha + age +N + TX.I + PrEP.I + PrEP.I*TX.I,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0 & ds.out.o$age>30),])
bns2d = lm(SuccNeg~alpha + age +N + TX.I +  PrEP.I + TX.I*PrEP.I,data=ds.out[which(ds.out$STIs.I==0 & ds.out$LATE.I==0 & ds.out.o$age>30),])
bpns.o2d = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I +TX.I*PrEP.I,data=ds.out.o[which(ds.out.o$STIs.I==0 & ds.out.o$LATE.I==0 & ds.out.o$age>30),])
bpns2d = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I + TX.I*PrEP.I,data=ds.out[which(ds.out$STIs.I==0 & ds.out$LATE.I==0 & ds.out.o$age>30),])


bns.o3 = lm(SuccNeg~alpha + age +N + TX.I + PrEP.I + LATE.I + TX.I*PrEP.I,data=ds.out.o[which(ds.out.o$STIs.I==0),])
bns3 = lm(SuccNeg~alpha + age +N + TX.I +  PrEP.I + LATE.I,data=ds.out[which(ds.out$STIs.I==0 ),])
bpns.o3 = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I +LATE.I,data=ds.out.o[which(ds.out.o$STIs.I==0 ),])
bpns3 = lm(NonSuccPos~alpha + age +N + TX.I + PrEP.I + LATE.I,data=ds.out[which(ds.out$STIs.I==0 ),])
##########################
# Ave Probability for TX #
##########################
ds.out$outcome1 = NA
ds.out[which(ds.out$TX.I==1 & ds.out$STIs.I==0),"outcome1"] = "Treatment"
ds.out[which(ds.out$PrEP.I==1 & ds.out$STIs.I==0),"outcome1"] = "PrEP"
ds.out[which(ds.out$TX.I==1 & ds.out$PrEP.I==1 & ds.out$STIs.I==0),"outcome1"] = "Treatment+PrEP"
ds.out[which(ds.out$TX.I==0 & ds.out$PrEP.I==0 & ds.out$STIs.I==0),"outcome1"] = "No Treatment or PrEP"

ds.out.o$outcome1 = NA
ds.out.o[which(ds.out.o$TX.I==1 & ds.out.o$STIs.I==0),"outcome1"] = "Treatment"
ds.out.o[which(ds.out.o$PrEP.I==1 & ds.out.o$STIs.I==0),"outcome1"] = "PrEP"
ds.out.o[which(ds.out.o$TX.I==1 & ds.out.o$PrEP.I==1 & ds.out.o$STIs.I==0),"outcome1"] = "Treatment+PrEP"
ds.out.o[which(ds.out.o$TX.I==0 & ds.out.o$PrEP.I==0 & ds.out.o$STIs.I==0),"outcome1"] = "No Treatment or PrEP"

ds.out$outcome2 = NA
ds.out[which(ds.out$TX.I==1 ),"outcome2"] = "Treatment"
ds.out[which(ds.out$PrEP.I==1),"outcome2"] = "PrEP"
ds.out[which(ds.out$TX.I==1 & ds.out$PrEP.I==1),"outcome2"] = "Treatment+PrEP"
ds.out[which(ds.out$TX.I==0 & ds.out$PrEP.I==0),"outcome2"] = "No Treatment or PrEP"

ds.out.o$outcome2 = NA
ds.out.o[which(ds.out.o$TX.I==1 ),"outcome2"] = "Treatment"
ds.out.o[which(ds.out.o$PrEP.I==1),"outcome2"] = "PrEP"
ds.out.o[which(ds.out.o$TX.I==1 & ds.out.o$PrEP.I==1),"outcome2"] = "Treatment+PrEP"
ds.out.o[which(ds.out.o$TX.I==0 & ds.out.o$PrEP.I==0),"outcome2"] = "No Treatment or PrEP"

#test for significance
w = ds.out.o[which(ds.out.o$outcome1=="No Treatment or PrEP"),"SuccNeg"]
x = ds.out.o[which(ds.out.o$outcome1=="Treatment"),"SuccNeg"]
y = ds.out.o[which(ds.out.o$outcome1=="Treatment+PrEP"),"SuccNeg"]
z = ds.out.o[which(ds.out.o$outcome1=="PrEP"),"SuccNeg"]
t.test(w,z)
t.test(w,x)
t.test(w,y)
t.test(x,y )
t.test(x,z)

ww = ds.out[which(ds.out$outcome1=="No Treatment or PrEP"),"SuccNeg"]
xx = ds.out[which(ds.out$outcome1=="Treatment"),"SuccNeg"]
yy = ds.out[which(ds.out$outcome1=="Treatment+PrEP"),"SuccNeg"]
zz = ds.out[which(ds.out$outcome1=="PrEP"),"SuccNeg"]
t.test(ww,zz)
t.test(ww,xx)
t.test(ww,yy)
t.test(xx,yy)
t.test(xx,zz)

t.test(xx,yy )
t.test(xx,zz)

a = ds.out.o[which(ds.out.o$outcome1=="Treatment"),"NonSuccPos"]
b = ds.out.o[which(ds.out.o$outcome1=="Treatment+PrEP"),"NonSuccPos"]
c = ds.out.o[which(ds.out.o$outcome1=="PrEP"),"NonSuccPos"]
t.test(a,b)
t.test(a,c)

aa = ds.out[which(ds.out$outcome1=="Treatment"),"NonSuccPos"]
bb = ds.out[which(ds.out$outcome1=="Treatment+PrEP"),"NonSuccPos"]
cc = ds.out[which(ds.out$outcome1=="PrEP"),"NonSuccPos"]
t.test(aa,bb)
t.test(aa,cc)

pdf("Figures/Optimal_Sub_1.pdf",height=5,width=7)
par(mfrow=c(1,1))
boxplot(Neg_Succ~outcome2,data=ds.out, main="Suboptimal: HIV-, Child", xlab="",ylab="Outcome Probability")#cex.axis=.8,cex.lab=.8,cex.main=.8)
boxplot(Neg_Succ~outcome2,data=ds.out.o, main="Optimal: HIV-, Child", xlab="",ylab="Outcome Probability")#cex.axis=.75)

boxplot(Neg_NonSucc~outcome2,data=ds.out, main="Suboptimal: HIV-, No Child", xlab="",ylab="Outcome Probability")#cex.lab=.8,cex.axis=.8,cex.main=.8)
boxplot(Neg_NonSucc~outcome2,data=ds.out.o, main="Optimal: HIV-, No Child", xlab="",ylab="Outcome Probability")#cex.axis=.6)

boxplot(Pos_NonSucc~outcome2,data=ds.out, main="Suboptimal: HIV+, No Child", xlab="",ylab="Outcome Probability")#cex.axis=.9)
boxplot(Pos_NonSucc~outcome2,data=ds.out.o, main="Optimal: HIV+, No Child", xlab="",ylab="Outcome Probability")#cex.axis=.85)
dev.off()
pdf("/Users/jaycocks/Documents/20110830_HIV_Wagner/R_Model/Optimal_Sub_1_annual.pdf",height=5,width=7)
par(mfrow=c(1,1))
boxplot(SuccNeg~outcome2,data=ds.out, main="Suboptimal: HIV-, Child", xlab="",ylab="Outcome Probability")#cex.axis=.8,cex.lab=.8,cex.main=.8)
boxplot(SuccNeg~outcome2,data=ds.out.o, main="Optimal: HIV-, Child", xlab="",ylab="Outcome Probability")#cex.axis=.75)

boxplot(NonSuccNeg~outcome2,data=ds.out, main="Suboptimal: HIV-, No Child", xlab="",ylab="Outcome Probability")#cex.lab=.8,cex.axis=.8,cex.main=.8)
boxplot(NonSuccNeg~outcome2,data=ds.out.o, main="Optimal: HIV-, No Child", xlab="",ylab="Outcome Probability")#cex.axis=.6)

boxplot(NonSuccPos~outcome2,data=ds.out, main="Suboptimal: HIV+, No Child", xlab="",ylab="Outcome Probability")#cex.axis=.9)
boxplot(NonSuccPos~outcome2,data=ds.out.o, main="Optimal: HIV+, No Child", xlab="",ylab="Outcome Probability")#cex.axis=.85)
dev.off()

#create a data frame with averages and standard deviations

ds.out$run.type = "suboptimal" 
ds.out.o$run.type= "optimal"
ds.all = rbind(ds.out,ds.out.o)
sn.avg<-ddply(ds.all, c("outcome1", "run.type"), function(df)
  return(c(sn.avg=mean(df$SuccNeg), sn.sd=sd(df$SuccNeg))))
nsn.avg <-ddply(ds.all, c("outcome1", "run.type"), function(df)
  return(c(nsn.avg=mean(df$NonSuccNeg), nsn.sd=sd(df$NonSuccNeg)))) 
nsp.avg <-ddply(ds.all, c("outcome1", "run.type"), function(df)
  return(c(nsp.avg=mean(df$NonSuccPos), nsp.sd=sd(df$NonSuccPos))))
sn.avg$diff <- NA
sn.avg[3,"diff"] = sn.avg[3,"sn.avg"]- sn.avg[1,"sn.avg"]
sn.avg[5,"diff"] = sn.avg[5,"sn.avg"]- sn.avg[3,"sn.avg"] 
sn.avg[7,"diff"] = sn.avg[7,"sn.avg"]- sn.avg[5,"sn.avg"] 

sn.avg[4,"diff"] = sn.avg[4,"sn.avg"]- sn.avg[2,"sn.avg"]
sn.avg[6,"diff"] = sn.avg[6,"sn.avg"]- sn.avg[4,"sn.avg"] 
sn.avg[8,"diff"] = sn.avg[8,"sn.avg"]- sn.avg[6,"sn.avg"] 
sn.avg <- sn.avg[1:8,]

nsp.avg$diff <- NA
nsp.avg[3,"diff"] = nsp.avg[1,"nsp.avg"]- nsp.avg[3,"nsp.avg"]
nsp.avg[5,"diff"] = nsp.avg[3,"nsp.avg"]- nsp.avg[5,"nsp.avg"] 
nsp.avg[7,"diff"] = nsp.avg[5,"nsp.avg"]- nsp.avg[7,"nsp.avg"] 

nsp.avg[4,"diff"] = nsp.avg[2,"nsp.avg"]- nsp.avg[4,"nsp.avg"]
nsp.avg[6,"diff"] = nsp.avg[4,"nsp.avg"]- nsp.avg[6,"nsp.avg"] 
nsp.avg[8,"diff"] = nsp.avg[6,"nsp.avg"]- nsp.avg[8,"nsp.avg"] 
nsp.avg <- nsp.avg[1:8,]

sn.avg2<-ddply(ds.all, c("outcome2", "run.type"), function(df)
  return(c(sn.avg2=mean(df$SuccNeg), sn.sd2=sd(df$SuccNeg))))
nsn.avg2 <-ddply(ds.all, c("outcome2", "run.type"), function(df)
  return(c(nsn.avg2=mean(df$NonSuccNeg), nsn.sd2=sd(df$NonSuccNeg)))) 
nsp.avg2 <-ddply(ds.all, c("outcome2", "run.type"), function(df)
  return(c(nsp.avg2=mean(df$NonSuccPos), nsp.sd2=sd(df$NonSuccPos))))

sn.avg2$diff <- NA
sn.avg2[3,"diff"] = sn.avg2[3,"sn.avg2"]- sn.avg2[1,"sn.avg2"]
sn.avg2[5,"diff"] = sn.avg2[5,"sn.avg2"]- sn.avg2[3,"sn.avg2"] 
sn.avg2[7,"diff"] = sn.avg2[7,"sn.avg2"]- sn.avg2[5,"sn.avg2"] 

sn.avg2[4,"diff"] = sn.avg2[4,"sn.avg2"]- sn.avg2[2,"sn.avg2"]
sn.avg2[6,"diff"] = sn.avg2[6,"sn.avg2"]- sn.avg2[4,"sn.avg2"] 
sn.avg2[8,"diff"] = sn.avg2[8,"sn.avg2"]- sn.avg2[6,"sn.avg2"] 
sn.avg2 <- sn.avg2[1:8,]

nsp.avg2$diff <- NA
nsp.avg2[3,"diff"] = nsp.avg2[1,"nsp.avg2"]- nsp.avg2[3,"nsp.avg2"]
nsp.avg2[5,"diff"] = nsp.avg2[3,"nsp.avg2"]- nsp.avg2[5,"nsp.avg2"] 
nsp.avg2[7,"diff"] = nsp.avg2[5,"nsp.avg2"]- nsp.avg2[7,"nsp.avg2"] 

nsp.avg2[4,"diff"] = nsp.avg2[2,"nsp.avg2"]- nsp.avg2[4,"nsp.avg2"]
nsp.avg2[6,"diff"] = nsp.avg2[4,"nsp.avg2"]- nsp.avg2[6,"nsp.avg2"] 
nsp.avg2[8,"diff"] = nsp.avg2[6,"nsp.avg2"]- nsp.avg2[8,"nsp.avg2"] 
nsp.avg2 <- nsp.avg2[1:8,]

#take diffs from baseline
diff.all.a <- cbind(sn.avg[3:8,],rep("Female HIV- Child",6))
diff.all.b <- cbind(nsp.avg[3:8,], rep("Female HIV+ No Child",6))
colnames(diff.all.b) <- colnames(diff.all.a)
diff.all <- rbind(diff.all.a,diff.all.b)
colnames(diff.all)[ncol(diff.all)] <- "outcome"
diff.all$run.type2 <- paste(diff.all$outcome,", All Ages, No STIs", sep="")
diff.all <- diff.all[,c(1,2,5,7)]
diff.all$run.type2 <- factor(diff.all$run.type2, levels = c("Female HIV- Child, All Ages, No STIs","Female HIV+ No Child, All Ages, No STIs" ))
diff.all$outcome1 <- factor(diff.all$outcome1, levels= c("PrEP","Treatment", "Treatment+PrEP"), labels = c("PrEP", "ART", "ART+PrEP"))

pdf("/Users/jaycocks/Documents/20110830_HIV_Wagner/R_Model/Figures/TreatPrEP_b_22Jan2103.pdf",height=5,width=8)

ggplot(diff.all, aes(x = factor(run.type), y=diff, fill = outcome1))+
  geom_bar() + facet_grid(.~run.type2) +
  ylab('average percent change in outcome') + xlab('') + opts(title = 'ART and PrEP Impact on Annual Outcome Probabilities') + scale_fill_discrete('') +
  theme_bw() +
  opts(axis.text.x=theme_text(vjust = 1), legend.position="top", legend.direction="horizontal") +
  scale_fill_manual(values = c("PrEP" = "lightgrey", "ART" = "darkgrey","ART+PrEP"="black"),name="")

diff.all.a <- cbind(sn.avg2[3:8,],rep("Female HIV- Child",6))
diff.all.b <- cbind(nsp.avg2[3:8,], rep("Female HIV+ No Child",6))
colnames(diff.all.b) <- colnames(diff.all.a)
diff.all2 <- rbind(diff.all.a,diff.all.b)
colnames(diff.all2)[ncol(diff.all2)] <- "outcome"
diff.all2$run.type2 <- paste(diff.all2$outcome,", All Ages", sep="")
diff.all2 <- diff.all2[,c(1,2,5,7)]
diff.all2$run.type2 <- factor(diff.all2$run.type2, levels = c("Female HIV- Child, All Ages","Female HIV+ No Child, All Ages" ))
diff.all2$outcome2 <- factor(diff.all2$outcome2, levels= c("PrEP","Treatment", "Treatment+PrEP"), labels = c("PrEP", "ART", "ART+PrEP"))
ggplot(diff.all2, aes(x = factor(run.type), y=diff, fill = outcome2))+
  geom_bar() + facet_grid(.~run.type2) +
  ylab('average percent change in outcome') + xlab('') + opts(title = 'ART and PrEP Impact on Annual Outcome Probabilities') + scale_fill_discrete('') +
  theme_bw() +
  opts(axis.text.x=theme_text(vjust = 1), legend.position="top", legend.direction="horizontal") +
  scale_fill_manual(values = c("PrEP" = "lightgrey", "ART" = "darkgrey","ART+PrEP"="black"),name="") +
  scale_y_continuous(breaks = round(seq(0,.5, by = 0.05),2))
dev.off()

#under and over 30
sn.avg.30<-ddply(ds.all[which(ds.all$age<=30),], c("outcome1", "run.type"), function(df)
  return(c(sn.avg=mean(df$SuccNeg), sn.sd=sd(df$SuccNeg), sn.med=median(ds$SuccNeg), sn.low=quantile(ds$SuccNeg,.25), sn.high=quantil(ds$SuccNeg,.75))))
nsp.avg.30 <-ddply(ds.all[which(ds.all$age<=30),], c("outcome1", "run.type"), function(df)
  return(c(nsp.avg=mean(df$NonSuccPos), nsp.sd=sd(df$NonSuccPos))))
sn.avg.31<-ddply(ds.all[which(ds.all$age>30),], c("outcome1", "run.type"), function(df)
  return(c(sn.avg=mean(df$SuccNeg), sn.sd=sd(df$SuccNeg))))
nsp.avg.31 <-ddply(ds.all[which(ds.all$age>30),], c("outcome1", "run.type"), function(df)
  return(c(nsp.avg=mean(df$NonSuccPos), nsp.sd=sd(df$NonSuccPos))))

pdf("/Users/jaycocks/Documents/20110830_HIV_Wagner/R_Model/Figures/TreatPrEP_AllAges_23Jan2103.pdf",height=5,width=7)
par(mfrow=c(1,1))

#create the barplot component
avg.plot<-qplot(outcome1, sn.avg, fill=factor(run.type), data=sn.avg[which(!is.na(sn.avg[,1])),], geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=sn.avg+sn.sd, ymin=sn.avg-sn.sd+sn.sd*.1), position=dodge)+
  theme_bw()+
   opts(title = "Average Probability of Female Remaining HIV-1-uninfected with Child",legend.position="top", legend.direction="horizontal") +
  labs(x = NULL, y = "Annual Probability") +
  #geom_text(aes(label = paste(sprintf("%.1f", sn.avg*100), "%", sep=""),
  #              y = sn.avg+0.015, x=outcome1),
  #          size = 4, position = position_dodge(width=1.34))+
  scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")+
scale_y_continuous(formatter = "percent",breaks = round(seq(0,.6, by = 0.05),1))

avg.plot<-qplot(outcome1, sn.avg, fill=factor(run.type), data=sn.avg[which(!is.na(sn.avg[,1])),], geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=sn.avg+sn.sd, ymin=sn.avg-sn.sd+sn.sd*.1), position=dodge)+
  theme_bw()+
  opts(title = "Average Probability of Female Remaining HIV-1-uninfected with Child",legend.position="top", legend.direction="horizontal") +
  labs(x = NULL, y = "Annual Probability") +
  geom_text(aes(label = paste(sprintf("%.1f", sn.avg*100), "%", sep=""),
                y = sn.avg+0.015, x=outcome1),
            size = 3, position = position_dodge(width=1.5))+
  scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")+
  scale_y_continuous(formatter = "percent",breaks = round(seq(0,.6, by = 0.05),1))

avg.plot<-qplot(outcome1, sn.avg, fill=factor(run.type), data=sn.avg[which(!is.na(sn.avg[,1])),], geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=sn.avg+sn.sd, ymin=sn.avg-sn.sd+sn.sd*.1), position=dodge)+
  theme_bw()+
  opts(title = "Average Probability of Female Remaining HIV-1-uninfected with Child",legend.position="top", legend.direction="horizontal") +
  labs(x = NULL, y = "Annual Probability") +
  geom_text(aes(label = paste( signif(sn.avg*100,0), "%", sep=""),
                y = sn.avg+0.015, x=outcome1),
            size = 3, position = position_dodge(width=1.34))+
              scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")+
              scale_y_continuous(formatter = "percent",breaks = round(seq(0,.6, by = 0.05),1))

dev.off()

pdf("/Users/jaycocks/Documents/20110830_HIV_Wagner/R_Model/Figures/TreatPrEP_22Jan2103.pdf",height=5,width=7)
par(mfrow=c(1,1))

#create the barplot component
avg.plot<-qplot(outcome1, sn.avg, fill=factor(run.type), data=sn.avg[which(!is.na(sn.avg[,1])),], geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=sn.avg+sn.sd, ymin=sn.avg-sn.sd+sn.sd*.1), position=dodge)+theme_bw()+
  opts(title = "Outcome (no STIs): Female HIV- with Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")
                                                                                                                                     
avg.plot<-qplot(outcome1, nsn.avg, fill=factor(run.type), data=nsn.avg[which(!is.na(nsn.avg[,1])),], geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=nsn.avg+nsn.sd, ymin=nsn.avg-nsn.sd+nsn.sd*.1), position=dodge)+theme_bw()+opts(title = "Outcome (no STIs): Female HIV- without Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")

avg.plot<-qplot(outcome1, nsp.avg, fill=factor(run.type), data=nsp.avg[which(!is.na(nsp.avg[,1])),], geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=nsp.avg+nsp.sd, ymin=nsp.avg-nsp.sd+nsp.sd*.1), position=dodge)+theme_bw()+opts(title = "Outcome (no STIs): Female HIV+ without Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")



######### with STIs
avg.plot<-qplot(outcome2, sn.avg2, fill=factor(run.type), data=sn.avg2, geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=sn.avg2+sn.sd2, ymin=sn.avg2-sn.sd2+sn.sd2*.1), position=dodge)+theme_bw()+opts(title = "Outcome: Female HIV- with Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")

avg.plot<-qplot(outcome2, nsn.avg2, fill=factor(run.type), data=nsn.avg2, geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=nsn.avg2+nsn.sd2, ymin=nsn.avg2-nsn.sd2+nsn.sd2*.1), position=dodge)+theme_bw()+opts(title = "Outcome: Female HIV- without Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")

avg.plot<-qplot(outcome2, nsp.avg2, fill=factor(run.type), data=nsp.avg2, geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=nsp.avg2+nsp.sd2, ymin=nsp.avg2-nsp.sd2+nsp.sd2*.1), position=dodge)+theme_bw()+opts(title = "Outcome: Female HIV+ without Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")

dev.off()

ds.out$outcome3 = NA
ds.out[which(ds.out$STIs.I==1 & ds.out$TX.I==0 & ds.out$PrEP.I==0),"outcome3"] = "STIs"
ds.out[which(ds.out$TX.I==0 & ds.out$PrEP.I==0 & ds.out$STIs.I == 0),"outcome3"] = "No STIs"
ds.out[which(ds.out$TX.I==1 & ds.out$PrEP.I==0 & ds.out$STIs.I == 1),"outcome3"] = "STIs+Treatment"
ds.out[which(ds.out$TX.I==0 & ds.out$PrEP.I==1 & ds.out$STIs.I == 1),"outcome3"] = "STIs+PrEP"

ds.out.o$outcome3 = NA
ds.out.o[which(ds.out.o$STIs.I==1 & ds.out.o$TX.I==0 & ds.out.o$PrEP.I==0),"outcome3"] = "STIs"
ds.out.o[which(ds.out.o$TX.I==0 & ds.out.o$PrEP.I==0 & ds.out.o$STIs.I == 0),"outcome3"] = "No STIs"
ds.out.o[which(ds.out.o$TX.I==1 & ds.out.o$PrEP.I==0 & ds.out.o$STIs.I == 1),"outcome3"] = "STIs+Treatment"
ds.out.o[which(ds.out.o$TX.I==0 & ds.out.o$PrEP.I==1 & ds.out.o$STIs.I == 1),"outcome3"] = "STIs+PrEP"

sn.avg3 <-ddply(ds.all[which(!is.na(ds.all$outcome3)),], c("outcome3", "run.type"), function(df)
  return(c(sn.avg3=mean(df$SuccNeg), sn.sd3=sd(df$SuccNeg))))
nsn.avg3 <-ddply(ds.all[which(!is.na(ds.all$outcome3)),], c("outcome3", "run.type"), function(df)
  return(c(nsn.avg3=mean(df$NonSuccNeg), nsn.sd3=sd(df$NonSuccNeg)))) 
nsp.avg3 <-ddply(ds.all[which(!is.na(ds.all$outcome3)),], c("outcome3", "run.type"), function(df)
  return(c(nsp.avg3=mean(df$NonSuccPos), nsp.sd3=sd(df$NonSuccPos))))

avg.plot<-qplot(outcome3, sn.avg3, fill=factor(run.type), data=sn.avg3, geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=sn.avg3+sn.sd3, ymin=sn.avg3-sn.sd3+sn.sd3*.1), position=dodge)+theme_bw()+opts(title = "Outcome: Female HIV- with Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")

avg.plot<-qplot(outcome3, nsn.avg3, fill=factor(run.type), data=nsn.avg3, geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=nsn.avg3+nsn.sd3, ymin=nsn.avg3-nsn.sd3+nsn.sd3*.1), position=dodge)+theme_bw()+opts(title = "Outcome: Female HIV- without Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")


avg.plot<-qplot(outcome3, nsp.avg3, fill=factor(run.type), data=nsp.avg3, geom="bar", position="dodge")
#first, define the width of the dodge
dodge <- position_dodge(width=0.9) 
#now add the error bars to the plot
avg.plot+geom_linerange(aes(ymax=nsp.avg3+nsp.sd3, ymin=nsp.avg3-nsp.sd3+nsp.sd3*.1), position=dodge)+theme_bw()+opts(title = "Outcome: Female HIV+ without Child",legend.position="top", legend.direction="horizontal") +labs(x = NULL, y = "Annual Probability") + scale_fill_manual(values = c("optimal" = "darkgrey", "suboptimal" = "lightgrey"),name="")


dev.off()
#p1 <- ggplot(ds.out[which(!is.na(ds.out$outcome2)),], aes(x=outcome2,y=Neg_Succ)) + geom_boxplot() +
#  stat_summary(fun.y=mean, geom="point", shape=5, size=4) + opts(title="Suboptimal: Female HIV-, Child")
#p1 <- p1 + scale_x_discrete(name="") + scale_y_continuous(name="P(outcome)")
#p2 <- ggplot(ds.out, aes(x=outcome2,y=Neg_NonSucc)) + geom_boxplot() +
#  stat_summary(fun.y=mean, geom="point", shape=5, size=4) + opts(title="Suboptimal: Female HIV-, No Child")
#p2 <- p2 + scale_x_discrete(name="") + scale_y_continuous(name="P(outcome)")
#p3 <- ggplot(ds.out, aes(x=outcome2,y=Pos_NonSucc)) + geom_boxplot() +
#  stat_summary(fun.y=mean, geom="point", shape=5, size=4) + opts(title="Suboptimal: Female HIV+, No Child")
#p3 <- p3 + scale_x_discrete(name="") + scale_y_continuous(name="P(outcome)")

ds.out$outcome3 = NA
ds.out[which(ds.out$STIs.I==1 & ds.out$TX.I==0 & ds.out$PrEP.I==0),"outcome3"] = "STIs"
ds.out[which(ds.out$TX.I==0 & ds.out$PrEP.I==0 & ds.out$STIs.I == 0),"outcome3"] = "No STIs"
ds.out[which(ds.out2$TX.I==1 & ds.out2$PrEP.I==0 & ds.out2$STIs.I == 1),"outcome3"] = "STIs+Treatment"
ds.out[which(ds.out2$TX.I==0 & ds.out2$PrEP.I==1 & ds.out2$STIs.I == 1),"outcome3"] = "STIs+PrEP"

pdf("Figures/Suboptimal_2.pdf",height=4,width=7)
par(mfrow=c(1,3))
boxplot(Neg_Succ~outcome3,data=ds.out, main="Suboptimal: HIV-, Child", xlab="",ylab="Outcome Probability")
boxplot(Neg_NonSucc~outcome3,data=ds.out, main="Suboptimal: HIV-, No Child", xlab="",ylab="Outcome Probability")
boxplot(Pos_NonSucc~outcome3,data=ds.out, main="Suboptimal: HIV+, No Child", xlab="",ylab="Outcome Probability")
dev.off()

ds.out.o$outcome3 = NA
ds.out.o[which(ds.out.o$STIs.I==1 & ds.out.o$TX.I==0 & ds.out.o$PrEP.I==0),"outcome3"] = "STIs"
ds.out.o[which(ds.out.o$TX.I==1 & ds.out.o$PrEP.I==1 & ds.out.o$STIs.I == 0),"outcome3"] = "No STIs"
#ds.out[which(ds.out2$TX.I==1 & ds.out2$PrEP.I==1 & ds.out2$STIs.I == 1),"outcome3"] = "Treatment+PrEP+STIs"
pdf("Figures/Optimal_Sub_2.pdf",height=4,width=7)
par(mfrow=c(1,3))
boxplot(Neg_Succ~outcome3,data=ds.out, main="Suboptimal: HIV-, Child", xlab="",ylab="Outcome Probability")
boxplot(Neg_NonSucc~outcome3,data=ds.out, main="Suboptimal: HIV-, No Child", xlab="",ylab="Outcome Probability")
boxplot(Pos_NonSucc~outcome3,data=ds.out, main="Suboptimal: HIV+, No Child", xlab="",ylab="Outcome Probability")
boxplot(Neg_Succ~outcome3,data=ds.out.o, main="Optimal: HIV-, Child", xlab="",ylab="Outcome Probability")
boxplot(Neg_NonSucc~outcome3,data=ds.out.o, main="Optimal: HIV-, No Child", xlab="",ylab="Outcome Probability")
boxplot(Pos_NonSucc~outcome3,data=ds.out.o, main="Optimal: HIV+, No Child", xlab="",ylab="Outcome Probability")
dev.off()


#####################
# Coefficient Plots #
#####################
coeff <- cbind(b1.o$coefficients,b1.s$coefficients, b2.o$coefficients, b2.s$coefficients, b3.o$coefficients, b3.s$coefficients) 
coeff <- coeff[c(6,8,7,5),]
rel.impact <- c(.33,.80,.32)
rel.impact.s <- c(.29,.75,.27)
ri <- rbind(rel.impact,rel.impact.s)
names(rel.impact) <- c("Optimal.PrEP","Suboptimal.PrEP","Optimal.STIs","Suboptimal.STIs","Optimal.Late","Suboptimal.Late")
rel.impact <- rbind(rel.imp)
pdf("Figures/coeff_opt_sub.pdf")
barplot(coeff[,c(5,3,1)], beside=TRUE, space=c(0,2), ylab="outcome probability impact",names.arg=c("Female HIV+, No Child","Female HIV-, No Child","Female HIV-, Child") ,legend.text=c("Treatment","PrEP","STIs","Late Stage"),main="Optimal Regression Coefficients",args.legend=list(c(location="topleft")))
barplot(coeff[,c(6,4,2)], beside=TRUE, space=c(0,2), ylab="outcome probability impact",names.arg=c("Female HIV+, No Child","Female HIV-, No Child","Female HIV-, Child") ,legend.text=c("Treatment","PrEP","STIs","Late Stage"),main="Suboptimal Regression Coefficients",args.legend=list(c(location="topleft")))
#barX <- barplot(ri, beside=TRUE, space=c(0,2), ylab="relative impact",names.arg=c("PrEP","STIs","Late Stage"),main="Impact of Options Relative To Treatment",legend.text=c("Optimal","Suboptimal"),args.legend=list(c(location="topleft")))
#text(cex=.5, x=barX, y=ri+par("cxy")[2]/2, round(ri,2), xpd=TRUE)
dev.off()


###############################
# Treatment and PrEp with Age #
###############################
#TX
opt.1 <- predict(b1.o, data.frame(alpha = 0.002161812, N = 6, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.1.h <- predict(b1.o, data.frame(alpha = 0.001021951, N = 3, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.1.l <- predict(b1.o, data.frame(alpha = 0.003071536, N = 12, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.2 <- predict(b2.o, data.frame(alpha = 0.002161812, N = 6, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.2.h <- predict(b2.o, data.frame(alpha = 0.001021951, N = 3, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.2.l <- predict(b2.o, data.frame(alpha = 0.003071536, N = 12, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.3 <- predict(b3.o, data.frame(alpha = 0.002161812, N = 6, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.3.l <- predict(b3.o, data.frame(alpha = 0.001021951, N = 3, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.3.h <- predict(b3.o, data.frame(alpha = 0.003071536, N = 12, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))

opt.1s <- predict(b1.s, data.frame(alpha = 0.002161812, N = 30, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.1s.h <- predict(b1.s, data.frame(alpha = 0.001021951, N = 1, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.1s.l <- predict(b1.s, data.frame(alpha = 0.003071536, N = 60, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.2s <- predict(b2.s, data.frame(alpha = 0.002161812, N = 3, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.2s.h <- predict(b2.s, data.frame(alpha = 0.001021951, N = 1, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.2s.l <- predict(b2.s, data.frame(alpha = 0.003071536, N = 60, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.3s <- predict(b3.s, data.frame(alpha = 0.002161812, N = 30, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.3s.l <- predict(b3.s, data.frame(alpha = 0.001021951, N = 1, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))
opt.3s.h <- predict(b3.s, data.frame(alpha = 0.003071536, N = 60, age=c(18:49),LATE.I=0,TX.I=1,STIs.I=0,PrEP.I=0))

#PrEP
opt.1p <- predict(b1.o, data.frame(alpha = 0.002161812, N = 6, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.1.hp <- predict(b1.o, data.frame(alpha = 0.001021951, N = 3, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.1.lp <- predict(b1.o, data.frame(alpha = 0.003071536, N = 12, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.2p <- predict(b2.o, data.frame(alpha = 0.002161812, N = 6, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.2.hp <- predict(b2.o, data.frame(alpha = 0.001021951, N = 3, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.2.lp <- predict(b2.o, data.frame(alpha = 0.003071536, N = 12, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.3p <- predict(b3.o, data.frame(alpha = 0.002161812, N = 6, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.3.lp <- predict(b3.o, data.frame(alpha = 0.001021951, N = 3, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.3.hp <- predict(b3.o, data.frame(alpha = 0.003071536, N = 12, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))

opt.1sp <- predict(b1.s, data.frame(alpha = 0.002161812, N = 30, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.1s.hp <- predict(b1.s, data.frame(alpha = 0.001021951, N = 1, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.1s.lp <- predict(b1.s, data.frame(alpha = 0.003071536, N = 60, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.2sp <- predict(b2.s, data.frame(alpha = 0.002161812, N = 30, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.2s.hp <- predict(b2.s, data.frame(alpha = 0.001021951, N = 1, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.2s.lp <- predict(b2.s, data.frame(alpha = 0.003071536, N = 60, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.3sp <- predict(b3.s, data.frame(alpha = 0.002161812, N = 30, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.3s.lp <- predict(b3.s, data.frame(alpha = 0.001021951, N = 1, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))
opt.3s.hp <- predict(b3.s, data.frame(alpha = 0.003071536, N = 60, age=c(18:49),LATE.I=0,TX.I=0,STIs.I=0,PrEP.I=1))


par (fig=c(0,1,0,1), omi=c(.1,0.1,0,0.1), mai=c(0.5,0.5,0.2,0.1),cex.lab=1.2)
layout(matrix(c(1,2,3,4,5,6,7,7,7), 3, 3, byrow = TRUE),heights=c(1,1,.1))
pdf("Figures/agetemp.pdf")
#par(mfrow=c(2,3))
#HIV- Child
plot(opt.1s,type="l",lwd=3,col="grey",main="Treatment: HIV-, Child",xlab="age",ylab="outcome probablity",ylim=c(min(opt.1.l,opt.1s.l),max(opt.1.h,opt.1s.h)))
#points(opt.1s.h, type="l", lty=2, col="grey")
#points(opt.1s.l, type="l",lty=2, col="grey")
points(opt.1,type="l",lwd=3)
#points(opt.1.h, type="l", lty=3)
#points(opt.1.l, type="l",lty=3)

#HIV- No Child
plot(opt.2s,type="l",lwd=3,col="grey",main="Treatment: HIV-, No Child",xlab="age",ylab="outcome probablity",ylim=c(min(opt.2.l,opt.2s.l),max(opt.2.h,opt.2s.h)))
#points(opt.2s.h, type="l", lty=2, col="grey")
#points(opt.2s.l, type="l",lty=2, col="grey")
points(opt.2,type="l",lwd=3)
#points(opt.2.h, type="l", lty=3)
#points(opt.2.l, type="l",lty=3)

#HIV+ No Child
plot(opt.3s,type="l",lwd=3,col="grey",main="Treatment: HIV+, No Child",xlab="age",ylab="outcome probablity",ylim=c(min(opt.3.l,opt.3s.l),max(opt.3.h,opt.3s.h)))
#points(opt.3s.h, type="l", lty=2, col="grey")
#points(opt.3s.l, type="l",lty=2, col="grey")
points(opt.3,type="l",lwd=3)
#points(opt.3.h, type="l", lty=3)
#points(opt.3.l, type="l",lty=3)
#HIV- Child
plot(opt.1sp,type="l",lwd=3,col="grey",main="PrEP: HIV-, Child",xlab="age",ylab="outcome probablity",ylim=c(min(opt.1.lp,opt.1s.lp),max(opt.1.hp,opt.1s.hp)))
#points(opt.1s.hp, type="l", lty=2, col="grey")
#points(opt.1s.lp, type="l",lty=2, col="grey")
points(opt.1p,type="l",lwd=3)
#points(opt.1.hp, type="l", lty=3)
#points(opt.1.lp, type="l",lty=3)

#HIV- No Child
plot(opt.2sp,type="l",lwd=3,col="grey",main="PrEP: HIV-, No Child",xlab="age",ylab="outcome probablity",ylim=c(min(opt.2.lp,opt.2s.lp),max(opt.2.hp,opt.2s.hp)))
#points(opt.2s.hp, type="l", lty=2, col="grey")
#points(opt.2s.lp, type="l",lty=2, col="grey")
points(opt.2p,type="l",lwd=3)
#points(opt.2.hp, type="l", lty=3)
#points(opt.2.lp, type="l",lty=3)

#HIV+ No Child
plot(opt.3sp,type="l",lwd=3,col="grey",main="PrEP: HIV+, No Child",xlab="age",ylab="outcome probablity",ylim=c(min(opt.3.lp,opt.3s.lp),max(opt.3.hp,opt.3s.hp)))
#points(opt.3s.hp, type="l", lty=2, col="grey")
#points(opt.3s.lp, type="l",lty=2, col="grey")
points(opt.3p,type="l",lwd=3)
#points(opt.3.hp, type="l", lty=3)
#points(opt.3.lp, type="l",lty=3)

par(mai=c(0, .35, 0, 0))
# c(bottom, left, top, right)
plot.new()
legend("center",c("Optimal", "Suboptimal"),lty=c(1,1), col=c("black","grey"),ncol=2,cex=1.2,pt.cex=1.2,seg.len=3,box.lwd =1)
dev.off()

########
# jsut plot
dst.s = ds.out[which(ds.out$TX.I==1),]
plot(dst.s$age,dst.s$Neg_Succ)
dst.s2 = ds.out[which(ds.out$TX.I==1 & ds.out$LATE.I==0 & ds.out$PrEP.I==0 & ds.out$STIs.I==0 & ds.out$N==30),]
dst.s2.o = ds.out.o[which(ds.out.o$TX.I==1 & ds.out.o$LATE.I==0 & ds.out.o$PrEP.I==0 & ds.out.o$STIs.I==0 & ds.out.o$N==6),]
plot(dst.s2$age,dst.s2$Neg_Succ)
plot(dst.s2.o$age,dst.s2.o$Neg_Succ,col="grey")

plot(dst.s2$age,dst.s2$Neg_NonSucc)
plot(dst.s2.o$age,dst.s2.o$Neg_NonSucc,col="grey",type="l")

plot(dst.s2$age,dst.s2$Pos_NonSucc)
plot(dst.s2.o$age,dst.s2.o$Pos_NonSucc,col="grey")

########
# PrEP #
########c





dotchart(x$mpg,labels=row.names(x),cex=.7,groups= x$cyl,
         main="Gas Milage for Car Models\ngrouped by cylinder",
         xlab="Miles Per Gallon", gcolor="black", color=x$color)
ds.out2 = ds.out
ds.t = ds.out2[which(ds.out2$TX.I==0),]
#ds.t = ds.out2[which(ds.out2$TX.I==1 & ds.out2$LATE.I==0 & ds.out2$STIs.I==0 & ds.out2$PrEP.I==0 & ds.out2$HAART.I==0),]
###New
b1t = lm(Neg_Succ ~ alpha + N + age + LATE.I +  STIs.I + PrEP.I,data=ds.t)
summary(b1t)
b2t = lm(Pos_Succ_Neg ~ alpha + p.MTCT + N + age + LATE.I +  STIs.I + PrEP.I + HAART.I ,data=ds.t)
summary(b2t)
b3t = lm(Pos_Succ_Pos ~ alpha + p.MTCT + N + age + LATE.I + STIs.I + PrEP.I + HAART.I,data=ds.t)
summary(b3t)
b4t = lm(Neg_NonSucc ~ alpha + N + age + LATE.I +  STIs.I + PrEP.I ,data=ds.t)
summary(b4t)
b5t = lm(Pos_NonSucc ~ alpha + N + age + LATE.I + STIs.I + PrEP.I, data=ds.t)
summary(b5t)

ds.out2 = ds.out
ds.t2 = ds.out2[which(ds.out2$TX.I==1),]
#ds.t = ds.out2[which(ds.out2$TX.I==1 & ds.out2$LATE.I==0 & ds.out2$STIs.I==0 & ds.out2$PrEP.I==0 & ds.out2$HAART.I==0),]
###New
b1t2 = lm(Neg_Succ ~ alpha + N + age + LATE.I +  STIs.I + PrEP.I,data=ds.t2)
summary(b1t2)
b2t2 = lm(Pos_Succ_Neg ~ alpha + p.MTCT + N + age + LATE.I +  STIs.I + PrEP.I + HAART.I ,data=ds.t2)
summary(b2t)
b3t2 = lm(Pos_Succ_Pos ~ alpha + p.MTCT + N + age + LATE.I + STIs.I + PrEP.I + HAART.I,data=ds.t2)
summary(b3t2)
b4t2 = lm(Neg_NonSucc ~ alpha + N + age + LATE.I +  STIs.I + PrEP.I ,data=ds.t2)
summary(b4t)
b5t2 = lm(Pos_NonSucc ~ alpha + N + age + LATE.I + STIs.I + PrEP.I, data=ds.t2)
summary(b5t2)
