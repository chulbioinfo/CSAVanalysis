#install.packages("dplyr")
library(dplyr)
#install.packages("corrplot")
library(corrplot)
#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")



# Random control sets #########################################################
setwd("E:\\1004_ctrl\\1004ctrl_TV_phy_20190202\\out\\wo_regional_variants\\1000random_sets_62core_sets\\1000ctrl")

# Loading data
data <- read.table("0.summary_Phy_TAAS_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_phylo_info = cbind.data.frame(ProTerBL,DistanceTN,DistanceTB,ProOriBL)
colnames(tmp_phylo_info)<-c("PTB","DTN","DTB","POB")
tmp_data_CSAV = cbind.data.frame(nTargets,ctrl,TV_A, TV_C, TV_D)
colnames(tmp_data_CSAV)<-c("nTargets","ctrl","CSAV","iCSAV","dCSAV")
detach(data)

data <- read.table("0.summary_Phy_TCC_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_data_CSCV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_CSCV)<-c("iCSCV","dCSCV","CSCV")
detach(data)

data <- read.table("0.summary_Phy_TSNV_20190321.marked_targets.info",sep='\t',header = T)
attach(data)
tmp_data_CSNV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_CSNV)<-c("iCSNV","dCSNV","CSNV")
detach(data)

data <- cbind.data.frame(tmp_data_CSAV$nTargets,
                                       tmp_data_CSAV$ctrl,
                                       tmp_phylo_info, 
                                       tmp_data_CSAV$iCSAV,
                                       tmp_data_CSCV$iCSCV,
                                       tmp_data_CSNV$iCSNV,
                                       tmp_data_CSAV$dCSAV,
                                       tmp_data_CSCV$dCSCV,
                                       tmp_data_CSNV$dCSNV,
                                       tmp_data_CSAV$CSAV,
                                       tmp_data_CSCV$CSCV,
                                       tmp_data_CSNV$CSNV
)
colnames(data)<-c("nTargets","ctrl","PTB","DTN","DTB","POB","iCSAV","iCSCV","iCSNV","dCSAV","dCSCV","dCSNV","CSAV","CSCV","CSNV")
rownames(data)<-data$nTargets
#####################################################################

Fig_1d <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  grid(lwd = 0.3)
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)
  
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line
  
  res = cor.test(tmp_ind_var, tmp_dep_var,  method = "spearman")
  print(res)
  
  xOrder  <- order(tmp_ind_var)
  x       <- tmp_ind_var[xOrder]  
  y       <- tmp_dep_var[xOrder]
  fit     <- lm(y ~ x, data=data.frame(x=x, y=y))
  newX    <- data.frame(x=jitter(x))
  #fitPred <- predict.lm(fit, newdata=newX, interval="prediction", level=0.95)
  #lines(newX$x, fitPred[,2], lty=2)
  #lines(newX$x, fitPred[,3], lty=2)
  
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[52][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[52][1])],
         col="black",pch=20)
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[353][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[353][1])],
         col="black",pch=20)
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[665][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[665][1])],
         col="black",pch=20)
  
  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)
  
  
  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = "1"
    if(tmp_p < 0.05){
      tmp_tmp_p = "0.05*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "0.01**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "0.001***"
        } 
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) < Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    "topleft",
    legend=rp, bty='n', text.col = "black")
  legend("bottomright",
         legend = c("6 vocal learners","Random control sets", "Outliers (Cook's Dist.)"),
         col=c("red","grey","black"),
         pch=c(15,20,20),
         text.col = c("red","grey","black")
  )
  
  
  plot(dep_ind_var.lm, pch = 18, col = "red", which=c(4))
  print(cooks.distance(dep_ind_var.lm))
}
par(mfrow=c(1,1),mar = c(2, 2, 2, 2))

# Figure 1C
pdf("Fig1d_iCSAV_dCSAV_randomctrl_v2020419.pdf",width=6,height=6)
layout(mat = matrix(c(1),1,1,byrow=T))

Fig_1d(data$dCSAV,data$iCSAV,"iCSAV/dCSAV of random controls (n=1,001)", "Different amino acid convergences (dCSAV)", "Identical amino acid convergences (iCSAV)")%%box()
data$nTargets[52]
data$nTargets[353]
data$nTargets[665]
dev.off()

Fig_1e <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  grid(lwd = 0.3)
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)
  
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line
  
  res = cor.test(tmp_ind_var, tmp_dep_var,  method = "spearman")
  print(res)
  
  xOrder  <- order(tmp_ind_var)
  x       <- tmp_ind_var[xOrder]  
  y       <- tmp_dep_var[xOrder]
  fit     <- lm(y ~ x, data=data.frame(x=x, y=y))
  newX    <- data.frame(x=jitter(x))
  #fitPred <- predict.lm(fit, newdata=newX, interval="prediction", level=0.95)
  #lines(newX$x, fitPred[,2], lty=2)
  #lines(newX$x, fitPred[,3], lty=2)
  
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[52][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[52][1])],
         col="black",pch=20)
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[353][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[353][1])],
         col="black",pch=20)
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[665][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[665][1])],
         col="black",pch=20)
  
  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)
  
  
  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = "1"
    if(tmp_p < 0.05){
      tmp_tmp_p = "0.05*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "0.01**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "0.001***"
        } 
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) < Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    "topleft",
    legend=rp, bty='n', text.col = "black")
  legend("bottomright",
         legend = c("6 vocal learners","Random control sets", "Outliers (Cook's Dist.)"),
         col=c("red","grey","black"),
         pch=c(15,20,20),
         text.col = c("red","grey","black")
  )
  
  # cook's distance plot
  plot(dep_ind_var.lm, pch = 18, col = "red", which=c(4))
  print(cooks.distance(dep_ind_var.lm))
  
  
}
par(mfrow=c(1,1),mar = c(2, 2, 2, 2))

pdf("Fig1e_iCSAV_POB_randomctrl_v2020419.pdf",width=6,height=6)
layout(mat = matrix(c(1),1,1,byrow=T))
Fig_1e(data$POB,data$iCSAV, "iCSAV/POB of random controls (n=1,001)", "Product of origin branch lengths (POB)", "Identical amino acid convergences (iCSAV)")%%box()
data$nTargets[52]
data$nTargets[353]
data$nTargets[665]
dev.off()

# Core control sets
############ Data loading #######################################

setwd("E:\\1004_ctrl\\1004ctrl_TV_phy_20190202\\out\\wo_regional_variants\\1000random_sets_62core_sets\\61ctrl\\")
#

data <- read.table("0.summary_Phy_TAAS_20181223.txt",sep='\t',header = T)
attach(data)
tmp_phylo_info = cbind.data.frame(ProTerBL,DistanceTN,DistanceTB,ProOriBL)
colnames(tmp_phylo_info)<-c("PTB","DTN","DTB","POB")
tmp_data_CSAV = cbind.data.frame(nTargets,ctrl,TV_A, TV_C, TV_D)
colnames(tmp_data_CSAV)<-c("nTargets","ctrl","CSAV","iCSAV","dCSAV")
detach(data)

data <- read.table("0.summary_Phy_TCC_20181223.txt",sep='\t',header = T)
attach(data)
tmp_data_CSCV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_CSCV)<-c("iCSCV","dCSCV","CSCV")
detach(data)

data <- read.table("0.summary_Phy_TSNV_20181223.txt",sep='\t',header = T)
attach(data)
tmp_data_CSNV= cbind.data.frame(TV_C,TV_D,TV_A)
colnames(tmp_data_CSNV)<-c("iCSNV","dCSNV","CSNV")
detach(data)

data <- cbind.data.frame(tmp_data_CSAV$nTargets,
                         tmp_data_CSAV$ctrl,
                         tmp_phylo_info, 
                         tmp_data_CSAV$iCSAV,
                         tmp_data_CSCV$iCSCV,
                         tmp_data_CSNV$iCSNV,
                         tmp_data_CSAV$dCSAV,
                         tmp_data_CSCV$dCSCV,
                         tmp_data_CSNV$dCSNV,
                         tmp_data_CSAV$CSAV,
                         tmp_data_CSCV$CSCV,
                         tmp_data_CSNV$CSNV
)
colnames(data)<-c("nTargets","ctrl","PTB","DTN","DTB","POB","iCSAV","iCSCV","iCSNV","dCSAV","dCSCV","dCSNV","CSAV","CSCV","CSNV")
rownames(data)<-data$nTargets
#################################################################

# Fig.1f
Fig_1f <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  grid(lwd = 0.3)
  # control sets
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)
  
  # regression analysis
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line
  
  
  res = cor.test(tmp_ind_var, tmp_dep_var,  method = "spearman")
  print(res)
  
  xOrder  <- order(tmp_ind_var)
  x       <- tmp_ind_var[xOrder]  
  y       <- tmp_dep_var[xOrder]
  fit     <- lm(y ~ x, data=data.frame(x=x, y=y))
  newX    <- data.frame(x=jitter(x))
  #fitPred <- predict.lm(fit, newdata=newX, interval="prediction", level=0.95)
  #lines(newX$x, fitPred[,2], lty=2)
  #lines(newX$x, fitPred[,3], lty=2)
  
  
  # outlier
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[12][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[12][1])],
         col="black",pch=20)
  
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[15][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[15][1])],
         col="black",pch=20)
  
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[35][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[35][1])],
         col="black",pch=20)
  
  # set of birds of prey
  points(tmp_ind_var[data$nTargets=="CARCR,CATAU,FALPE,HALAL,HALLE,LEPDI,TYTAL"],
         tmp_dep_var[data$nTargets=="CARCR,CATAU,FALPE,HALAL,HALLE,LEPDI,TYTAL"],
         col="purple",pch=17)
  
  # set of waterbirds
  points(tmp_ind_var[data$nTargets=="APTFO,BALRE,CHAVO,EGRGA,FULGL,GAVST,NIPNI,PELCR,PHACA,PHALE,PHORU,PODCR,PYGAD"],
         tmp_dep_var[data$nTargets=="APTFO,BALRE,CHAVO,EGRGA,FULGL,GAVST,NIPNI,PELCR,PHACA,PHALE,PHORU,PODCR,PYGAD"],
         col="blue",pch=16)
  
  # set of closest control set
  points(tmp_ind_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="orange",pch=18)
  
  # set of vocal learners
  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)
  
  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = "1"
    if(tmp_p < 0.05){
      tmp_tmp_p = "0.05*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "0.01**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "0.001***"
        } 
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) < Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    "bottomright",
    legend=rp, bty='n', text.col = "black")
  legend("topleft",
         legend = c("6 vocal learners","5 vocal learners+Swift","6 birds of prey","15 waterbirds","Core control sets", "Outlier (Cook's Dist.)"),
         col=c("red","orange","purple","blue","grey","black"),
         pch=c(15,18,17,16,20,20),
         text.col = c("red","orange","purple","blue","grey","black")
  )
  
  
  
  plot(dep_ind_var.lm, pch = 18, col = "red", which=c(4))
  print(sort(cooks.distance(dep_ind_var.lm)))
}
par(mfrow=c(1,1),mar = c(2, 2, 2, 2))

pdf("Fig1f_iCSAV_dCSAV_corectrl_v2020419.pdf",width=6,height=6)
layout(mat = matrix(c(1),1,1,byrow=T))
Fig_1f(data$dCSAV,data$iCSAV,"iCSAV/dCSAV of core controls (n=62)", "Different amino acid convergences (dCSAV)", "Identical amino acid convergences (iCSAV)")%%box()
data$nTargets[12]
data$nTargets[15]
data$nTargets[35]
dev.off()
?lm
?abline
# Fig.1g
Fig_1g <- function(tmp_ind_var, tmp_dep_var, nTitle, name_ind_var, name_dep_var){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, main = nTitle, col="grey",type="n")
  grid(lwd = 0.3)
  # control sets
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)
  
  # regression analysis
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line
  
  
  res = cor.test(tmp_ind_var, tmp_dep_var,  method = "spearman", exact=F)
  print(res)
  
  xOrder  <- order(tmp_ind_var)
  x       <- tmp_ind_var[xOrder]  
  y       <- tmp_dep_var[xOrder]
  fit     <- lm(y ~ x, data=data.frame(x=x, y=y))
  newX    <- data.frame(x=jitter(x))
  #fitPred <- predict.lm(fit, newdata=newX, interval="prediction", level=0.95)
  #lines(newX$x, fitPred[,2], lty=2)
  #lines(newX$x, fitPred[,3], lty=2)
  
  
  # outlier
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[12][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[12][1])],
         col="black",pch=20)
  
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[57][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[57][1])],
         col="black",pch=20)
  
  points(tmp_ind_var[data$nTargets==as.vector(data$nTargets[59][1])],
         tmp_dep_var[data$nTargets==as.vector(data$nTargets[59][1])],
         col="black",pch=20)
  
  # set of birds of prey
  points(tmp_ind_var[data$nTargets=="CARCR,CATAU,FALPE,HALAL,HALLE,LEPDI,TYTAL"],
         tmp_dep_var[data$nTargets=="CARCR,CATAU,FALPE,HALAL,HALLE,LEPDI,TYTAL"],
         col="purple",pch=17)
  
  # set of waterbirds
  points(tmp_ind_var[data$nTargets=="APTFO,BALRE,CHAVO,EGRGA,FULGL,GAVST,NIPNI,PELCR,PHACA,PHALE,PHORU,PODCR,PYGAD"],
         tmp_dep_var[data$nTargets=="APTFO,BALRE,CHAVO,EGRGA,FULGL,GAVST,NIPNI,PELCR,PHACA,PHALE,PHORU,PODCR,PYGAD"],
         col="blue",pch=16)
  
  # set of closest control set
  points(tmp_ind_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="orange",pch=18)
  
  # set of vocal learners
  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)
  
  
  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  if(tmp_p <= 1){
    tmp_tmp_p = "1"
    if(tmp_p < 0.05){
      tmp_tmp_p = "0.05*"
      if(tmp_p < 0.01){
        tmp_tmp_p = "0.01**"
        if(tmp_p < 0.001){
          tmp_tmp_p = "0.001***"
        } 
      }
    }
  }
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) < Tmp_Pvalue),
                     list(Tmp_Pvalue = tmp_tmp_p))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    "bottomright",
    legend=rp, bty='n', text.col = "black")
  legend("topleft",
         legend = c("6 vocal learners","5 vocal learners+Swift","6 birds of prey","15 waterbirds","Core control sets", "Outlier (Cook's Dist.)"),
         col=c("red","orange","purple","blue","grey","black"),
         pch=c(15,18,17,16,20,20),
         text.col = c("red","orange","purple","blue","grey","black")
  )
  
  
  
  plot(dep_ind_var.lm, pch = 18, col = "red", which=c(4))
  print(sort(cooks.distance(dep_ind_var.lm)))
  
  
  
}

pdf("Fig1g_iCSAV_POB_corectrl_v2020419.pdf",width=6,height=6)
layout(mat = matrix(c(1),1,1,byrow=T))
Fig_1g(data$POB,data$iCSAV,"iCSAV/POB of core controls (n=62)", "Product of origin branch lengths (POB)", "Identical amino acid convergences (iCSAV)")%%box()
data$nTargets[12]
data$nTargets[57]
data$nTargets[59]
dev.off()


