#install.packages("dplyr")
library(dplyr)
#install.packages("corrplot")
library(corrplot)
#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")



# Random control sets #########################################################
setwd("E:\\1004_ctrl\\1004ctrl_TV_phy_20190202\\out\\wo_regional_variants\\1000random_sets_62core_sets\\1000ctrl")

#####################################################################

# Figure 2
cor.plot.wo.legend <- function(tmp_ind_var, tmp_dep_var,name_ind_var,name_dep_var,legend_pos_text, nTitle, data){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, col="grey",main=nTitle,type="n")
  grid(lwd = 1)
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)
  
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line
  
  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)
  #points(tmp_ind_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
  #       tmp_dep_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
  #       col="orange",pch=18)
  
  
  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == Tmp_Pvalue),
                     list(Tmp_Pvalue = format(tmp_p, digits = 2)))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    legend_pos_text, 
    legend=rp, bty='n', text.col = "black")
}
fNAME = "0.summary_Phy_TAAS_20190321.marked_targets.info"
nTitle = ""
phylo_TV_plot<- function(fNAME,nTitle){
  tmp_data <- read.table(fNAME,sep='\t',header = T)
  summary(tmp_data)
  attach(tmp_data)
  cor.plot.wo.legend(ProOriBL, TV_C,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProTerBL, TV_C,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTB, TV_C,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTN, TV_C,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProOriBL, TV_D,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProTerBL, TV_D,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTB, TV_D,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTN, TV_D,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProOriBL, TV_A,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProTerBL, TV_A,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTB, TV_A,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTN, TV_A,"","","topleft",nTitle, tmp_data)
  detach(tmp_data)
}




#jpeg(filename = "Fig2_test.jpeg",
#     width = 17, height = 17, units = "cm", pointsize = 12,
#     res = 300, bg = "white")
pdf("Fig2_1000randomctrl_v20200420.pdf",width=8,height=6)
par(mfrow=c(3,4), mar = c(3, 3, 0, 0))
#layout(mat = matrix(c(1:12),3,4,byrow=T))
phylo_TV_plot("0.summary_Phy_TAAS_20190321.marked_targets.info","")
dev.off()



# Figure S2

setwd("E:\\1004_ctrl\\1004ctrl_TV_phy_20190202\\out\\wo_regional_variants\\1000random_sets_62core_sets\\61ctrl")



cor.plot.wo.legend <- function(tmp_ind_var, tmp_dep_var,name_ind_var,name_dep_var,legend_pos_text, nTitle, data){
  plot(tmp_ind_var, tmp_dep_var, xlab = name_ind_var, ylab = name_dep_var, col="grey",main=nTitle,type="n")
  grid(lwd = 1)
  points(tmp_ind_var[data$ctrl=="1"],
         tmp_dep_var[data$ctrl=="1"],
         col="grey",pch=20)
  
  dep_ind_var.lm = lm(tmp_dep_var ~ tmp_ind_var)
  abline(dep_ind_var.lm,col="blue",lwd = 1,lty=1)# lty=2 # dashed line
  
  points(tmp_ind_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CALAN,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="red",pch=15)
  points(tmp_ind_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         tmp_dep_var[data$nTargets=="CHAPE,CORBR,GEOFO,MELUN,NESNO,TAEGU"],
         col="orange",pch=18)
  points(tmp_ind_var[data$nTargets=="CARCR,CATAU,FALPE,HALAL,HALLE,LEPDI,TYTAL"],
         tmp_dep_var[data$nTargets=="CARCR,CATAU,FALPE,HALAL,HALLE,LEPDI,TYTAL"],
         col="purple",pch=17)
  points(tmp_ind_var[data$nTargets=="APTFO,BALRE,CHAVO,EGRGA,FULGL,GAVST,NIPNI,PELCR,PHACA,PHALE,PHORU,PODCR,PYGAD"],
         tmp_dep_var[data$nTargets=="APTFO,BALRE,CHAVO,EGRGA,FULGL,GAVST,NIPNI,PELCR,PHACA,PHALE,PHORU,PODCR,PYGAD"],
         col="blue",pch=16)
  
  
  lm.sum = summary(dep_ind_var.lm)
  tmp_r2 = lm.sum$adj.r.squared
  tmp_p = lm.sum$coefficients[2,4]
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == Tmp_R2),
                     list(Tmp_R2 = format(tmp_r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == Tmp_Pvalue),
                     list(Tmp_Pvalue = format(tmp_p, digits = 2)))[2]
  legend(#min(tmp_ind_var)+(max(tmp_ind_var)-min(tmp_ind_var))*0.4,
    #min(tmp_dep_var)+(max(tmp_dep_var)-min(tmp_dep_var))*0.95,
    legend_pos_text, 
    legend=rp, bty='n', text.col = "black")
}
fNAME = "0.summary_Phy_TAAS_20181223.txt"
nTitle = ""
phylo_TV_plot<- function(fNAME,nTitle){
  tmp_data <- read.table(fNAME,sep='\t',header = T)
  summary(tmp_data)
  attach(tmp_data)
  cor.plot.wo.legend(ProOriBL, TV_A,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProTerBL, TV_A,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTB, TV_A,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTN, TV_A,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProOriBL, TV_C,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProTerBL, TV_C,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTB, TV_C,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTN, TV_C,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProOriBL, TV_D,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(ProTerBL, TV_D,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTB, TV_D,"","","topleft",nTitle, tmp_data)
  cor.plot.wo.legend(DistanceTN, TV_D,"","","topleft",nTitle, tmp_data)

  detach(tmp_data)
}



#jpeg(filename = "Fig2_test.jpeg",
#     width = 17, height = 17, units = "cm", pointsize = 12,
#     res = 300, bg = "white")
pdf("FigS2_61corectrl_v20200420.pdf",width=8,height=6)
par(mfrow=c(3,4), mar = c(3, 3, 0, 0))
#layout(mat = matrix(c(1:12),3,4,byrow=T))
phylo_TV_plot("0.summary_Phy_TAAS_20181223.txt","")
dev.off()



