library(tidyverse)
library(caret)#Sample partition, Support Vector Machine 
library(glmnet)#LASSO logistic regression
library(randomForest)#RF
library(ROCR)#ROC
library(pROC)#ROC
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time

ROC<-list.files("../",pattern = "ROC.rds",recursive = TRUE)
ROC_list<-list()
for (file in ROC) {
  temp<-readRDS(paste0("../",file))
  name<-word(file,start = 1,sep = "/")
  name<-word(name,start = 1,end = 4,sep = "_")
  if (grepl("training",file)) {
    ROC_list[["training"]][[name]]<-temp
  }
  if (grepl("testing",file)) {
    ROC_list[["testing"]][[name]]<-temp
  }
}
AUC<-list.files("../",pattern = "AUC.rds",recursive = TRUE)
AUC_list<-list()
for (file in AUC) {
  temp<-readRDS(paste0("../",file))
  name<-word(file,start = 1,sep = "/")
  name<-word(name,start = 1,end = 4,sep = "_")
  if (grepl("training",file)) {
    AUC_list[["training"]][[name]]<-temp
  }
  if (grepl("testing",file)) {
    AUC_list[["testing"]][[name]]<-temp
  }
}

tests<-names(AUC_list[["testing"]])
#####PLOT----
pdf("ROC-Braintumor_vs_Healthy.pdf", width = 7, height = 7)
#png("ROC-Braintumor_vs_Healthy.png", width = 7, height = 7, units = "in", res=500)
plot(ROC_list[["training"]][["CSF_Braintumor_vs_Healthy"]], col="#A6CEE3",lwd=2,lty=1,cex.lab=1.5, cex.axis=1.5, main="Brain cancers vs. Cancer-free")
plot(ROC_list[["testing"]][["CSF_Braintumor_vs_Healthy"]], add = TRUE, col="#1F78B4",lwd=3,lty=1)
plot(ROC_list[["training"]][["plasma_Braintumor_vs_Healthy"]], add = TRUE, col="#B2DF8A",lwd=2,lty=1)
plot(ROC_list[["testing"]][["plasma_Braintumor_vs_Healthy"]], add = TRUE, col="#33A02C",lwd=3,lty=1)
plot(ROC_list[["training"]][["CSF&plasma_Braintumor_vs_Healthy"]], add = TRUE, col="#FB9A99",lwd=2,lty=1)
plot(ROC_list[["testing"]][["CSF&plasma_Braintumor_vs_Healthy"]], add = TRUE, col="#E31A1C",lwd=3,lty=1)
legend(x=0.7, y=0.35, legend=paste0("CSF training AUC=",AUC_list[["training"]][["CSF_Braintumor_vs_Healthy"]]),col="#A6CEE3", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.3, legend=paste0("CSF testing AUC=",AUC_list[["testing"]][["CSF_Braintumor_vs_Healthy"]]),col="#1F78B4", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.25, legend=paste0("plasma training AUC=",AUC_list[["training"]][["plasma_Braintumor_vs_Healthy"]]),col="#B2DF8A", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.2, legend=paste0("plasma testing AUC=",AUC_list[["testing"]][["plasma_Braintumor_vs_Healthy"]]),col="#33A02C", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.15, legend=paste0("CSF&plasma training AUC=",AUC_list[["training"]][["CSF&plasma_Braintumor_vs_Healthy"]]),col="#FB9A99", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.1, legend=paste0("CSF&plasma testing AUC=",AUC_list[["testing"]][["CSF&plasma_Braintumor_vs_Healthy"]]),col="#E31A1C", lty=1,lwd=2, bty="n",cex=1.2)
dev.off()

pdf("ROC-Gliomas_vs_others.pdf", width = 7, height = 7)
#png("ROC-Gliomas_vs_others.png", width = 7, height = 7, units = "in", res=500)
plot(ROC_list[["training"]][["CSF_Gliomas_vs_others"]], col="#A6CEE3",lwd=2,lty=1,cex.lab=1.5, cex.axis=1.5, main="Gliomas vs. Others")
plot(ROC_list[["testing"]][["CSF_Gliomas_vs_others"]], add = TRUE, col="#1F78B4",lwd=3,lty=1)
plot(ROC_list[["training"]][["plasma_Gliomas_vs_others"]], add = TRUE, col="#B2DF8A",lwd=2,lty=1)
plot(ROC_list[["testing"]][["plasma_Gliomas_vs_others"]], add = TRUE, col="#33A02C",lwd=3,lty=1)
plot(ROC_list[["training"]][["CSF&plasma_Gliomas_vs_others"]], add = TRUE, col="#FB9A99",lwd=2,lty=1)
plot(ROC_list[["testing"]][["CSF&plasma_Gliomas_vs_others"]], add = TRUE, col="#E31A1C",lwd=3,lty=1)
legend(x=0.7, y=0.35, legend=paste0("CSF training AUC=",AUC_list[["training"]][["CSF_Gliomas_vs_others"]]),col="#A6CEE3", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.3, legend=paste0("CSF testing AUC=",AUC_list[["testing"]][["CSF_Gliomas_vs_others"]]),col="#1F78B4", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.25, legend=paste0("plasma training AUC=",AUC_list[["training"]][["plasma_Gliomas_vs_others"]]),col="#B2DF8A", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.2, legend=paste0("plasma testing AUC=",AUC_list[["testing"]][["plasma_Gliomas_vs_others"]]),col="#33A02C", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.15, legend=paste0("CSF&plasma training AUC=",AUC_list[["training"]][["CSF&plasma_Gliomas_vs_others"]]),col="#FB9A99", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.1, legend=paste0("CSF&plasma testing AUC=",AUC_list[["testing"]][["CSF&plasma_Gliomas_vs_others"]]),col="#E31A1C", lty=1,lwd=2, bty="n",cex=1.2)
dev.off()

pdf("ROC-Meningiomas_vs_others.pdf", width = 7, height = 7)
#png("ROC-Meningiomas_vs_others.png", width = 7, height = 7, units = "in", res=500)
plot(ROC_list[["training"]][["CSF_Meningiomas_vs_others"]], col="#A6CEE3",lwd=2,lty=1,cex.lab=1.5, cex.axis=1.5, main="Meningiomas vs. Others")
plot(ROC_list[["testing"]][["CSF_Meningiomas_vs_others"]], add = TRUE, col="#1F78B4",lwd=3,lty=1)
plot(ROC_list[["training"]][["plasma_Meningiomas_vs_others"]], add = TRUE, col="#B2DF8A",lwd=2,lty=1)
plot(ROC_list[["testing"]][["plasma_Meningiomas_vs_others"]], add = TRUE, col="#33A02C",lwd=3,lty=1)
plot(ROC_list[["training"]][["CSF&plasma_Meningiomas_vs_others"]], add = TRUE, col="#FB9A99",lwd=2,lty=1)
plot(ROC_list[["testing"]][["CSF&plasma_Meningiomas_vs_others"]], add = TRUE, col="#E31A1C",lwd=3,lty=1)
legend(x=0.7, y=0.35, legend=paste0("CSF training AUC=",AUC_list[["training"]][["CSF_Meningiomas_vs_others"]]),col="#A6CEE3", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.3, legend=paste0("CSF testing AUC=",AUC_list[["testing"]][["CSF_Meningiomas_vs_others"]]),col="#1F78B4", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.25, legend=paste0("plasma training AUC=",AUC_list[["training"]][["plasma_Meningiomas_vs_others"]]),col="#B2DF8A", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.2, legend=paste0("plasma testing AUC=",AUC_list[["testing"]][["plasma_Meningiomas_vs_others"]]),col="#33A02C", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.15, legend=paste0("CSF&plasma training AUC=",AUC_list[["training"]][["CSF&plasma_Meningiomas_vs_others"]]),col="#FB9A99", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.1, legend=paste0("CSF&plasma testing AUC=",AUC_list[["testing"]][["CSF&plasma_Meningiomas_vs_others"]]),col="#E31A1C", lty=1,lwd=2, bty="n",cex=1.2)
dev.off()
