library(tidyverse)
library(DescTools)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
library(RColorBrewer)
library(rstatix)
library(reshape2)
library(patchwork)#Combine ggplot
library(ggvenn)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time

DE_results<-readRDS("../DE_results_distribution.rds")
names(DE_results)
Gli_Nor<-DE_results[["plasma - Gliomas vs. Healthy"]][["candi_list"]][["cfRNA"]]
Men_Nor<-DE_results[["plasma - Meningiomas vs. Healthy"]][["candi_list"]][["cfRNA"]]
Gli_Men<-DE_results[["plasma - Gliomas vs. Meningiomas"]][["candi_list"]][["cfRNA"]]
Men_Gli<-DE_results[["plasma - Meningiomas vs. Gliomas"]][["candi_list"]][["cfRNA"]]
##gain
Gli_Nor_gain<-rownames(Gli_Nor[Gli_Nor$log2FoldChange>0,])
Men_Nor_gain<-rownames(Men_Nor[Men_Nor$log2FoldChange>0,])
Gli_Men_gain<-rownames(Gli_Men[Gli_Men$log2FoldChange>0,])
Men_Gli_gain<-rownames(Men_Gli[Men_Gli$log2FoldChange>0,])
##loss
Gli_Nor_loss<-rownames(Gli_Nor[Gli_Nor$log2FoldChange<0,])
Men_Nor_loss<-rownames(Men_Nor[Men_Nor$log2FoldChange<0,])
Gli_Men_loss<-rownames(Gli_Men[Gli_Men$log2FoldChange<0,])
Men_Gli_loss<-rownames(Men_Gli[Men_Gli$log2FoldChange<0,])

###Venn diagram----
dir.create(paste0(getwd(),"/Venn diagram"), showWarnings = FALSE)
###Gli_gain
Venn_list<-list(
  `Gli vs. Nor`=Gli_Nor_gain,
  `Men vs. Nor`=Men_Nor_gain,
  `Gli vs. Men`=Gli_Men_gain)
Venn <- ggvenn(Venn_list, show_percentage = T, set_name_size = 6, text_size = 5,fill_color = c("red", "orange","yellow"),fill_alpha = 0.5)
Venn+labs(title = "plasma - gain")+theme(plot.title = element_text(size=16,face = "bold",hjust=0.5))
ggsave("Venn diagram/plasma_Gli_gain.pdf", width=7, height = 7)
ggsave("Venn diagram/plasma_Gli_gain.png", width=7, height = 7, bg = "white")
###Gli_loss
Venn_list<-list(
  `Gli vs. Nor`=Gli_Nor_loss,
  `Men vs. Nor`=Men_Nor_loss,
  `Gli vs. Men`=Gli_Men_loss)
Venn <- ggvenn(Venn_list, show_percentage = T, set_name_size = 6, text_size = 5,fill_color = c("red", "orange","yellow"),fill_alpha = 0.5)
Venn+labs(title = "plasma - loss")+theme(plot.title = element_text(size=16,face = "bold",hjust=0.5))
ggsave("Venn diagram/plasma_Gli_loss.pdf", width=7, height = 7)
ggsave("Venn diagram/plasma_Gli_loss.png", width=7, height = 7, bg = "white")

###Men_gain
Venn_list<-list(
  `Gli vs. Nor`=Gli_Nor_gain,
  `Men vs. Nor`=Men_Nor_gain,
  `Men vs. Gli`=Men_Gli_gain)
Venn <- ggvenn(Venn_list, show_percentage = T, set_name_size = 6, text_size = 5,fill_color = c("red", "orange","yellow"),fill_alpha = 0.5)
Venn+labs(title = "plasma - gain")+theme(plot.title = element_text(size=16,face = "bold",hjust=0.5))
ggsave("Venn diagram/plasma_Men_gain.pdf", width=7, height = 7)
ggsave("Venn diagram/plasma_Men_gain.png", width=7, height = 7, bg = "white")
###Men_loss
Venn_list<-list(
  `Gli vs. Nor`=Gli_Nor_loss,
  `Men vs. Nor`=Men_Nor_loss,
  `Men vs. Gli`=Men_Gli_loss)
Venn <- ggvenn(Venn_list, show_percentage = T, set_name_size = 6, text_size = 5,fill_color = c("red", "orange","yellow"),fill_alpha = 0.5)
Venn+labs(title = "plasma - loss")+theme(plot.title = element_text(size=16,face = "bold",hjust=0.5))
ggsave("Venn diagram/plasma_Men_loss.pdf", width=7, height = 7)
ggsave("Venn diagram/plasma_Men_loss.png", width=7, height = 7, bg = "white")

###goi----
##three groups overlapping
RNA_names<-readRDS("../../Quality controls/RNA_names.rds")
BrainTumor <- c(c(Gli_Nor_gain,Men_Nor_gain)[duplicated(c(Gli_Nor_gain,Men_Nor_gain))],
                c(Gli_Nor_loss,Men_Nor_loss)[duplicated(c(Gli_Nor_loss,Men_Nor_loss))])
BrainTumor_RNA<-list()
for (RNA in names(RNA_names)[-8]) {
  BrainTumor_RNA[[RNA]]<-BrainTumor[BrainTumor %in% RNA_names[[RNA]]]
}
Gliomas <-c(c(Gli_Nor_gain,Gli_Men_gain)[duplicated(c(Gli_Nor_gain,Gli_Men_gain))],
            c(Gli_Nor_loss,Gli_Men_loss)[duplicated(c(Gli_Nor_loss,Gli_Men_loss))])
Gliomas_RNA<-list()
for (RNA in names(RNA_names)[-8]) {
  Gliomas_RNA[[RNA]]<-Gliomas[Gliomas %in% RNA_names[[RNA]]]
}
Meningiomas<-c(c(Men_Nor_gain,Men_Gli_gain)[duplicated(c(Men_Nor_gain,Men_Gli_gain))],
               c(Men_Nor_loss,Men_Gli_loss)[duplicated(c(Men_Nor_loss,Men_Gli_loss))])
Meningiomas_RNA<-list()
for (RNA in names(RNA_names)[-8]) {
  Meningiomas_RNA[[RNA]]<-Meningiomas[Meningiomas %in% RNA_names[[RNA]]]
}

goi<-list(BrainTumor=BrainTumor_RNA,Gliomas=Gliomas_RNA,Meningiomas=Meningiomas_RNA)
saveRDS(goi,"plasma_GOI.rds")

####boxplot----
dir.create(paste0(getwd(),"/Boxplot"), showWarnings = FALSE)
meta<-read.csv("../../Quality controls/BrainTumor_Summary_good_samples.csv", row.names = 1)
meta<-meta[OrderMixed(meta$seqID),]
table(meta$Group)
cfRNA<-read.csv("../../Quality controls/BrainTumor_cfRNA_log2RPM_good_samples.csv", row.names = 1)
cfRNA<-cfRNA[,OrderMixed(colnames(cfRNA))]
table(meta$seqID==colnames(cfRNA))
for (name in names(goi)) {
  dir.create(paste0(getwd(),"/Boxplot/",name), showWarnings = FALSE)
  for (RNA in names(RNA_names)[-8]) {
    dir.create(paste0(getwd(),"/Boxplot/",name,"/",RNA), showWarnings = FALSE)
    subset<-cfRNA[goi[[name]][[RNA]],]
    if (nrow(subset)==0) next
    subset<-data.frame(t(subset),check.names = FALSE)
    subset$Group<-meta$Group
    N=1:(ncol(subset)-1)
    registerDoParallel(cl<-makeCluster(10))
    result_seed<-foreach(i=N, .combine="c",.packages=c("tidyverse","ggplot2","ggpubr","rstatix","reshape2")) %dopar% {
      df<-subset[,c(i,ncol(subset))]
      df$Group<-factor(df$Group,levels = c("GLI_CSF","MEN_CSF","NOR_CSF","GLI_plasma","MEN_plasma","NOR_plasma"))
      colnames(df)[1]<-"RNA"
      variance<-aggregate(df$RNA,by=list(df$Group),var)
      if (table(variance$x>5)['FALSE']==6) {
        sum<-aggregate(df$RNA,by=list(df$Group),sum)
        if (table(sum$x==0)['FALSE']==6) {
          stat.test <- df %>%
            wilcox_test(RNA ~ Group, paired = FALSE, p.adjust.method = "none")%>%
            add_significance("p")%>%
            add_xy_position(fun="max")
          stat.test<-stat.test[c(1,2,6,13,14,15),]
          max<-max(df$RNA)
          stat.test$y.position[c(1,4)]<-max+0.5
          stat.test$y.position[c(2,5)]<-max+1.5
          stat.test$y.position[c(3,6)]<-max+2.5
          if (name=="BrainTumor" & stat.test$p[5]<0.05 & stat.test$p[6]<0.05) {
            plot <- ggplot(df, aes(x=Group, y=RNA)) + geom_boxplot(aes(fill = Group),outlier.shape=NA)+
              scale_fill_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
              geom_jitter(shape=16, position=position_jitter(0.2))+
              labs(title=paste0(colnames(subset)[i]), x="", y = expression('Log'[2]*' (RPM)'))+ theme_classic()+
              theme(text = element_text(size = 14),legend.position = "none",
                    axis.text.x = element_text(colour="black",angle = 45,hjust=1),axis.text.y = element_text(colour="black"))
            plot + stat_pvalue_manual(stat.test, label = "p",tip.length = 0)
            ggsave(paste0("Boxplot/",name,"/",RNA,"/",colnames(subset)[i],".png"),device = "png",bg="white",width = 6, height = 6)
          }
          if (name=="Gliomas" & stat.test$p[4]<0.05 & stat.test$p[5]<0.05) {
            plot <- ggplot(df, aes(x=Group, y=RNA)) + geom_boxplot(aes(fill = Group),outlier.shape=NA)+
              scale_fill_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
              geom_jitter(shape=16, position=position_jitter(0.2))+
              labs(title=paste0(colnames(subset)[i]), x="", y = expression('Log'[2]*' (RPM)'))+ theme_classic()+
              theme(text = element_text(size = 14),legend.position = "none",
                    axis.text.x = element_text(colour="black",angle = 45,hjust=1),axis.text.y = element_text(colour="black"))
            plot + stat_pvalue_manual(stat.test, label = "p",tip.length = 0)
            ggsave(paste0("Boxplot/",name,"/",RNA,"/",colnames(subset)[i],".png"),device = "png",bg="white",width = 6, height = 6)
          }
          if (name=="Meningiomas" & stat.test$p[4]<0.05 & stat.test$p[6]<0.05) {
            plot <- ggplot(df, aes(x=Group, y=RNA)) + geom_boxplot(aes(fill = Group),outlier.shape=NA)+
              scale_fill_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
              geom_jitter(shape=16, position=position_jitter(0.2))+
              labs(title=paste0(colnames(subset)[i]), x="", y = expression('Log'[2]*' (RPM)'))+ theme_classic()+
              theme(text = element_text(size = 14),legend.position = "none",
                    axis.text.x = element_text(colour="black",angle = 45,hjust=1),axis.text.y = element_text(colour="black"))
            plot + stat_pvalue_manual(stat.test, label = "p",tip.length = 0)
            ggsave(paste0("Boxplot/",name,"/",RNA,"/",colnames(subset)[i],".png"),device = "png",bg="white",width = 6, height = 6)
          }
        }
      }
    }
    stopCluster(cl)
  }
}

ploted_goi<-list.files("Boxplot/Gliomas/",pattern = ".png",recursive = TRUE)
ploted_goi<-word(ploted_goi,start = 2,sep = "/")
ploted_goi<-gsub(".png","",ploted_goi)
saveRDS(ploted_goi,"ploted_goi_plasma_Gliomas.rds")
writeLines(ploted_goi,"ploted_goi_plasma_Gliomas.txt")
ploted_goi<-list.files("Boxplot/Meningiomas/",pattern = ".png",recursive = TRUE)
ploted_goi<-word(ploted_goi,start = 2,sep = "/")
ploted_goi<-gsub(".png","",ploted_goi)
saveRDS(ploted_goi,"ploted_goi_plasma_Meningiomas.rds")
writeLines(ploted_goi,"ploted_goi_plasma_Meningiomas.txt")
ploted_goi<-list.files("Boxplot/BrainTumor/",pattern = ".png",recursive = TRUE)
ploted_goi<-word(ploted_goi,start = 2,sep = "/")
ploted_goi<-gsub(".png","",ploted_goi)
saveRDS(ploted_goi,"ploted_goi_plasma_BrainTumor.rds")
writeLines(ploted_goi,"ploted_goi_plasma_BrainTumor.txt")
