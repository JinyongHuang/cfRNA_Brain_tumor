library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
library(RColorBrewer)
library(patchwork)#Combine ggplot
library(DESeq2)#Differential analysis
library(EnhancedVolcano)#Volcano
library(rstatix)
library(reshape2)
library(DescTools)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time
library(ComplexHeatmap)#Heatmap
library(circlize)#colorRamp2
library(Rtsne)

CSF_Braintumor_goi<-readLines("../CSF_Braintumor_vs_Healthy_p=0.6/goi/5-name.txt")
CSF_Gliomas_goi<-readLines("../CSF_Gliomas_vs_others_p=0.6/goi/5-name.txt")
CSF_Meningiomas_goi<-readLines("../CSF_Meningiomas_vs_others_p=0.6/goi/5-name.txt")
plasma_Braintumor_goi<-readLines("../plasma_Braintumor_vs_Healthy_p=0.6/goi/1-name.txt")
plasma_Gliomas_goi<-readLines("../plasma_Gliomas_vs_others_p=0.6/goi/1-name.txt")
plasma_Meningiomas_goi<-readLines("../plasma_Meningiomas_vs_others_p=0.6/goi/5-name.txt")

goi<-list(CSF_Braintumor=CSF_Braintumor_goi,CSF_Gliomas=CSF_Gliomas_goi,CSF_Meningiomas=CSF_Meningiomas_goi,
          plasma_Braintumor=plasma_Braintumor_goi,plasma_Gliomas=plasma_Gliomas_goi,plasma_Meningiomas=plasma_Meningiomas_goi)
####boxplot all goi----
dir.create(paste0(getwd(),"/Boxplot"), showWarnings = FALSE)
meta<-read.csv("../../Quality controls/BrainTumor_Summary_good_samples.csv", row.names = 1)
meta<-meta[OrderMixed(meta$seqID),]
table(meta$Group)
cfRNA<-read.csv("../../Quality controls/BrainTumor_cfRNA_log2RPM_good_samples.csv", row.names = 1)
cfRNA<-cfRNA[,OrderMixed(colnames(cfRNA))]
table(meta$seqID==colnames(cfRNA))
for (name in names(goi)) {
  dir.create(paste0(getwd(),"/Boxplot/",name), showWarnings = FALSE)
  subset<-cfRNA[goi[[name]],]
  if (nrow(subset)==0) next
  subset<-data.frame(t(subset),check.names = FALSE)
  subset$Group<-meta$Group
  N=1:(ncol(subset)-1)
  registerDoParallel(cl<-makeCluster(10))
  result_seed<-foreach(i=N, .combine="c",.packages=c("tidyverse","ggplot2","ggpubr","rstatix","reshape2")) %dopar% {
    df<-subset[,c(i,ncol(subset))]
    df$Group<-factor(df$Group,levels = c("GLI_CSF","MEN_CSF","NOR_CSF","GLI_plasma","MEN_plasma","NOR_plasma"))
    colnames(df)[1]<-"RNA"
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
      plot <- ggplot(df, aes(x=Group, y=RNA)) + geom_boxplot(aes(fill = Group),outlier.shape=NA)+
        scale_fill_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
        geom_jitter(shape=16, position=position_jitter(0.2))+
        labs(title=paste0(colnames(subset)[i]), x="", y = expression('Log'[2]*' (RPM)'))+ theme_classic()+
        theme(text = element_text(size = 14),legend.position = "none",
              axis.text.x = element_text(colour="black",angle = 45,hjust=1),axis.text.y = element_text(colour="black"))
      plot + stat_pvalue_manual(stat.test, label = "p",tip.length = 0)
      ggsave(paste0("Boxplot/",name,"/",colnames(subset)[i],".png"),device = "png",bg="white",width = 6, height = 6)
      ggsave(paste0("Boxplot/",name,"/",colnames(subset)[i],".pdf"),width = 5, height = 5)
    }
  }
  stopCluster(cl)
}

# ####Volcano plot----
DE_results_csf<-readRDS("../../DESeq2_One_vs_One_csf/DE_results_distribution.rds")
DE_results_plasma<-readRDS("../../DESeq2_One_vs_One_plasma//DE_results_distribution.rds")
DE_results<-c(DE_results_csf[1:2],DE_results_plasma[1:2])
dir.create(paste0(getwd(),"/Volcano plot"), showWarnings = FALSE)
for (i in 1:length(DE_results)) {
  temp<-DE_results[[i]][["all_list"]][["cfRNA"]]
  increase<-nrow(subset(temp, padj<0.1 & log2FoldChange>0))
  decrease<-nrow(subset(temp, padj<0.1 & log2FoldChange<0))
  temp<-temp[temp$log2FoldChange>(-10),]
  x1=min(temp$log2FoldChange)
  if (x1>(-2)) {x1=(-2)}
  x2=max(temp$log2FoldChange)
  if (x2<=2) {x2=2}
  y=-log10(min(temp$padj))
  if (y<2) {y=2}
  if (i==1) {label<-rownames(temp[rownames(temp) %in% goi[[2]],])}
  if (i==2) {label<-rownames(temp[rownames(temp) %in% goi[[3]],])}
  if (i==3) {label<-rownames(temp[rownames(temp) %in% goi[[5]],])}
  if (i==4) {label<-rownames(temp[rownames(temp) %in% goi[[6]],])}
  Volcano<-EnhancedVolcano(temp, lab=rownames(temp), selectLab=label,labSize = 3, boxedLabels = FALSE, drawConnectors = TRUE, arrowheads = FALSE, widthConnectors = 0.1,
                           x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' adjusted p value'),
                           title = paste0(names(DE_results)[i]), subtitle = paste0("loss ",decrease,"  |  gain ",increase), caption = NULL,
                           pCutoff = 0.1, FCcutoff = 0.5, xlim = c(x1,x2), ylim = c(0,y), pointSize = 0.8,
                           col = c("grey","grey","#4E62AB","#D6404E"), colAlpha = 0.9, legendLabels = c("NS", "FC", "FDR", "FDR & FC"))
  Volcano<-Volcano + theme_classic() + theme(axis.text = element_text(color = "black"), legend.position = "none",
                                             panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25),
                                                       plot.title = element_text(face="bold",hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
  ggsave(paste0("Volcano plot/",names(DE_results)[i],".pdf"), width=4, height=5, units="in")
  ggsave(paste0("Volcano plot/",names(DE_results)[i],".png"), width=4, height=5, units="in", bg="white")
}

##Heatmap----
##csf
GOI<-c(goi[[2]],goi[[3]],goi[[1]])
table(meta$Group)
data<-cfRNA[GOI,grep("csf",colnames(cfRNA))]
z.score <- t(apply(data,1,scale)) #z-score
colnames(z.score) <- colnames(data)
reorder<-c(meta[meta$Group=="GLI_CSF",][,"seqID"],meta[meta$Group=="MEN_CSF",][,"seqID"],meta[meta$Group=="NOR_CSF",][,"seqID"])
z.score <- z.score[,reorder]
type <- c(rep("GLI_CSF",16),rep("MEN_CSF",43),rep("NOR_CSF",18))
type<-factor(type,levels = c("GLI_CSF","MEN_CSF","NOR_CSF"))
col_type <- c("GLI_CSF" = "#E31A1C", "MEN_CSF" = "#FF7F00", "NOR_CSF" = "#1F78B4")
top<- HeatmapAnnotation(`Sample group` = type, col = list(`Sample group` = col_type), show_annotation_name = F, which="column")
GOItype<-c(rep("Gliomas",13),rep("Meningiomas",3),rep("Brain tumor",14))
GOItype<-factor(GOItype,levels = c("Gliomas","Meningiomas","Brain tumor"))
col_GOItype<-c("Gliomas" = "#FB9A99", "Meningiomas" = "#FDBF6F", "Brain tumor" = "#A6CEE3")
left<- HeatmapAnnotation(`GOI cfRNAs` = GOItype, col = list(`GOI cfRNAs` = col_GOItype), show_annotation_name = F, which="row")
h<-Heatmap(z.score, cluster_rows = T, show_row_dend = F, cluster_columns = F, show_column_dend = F, 
           clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
           # clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
           col = colorRamp2(c(quantile(z.score,0.05), 0, quantile(z.score,0.95)), c("#4E62AB", "white", "#9E0142")), 
           show_column_names = FALSE, show_row_names = TRUE, row_title = NULL,column_title = NULL, column_split = type, row_split = GOItype,
           name="z-score (log2 RPM)",top_annotation = top, left_annotation = left, rect_gp = gpar(col = "grey", lwd = 0.5))
print(h)

##plasma
GOI<-c(goi[[5]],goi[[6]],goi[[4]])
table(meta$Group)
data<-cfRNA[GOI,]
z.score <- t(apply(data,1,scale)) #z-score
colnames(z.score) <- colnames(cfRNA)
reorder<-c(meta[meta$Group=="GLI_plasma",][,"seqID"],meta[meta$Group=="MEN_plasma",][,"seqID"],meta[meta$Group=="NOR_plasma",][,"seqID"])
z.score <- z.score[,reorder]
type <- c(rep("GLI_plasma",18),rep("MEN_plasma",46),rep("NOR_plasma",18))
type<-factor(type,levels = c("GLI_plasma","MEN_plasma","NOR_plasma"))
col_type <- c("GLI_plasma" = "#FB9A99", "MEN_plasma" = "#FDBF6F", "NOR_plasma" = "#A6CEE3")
top<- HeatmapAnnotation(Group = type, col = list(Group = col_type), show_annotation_name = F, which="column")
GOItype<-c(rep("Gliomas",2),rep("Meningiomas",13),rep("Brain tumor",5))
GOItype<-factor(GOItype,levels = c("Gliomas","Meningiomas","Brain tumor"))
col_GOItype<-c("Gliomas" = "#E31A1C", "Meningiomas" = "#FF7F00", "Brain tumor" = "#1F78B4")
left<- HeatmapAnnotation(`GOI cfRNAs` = GOItype, col = list(`GOI cfRNAs` = col_GOItype), show_annotation_name = F, which="row")

h<-Heatmap(z.score, cluster_rows = T, show_row_dend = F, cluster_columns = F, show_column_dend = F, 
           clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
           # clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
           col = colorRamp2(c(quantile(z.score,0.05), 0, quantile(z.score,0.95)), c("#4E62AB", "white", "#9E0142")), 
           show_column_names = FALSE, show_row_names = TRUE, row_title = NULL,column_title = NULL, column_split = type, row_split = GOItype,
           name="z-score (log2 RPM)",top_annotation = top, left_annotation = left, rect_gp = gpar(col = "grey", lwd = 0.5))
print(h)

####Tsne----
##csf
GOI<-c(goi[[2]],goi[[3]],goi[[1]])
temp<-meta[grep("csf",meta$seqID),]
data<-cfRNA[GOI,temp$seqID]
data<-data.frame(t(data),check.names = FALSE)
for (i in 1:10){
  set.seed(i)
  tsne <- Rtsne(data, dims = 2, perplexity=10, check_duplicates = FALSE)
  tsne <-as.data.frame(tsne$Y)
  rownames(tsne)<-temp$seqID
  colnames(tsne)<-c("tSNE1","tSNE2")
  table(rownames(tsne)==temp$seqID)
  tsne$Group<-factor(temp$Group)
  ggplot(tsne, aes(x=tSNE1, y=tSNE2,col= Group,fill=Group)) + 
    geom_point() + theme_classic()+ labs(x="t-SNE 1",y="t-SNE 2", title="csf cfRNAs of interest")+
    stat_ellipse(level=0.5, geom="polygon", alpha=1/3)+
    scale_fill_manual(values =c("#E31A1C","#FF7F00","#1F78B4"))+
    scale_color_manual(values =c("#E31A1C","#FF7F00","#1F78B4"))+
    theme(legend.title = element_blank(),legend.position = c(0.15,0.15),legend.background = element_blank(),
          axis.text = element_text(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),
          axis.line = element_line(linewidth=0.15), text = element_text(size=10))
  ggsave(paste0("tsne/csf/csf-2D tsne_",i,".png"),width = 5, height = 5, bg="white")
  ggsave(paste0("tsne/csf/csf-2D tsne_",i,".pdf"),width = 5, height = 5)
}
##plasma
GOI<-c(goi[[5]],goi[[6]],goi[[4]])
temp<-meta[grep("csf",meta$seqID, invert = TRUE),]
data<-cfRNA[GOI,temp$seqID]
data<-data.frame(t(data),check.names = FALSE)
for (i in 1:10){
  set.seed(i)
  tsne <- Rtsne(data, dims = 2, perplexity=10, check_duplicates = FALSE)
  tsne <-as.data.frame(tsne$Y)
  rownames(tsne)<-temp$seqID
  colnames(tsne)<-c("tSNE1","tSNE2")
  table(rownames(tsne)==temp$seqID)
  tsne$Group<-factor(temp$Group)
  ggplot(tsne, aes(x=tSNE1, y=tSNE2,col= Group,fill=Group)) + 
    geom_point() + theme_classic()+ labs(x="t-SNE 1",y="t-SNE 2", title="plasma cfRNAs of interest")+
    stat_ellipse(level=0.5, geom="polygon", alpha=1/3)+
    scale_fill_manual(values =c("#E31A1C","#FF7F00","#1F78B4"))+
    scale_color_manual(values =c("#E31A1C","#FF7F00","#1F78B4"))+
    theme(legend.title = element_blank(),legend.position = c(0.15,0.15),legend.background = element_blank(),
          axis.text = element_text(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),
          axis.line = element_line(linewidth=0.15), text = element_text(size=10))
  ggsave(paste0("tsne/plasma/plasma-2D tsne_",i,".png"),width = 5, height = 5, bg="white")
  ggsave(paste0("tsne/plasma/plasma-2D tsne_",i,".pdf"),width = 5, height = 5)
}