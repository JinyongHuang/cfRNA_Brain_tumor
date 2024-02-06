library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(DescTools)
library(openxlsx)

sample<-read.csv("../CSXY_StudyDesign.csv")
meta<-read.csv("../../../Raw read counts and Summary by ID/Summary.csv", row.names = 1)
meta<-meta[rownames(meta) %in% sample$ID,]
sample$ID[!sample$ID %in% rownames(meta)]
meta<-merge(sample,meta,by.x = "ID",by.y = 0)
table(meta$Group)
colnames(meta)<-gsub("tRFs","tsRNA",colnames(meta))
meta<-meta[OrderMixed(meta$ID),]
write.csv(meta,"BrainTumor_Summary_all_libraries.csv")

cfRNA<-readRDS("../../../Raw read counts and Summary by ID/cfRNA.rds")
cfRNA<-cfRNA[,meta$ID]
cfRNA<-cfRNA[,OrderMixed(colnames(cfRNA))]
table(meta$ID==colnames(cfRNA))
##remove rsRNA and ysRNA
files = list.files(path = "../../../Raw read counts and Summary by ID/", pattern = "RNA.csv$", full.names = TRUE) 
RNA_names<-list()
for (file in files) {
  temp<-read.csv(file, row.names=1, check.names=FALSE)
  name<-word(file, start = -1, sep = "/")
  name<-word(name, start = 1, sep = "\\.")
  RNA_names[[name]]<-rownames(temp)
}
RNA_names<-RNA_names[-c(5,9)]
RNA_names[["cfRNA"]]<-c(RNA_names[["lncRNA"]],RNA_names[["miRNA"]],RNA_names[["mRNA"]],RNA_names[["piRNA"]],RNA_names[["snoRNA"]],RNA_names[["snRNA"]],RNA_names[["tsRNA"]])
saveRDS(RNA_names,"RNA_names.rds")
cfRNA<-cfRNA[rownames(cfRNA) %in% RNA_names[["cfRNA"]],]
##combine redo samples
temp<-data.frame(t(cfRNA),check.names = FALSE)
table(meta$ID==rownames(temp))
temp<-merge(temp,meta[,c(1,2,5)],by.x=0,by.y="ID")
rownames(temp)<-temp[,1];temp<-temp[,-1]
cfRNA_unique<-aggregate(temp[,-((ncol(temp)-1):ncol(temp))],by=list(temp$seqID),sum) 
rownames(cfRNA_unique)<-cfRNA_unique[,1];cfRNA_unique<-cfRNA_unique[,-1]
cfRNA_unique<-data.frame(t(cfRNA_unique),check.names = FALSE)
write.csv(cfRNA_unique,"BrainTumor_cfRNA_rawreads_all_samples.csv")

##unique summary----
Summary<-meta
Summary<-Summary[,-grep("_ratio",colnames(Summary))]
Summary<-Summary[,-grep("_detected",colnames(Summary))]
Summary<-Summary[,-grep("CleanRatio",colnames(Summary))]
Summary<-Summary[,-grep("MapRate",colnames(Summary))]
rownames(Summary)<-Summary[,1];Summary<-Summary[,-1]

Summary_unique<-aggregate(Summary[,-(1:10)],by=list(Summary$seqID),sum)
rownames(Summary_unique)<-Summary_unique[,1];Summary_unique<-Summary_unique[,-1]
Summary_unique <- Summary_unique %>%add_column(LibraryID=tapply(rownames(Summary), Summary$seqID, FUN = function(x) paste(x, collapse = ",")),.before="Total")
Summary_unique <- Summary_unique %>%add_column(CleanRatio = round(Summary_unique$Trimmed/Summary_unique$Total*100,2), .after = "Unmapped")
Summary_unique <- Summary_unique %>%add_column(MapRate = round((Summary_unique$Total-Summary_unique$Unmapped)/Summary_unique$Total*100,2), .after = "CleanRatio")
Summary_unique<-Summary_unique[,-grep("rsRNA",colnames(Summary_unique))]
Summary_unique<-Summary_unique[,-grep("ysRNA",colnames(Summary_unique))]
ratio<-round(Summary_unique[,grep("sum",colnames(Summary_unique))]/Summary_unique$Trimmed*100,2)
colnames(ratio)<-gsub("sum","ratio",colnames(ratio))
ratio$Unspecified_ratio<-100-apply(ratio[,-c(1:2)],1,sum)
Summary_unique<-Summary_unique[OrderMixed(rownames(Summary_unique)),]
ratio<-ratio[OrderMixed(rownames(ratio)),]
Summary_unique<-cbind(Summary_unique,ratio)

temp<-meta[,c(2,3,5:11)]
temp<-temp[!duplicated(temp$seqID),]
Summary_unique<-merge(temp,Summary_unique,by.x = "seqID",by.y = 0)
Summary_unique$Group<-gsub("_CSXY","",Summary_unique$Group)
Summary_unique$Group<-paste0(Summary_unique$Group,"_plasma")
Summary_unique$Group<-gsub("_CSF_plasma","_CSF",Summary_unique$Group)
#survival
survival<-readRDS("../survival/survival.rds")
Summary_unique<-merge(Summary_unique,survival,by.x = "clinicalID", by.y = "ID", all.x = TRUE)
write.csv(Summary_unique,"BrainTumor_Summary_all_samples.csv")

#failed sample----
cleanratiofail<-Summary_unique[Summary_unique$CleanRatio<20,]
table(cleanratiofail$Group)
cleanratiofail<-cleanratiofail$seqID

cleanreadfail<-Summary_unique[Summary_unique$Trimmed<3000000,]
table(cleanreadfail$Group)
cleanreadfail<-cleanreadfail$seqID

fail<-unique(c(cleanratiofail,cleanreadfail))


meta_good<-Summary_unique[!(Summary_unique$seqID %in% fail),]
table(meta_good$Group)
meta_good$specimen<-word(meta_good$Group,start = 2,sep = "_")
write.csv(meta_good,"BrainTumor_Summary_good_samples.csv")
cfRNA_good<-cfRNA_unique[,meta_good$seqID]
write.csv(cfRNA_good,"BrainTumor_cfRNA_rawreads_good_samples.csv")

##summary
summary(meta_good$Total)
summary(meta_good$Trimmed)
summary(meta_good$CleanRatio)
aggregate(meta_good[,"tsRNA_ratio"],by=list(meta_good$specimen),summary)
aggregate(meta_good[,"piRNA_ratio"],by=list(meta_good$specimen),summary)
aggregate(meta_good[,"miRNA_ratio"],by=list(meta_good$specimen),summary)
aggregate(meta_good[,"lncRNA_ratio"],by=list(meta_good$specimen),summary)
aggregate(meta_good[,"mRNA_ratio"],by=list(meta_good$specimen),summary)

##RPM normalization
table(meta_good$seqID==colnames(cfRNA_good))
trimmedreads<-meta_good$Trimmed
CPM<-apply(cfRNA_good,1,function(x) x/trimmedreads*1e+6)
CPM<-as.data.frame(t(CPM))
write.csv(CPM,"BrainTumor_cfRNA_RPM_good_samples.csv")
log2CPM<-log2(CPM+1)
write.csv(log2CPM,"BrainTumor_cfRNA_log2RPM_good_samples.csv")

##export
meta_good<-meta_good[order(meta_good$Group),]
table(meta_good$Group)
meta_good$newID<-c(paste0("GLI_CSF",1:16),paste0("GLI_plasma",1:18),paste0("MEN_CSF",1:43),
                   paste0("MEN_plasma",1:46),paste0("NOR_CSF",1:18),paste0("NOR_plasma",1:18))
log2CPM<-log2CPM[,OrderMixed(colnames(log2CPM))]
meta_good<-meta_good[OrderMixed(meta_good$seqID),]
table(colnames(log2CPM)==meta_good$seqID)
colnames(log2CPM)=meta_good$newID
table(colnames(log2CPM)==meta_good$newID)
log2CPM<-log2CPM[,OrderMixed(colnames(log2CPM))]
data<-log2CPM[rowSums(log2CPM[])>0,]
write.xlsx(data,"Supplementary data.xlsx", rowNames=T)
###quality controls
meta_good$Group<-gsub("NOR","CON",meta_good$Group)
meta_good$Group<-factor(meta_good$Group,levels = c("GLI_CSF","MEN_CSF","CON_CSF","GLI_plasma","MEN_plasma","CON_plasma"))
ggplot(meta_good, aes(x = Group, y = CleanRatio, color = Group)) + geom_boxplot(outlier.shape = NA) + geom_jitter() +
    labs(x="", y = "Clean read ratio (%)") +theme_classic()+ scale_color_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
    theme(text = element_text(size=12, color = "black"),legend.position = "none",
          axis.text.x = element_text(angle = 30, vjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15))
ggsave(paste0("Clean read ratio boxplot.pdf"), width=6, height = 6,units = "in")
ggsave(paste0("Clean read ratio boxplot.png"), width=6, height = 6,units = "in", bg="white")
ggplot(meta_good, aes(x = Group, y = Trimmed, color = Group)) + geom_boxplot(outlier.shape = NA) + geom_jitter() +
  labs(x="", y = "Clean reads") +theme_classic()+ scale_color_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
  theme(text = element_text(size=12, color = "black"),legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15))
ggsave(paste0("Clean reads boxplot.pdf"), width=6, height = 6,units = "in")
ggsave(paste0("Clean reads boxplot.png"), width=6, height = 6,units = "in", bg="white")
ggplot(meta_good, aes(x = Group, y = Total, color = Group)) + geom_boxplot(outlier.shape = NA) + geom_jitter() +
  labs(x="", y = "Total reads") +theme_classic()+ scale_color_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
  theme(text = element_text(size=12, color = "black"),legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15))
ggsave(paste0("Total reads boxplot.pdf"), width=6, height = 6,units = "in")
ggsave(paste0("Total reads boxplot.png"), width=6, height = 6,units = "in", bg="white")
ggplot(meta_good, aes(x = Group, y = tsRNA_ratio, color = Group)) + geom_boxplot(outlier.shape = NA) + geom_jitter() +
  labs(x="", y = "tsRNA ratio (%)") +theme_classic()+ scale_color_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
  theme(text = element_text(size=12, color = "black"),legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15))
ggsave(paste0("tsRNA ratio boxplot.pdf"), width=6, height = 3,units = "in")
ggsave(paste0("tsRNA ratio boxplot.png"), width=6, height = 3,units = "in", bg="white")
ggplot(meta_good, aes(x = Group, y = piRNA_ratio, color = Group)) + geom_boxplot(outlier.shape = NA) + geom_jitter() +
  labs(x="", y = "piRNA ratio (%)") +theme_classic()+ scale_color_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
  theme(text = element_text(size=12, color = "black"),legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15))
ggsave(paste0("piRNA ratio boxplot.pdf"), width=6, height = 6,units = "in")
ggsave(paste0("piRNA ratio boxplot.png"), width=6, height = 6,units = "in", bg="white")
ggplot(meta_good, aes(x = Group, y = miRNA_ratio, color = Group)) + geom_boxplot(outlier.shape = NA) + geom_jitter() +
  labs(x="", y = "miRNA ratio (%)") +theme_classic()+ scale_color_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
  theme(text = element_text(size=12, color = "black"),legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15))
ggsave(paste0("miRNA ratio boxplot.pdf"), width=6, height = 3,units = "in")
ggsave(paste0("miRNA ratio boxplot.png"), width=6, height = 3,units = "in", bg="white")
ggplot(meta_good, aes(x = Group, y = mRNA_ratio, color = Group)) + geom_boxplot(outlier.shape = NA) + geom_jitter() +
  labs(x="", y = "mRNA ratio (%)") +theme_classic()+ scale_color_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
  theme(text = element_text(size=12, color = "black"),legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15))
ggsave(paste0("mRNA ratio boxplot.pdf"), width=6, height = 6,units = "in")
ggsave(paste0("mRNA ratio boxplot.png"), width=6, height = 6,units = "in", bg="white")
ggplot(meta_good, aes(x = Group, y = lncRNA_ratio, color = Group)) + geom_boxplot(outlier.shape = NA) + geom_jitter() +
  labs(x="", y = "lncRNA ratio (%)") +theme_classic()+ scale_color_manual(values=c("#E31A1C","#FF7F00","#1F78B4","#FB9A99","#FDBF6F","#A6CEE3"))+
  theme(text = element_text(size=12, color = "black"),legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15))
ggsave(paste0("lncRNA ratio boxplot.pdf"), width=6, height = 6,units = "in")
ggsave(paste0("lncRNA ratio boxplot.png"), width=6, height = 6,units = "in", bg="white")


####read count----
Sum<-meta_good[,c(3,grep("sum",colnames(meta_good)))]
Sum<-apply(Sum[-1], 1, sum)
reads<-meta_good[,c(2,3,6,7)]
reads$Sum<-Sum
reads$Removed<-reads$Total-reads$Trimmed
reads$unspecified<-reads$Trimmed-reads$Sum

##individual total read count
temp<-reads[,c(1,5:7)]
colnames(temp)[2]<-"High quality annotated reads"
colnames(temp)[3]<-"Low quality reads"
colnames(temp)[4]<-"High quality unspecified reads"
plot_df<-melt(temp,"seqID")
table(plot_df$seqID==rep(reads$seqID,3))
plot_df$Group<-rep(reads$Group,3)
plot_df$Group<-factor(plot_df$Group,levels = c("GLI_CSF","MEN_CSF","CON_CSF","GLI_plasma","MEN_plasma","CON_plasma"))
temp<-temp[OrderMixed(temp$seqID),]
plot_df$seqID<-factor(plot_df$seqID,levels=temp$seqID)
plot_df$variable<-factor(plot_df$variable,levels=c("High quality annotated reads","High quality unspecified reads","Low quality reads"))
ggplot(plot_df, aes(x = seqID, y = value, fill = variable)) +geom_bar(stat = "identity") +
  labs(x="", y = "Total reads", title="") +theme_classic()+ scale_fill_brewer(palette = "Set1")+facet_wrap(~Group, scales = "free_x", nrow = 2)+
  theme(text = element_text(size=12, color = "black"),legend.title=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, margin = margin(t = -10)),
        axis.line = element_blank(),axis.ticks.x = element_blank(),strip.text = element_text(size=12,face="bold"))
ggsave(paste0("Stacked Bar Graph-Reads-all.pdf"), width=24, height = 16,units = "in")
ggsave(paste0("Stacked Bar Graph-Reads-all.png"), width=24, height = 16,units = "in", bg="white")
