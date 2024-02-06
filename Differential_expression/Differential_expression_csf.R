library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
library(RColorBrewer)
library(patchwork)#Combine ggplot
library(DESeq2)#Differential analysis
library(EnhancedVolcano)#Volcano

meta<-read.csv("../Quality controls/BrainTumor_Summary_good_samples.csv", row.names = 1)
meta<-meta[grep("CSF",meta$Group),]
table(meta$Group)
rownames(meta)<-meta[,"seqID"]
meta$Group<-as.factor(meta$Group)

cfRNA<-read.csv("../Quality controls/BrainTumor_cfRNA_rawreads_good_samples.csv", row.names = 1)
cfRNA<-cfRNA[,meta$seqID]
RNA_names<-readRDS("../Quality controls/RNA_names.rds")


####DESeq2----
df<-data.frame(specimen=c("CSF","CSF","CSF","CSF"),
               case=c("GLI_CSF","MEN_CSF","GLI_CSF","MEN_CSF"),
               control=c("NOR_CSF","NOR_CSF","MEN_CSF","GLI_CSF"),
               case_name=c("Gliomas","Meningiomas","Gliomas","Meningiomas"),
               control_name=c("Healthy","Healthy","Meningiomas","Gliomas"))
dds<-DESeqDataSetFromMatrix(cfRNA, meta, ~ Group)
dds<-DESeq(dds)  

DE_results<-list()
for (i in 1:nrow(df)) {
  test_name<-paste0(df$specimen[i]," - ",df$case_name[i]," vs. ",df$control_name[i])
  print(test_name)
  results<-results(dds, contrast = c("Group", df$case[i], df$control[i])) 
  #summary(results)
  results<-data.frame(results)
  results<-results[!(is.na(results$padj)),]
  results<-results[order(results$pvalue),]
  DE_results[[test_name]]<-results
}
saveRDS(DE_results,"DE_results.rds")

####Type distribution----
dir.create(paste0(getwd(),"/DESeq2 results"), showWarnings = FALSE)
dir.create(paste0(getwd(),"/DESeq2 significant results"), showWarnings = FALSE)
dir.create(paste0(getwd(),"/DESeq2 candidate results"), showWarnings = FALSE)
DE_results_distribution<-list()
for (i in 1:nrow(df)) {
  test_name<-paste0(df$specimen[i]," - ",df$case_name[i]," vs. ",df$control_name[i])
  results<-DE_results[[test_name]]
  all_list<-list()
  for (RNA in names(RNA_names)) {
    table<-results[(rownames(results) %in% RNA_names[[RNA]]),]
    all_list[[RNA]]<-table
    write.csv(table, paste0("DESeq2 results/",df$specimen[i],"-",RNA,"-",df$case_name[i]," vs. ",df$control_name[i],".csv"))
  }
  temp<-DE_results[[test_name]][["sig_results"]]
  sig_results<-results[results$padj < 0.1,]
  sig_results<-sig_results[order(sig_results$pvalue),]
  sig_list<-list()
  for (RNA in names(RNA_names)) {
    table<-sig_results[(rownames(sig_results) %in% RNA_names[[RNA]]),]
    sig_list[[RNA]]<-table
    write.csv(table, paste0("DESeq2 significant results/",df$specimen[i],"-",RNA,"-",df$case_name[i]," vs. ",df$control_name[i],".csv"))
  }
  temp<-DE_results[[test_name]][["candi_results"]]
  candi_results<-results[results$padj < 0.1 & results$baseMean>=10 & abs(results$log2FoldChange)>=0.5,]
  candi_results<-candi_results[order(candi_results$pvalue),]
  candi_list<-list()
  for (RNA in names(RNA_names)) {
    table<-candi_results[(rownames(candi_results) %in% RNA_names[[RNA]]),]
    candi_list[[RNA]]<-table
    write.csv(table, paste0("DESeq2 candidate results/",df$specimen[i],"-",RNA,"-",df$case_name[i]," vs. ",df$control_name[i],".csv"))
  }
  DE_results_distribution[[test_name]]<-list(all_list=all_list,sig_list=sig_list,candi_list=candi_list)
}
saveRDS(DE_results_distribution,"DE_results_distribution.rds")

####Summary----
Summary<-matrix(ncol=nrow(df)*2,nrow=length(RNA_names)+2)
Summary[1,]<-c(rep(names(DE_results),each = 2))
Summary[2,]<-c("gain","loss")
colnames(Summary)<-paste0(Summary[1,],"-",Summary[2,])
rownames(Summary)<-c("Comparison","GainLoss",names(RNA_names))
for (test in names(DE_results)) {
  name<-grep(test,colnames(Summary))
  for (RNA in names(RNA_names)) {
    temp<-DE_results_distribution[[test]][["sig_list"]][[RNA]]
    Summary[RNA,name[1]]<-nrow(subset(temp, log2FoldChange>0))
    Summary[RNA,name[2]]<-nrow(subset(temp, log2FoldChange<0))
  }
} 
Summary<-data.frame(Summary)
write.csv(Summary,"gain and loss.csv")

Summary_filter<-matrix(ncol=nrow(df)*2,nrow=length(RNA_names)+2)
Summary_filter[1,]<-c(rep(names(DE_results),each = 2))
Summary_filter[2,]<-c("gain_filter","loss_filter")
colnames(Summary_filter)<-paste0(Summary_filter[1,],"-",Summary_filter[2,])
rownames(Summary_filter)<-c("Comparison","GainLoss",names(RNA_names))
for (test in names(DE_results)) {
  name<-grep(test,colnames(Summary_filter))
  for (RNA in names(RNA_names)) {
    temp<-DE_results_distribution[[test]][["candi_list"]][[RNA]]
    Summary_filter[RNA,name[1]]<-nrow(subset(temp, log2FoldChange>0))
    Summary_filter[RNA,name[2]]<-nrow(subset(temp, log2FoldChange<0))
  }
} 
Summary_filter<-data.frame(Summary_filter)
write.csv(Summary_filter,"gain and loss filtered.csv")

####Volcano plot----
dir.create(paste0(getwd(),"/Volcano plot"), showWarnings = FALSE)
###cfRNA
for (i in 1:nrow(df)) {
  temp<-DE_results[[i]]
  x1=min(temp$log2FoldChange)
  if (x1>(-2)) {x1=(-2)}
  x2=max(temp$log2FoldChange)
  if (x2<=2) {x2=2}
  y=-log10(min(temp$padj))
  if (y<2) {y=2}
  increase<-nrow(subset(temp, padj<0.1 & log2FoldChange>0.5))
  decrease<-nrow(subset(temp, padj<0.1 & log2FoldChange<(-0.5)))
  Volcano<-EnhancedVolcano(temp, lab=NA,
                           x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' adjusted p value'), 
                           title = paste0(names(DE_results)[i]," - cfRNA"), subtitle = paste0("loss ",decrease,"  |  gain ",increase), caption = NULL,
                           pCutoff = 0.1, FCcutoff = 0.5, xlim = c(x1,x2), ylim = c(0,y), pointSize = 0.8, 
                           col = c("black","black","#4E62AB","#D6404E"), colAlpha = 0.9, legendLabels = c("NS", "FC", "FDR", "FDR & FC"))
  Volcano<-Volcano + theme_classic() + theme(axis.text = element_text(color = "black"), legend.position = "none",
                                             panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25),
                                                       plot.title = element_text(face="bold",hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
  ggsave(paste0("Volcano plot/all cfRNA-",df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i],".pdf"), width=6, height=8, units="in")
  ggsave(paste0("Volcano plot/all cfRNA-",df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i],".png"), width=6, height=8, units="in", bg="white")
}

###9 RNA types
for (i in 1:nrow(df)) {
  subset<-DE_results_distribution[[i]][["all_list"]]
  plot<-list()
  for (RNA in names(RNA_names)) {
    temp<-subset[[RNA]]
    x1=min(temp$log2FoldChange)
    if (x1>(-2)) {x1=(-2)}
    x2=max(temp$log2FoldChange)
    if (x2<=2) {x2=2}
    y=-log10(min(temp$padj))
    if (y<2) {y=2}
    increase<-nrow(subset(temp, padj<0.1 & log2FoldChange>0.5))
    decrease<-nrow(subset(temp, padj<0.1 & log2FoldChange<(-0.5)))
    Volcano<-EnhancedVolcano(temp, lab=NA,
                             x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' adjusted p value'), 
                             title = paste0(RNA), subtitle = paste0("loss ",decrease,"  |  gain ",increase), caption = NULL,
                             pCutoff = 0.1, FCcutoff = 0.5, xlim = c(x1,x2), ylim = c(0,y), pointSize = 0.8, 
                             col = c("black","black","#4E62AB","#D6404E"), colAlpha = 0.9, legendLabels = c("NS", "FC", "FDR", "FDR & FC"))
    Volcano<-Volcano + theme_classic() + theme(axis.title = element_blank(),axis.text = element_text(color = "black"), legend.position = "none",
                                               panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25),
                                               plot.title = element_text(face="bold",hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
    plot[[RNA]]<-Volcano
  }
  combine<-ggarrange(plot[["mRNA"]],plot[["lncRNA"]],plot[["miRNA"]],
                     plot[["piRNA"]],plot[["snRNA"]],plot[["snoRNA"]],
                     plot[["tsRNA"]],
                     ncol = 3, nrow = 3)
  annotate_figure(combine,
                  top = text_grob(paste0(names(DE_results)[i]),face="bold", size=15),
                  left = text_grob(expression('-Log'[10]*' adjusted p value'),rot = 90,size=12),
                  bottom= text_grob(expression('Log'[2]*' fold change'),size=12))
  ggsave(paste0("Volcano plot/combined RNA-",df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i],".pdf"), width=12, height=16, units="in")
  ggsave(paste0("Volcano plot/combined RNA-",df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i],".png"), width=12, height=16, units="in", bg="white")
}

##MA plot----
dir.create(paste0(getwd(),"/MA plot"), showWarnings = FALSE)
for (i in 1:nrow(df)) {
  temp<-DE_results[[i]]
  MA<-ggmaplot(temp, fdr=0.1, fc=(2^0.5), size=0.8, top = 0, 
               main=paste0(df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i]), legend = "top",
               xlab = "Log2 mean of normalized counts", ylab=expression('Log'[2]*' fold change'), ggtheme = ggplot2::theme_minimal())
  ggsave(paste0("MA plot/all cfRNA-",df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i],".pdf"), height = 5, width = 5, units="in")
  ggsave(paste0("MA plot/all cfRNA-",df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i],".png"), height = 5, width = 5, units="in", bg="white")
}

for (i in 1:nrow(df)) {
  subset<-DE_results_distribution[[i]][["all_list"]]
  plot<-list()
  for (RNA in names(RNA_names)) {
    temp<-subset[[RNA]]
    MA<-ggmaplot(temp, fdr=0.1, fc=(2^0.5), size=0.8, top = 0, 
                 main=paste0(RNA), legend = "top",
                 xlab = "", ylab="", ggtheme = ggplot2::theme_minimal())
    plot[[RNA]]<-MA
  }
  combine<-ggarrange(plot[["mRNA"]],plot[["lncRNA"]],plot[["miRNA"]],
                     plot[["piRNA"]],plot[["snRNA"]],plot[["snoRNA"]],
                     plot[["tsRNA"]],
                     ncol = 3, nrow = 3)
  annotate_figure(combine,
                  top = text_grob(paste0(names(DE_results)[i]),face="bold", size=15),
                  left = text_grob(expression('-Log'[2]*' fold change'),rot = 90,size=12),
                  bottom= text_grob(expression('Log'[2]*' mean of normalized counts'),size=12))
  ggsave(paste0("MA plot/combined RNA-",df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i],".pdf"), width=12, height=12, units="in")
  ggsave(paste0("MA plot/combined RNA-",df$specimen[i],"-",df$case_name[i]," vs. ",df$control_name[i],".png"), width=12, height=12, units="in", bg="white")
}

###Distribution----
#"CSF - Gliomas vs. Healthy"
sig_results<-DE_results_distribution[["CSF - Gliomas vs. Healthy"]][["candi_list"]]
Type_number<-data.frame(group=c("lncRNA","mRNA","miRNA","piRNA","snRNA","snoRNA","tsRNA"),
                        value=c(nrow(sig_results[["lncRNA"]]),nrow(sig_results[["mRNA"]]),nrow(sig_results[["miRNA"]]),nrow(sig_results[["piRNA"]]),
                                nrow(sig_results[["snRNA"]]),nrow(sig_results[["snoRNA"]]),nrow(sig_results[["tsRNA"]])))
Type_number <- Type_number %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(Type_number$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
Type_number$prop2<-paste0(round(Type_number$prop,1),"%")
Type_number$group<-factor(Type_number$group,levels=c("lncRNA","mRNA","miRNA","piRNA","snRNA","snoRNA","tsRNA"))
Type_number_out<-Type_number[Type_number$group== "miRNA" | Type_number$group== "piRNA" ,]
Type_number_in<-Type_number[Type_number$group== "lncRNA" | Type_number$group== "tsRNA" | Type_number$group== "snoRNA" |Type_number$group== "mRNA" | Type_number$group== "snRNA",]
ggplot(Type_number, aes(x="", y=prop, fill=group)) + geom_bar(stat="identity", width=0.5, color="black") + coord_polar("y", start=0) + 
  geom_text(data = Type_number_in, aes(y = ypos, label = paste(prop2)), color = "black", size=5, nudge_x = 0.1) +
  geom_label_repel(data = Type_number_out, aes(y = ypos, label = paste0(prop2)), size = 5, nudge_x = 0.6, show.legend = FALSE) +
  scale_fill_manual(values=c("#D6404E","#f78e26","#f7afb9","#F0D43A","#87CFA4","#4E62AB","#8b66b8"))+
  labs(title="CSF - Gliomas vs. Healthy") +theme_void() + 
  theme(plot.title = element_text(hjust=0.5,face="bold"),legend.text=element_text(size=14), legend.title=element_blank(), legend.position = "bottom")
ggsave(paste0("Pie plot-CSF - Gliomas vs. Healthy-Candidate results.pdf"), width=6, height = 8,units = "in")
#"CSF - Meningiomas vs. Healthy"
sig_results<-DE_results_distribution[["CSF - Meningiomas vs. Healthy"]][["candi_list"]]
Type_number<-data.frame(group=c("lncRNA","mRNA","miRNA","piRNA","snRNA","snoRNA","tsRNA"),
                        value=c(nrow(sig_results[["lncRNA"]]),nrow(sig_results[["mRNA"]]),nrow(sig_results[["miRNA"]]),nrow(sig_results[["piRNA"]]),
                                nrow(sig_results[["snRNA"]]),nrow(sig_results[["snoRNA"]]),nrow(sig_results[["tsRNA"]])))
Type_number <- Type_number %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(Type_number$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
Type_number$prop2<-paste0(round(Type_number$prop,1),"%")
Type_number$group<-factor(Type_number$group,levels=c("lncRNA","mRNA","miRNA","piRNA","snRNA","snoRNA","tsRNA"))
Type_number_out<-Type_number[Type_number$group== "snRNA" | Type_number$group== "snoRNA" ,]
Type_number_in<-Type_number[Type_number$group== "miRNA" | Type_number$group== "piRNA" | Type_number$group== "lncRNA" | Type_number$group== "tsRNA" |Type_number$group== "mRNA" ,]
ggplot(Type_number, aes(x="", y=prop, fill=group)) + geom_bar(stat="identity", width=0.5, color="black") + coord_polar("y", start=0) + 
  geom_text(data = Type_number_in, aes(y = ypos, label = paste(prop2)), color = "black", size=5, nudge_x = 0.1) +
  geom_label_repel(data = Type_number_out, aes(y = ypos, label = paste0(prop2)), size = 5, nudge_x = 0.6, show.legend = FALSE) +
  scale_fill_manual(values=c("#D6404E","#f78e26","#f7afb9","#F0D43A","#87CFA4","#4E62AB","#8b66b8"))+
  labs(title="CSF - Meningiomas vs. Healthy") +theme_void() + 
  theme(plot.title = element_text(hjust=0.5,face="bold"),legend.text=element_text(size=14), legend.title=element_blank(), legend.position = "bottom")
ggsave(paste0("Pie plot-CSF - Meningiomas vs. Healthy-Candidate results.pdf"), width=6, height = 8,units = "in")