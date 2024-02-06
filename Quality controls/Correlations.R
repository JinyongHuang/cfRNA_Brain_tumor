library(tidyverse)
library(DescTools)
#library(edgeR) #cpm()
library(gtools)
library(ggplot2)
library(ggpubr)

####meta----
meta<-read.csv("../Quality controls/BrainTumor_Summary_good_samples.csv", row.names = 1)
table(meta$Group)
cfRNA_rpm<-read.csv("../Quality controls/BrainTumor_cfRNA_log2RPM_good_samples.csv", row.names = 1)

####RNA names
RNA_names<-readRDS("../Quality controls/RNA_names.rds")

RNA_rpm<-list()
for (RNA in names(RNA_names)) {
  RNA_rpm[[RNA]]<-cfRNA_rpm[rownames(cfRNA_rpm) %in% RNA_names[[RNA]],]
}

####Pearson's correlation coefficient (R) & Coefficient of Determination (R-Squared)----
combination<-data.frame(X1=c("GLI_CSF","GLI_CSF","MEN_CSF","GLI_plasma","GLI_plasma","MEN_plasma"),
                        X2=c("MEN_CSF","NOR_CSF","NOR_CSF","MEN_plasma","NOR_plasma","NOR_plasma"))
combination$X3=paste(combination$X1,"vs.",combination$X2)

corcoef <-matrix(nrow=nrow(combination),ncol=length(RNA_rpm))
rownames(corcoef)<-str_c(combination[,1]," vs. ",combination[,2])
colnames(corcoef)<-names(RNA_rpm)
RSquared<-matrix(nrow=nrow(combination),ncol=length(RNA_rpm))
rownames(RSquared)<-str_c(combination[,1]," vs. ",combination[,2])
colnames(RSquared)<-names(RNA_rpm)
for (i in 1:nrow(combination)) {
  print(combination[i,3])
  group1<-combination[i,1]
  id1<-meta[meta$Group==group1,][,"seqID"]
  group2<-combination[i,2]
  id2<-meta[meta$Group==group2,][,"seqID"]
  for (r in 1:length(RNA_rpm)){
    temp1<-RNA_rpm[[r]][,id1]
    temp2<-RNA_rpm[[r]][,id2]
    temp <- merge(apply(temp1,1,mean),apply(temp2,1,mean),by=0)
    corcoef[i,r]=cor(temp$x,temp$y)
    model <- lm(y~x, data=temp)
    RSquared[i,r]=summary(model)$r.squared
  }
}
write.csv(corcoef,"Correlation coefficient (R) by cohort.csv")
write.csv(RSquared,"Coefficient of Determination (R-Squared) by cohort.csv")

# ####Scatterplot & Pearson's correlation----
dir.create(paste0(getwd(),"/Scatterplot"), showWarnings = FALSE)
plot_list<-list()
plot_cfRNA_list<-list()
for (i in 1:nrow(combination)) {
  test<-combination[i,3]
  print(test)
  group1<-combination[i,1]
  id1<-meta[meta$Group==group1,][,"seqID"]
  group2<-combination[i,2]
  id2<-meta[meta$Group==group2,][,"seqID"]
  for (RNA in names(RNA_rpm)){
    temp1<-RNA_rpm[[RNA]][,id1]
    temp2<-RNA_rpm[[RNA]][,id2]
    temp <- merge(apply(temp1,1,mean),apply(temp2,1,mean),by=0)
    plot <- ggplot(temp, aes(x=x, y=y)) + geom_point(size=0.75,alpha=0.8,color="#D6404E") + geom_rug() + theme_minimal() +
      stat_cor(aes(label = paste(after_stat(r.label))),digits = 3,method = "pearson", label.x = range(temp1[,1])[2]*0.5, label.y = range(temp1[,1])[2]*0.15, size=4)+
      labs(title=paste0(RNA),x="",y="")
    plot_list[[test]][[RNA]]<-plot
  }
  plot<-ggarrange(plot_list[[test]][["lncRNA"]],plot_list[[test]][["mRNA"]],plot_list[[test]][["miRNA"]],
                  plot_list[[test]][["piRNA"]],plot_list[[test]][["snRNA"]],plot_list[[test]][["snoRNA"]],
                  plot_list[[test]][["tsRNA"]],plot_list[[test]][["cfRNA"]],
                  ncol = 3, nrow = 3)
  plot<-annotate_figure(plot,
                  left = text_grob(paste0("Log2(RPM) of ",group1," (N=",length(id1),")"),rot = 90),
                  bottom=text_grob(paste0("Log2(RPM) of ",group2," (N=",length(id2),")")))
  ggsave(paste0("Scatterplot/",test,"-each RNAs.pdf"), width = 9, height = 9)
  ggsave(paste0("Scatterplot/",test,"-each RNAs.png"), width = 9, height = 9, units = "in", bg="white")
  plot_cfRNA<-plot_list[[test]][["cfRNA"]]+labs(x=paste0("Log2(RPM) of ",group1," (N=",length(id1),")"), y = paste0("Log2(RPM) of ",group2," (N=",length(id2),")"),title = paste0(test))
  plot_cfRNA_list[[test]]<-plot_cfRNA
  ggsave(paste0("Scatterplot/",test,"-cfRNA.pdf"), width = 3, height = 3, units="in")
  ggsave(paste0("Scatterplot/",test,"-cfRNA.png"), width = 3, height = 3, units = "in", bg="white")
}

plot_cfRNA<-ggarrange(plot_cfRNA_list[[1]],plot_cfRNA_list[[2]],plot_cfRNA_list[[3]],
                      plot_cfRNA_list[[4]],plot_cfRNA_list[[5]],plot_cfRNA_list[[6]],
                      ncol = 3, nrow = 2)
ggsave(paste0("Scatterplot/All test-cfRNAs.pdf"), width = 9, height = 6)
ggsave(paste0("Scatterplot/All test-cfRNAs.png"), width = 9, height = 6, units = "in", bg="white")