library(tidyverse)
library(survminer)
library(survival)
library(ggplot2)
library(ggpubr)
library(DescTools)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time

####Data preparation----
meta<-read.csv("../../Quality controls/BrainTumor_Summary_good_samples.csv", header = T)
table(meta$Group)
table(meta$death)

csf_change<-read.csv("../CSF_gliomas/csf_goi_change.csv", row.names = 1, check.names = FALSE)
csf_change<-merge(meta[,c("seqID","clinicalID","death","OS")],csf_change,by.x = "seqID",by.y = 0)
colnames(csf_change)[5:ncol(csf_change)]<-paste0("csf_",colnames(csf_change)[5:ncol(csf_change)])
csf_change<-csf_change[,-1]
plasma_change<-read.csv("../plasma_gliomas/plasma_goi_change.csv", row.names = 1, check.names = FALSE)
plasma_change<-merge(meta[,c("seqID","clinicalID","death","OS")],plasma_change,by.x = "seqID",by.y = 0)
colnames(plasma_change)[5:ncol(plasma_change)]<-paste0("plasma_",colnames(plasma_change)[5:ncol(plasma_change)])
plasma_change<-plasma_change[,-1]
merge_change<-merge(csf_change,plasma_change,by="clinicalID")
merge_change<-merge_change[,-grep(".y",colnames(merge_change))]
colnames(merge_change)<-gsub("\\.x","",colnames(merge_change))

####risk scores----
riskscores <- matrix(ncol = 3, nrow = nrow(merge_change))
rownames(riskscores)<-merge_change$clinicalID
colnames(riskscores)<-c("csf","plasma","csf&plasma")
merge_change<-merge_change[OrderMixed(rownames(merge_change)),]
riskscores<-riskscores[OrderMixed(rownames(riskscores)),]
##csf
change<-merge_change[,c(1,grep("csf",colnames(merge_change)))]
rownames(change)<-change[,1];change<-change[,-1]
weight <- 0
for (row in 1:nrow(change)) {
  for (col in 1:ncol(change))  {
    if (change[row,col]=="gain" | change[row,col]=="loss") {weight <- weight + 1}
    if (change[row,col]=="non")                                     {weight <- weight}
  }
  riskscores[row,"csf"] <- weight
  weight <- 0
}
##plasma
change<-merge_change[,c(1,grep("plasma",colnames(merge_change)))]
rownames(change)<-change[,1];change<-change[,-1]
weight <- 0
for (row in 1:nrow(change)) {
  for (col in 1:ncol(change))  {
    if (change[row,col]=="gain" | change[row,col]=="loss") {weight <- weight + 1}
    if (change[row,col]=="non")                                     {weight <- weight}
  }
  riskscores[row,"plasma"] <- weight
  weight <- 0
}
##csf&plasma
change<-merge_change[,c(1,grep("csf",colnames(merge_change)),grep("plasma",colnames(merge_change)))]
rownames(change)<-change[,1];change<-change[,-1]
weight <- 0
for (row in 1:nrow(change)) {
  for (col in 1:ncol(change))  {
    if (change[row,col]=="gain" | change[row,col]=="loss") {weight <- weight + 1}
    if (change[row,col]=="non")                            {weight <- weight}
  }
  riskscores[row,"csf&plasma"] <- weight
  weight <- 0
}

riskscores<-data.frame(riskscores)
riskscores<-merge(riskscores,merge_change[,c("clinicalID","death","OS")],by.x = 0, by.y = "clinicalID")
rownames(riskscores)<-riskscores[,1];riskscores<-riskscores[,-1]
write.csv(riskscores,"riskscores.csv")
riskscores$vital<-ifelse(riskscores$death==1,"Dead","Alive")

ggplot(riskscores, aes(x = reorder(rownames(riskscores), -csf.plasma), y = csf.plasma)) +
  geom_bar(aes(fill = 'CSF'), stat = 'identity') +
  geom_bar(aes(y = plasma, fill = 'Plasma'), stat = 'identity') +
  labs(x = '', y = 'Risk scores', title = "", fill = 'Type') +
  geom_tile(aes(y = 25, fill = vital), width = 1, height = 0.5) +
  scale_fill_manual(values = c("CSF" = "#D6404E", "Plasma" = "#FB9A99", 'Alive' = 'grey', 'Dead' = 'black'),
                    breaks = c("CSF", "Plasma", "Alive", "Dead")) +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.title = element_blank(), legend.position = "bottom",
        axis.line = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("Combined csf and plasma risk scores.png"),bg="white",width = 6, height = 6)
ggsave(paste0("Combined csf and plasma risk scores.pdf"),width = 6, height = 6)

####p value----
colnames(riskscores)[1]<-"CSF 23-cfRNA scores"
colnames(riskscores)[2]<-"plasma 4-cfRNA scores"
colnames(riskscores)[3]<-"CSF&plasma 27-cfRNA scores"
combinations<-colnames(riskscores)[1:3]
median_survival<-list()
for (test in combinations) {
  df<-riskscores[,c(test,"death","OS")]
  colnames(df)[1]<-"scores"
  df$risk <- ifelse (df$scores > quantile(df$scores,0.6) , "high" , "low")
  df$risk<- factor(df$risk, levels = c("low","high"))
  highdf=df[df$risk=="high",]
  nhigh<-nrow(highdf)
  nhighevent<-nrow(highdf[highdf$death==1,])
  lowdf=df[df$risk=="low",]
  nlow<-nrow(lowdf)
  nlowevent<-nrow(lowdf[lowdf$death==1,])
  cox <- coxph(Surv(OS, death) ~ risk, data=df)
  coef<-round(cox$coefficients,2)
  HR <- round(exp(cox$coefficients),2)
  pvalue <- summary(cox)$sctest[3]
  OSfit <- survfit(Surv(OS, death) ~ risk, data=df)
  median_survival[[test]]<-data.frame(summary(OSfit)$table)
  OS_plot <- ggsurvplot(OSfit, pval=TRUE, censor = TRUE,
                        xlab = "Time in months",title=paste(test),
                        legend="bottom", legend.title = "",
                        #risk.table = T, risk.table.fontsize=5, risk.table.y.text.col=FALSE,
                        font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
  OS_plot$plot <- OS_plot$plot+scale_color_manual(values = c("#4E62AB","#D6404E"))+
    annotate("text", x = 7, y = 0.1, label = paste0("high = ",nhighevent," death in ",nhigh," patients"),colour = "#D6404E",size=5)+annotate("text", x = 7, y = 0, label =  paste0("low = ",nlowevent," death in ",nlow," patients"),colour = "#4E62AB",size=5)+
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25))
  if (HR<100){
    OS_plot$plot <- OS_plot$plot+annotate("text", x = 8, y = 0.2, label = paste("Hazard ratio = ",HR,sep = ""),size=5,hjust = 0)
  }
  ggsave(filename = paste(test,"Multiple genes Kaplan-Meier.png"),
         plot = OS_plot$plot,
         width = 18,
         height = 18,
         units = 'cm')
  ggsave(filename = paste(test,"Multiple genes Kaplan-Meier.pdf"),
         plot = OS_plot$plot,
         width = 18,
         height = 18,
         units = 'cm')
}

summary(df[df$risk=="high","scores"])
summary(df[df$risk=="low","scores"])
