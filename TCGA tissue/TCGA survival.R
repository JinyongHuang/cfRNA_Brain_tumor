library(tidyverse)
library(survminer)
library(survival)
library(ggplot2)
library(ggpubr)
library(openxlsx)

#loop---
genelist_names = list.files("gene_data/")
genelist_names<-gsub(".csv","",genelist_names)
genelist = list()
for (genelist_name in genelist_names) {
  temp = read.csv(paste0("gene_data/",genelist_name,".csv"))
  temp <- separate(temp, Description, into = c("Age", "Sex","Race","Vital","Survival"), sep = ",",fill = "left")
  temp$Death<-ifelse(temp$Vital==" dead",1,0)
  temp$Survival<-gsub(" days","",temp$Survival)
  temp$Survival<-as.numeric(temp$Survival)
  temp$Survival<-temp$Survival/30
  genelist[[genelist_name]]<-temp
}


#surival 
dir.create(paste0(getwd(),'/surival_plot'), showWarnings = FALSE)
for (genelist_name in genelist_names) {
  cutpoint <- surv_cutpoint(genelist[[genelist_name]], time = "Survival", event = "Death", variables = c("FPKM"))#cutoff from maxstat (maximally selected rank statistics)
  cutoff<-cutpoint[["cutpoint"]]$cutpoint
  df<-surv_categorize(cutpoint)
  df$FPKM<- factor(df$FPKM, levels = c("high","low"))
  
  highdf=df[df$FPKM=="high",]
  nhigh<-nrow(highdf)
  nhighevent<-nrow(highdf[highdf$Death==1,])
  lowdf=df[df$FPKM=="low",]
  nlow<-nrow(lowdf)
  nlowevent<-nrow(lowdf[lowdf$Death==1,])
  
  cox <- coxph(Surv(Survival, Death) ~ FPKM, data=df)
  HR <- round(exp(cox$coefficients),2)
  pvalue <- signif(summary(cox)$sctest[3],4)
  
  OSfit <- survfit(Surv(Survival, Death) ~ FPKM, data=df)
  OS_plot <- ggsurvplot(OSfit, pval=FALSE, censor = TRUE,
                        xlab = "Time in months",title=paste0(genelist_name," in TCGA glioma tissues"),
                        legend="bottom", legend.title = "",
                        font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
  OS_plot$plot <- OS_plot$plot+scale_color_manual(values = c("#D6404E","#4E62AB"))+
    annotate("text", x = 50, y = 0.6, label = paste0("p = ",pvalue),size=5,hjust = 0)+
    annotate("text", x = 50, y = 0.5, label = paste0("Hazard ratio = ",HR),size=5,hjust = 0)+
    annotate("text", x = 60, y = 0.4, label = paste0("high = ",nhighevent," death in ",nhigh," patients"),colour = "#4E62AB",size=5)+
    annotate("text", x = 60, y = 0.3, label =  paste0("low = ",nlowevent," death in ",nlow," patients"),colour = "#D6404E",size=5)+
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25))
  
  ggsave(filename = paste0("surival_plot/",genelist_name,'.png'),
         plot = OS_plot$plot,
         width = 18,
         height = 18,
         units = 'cm')
  ggsave(filename = paste0("surival_plot/",genelist_name,'.pdf'),
         plot = OS_plot$plot,
         width = 18,
         height = 18,
         units = 'cm')
}

###########
list_of_dfs <- lapply(names(genelist), function(name) {
  df <- genelist[[name]]
  names(df)[7] <- name  # Rename the 7th column
  return(df)
})
merged_df <- Reduce(function(x, y) merge(x, y[,c(1,7)], by = names(x)[1], all = TRUE), list_of_dfs)
write.xlsx(merged_df,"Supplementary Data 3.xlsx")
