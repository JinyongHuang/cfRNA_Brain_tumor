library(tidyverse)
library(survminer)
library(survival)
library(ggplot2)
library(ggpubr)
library(DescTools)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time
library(ComplexHeatmap)

####Data preparation----
meta<-read.csv("../../Quality controls/BrainTumor_Summary_good_samples.csv", header = T)
meta <- meta[grep("plasma",meta$Group),]##
table(meta$Group)
table(meta$death)
# case_clinical<-readLines("../CSF_gliomas/case_clinical.txt")
# case<-meta[meta$clinicalID %in% case_clinical,"seqID"]
case<-meta[meta$Group=="GLI_plasma",]$seqID
control<-meta[meta$Group!="GLI_plasma",]$seqID

sig_cfRNA<-readRDS("../../DESeq2_One_vs_One_plasma//DE_results_distribution.rds")
sig_cfRNA<-sig_cfRNA[["plasma - Gliomas vs. Healthy"]][["sig_list"]][["cfRNA"]]##
sig_cfRNA<-sig_cfRNA[sig_cfRNA$baseMean>10 & abs(sig_cfRNA$log2FoldChange)>0.5,]
goi<-rownames(sig_cfRNA)
sig_cfRNA$change<-ifelse(sig_cfRNA$log2FoldChange>0,"gain","loss")
gain<-rownames(sig_cfRNA[sig_cfRNA$change=="gain",])
loss<-rownames(sig_cfRNA[sig_cfRNA$change=="loss",])


cfRNA<-read.csv("../../Quality controls/BrainTumor_cfRNA_log2RPM_good_samples.csv", row.names = 1)
cfRNA<-cfRNA[goi,meta$seqID]
cfRNA<-data.frame(t(cfRNA),check.names = FALSE)

####Condition----
cfRNA_control<-cfRNA[control,]
summary_list<-list()
for (i in 1:length(goi)) {
  summary_list[[colnames(cfRNA_control)[i]]]<-summary(cfRNA_control[,i])
}
summary_df<-as.data.frame(do.call(rbind, summary_list))
upper_cutoff<-summary_df$`3rd Qu.`
names(upper_cutoff)<-rownames(summary_df)
lower_cutoff<-summary_df$`1st Qu.`
names(lower_cutoff)<-rownames(summary_df)

####Change in cases----
cfRNA_case<-cfRNA[case,]
change_case<-matrix(ncol=ncol(cfRNA_case),nrow=nrow(cfRNA_case))
rownames(change_case)<-rownames(cfRNA_case)
colnames(change_case)<-colnames(cfRNA_case)
for (i in 1:length(goi)) {
  for (p in 1:length(case)) {
    if (goi[i] %in% gain) {
      change_case[p,i]<-ifelse(cfRNA_case[p,i]>upper_cutoff[i],"gain","non")
    }
    if (goi[i] %in% loss) {
      change_case[p,i]<-ifelse(cfRNA_case[p,i]<lower_cutoff[i],"loss","non")
    }
  }
}
change_case<-as.data.frame(change_case)
change_case<-merge(change_case,meta[,c("seqID","death","OS")],by.x = 0, by.y = "seqID")
rownames(change_case)<-change_case[,1];change_case<-change_case[,-1]
write.csv(change_case,"plasma_goi_change.csv")
####Single RNA survival----
registerDoParallel(cl<-makeCluster(10))
goi_good<-foreach(g=1:length(goi), .combine="c", .packages=c("tidyverse","survminer","survival")) %dopar% {
  df <- data.frame(OS=change_case$OS,Death=change_case$death, change=change_case[,g])
  if (goi[g] %in% gain) {df$change<- factor(df$change, levels = c("non","gain"))}
  if (goi[g] %in% loss) {df$change<- factor(df$change, levels = c("non","loss"))}
  cox <- coxph(Surv(OS, Death) ~ change, data=df)
  HR <- round(exp(cox$coefficients),2)
  pvalue <- summary(cox)$sctest[3]
  if (HR>1 & pvalue<0.05) {
    goi[g]
  }
}
change_case_top<-change_case[,goi_good]
write.csv(change_case_top,"plasma_goi_change.csv")
writeLines(goi_good,"goi.txt")

HR_pvalue<-foreach(g=goi_good, .combine="rbind", .packages=c("tidyverse","survminer","survival")) %dopar% {
  df <- data.frame(OS=change_case$OS,Death=change_case$death, change=change_case[,g])
  if (g %in% gain) {df$change<- factor(df$change, levels = c("non","gain"))}
  if (g %in% loss) {df$change<- factor(df$change, levels = c("non","loss"))}
  cox <- coxph(Surv(OS, Death) ~ change, data=df)
  coef<-round(cox$coefficients,2)
  HR <- round(exp(cox$coefficients),2)
  lowerCI<-round(summary(cox)$conf.int[,"lower .95"],2)
  upperCI<-round(summary(cox)$conf.int[,"upper .95"],2)
  pvalue <- round(summary(cox)$sctest[3],4)
  c(coef,HR,lowerCI,upperCI,pvalue)
}
rownames(HR_pvalue)<-goi_good
colnames(HR_pvalue)[1]<-"coefficients"
colnames(HR_pvalue)[2]<-"HR"
colnames(HR_pvalue)[3]<-"lowerCI"
colnames(HR_pvalue)[4]<-"upperCI"
HR_pvalue<-data.frame(HR_pvalue)
HR_pvalue$goi<-rownames(HR_pvalue)

median_survival<-foreach(g=goi_good, .packages=c("tidyverse","survminer","survival")) %dopar% {
  df <- data.frame(OS=change_case$OS,Death=change_case$death, change=change_case[,g])
  if (g %in% gain) {df$change<- factor(df$change, levels = c("non","gain"))}
  if (g %in% loss) {df$change<- factor(df$change, levels = c("non","loss"))}
  OSfit <- survfit(Surv(OS, Death) ~ change, data=df)

  table<-data.frame(summary(OSfit)$table)
  table
}
stopCluster(cl)
names(median_survival)<-goi_good
median_survival<- lapply(median_survival, function(x) data.frame(x$rmean, row.names = rownames(x)))


dir.create(paste0(getwd(),"/Single gene Kaplan-Meier"), showWarnings = FALSE)
for (g in goi_good) {
  df <- data.frame(OS=change_case$OS,Death=change_case$death, change=change_case[,g])
  if (g %in% gain) {
    df$change<- factor(df$change, levels = c("non","gain"))
    gaindf=df[df$change=="gain",]
    ngain<-nrow(gaindf)
    ngainevent<-nrow(gaindf[gaindf$Death==1,])
    nondf=df[df$change=="non",]
    nnon<-nrow(nondf)
    nnonevent<-nrow(nondf[nondf$Death==1,])
    cox <- coxph(Surv(OS, Death) ~ change, data=df)
    HR <- round(exp(cox$coefficients),2)
    OSfit <- survfit(Surv(OS, Death) ~ change, data=df)
    OS_plot <- ggsurvplot(OSfit, pval=TRUE,conf.int = FALSE,
                          xlab = "Time in months",title=paste(g),
                          legend="bottom", legend.title = "",
                          font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
    OS_plot$plot <- OS_plot$plot+scale_color_manual(values = c("#4E62AB","#D6404E"))+
      annotate("text", x = 9, y = 0.12, label = paste0("gain = ",ngainevent," death in ",ngain," patients"),colour = "#D6404E",size=5)+annotate("text", x = 9, y = 0.04, label =  paste0("non = ",nnonevent," death in ",nnon," patients"),colour = "#4E62AB",size=5)+
      theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25))
    if (HR<100){
      OS_plot$plot <- OS_plot$plot+annotate("text", x = 0.5, y = 0.28, label = paste("Hazard ratio = ",HR,sep = ""),size=5,hjust = 0)
    }
    png(file=paste("Single gene Kaplan-Meier/",g,"-OS.png",sep = ""),width = 400, height = 400, units = "px")
    print(OS_plot)
    dev.off()
    pdf(file=paste("Single gene Kaplan-Meier/",g,"-OS.pdf",sep = ""))
    print(OS_plot)
    dev.off()
    }
  if (g %in% loss) {
    df$change<- factor(df$change, levels = c("non","loss"))
    lossdf=df[df$change=="loss",]
    nloss<-nrow(lossdf)
    nlossevent<-nrow(lossdf[lossdf$Death==1,])
    nondf=df[df$change=="non",]
    nnon<-nrow(nondf)
    nnonevent<-nrow(nondf[nondf$Death==1,])
    cox <- coxph(Surv(OS, Death) ~ change, data=df)
    HR <- round(exp(cox$coefficients),2)
    OSfit <- survfit(Surv(OS, Death) ~ change, data=df)
    OS_plot <- ggsurvplot(OSfit, pval=TRUE,conf.int = FALSE,
                          xlab = "Time in months",title=paste(g),
                          legend="bottom", legend.title = "",
                          font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
    OS_plot$plot <- OS_plot$plot+scale_color_manual(values = c("#4E62AB","#D6404E"))+
      annotate("text", x = 9, y = 0.12, label = paste0("loss = ",nlossevent," death in ",nloss," patients"),colour = "#D6404E",size=5)+annotate("text", x = 9, y = 0.04, label =  paste0("non = ",nnonevent," death in ",nnon," patients"),colour = "#4E62AB",size=5)+
      theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25))
    if (HR<100){
      OS_plot$plot <- OS_plot$plot+annotate("text", x = 0.5, y = 0.28, label = paste("Hazard ratio = ",HR,sep = ""),size=5,hjust = 0)
    }
    png(file=paste("Single gene Kaplan-Meier/",g,"-OS.png",sep = ""),width = 400, height = 400, units = "px")
    print(OS_plot)
    dev.off()
    pdf(file=paste("Single gene Kaplan-Meier/",g,"-OS.pdf",sep = ""),width = 6, height = 6)
    print(OS_plot)
    dev.off()
    }
}

####forest plot----
df<-HR_pvalue[grep("Inf",HR_pvalue$upperCI,invert = TRUE),]#remove inf sample
df$variance <- (df$upperCI - df$lowerCI)^2
df$logHR<-log2(df$HR)
df$logupperCI<-log2(df$upperCI)
df$loglowerCI<-log2(df$lowerCI)
df<-df[order(df$pvalue),] 

forestplot <- ggplot(data=df, aes(y=reorder(goi,-pvalue), x=logHR, xmin=loglowerCI, xmax=logupperCI, size = variance)) +
  geom_pointrange() + scale_size_continuous(range = c(0.75, 0.25))+
  geom_vline(xintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  xlim(-1,12) + xlab("log2 Hazard ratio") + ylab("") +  ggtitle("")+
  geom_text(data = df, aes(y = reorder(goi,-pvalue), label = paste0(HR, " [", lowerCI, "-", upperCI, "]")), x = 7.2, hjust = 0, size = 3.5) +
  geom_text(data = df, aes(y = reorder(goi,-pvalue), label = pvalue), x = 11, hjust = 0, size = 3.5) +
  theme_classic()+theme(legend.position = "none",axis.text=element_text(size=10,colour="black"))
print(forestplot)
ggsave("plasma forest plot.png", width=6.5, height=3, units="in", bg="white")
ggsave("plasma forest plot.pdf", width=6.5, height=3, units="in")

####Oncoprint----
x=1
y=1
w=1
h=1
col = c(gain="#D6404E", loss="#4E62AB")
alter_fun = list(
  gain = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                        gp = gpar(fill = col["gain"], col = "black")),
  loss = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                        gp = gpar(fill = col["loss"], col = "black"))
)
temp<-change_case_top
temp[]<-sapply(temp, function(x) gsub("non","",x))
temp <-t(temp)
sample_order <- names(sort(apply(temp, 2, function(x) sum(x==""))))
death<-meta[meta$seqID %in% case,c("seqID","death")]
table(death$seqID==colnames(temp))
death<-death$death
death<-gsub("1","dead",death)
death<-gsub("0","alive",death)
col_death <- c("dead" = "black", "alive" = "grey")
top1<- HeatmapAnnotation(`Vital status` = death, col = list(`Vital status` = col_death), show_annotation_name = T, which="column")
grade<-meta[meta$seqID %in% case,c("seqID","Grade")]
table(grade$seqID==colnames(temp))
grade<-grade$Grade
table(grade)
col_grade <- c("2" = "#FCAE91", "3" = "#FB6A4A", "4" = "#CB181D")
top2<- HeatmapAnnotation(`WHO grade` = grade, col = list(`WHO grade` = col_grade), show_annotation_name = T, which="column")

oncoPrint(temp, alter_fun = alter_fun, col = col, top_annotation = c(top1,top2), right_annotation =NULL,
          alter_fun_is_vectorized =FALSE,
          column_title = "Prognostic cfRNAs in gliomas plasma",column_order = sample_order)

####Multiple RNA survival ----
HR_pvalue<-HR_pvalue[order(HR_pvalue$pvalue),]
riskscores <- matrix(ncol = length(goi_good)-1, nrow = length(case))
rownames(riskscores)<-case
colnames(riskscores)<-paste0("Riskscores_topN_",2:length(goi_good))
change_case<-change_case[OrderMixed(rownames(change_case)),]
riskscores<-riskscores[OrderMixed(rownames(riskscores)),]

weight <- 0
for (top in 2:length(goi_good)) {
  goi_top<-rownames(HR_pvalue)[1:top]
  change_case_top<-change_case[,goi_top]
  for (row in 1:nrow(change_case_top)) {
    for (col in 1:ncol(change_case_top))  {
      if (change_case_top[row,col]=="gain" | change_case_top[row,col]=="loss") {weight <- weight + 1}
      if (change_case_top[row,col]=="non")                                     {weight <- weight}
    }
    riskscores[row,top-1] <- weight
    weight <- 0
  }
}
riskscores<-data.frame(riskscores)
riskscores<-merge(riskscores,meta[,c("seqID","death","OS")],by.x = 0, by.y = "seqID")
rownames(riskscores)<-riskscores[,1];riskscores<-riskscores[,-1]

dir.create(paste0(getwd(),"/Multiple genes Kaplan-Meier"), showWarnings = FALSE)
multiple_HR_pvalue<-data.frame(coef=rep(NA,ncol(riskscores)),HR=rep(NA,ncol(riskscores)),pvalue=rep(NA,ncol(riskscores)),row.names = colnames(riskscores))
multiple_HR_pvalue<-multiple_HR_pvalue[1:length(goi_good)-1,]
for (top in 2:length(goi_good)) {
  df<-riskscores[,c(top-1,ncol(riskscores)-1,ncol(riskscores))]
  colnames(df)[1]<-"riskscores"
  df$risk <- ifelse (df$riskscores > quantile(df$riskscores,0.5) , "high" , "low")
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
  multiple_HR_pvalue[top-1,"coef"]<-coef
  multiple_HR_pvalue[top-1,"HR"]<-HR
  multiple_HR_pvalue[top-1,"pvalue"]<-pvalue
  OSfit <- survfit(Surv(OS, death) ~ risk, data=df)
  OS_plot <- ggsurvplot(OSfit, pval=TRUE, censor = TRUE,
                        xlab = "Time in months",title=paste("Top ",top," cfRNA"),
                        legend="bottom", legend.title = "",
                        #risk.table = T, risk.table.fontsize=5, risk.table.y.text.col=FALSE,
                        font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
  OS_plot$plot <- OS_plot$plot+scale_color_manual(values = c("#4E62AB","#D6404E"))+
    annotate("text", x = 7, y = 0.1, label = paste0("high = ",nhighevent," death in ",nhigh," patients"),colour = "#D6404E",size=5)+annotate("text", x = 7, y = 0, label =  paste0("low = ",nlowevent," death in ",nlow," patients"),colour = "#4E62AB",size=5)+
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25))
  if (HR<100){
    OS_plot$plot <- OS_plot$plot+annotate("text", x = 8, y = 0.2, label = paste("Hazard ratio = ",HR,sep = ""),size=5,hjust = 0)
  }
  png(file=paste("Multiple genes Kaplan-Meier/Top ",top," cfRNA.png",sep = ""),width = 480, height = 480, units = "px")
  print(OS_plot)
  dev.off()
}


