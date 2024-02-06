library(tidyverse)
library(survminer)
library(survival)
library(ggplot2)
library(ggpubr)
library(DescTools)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time
library(ComplexHeatmap)
#library(forestplot)

####Data preparation----
meta<-read.csv("../../Quality controls/BrainTumor_Summary_good_samples.csv", header = T)
meta <- meta[grep("CSF",meta$Group),]##
table(meta$Group)
table(meta$death)
case<-meta[meta$Group=="GLI_CSF",]
riskscore<-read.csv("../CSF&plasma/riskscores.csv",row.names = 1)

##Grade----
case$WHO<-ifelse(case$Grade>2,"high","low")
case$WHO<- factor(case$WHO, levels = c("low","high"))
highdf=case[case$WHO=="high",]
nhigh<-nrow(highdf)
nhighevent<-nrow(highdf[highdf$death==1,])
lowdf=case[case$WHO=="low",]
nlow<-nrow(lowdf)
nlowevent<-nrow(lowdf[lowdf$death==1,])
cox <- coxph(Surv(OS, death) ~ WHO, data=case)
coef<-round(cox$coefficients,2)
HR <- round(exp(cox$coefficients),2)
pvalue <- summary(cox)$sctest[3]
OSfit <- survfit(Surv(OS, death) ~ WHO, data=case)
OS_plot <- ggsurvplot(OSfit, pval=TRUE, censor = TRUE,
                      xlab = "Time in months",title=paste("WHO grade"),
                      legend="bottom", legend.title = "",
                      #risk.table = T, risk.table.fontsize=5, risk.table.y.text.col=FALSE,
                      font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
OS_plot$plot <- OS_plot$plot+scale_color_manual(values = c("#4E62AB","#D6404E"))+
  annotate("text", x = 7, y = 0.1, label = paste0("high = ",nhighevent," death in ",nhigh," patients"),colour = "#D6404E",size=5)+annotate("text", x = 7, y = 0, label =  paste0("low = ",nlowevent," death in ",nlow," patients"),colour = "#4E62AB",size=5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25))
OS_plot$plot <- OS_plot$plot+annotate("text", x = 8, y = 0.2, label = paste("Hazard ratio = ",HR,sep = ""),size=5,hjust = 0)
png(file=paste("WHO grade.png",sep = ""),width = 480, height = 480, units = "px")
print(OS_plot)
dev.off()

##Age----
median(case$age)
case$Age<-ifelse(case$age>45.5,">45.5","<45.5")
case$Age<- factor(case$Age, levels = c("<45.5",">45.5"))
highdf=case[case$Age==">45.5",]
nhigh<-nrow(highdf)
nhighevent<-nrow(highdf[highdf$death==1,])
lowdf=case[case$Age=="<45.5",]
nlow<-nrow(lowdf)
nlowevent<-nrow(lowdf[lowdf$death==1,])
cox <- coxph(Surv(OS, death) ~ Age, data=case)
coef<-round(cox$coefficients,2)
HR <- round(exp(cox$coefficients),2)
pvalue <- summary(cox)$sctest[3]
OSfit <- survfit(Surv(OS, death) ~ Age, data=case)
OS_plot <- ggsurvplot(OSfit, pval=TRUE, censor = TRUE,
                      xlab = "Time in months",title=paste("Age"),
                      legend="bottom", legend.title = "",
                      #risk.table = T, risk.table.fontsize=5, risk.table.y.text.col=FALSE,
                      font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
OS_plot$plot <- OS_plot$plot+scale_color_manual(values = c("#4E62AB","#D6404E"))+
  annotate("text", x = 7, y = 0.1, label = paste0("high = ",nhighevent," death in ",nhigh," patients"),colour = "#D6404E",size=5)+annotate("text", x = 7, y = 0, label =  paste0("low = ",nlowevent," death in ",nlow," patients"),colour = "#4E62AB",size=5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.25))
OS_plot$plot <- OS_plot$plot+annotate("text", x = 8, y = 0.2, label = paste("Hazard ratio = ",HR,sep = ""),size=5,hjust = 0)
png(file=paste("Age.png",sep = ""),width = 480, height = 480, units = "px")
print(OS_plot)
dev.off()

###multivariates----
df<-merge(riskscore,case[,c("clinicalID","WHO","Age")],by.x = 0,by.y = "clinicalID")
df$CSFrisk <- ifelse (df$csf > quantile(df$csf,0.6) , "high" , "low")
df$CSFrisk<- factor(df$CSFrisk, levels = c("low","high"))
df$plasmarisk <- ifelse (df$plasma > quantile(df$plasma,0.6) , "high" , "low")
df$plasmarisk<- factor(df$plasmarisk, levels = c("low","high"))
df<-df[,c(1,5:10)]

cox <- coxph(Surv(OS, death) ~ WHO + Age + CSFrisk, data=df)
cox
coef<-round(cox$coefficients,2)
HR <- round(exp(cox$coefficients),2)
lowerCI<-round(summary(cox)$conf.int[,"lower .95"],2)
upperCI<-round(summary(cox)$conf.int[,"upper .95"],2)
summary<-data.frame(coef=coef,HR=HR,lowerCI=lowerCI,upperCI=upperCI,pvalue=pvalue)
summary$variance <- (summary$upperCI - summary$lowerCI)^2
summary$logHR<-log2(summary$HR)
summary$logupperCI<-log2(summary$upperCI)
summary$loglowerCI<-log2(summary$lowerCI)
forestplot <- ggplot(data=summary, aes(y=reorder(rownames(summary),HR), x=logHR, xmin=loglowerCI, xmax=logupperCI, size = variance)) +
  geom_pointrange() + scale_size_continuous(range = c(0.75, 0.25))+
  geom_vline(xintercept=0, lty=2) + 
  xlim(-2.1,12) + xlab("log2 Hazard ratio") + ylab("") +  ggtitle("")+
  geom_text(data = summary, aes(y = reorder(rownames(summary),HR), label = paste0(HR, " [", lowerCI, "-", upperCI, "]")), x = 7, hjust = 0, size = 3.5) +
  theme_classic()+theme(legend.position = "none",axis.text=element_text(size=10,colour="black"))
print(forestplot)
ggsave("multivariates forest plot.png", width=6, height=3, units="in", bg="white")
ggsave("multivariates forest plot.pdf", width=6, height=3, units="in")
