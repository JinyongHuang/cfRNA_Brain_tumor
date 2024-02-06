library(tidyverse)
library(caret)#Sample partition, Support Vector Machine 
library(glmnet)#LASSO logistic regression
library(randomForest)#RF
library(MASS)#LDA
library(ROCR)#ROC
library(pROC)#ROC
library(reshape2)
library(DescTools)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time

####Sample ID----
meta<-read.csv("../../Quality controls/BrainTumor_Summary_good_samples.csv", header = T)
table(meta$Group)
csf_sample<-meta[grep("CSF",meta$Group),]
plasma_sample<-meta[grep("plasma",meta$Group),]
meta<-meta[duplicated(meta$clinicalID),]
meta$Group<-word(meta$Group,start=1,sep = "_")


####Data for machine learning----
cfRNA<-read.csv("../../Quality controls/BrainTumor_cfRNA_log2RPM_good_samples.csv", row.names = 1)
csf_goi<-readLines("../CSF_Gliomas_vs_others_p=0.6/goi/5-name.txt")##
ncsf<-length(csf_goi)
csf_cfRNA<-cfRNA[csf_goi,csf_sample$seqID]
rownames(csf_cfRNA)<-paste0("CSF_",rownames(csf_cfRNA))
csf_sample<-csf_sample[OrderMixed(csf_sample$seqID),]
csf_cfRNA<-csf_cfRNA[,OrderMixed(colnames(csf_cfRNA))]
colnames(csf_cfRNA)=csf_sample$clinicalID
csf_cfRNA<-t(csf_cfRNA)
plasma_goi<-readLines("../plasma_Gliomas_vs_others_p=0.6/goi/1-name.txt")##
nplasma<-length(plasma_goi)
plasma_cfRNA<-cfRNA[plasma_goi,plasma_sample$seqID]
rownames(plasma_cfRNA)<-paste0("plasma_",rownames(plasma_cfRNA))
plasma_sample<-plasma_sample[OrderMixed(plasma_sample$seqID),]
plasma_cfRNA<-plasma_cfRNA[,OrderMixed(colnames(plasma_cfRNA))]
colnames(plasma_cfRNA)=plasma_sample$clinicalID
plasma_cfRNA<-t(plasma_cfRNA)
cfRNA<-merge(csf_cfRNA,plasma_cfRNA,by=0)
rownames(cfRNA)<-cfRNA[,1];cfRNA<-cfRNA[,-1]
cfRNA<-as.matrix(cfRNA)

###Discovery cohort
meta$group<-ifelse(meta$Group=="GLI","Gliomas","others")##
dis_case<-subset(meta,Group=="GLI")##
dis_control<-subset(meta,Group!="GLI")##
dis_ID<-as.character(c(dis_case$clinicalID,dis_control$clinicalID))

####All in one----
dis_type<-factor(c(dis_case$group,dis_control$group),levels = c("Gliomas","others"))##
seeds<-1:99 #repeat 99 times
discovery_x <- cfRNA[dis_ID,]
discovery_y <- dis_type
###Ridge logistic regression
registerDoParallel(cl<-makeCluster(10))
result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
  set.seed(seed) #for loop
  splitSample <- createDataPartition(discovery_y, p = 0.6, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  ##Fit
  cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
  LR <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
  ##Training set
  pred.class <- predict(LR, training_x, type = "class")
  ACC_training<-round(100* mean(pred.class[,1] == training_y),2)
  prob <- predict(LR, training_x, type = "response")
  pred <- prediction(prob, training_y)
  AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  roc<-roc(training_y,as.numeric(prob))
  coord_list <- coords(roc, x = "all")
  SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
  SenAt1Spe_training<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
  ##Held-out testing set
  pred.class <- predict(LR, testing_x, type = "class")
  ACC_testing<-round(100* mean(pred.class[,1] == testing_y),2)
  prob <- predict(LR, testing_x, type = "response")
  pred <- prediction(prob,testing_y)
  AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  roc<-roc(testing_y,as.numeric(prob))
  coord_list <- coords(roc, x = "all")
  SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
  SenAt1Spe_testing<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
  c(ACC_training,ACC_testing,AUC_training,AUC_testing,SenAt1Spe_training,SenAt1Spe_testing)
}
stopCluster(cl)
colnames(result_seed)= c('ACC_training','ACC_testing','AUC_training','AUC_testing','SenAt1Spe_training','SenAt1Spe_testing')
LR_df<-result_seed
###Random forest
registerDoParallel(cl<-makeCluster(10))
result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","randomForest","ROCR","pROC")) %dopar% {
  set.seed(seed) #for loop
  splitSample <- createDataPartition(discovery_y, p = 0.6, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  ##Fit
  RF = randomForest(x = training_x, y = training_y)
  ##Training set
  pred.class <- predict(RF, training_x, type = "class")
  ACC_training<-round(100* mean(pred.class == training_y),2)
  prob <- predict(RF, training_x, type = "prob")
  pred <- prediction(prob[,2], training_y)
  AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  roc<-roc(training_y,prob[,2])
  coord_list <- coords(roc, x = "all")
  SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
  SenAt1Spe_training<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
  ##Held-out testing set
  pred.class <- predict(RF, testing_x, type = "class")
  ACC_testing<-round(100* mean(pred.class == testing_y),2)
  prob <- predict(RF, testing_x, type = "prob")
  pred <- prediction(prob[,2], testing_y)
  AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  roc<-roc(testing_y,prob[,2])
  coord_list <- coords(roc, x = "all")
  SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
  SenAt1Spe_testing<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
  c(ACC_training,ACC_testing,AUC_training,AUC_testing,SenAt1Spe_training,SenAt1Spe_testing)
}
stopCluster(cl)
colnames(result_seed)= c('ACC_training','ACC_testing','AUC_training','AUC_testing','SenAt1Spe_training','SenAt1Spe_testing')
RF_df<-result_seed

###Support vector machine
registerDoParallel(cl<-makeCluster(10))
result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","ROCR","pROC")) %dopar% {
  set.seed(seed) #for loop
  splitSample <- createDataPartition(discovery_y, p = 0.6, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  ##Fit
  trControl <- trainControl(method="cv", number=10, repeats=NA, p = 0.6, classProbs = TRUE)# Set up Repeated k-fold Cross Validation
  SVM <- train(training_x, training_y, method="svmLinear", trControl=trControl, preProcess=c("center","scale"))
  ##Training set
  pred.class <- predict(SVM, training_x, type = "raw")
  ACC_training<-round(100* mean(pred.class == training_y),2)
  prob <- predict(SVM, training_x, type = "prob")
  pred <- prediction(prob[,2], training_y)
  AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  roc<-roc(training_y,prob[,2])
  coord_list <- coords(roc, x = "all")
  SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
  SenAt1Spe_training<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
  ##Held-out testing set
  pred.class <- predict(SVM, testing_x, type = "raw")
  ACC_testing<-round(100* mean(pred.class == testing_y),2)
  prob <- predict(SVM, testing_x, type = "prob")
  pred <- prediction(prob[,2], testing_y)
  AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  roc<-roc(testing_y,prob[,2])
  coord_list <- coords(roc, x = "all")
  SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
  SenAt1Spe_testing<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
  c(ACC_training,ACC_testing,AUC_training,AUC_testing,SenAt1Spe_training,SenAt1Spe_testing)
}
stopCluster(cl)
colnames(result_seed)= c('ACC_training','ACC_testing','AUC_training','AUC_testing','SenAt1Spe_training','SenAt1Spe_testing')
SVM_df<-result_seed
All<-rbind(LR_df,RF_df,SVM_df)


####AUC box plot
df<-data.frame(All[,3:4])
df$method<-c(rep("LR",99),rep("RF",99),rep("SVM",99))
df<-melt(df)
df$variable<-factor(c(rep("Training",297),rep("Testing",297)),levels = c("Training","Testing"))
ggplot(df, aes(x=variable, y= value,colour=variable))+geom_boxplot(outlier.shape = NA)+  geom_jitter(shape=16, position=position_jitter(0.2)) + facet_wrap(~method, scale="free")+
  scale_colour_brewer(palette="Set1")+labs(title=paste("Combine",ncsf, "csf cfRNAs &",nplasma, "plasma cfRNAs") ,x="", y = "AUC in 100 iterations")+ theme_classic()+
  theme(axis.text.y = element_text(colour = "black",size = 14), axis.text.x = element_text(colour = "black",size = 14),
        strip.text.x = element_text(colour = "black",size = 14,face = "bold"), legend.position = "none", 
        axis.title = element_text(size = 16), plot.title = element_text(size = 16,hjust = 0.5,face = "bold"), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggsave(paste("AUC boxplots-Combine",ncsf, "csf cfRNAs &",nplasma, "plasma cfRNAs.pdf"), width = 7, height = 7)
ggsave(paste("AUC boxplots-Combine",ncsf, "csf cfRNAs &",nplasma, "plasma cfRNAs.png"), width = 7, height = 7, bg="white")


####median ROC
df<-data.frame(All[,3:4])
df<-df[1:99,]
df<-df[order(df$AUC_testing),]
seed<-as.numeric(word(rownames(df)[50],start=-1,sep = "\\."))
dis_type<-factor(c(dis_case$group,dis_control$group),levels = c("Gliomas","others"))##
discovery_x <- cfRNA[dis_ID,]
discovery_y <- dis_type
set.seed(seed)
splitSample <- createDataPartition(discovery_y, p = 0.6, list = FALSE)
training_x <- discovery_x[splitSample,]
training_y <- discovery_y[splitSample]
testing_x <- discovery_x[-splitSample,]
testing_y <- discovery_y[-splitSample]
cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
Ridge <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
#training
LR_prob <- predict(Ridge, training_x, type = "response")[,1]
LR_roc <- plot.roc(training_y,LR_prob)
LR_pred <- prediction(LR_prob, training_y)
LR_AUC <- round(attr(performance(LR_pred, "auc"), "y.values")[[1]],3)
saveRDS(LR_roc,"CSF&plasma_Gliomas_vs_others_p=0.6_training_ROC.rds")##
saveRDS(LR_AUC,"CSF&plasma_Gliomas_vs_others_p=0.6_training_AUC.rds")##
#testing
LR_prob <- predict(Ridge, testing_x, type = "response")[,1]
LR_roc <- plot.roc(testing_y,LR_prob)
LR_pred <- prediction(LR_prob, testing_y)
LR_AUC <- round(attr(performance(LR_pred, "auc"), "y.values")[[1]],3)
saveRDS(LR_roc,"CSF&plasma_Gliomas_vs_others_p=0.6_testing_ROC.rds")##
saveRDS(LR_AUC,"CSF&plasma_Gliomas_vs_others_p=0.6_testing_AUC.rds")##

##
summary<-as.data.frame(All)
summary$method<-c(rep("LR",99),rep("RF",99),rep("SVM",99))
summary<-summary[,c("AUC_testing","method")]
summary(summary[summary$method=="LR","AUC_testing"])
summary(summary[summary$method=="RF","AUC_testing"])
summary(summary[summary$method=="SVM","AUC_testing"])
