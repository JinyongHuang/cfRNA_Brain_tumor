library(tidyverse)
library(caret)#Sample partition, Support Vector Machine 
library(glmnet)#LASSO logistic regression
library(randomForest)#RF
library(MASS)#LDA
library(ROCR)#ROC
library(pROC)#ROC
library(reshape2)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time

####Sample ID----
meta<-read.csv("../../Quality controls/BrainTumor_Summary_good_samples.csv", header = T)
meta <- meta[grep("CSF",meta$Group),]##
table(meta$Group)
meta$group<-ifelse(meta$Group=="GLI_CSF","Gliomas","others")##
###Discovery cohort
dis_case<-subset(meta,Group=="GLI_CSF")##
dis_control<-subset(meta,Group!="GLI_CSF")##
dis_ID<-c(dis_case$seqID,dis_control$seqID)
dis_type<-factor(c(dis_case$group,dis_control$group),levels = c("Gliomas","others"))##
ndiscovery<-length(dis_ID)

####cfRNA of interest 
goi_cfRNA<-readRDS("../../DESeq2_One_vs_One_csf/CSF_intersection/ploted_goi_csf_Gliomas.rds")
####Data for machine learning----
cfRNA<-read.csv("../../Quality controls/BrainTumor_cfRNA_log2RPM_good_samples.csv", row.names = 1)
cfRNA<-cfRNA[goi_cfRNA,dis_ID]
cfRNA<-t(cfRNA)

####Filter by LASSO logistic regression---
goi_N_list<-list()
condition<-c(0,1,5,10,20,30,40,50,60,70,80,90,100)
discovery_x <- cfRNA
discovery_y <- dis_type
seeds<-1:100 #repeat 100 times
registerDoParallel(cl<-makeCluster(10))
coef<-foreach(seed=seeds, .combine="cbind",.packages=c("tidyverse","caret","glmnet")) %dopar% {
  set.seed(seed)
  print(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.6, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  ##Fit
  cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 1, nfolds = 10)
  LASSO <- glmnet(training_x, training_y, family = "binomial", alpha = 1, lambda = cvfit$lambda.1se)
  temp <- data.frame(matrix(coef(LASSO)))
  temp
}
stopCluster(cl)
rownames(coef)<-c("Intercept",colnames(cfRNA))
colnames(coef)<-1:100
for (N in condition) {
  coef_filter <-data.frame(apply(coef[-1,], 1, function(x) length(which(x!=0))))
  coef_filter$name<-rownames(coef_filter)
  coef_filter <-coef_filter[coef_filter$apply.coef..1.....1..function.x..length.which.x....0... >= N,]
  goi_N_list[[as.character(N)]] <-rownames(coef_filter)
}

dir.create(paste0(getwd(),"/goi"), showWarnings = FALSE)
for (N in condition) {
    goi=goi_N_list[[as.character(N)]]
    writeLines(goi,paste0("goi/",N,"-name.txt"))
}  
####All in one----
All<-list()
seeds<-1:99 #repeat 99 times
for (N in condition) {
  print(paste("Start working on",N))
  goi=goi_N_list[[as.character(N)]]
  if (length(goi)<=1) break 
  subset<-cfRNA[,goi]
  discovery_x <- subset[dis_ID,]
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
  
  All[[as.character(N)]]<-rbind(LR_df,RF_df,SVM_df)
}

####AUC box plot
dir.create(paste0(getwd(),"/AUC boxplots"), showWarnings = FALSE)
for (N in as.character(condition)) {
  n=length(goi_N_list[[N]])
  if (n<2) next
  df<-data.frame(All[[N]][,3:4])
  df$method<-c(rep("LR",99),rep("RF",99),rep("SVM",99))
  df<-melt(df)
  df$variable<-factor(c(rep("Training",297),rep("Testing",297)),levels = c("Training","Testing"))
  ggplot(df, aes(x=variable, y= value,colour=variable))+geom_boxplot(outlier.shape = NA)+  geom_jitter(shape=16, position=position_jitter(0.2)) + facet_wrap(~method, scale="free")+
    scale_colour_brewer(palette="Set1")+labs(title=paste0("CSF - Gliomas vs. Others - ",n, " cfRNA") ,x="", y = "AUC in 100 iterations")+ theme_classic()+
    theme(axis.text.y = element_text(colour = "black",size = 14), axis.text.x = element_text(colour = "black",size = 14),
          strip.text.x = element_text(colour = "black",size = 14,face = "bold"), legend.position = "none", 
          axis.title = element_text(size = 16), plot.title = element_text(size = 16,hjust = 0.5,face = "bold"), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  ggsave(paste0("AUC boxplots/Machine learning-LASSO selected-",N," in 100 repeats.pdf"), width = 7, height = 4)
  ggsave(paste0("AUC boxplots/Machine learning-LASSO selected-",N," in 100 repeats.png"), width = 7, height = 4, bg="white")
}

####median ROC&AUC----
##LASSO5in100
df<-data.frame(All[["5"]][,3:4])
df<-df[1:99,]
df<-df[order(df$AUC_testing),]
seed<-as.numeric(word(rownames(df)[50],start=-1,sep = "\\."))
goi=goi_N_list[["5"]]
dis_type<-factor(c(dis_case$group,dis_control$group),levels = c("Gliomas","others"))##
discovery_x <- cfRNA[dis_ID,goi]
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
saveRDS(LR_roc,"CSF_Gliomas_vs_others_p=0.6_LASSO5_training_ROC.rds")
saveRDS(LR_AUC,"CSF_Gliomas_vs_others_p=0.6_LASSO5_training_AUC.rds")
#testing
LR_prob <- predict(Ridge, testing_x, type = "response")[,1]
LR_roc <- plot.roc(testing_y,LR_prob)
LR_pred <- prediction(LR_prob, testing_y)
LR_AUC <- round(attr(performance(LR_pred, "auc"), "y.values")[[1]],3)
saveRDS(LR_roc,"CSF_Gliomas_vs_others_p=0.6_LASSO5_testing_ROC.rds")
saveRDS(LR_AUC,"CSF_Gliomas_vs_others_p=0.6_LASSO5_testing_AUC.rds")

##
summary<-as.data.frame(All[["5"]])
summary$method<-c(rep("LR",99),rep("RF",99),rep("SVM",99))
summary<-summary[,c("AUC_testing","method")]
summary(summary[summary$method=="LR","AUC_testing"])
summary(summary[summary$method=="RF","AUC_testing"])
summary(summary[summary$method=="SVM","AUC_testing"])
