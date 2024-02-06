#library(xlsx)
library(openxlsx)

##create the loop to read the excel data
csv_path  = paste0(getwd(),'/box_plot')
csv.files  = list.files(csv_path)

csv_datas = list()
setwd(paste0(getwd(),'/box_plot'))
for (csv.file in csv.files) {
  csv_datas[[csv.file]] = read.xlsx(csv.file,sheetIndex = 1)

  csv_datas[[csv.file]]$type <- sapply(strsplit(csv_datas[[csv.file]]$type, " "), function(x) {
    x[1] <- paste0(toupper(substr(x[1], 1, 1)), tolower(substr(x[1], 2, nchar(x[1]))))
    if(length(x) > 1) {
      x[2:length(x)] <- tolower(x[2:length(x)])
    }
    return(paste(x, collapse = " "))
  })
  
  
}

csv_datas[[csv.file]]$type <- sapply(strsplit(csv_datas[[csv.file]]$type, " "), function(x) {
  x[1] <- paste0(toupper(substr(x[1], 1, 1)), tolower(substr(x[1], 2, nchar(x[1]))))
  if(length(x) > 1) {
    x[2:length(x)] <- tolower(x[2:length(x)])
  }
  return(paste(x, collapse = " "))
})




##boxplot
library(ggplot2)
for (i in 1:length(csv_datas)) {
  

plottitle = gsub('.xlsx','',names(csv_datas[i]))
plot = ggplot(csv_datas[[i]],aes(x = type,y = FPKM))+
  geom_boxplot(outlier.fill = NA,
               outlier.alpha = 1,
               outlier.shape = 1,
               aes(fill = type))+
  
  
  
  labs(title = plottitle)+
  theme_classic()+
  theme(axis.text.x  = element_text(angle = 45,
                                 hjust = 1,
                                 ),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.length.y = unit(0.2,'cm'),
        plot.title  = element_text(hjust = 0.5)
  )

ggsave(plot,
       filename = paste0(plottitle,'_boxplot.png'),
       width = 6,
       height = 4,
       path = getwd())

ggsave(plot,
       filename = paste0(plottitle,'_boxplot.pdf'),
       width = 6,
       height = 4,
       path = getwd())
}
plot

###########
df1<-read.xlsx("C10orf90.xlsx")
colnames(df1)[3]<-"C10orf90"
df2<-read.xlsx("DPP10.xlsx")
colnames(df2)[3]<-"DPP10"
df3<-read.xlsx("PDE4D.xlsx")
colnames(df3)[3]<-"PDE4D"
df4<-read.xlsx("PRKN.xlsx")
colnames(df4)[3]<-"PRKN"
df<-merge(df1,df2,by="Sample")
df<-merge(df,df3,by="Sample")
df<-merge(df,df4,by="Sample")
df<-df[,c("Sample","type.x","C10orf90","DPP10","PDE4D","PRKN")]
colnames(df)[2]<-"Type"
df<-df[order(df$Type),]
write.xlsx(df,"Supplementary Data 2.xlsx")
