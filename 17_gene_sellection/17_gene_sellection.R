#加载包
library(tidyverse)
library(ggpubr)
library(survminer)
library(survival)
library(survivalROC)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(ReactomePA)
library(igraph)
library(ggraph)
library(ggradar)
library(tuneR)
library(limma)
library(stringr)
library(RColorBrewer)
library(forestplot)
library(fmsb)
library(circlize)
library(ggsci)
library(parallel)
library(maftools)
library(circlize)
library(ggsci)
library(parallel)
library(patchwork)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
rm(list = ls())

load("DATA.RDATA")

ls()

loaded_data <- load("DATA.RDATA")
loaded_data 
#write.csv(clin,"clin_output.csv",row.names = FALSE) 没有列名
#write.csv(outTab,"outTab_output.csv",row.names = FALSE)   没有列名
write.csv(clin,"clin_output.csv",row.names = TRUE,col.names = TRUE)
write.csv(outTab,"outTab_output.csv",row.names = TRUE,col.names = TRUE)



genename<-rownames(outTab)
aaa=as.data.frame(t(outTab))
aaa$ID=rownames(aaa)
aaa=dplyr::select(aaa,ID,everything())

OSdata=inner_join(clin,aaa,by="ID")
#OS
Index<-sort(unique(clin$cohort))

library(survival)   
library(tidyverse)
OSdata$OS.time=OSdata$OS.time/365  

for (i in 1:length(Index)) {
  rt<-OSdata%>%
    dplyr::select(ID,cohort,OS,OS.time,genename)%>%
    dplyr::filter(cohort == Index[i])%>%
    dplyr::select(-cohort)

  # rt1<-dplyr::select(rt,-"OS",-"OS.time")
  # rownames(rt1)<-rt1[,1]
  # rt1<-rt1[,-1]
  # rt1<-as.data.frame(t(rt1))
  # rt1<-rt1[-which(rowSums(rt1 - apply(rt1, 2,mean) == 0) == length(colnames(rt1))),]
  # rt1<-as.data.frame(t(rt1))
  # rt1$ID<-rownames(rt1)
  # rt<-rt[,1:3]
  # rt<-inner_join(rt,rt1,by = "ID")   
  
  outTab=data.frame()
  for(gene in    colnames(rt)[4:ncol(rt)]       ){
    print(which(gene== colnames(rt)[4:ncol(rt)]  ))
    cox=coxph(Surv(OS.time, OS) ~ rt[,gene], data = rt)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP) )
  }
  outTab$HR=as.numeric(outTab$HR)
  outTab$HR.95L=as.numeric(outTab$HR.95L)
  outTab$HR.95H=as.numeric(outTab$HR.95H)
  outTab$pvalue=as.numeric(outTab$pvalue)
  write.table(outTab,file=paste0(Index[i],"单因素回归分析结果.txt"),sep="\t",row.names=F,quote=F)
  badgene=dplyr::filter(outTab,pvalue < 0.1,HR>1)#如果无交集，改变条件
  assign(paste0(Index[i],".badgene"),badgene$gene)
  goodgene=dplyr::filter(outTab,pvalue < 0.1,HR<1)
  assign(paste0(Index[i],".goodgene"),goodgene$gene)
}

rt<-OSdata%>%
  dplyr::select(ID,cohort,OS,OS.time,genename)%>%
  dplyr::select(-cohort)

# rt1<-dplyr::select(rt,-"OS",-"OS.time")
# rownames(rt1)<-rt1[,1]
# rt1<-rt1[,-1]
# rt1<-as.data.frame(t(rt1))
# rt1<-rt1[-which(rowSums(rt1 - apply(rt1, 2,mean) == 0) == length(colnames(rt1))),]
# rt1<-as.data.frame(t(rt1))
# rt1$ID<-rownames(rt1)
# rt<-rt[,1:3]
# rt<-inner_join(rt,rt1,by = "ID")   

outTab=data.frame()
for(gene in    colnames(rt)[4:ncol(rt)]       ){
  print(which(gene== colnames(rt)[4:ncol(rt)]  ))
  cox=coxph(Surv(OS.time, OS) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(gene=gene,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxP) )
}
outTab$HR=as.numeric(outTab$HR)
outTab$HR.95L=as.numeric(outTab$HR.95L)
outTab$HR.95H=as.numeric(outTab$HR.95H)
outTab$pvalue=as.numeric(outTab$pvalue)
write.table(outTab,file=paste0("merge","单因素回归分析结果.txt"),sep="\t",row.names=F,quote=F)
badgene=dplyr::filter(outTab,pvalue < 0.001,HR>1)#如果无交集，改变条件
assign(paste0("merge",".badgene"),badgene$gene)
goodgene=dplyr::filter(outTab,pvalue < 0.001,HR<1)
assign(paste0("merge",".goodgene"),goodgene$gene)

# 
# 
# 
# # 创建向量列表
# vecs <- list(GSE53622.badgene, GSE53624.badgene, TCGA.badgene)
# 
# #选择至少在4个cohort中出现的交集基因
# # 将所有向量放入列表
# 
# # 初始化结果为空向量
# result <- c()
# aaa=c(GSE53622.badgene, GSE53624.badgene, TCGA.badgene)
# aaa=unique(aaa)
# # 循环遍历向量的每个元素
# for (elem in aaa) {
#   # 统计该元素在多少个向量中出现
#   count <- sum(sapply(vecs, function(x) elem %in% x))
#   
#   # 如果在至少3个向量中出现，则加入结果向量
#   if (count == 3) {
#     result <- c(result, elem)
#   }
# }
# 
# # 打印结果
# print(result)
# badgene=result
# 
# 
# # 创建向量列表
# vecs <- list(GSE53622.goodgene, GSE53624.goodgene, TCGA.goodgene)
# 
# #选择至少在4个cohort中出现的交集基因
# # 将所有向量放入列表
# 
# # 初始化结果为空向量
# result <- c()
# aaa=c(GSE53622.goodgene, GSE53624.goodgene, TCGA.goodgene)
# aaa=unique(aaa)
# # 循环遍历向量的每个元素
# for (elem in aaa) {
#   # 统计该元素在多少个向量中出现
#   count <- sum(sapply(vecs, function(x) elem %in% x))
#   
#   # 如果在至少3个向量中出现，则加入结果向量
#   if (count == 3) {
#     result <- c(result, elem)
#   }
# }


result=c(merge.badgene,merge.goodgene)
# 打印结果
print(result)
badgene=merge.badgene
goodgene=merge.goodgene
sss=Index[3]
sss="merge"
rt=read.table(paste0(sss,"单因素回归分析结果.txt"),header=T,sep="\t",row.names=1,check.names=F)
rt=rt[c(badgene,goodgene),]
rt=arrange(rt,desc(HR))
#把含有0或者NA的行全都去掉
rt[rt == 0] <- NA
rt<-na.omit(rt)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue") 
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.05, "<0.05", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
pdf(file=paste0(sss,"_OS_forest.pdf"),
    width = 10,            
    height = 10,           
)

forestplot(tabletext, 
           zero = 1,
           lwd.zero = 2,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4,
           xlab="Hazard ratio"
)
dev.off()
OSdata=OSdata[,c(colnames(OSdata)[1:7],badgene,goodgene)]
save(goodgene,badgene,OSdata,file =  "DATA2.RDATA")
