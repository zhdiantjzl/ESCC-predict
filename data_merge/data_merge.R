rm(list = ls())
#引用包
library(limma)
library(sva)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(tidyverse)
files=c("TCGA.txt","GSE53622.txt","GSE53624.txt")       
clin=read.table("clinicaldata.txt",header = T,sep = "\t")

#获取交集基因
geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile, header=T, sep="\t",check.names=F)
  ttt=intersect(clin$ID,colnames(rt))
  rt=rt[,c(colnames(rt)[1],ttt)]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)
table(clin$cohort)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  ttt=intersect(clin$ID,colnames(rt))
  rt=rt[,c(colnames(rt)[1],ttt)]
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))

  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  if(header[1] != "TCGA"){
    rt=normalizeBetweenArrays(rt)
  }
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}
dim(allTab)
#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)

library(FactoMineR)
library(factoextra)
library(pca3d) # 加载R包
pca <- prcomp(t(allTab),scale. = TRUE) # 使用R自带的主成分分析函数
pheno<-data.frame(ID=colnames(allTab))
pheno$cancer<-pheno$ID
pheno[1:81,2]<-"TCGA"
pheno[82:141,2]<-"GSE53622"
pheno[142:260,2]<-"GSE53624"

table(pheno$cancer)
ddb.pca <- PCA(t(allTab), graph = FALSE)

pdf(file = "去除批次效应前.pdf",height=8,width=8)
fviz_pca_ind(ddb.pca,
             geom.ind = "point", # 只显示点
             pointsize =1, # 点的大小
             pointshape = 25, # 点的形状
             fill.ind = pheno$cancer, # 分组颜色
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # 增加置信椭圆
             legend.title = "Groups", # 图例标题
             title="") +
  theme_bw() + # 和ggplot2对接进行美化
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )
dev.off()



library(FactoMineR)
library(factoextra)
ddb.pca <- PCA(t(outTab), graph = FALSE)
pdf(file = "去除批次效应后.pdf",height=8,width=8)
fviz_pca_ind(ddb.pca,
             geom.ind = "point", # 只显示点
             pointsize =1, # 点的大小
             pointshape = 25, # 点的形状
             fill.ind = pheno$cancer, # 分组颜色
             palette =  c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # 增加置信椭圆
             legend.title = "Groups", # 图例标题
             title="") +
  theme_bw() + # 和ggplot2对接进行美化
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )

dev.off()


outTab=as.data.frame(outTab)
# colnames(outTab)=str_replace_all(colnames(outTab),"TCGA_TCGA_","TCGA_")
# colnames(outTab)=str_replace_all(colnames(outTab),"CGGA_301_","")
# colnames(outTab)=str_replace_all(colnames(outTab),"CGGA_325_","")
# colnames(outTab)=str_replace_all(colnames(outTab),"CGGA_693_","")
# colnames(outTab)=str_replace_all(colnames(outTab),"Rembrandt__475_","")
# colnames(outTab)=str_replace_all(colnames(outTab),"GSE13041_","")
clin$ID=paste0(clin$cohort,"_",clin$ID)
clin=dplyr::filter(clin,ID %in% colnames(outTab))
#$clin=clin[,-1]
clin=dplyr::select(clin,ID,everything())
save(outTab,clin,file = "DATA.RDATA")
