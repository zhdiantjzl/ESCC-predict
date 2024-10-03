rm(list = ls())



load("All cohorts results.RDATA")
load("DATA.RDATA")

# clin=read.table("clinicaldata.txt",header = T,sep = "\t")

#加载包
library(tidyverse)
library(survminer)
library(survival)
library(survivalROC)
library(RColorBrewer)
library(forestplot)
colnames(clin)
clin$gender[clin$gender=="female"]=0
clin$gender[clin$gender=="male"]=1

clin$stage[clin$stage=="stage I"]=1
clin$stage[clin$stage=="stage II"]=2
clin$stage[clin$stage=="stage III"]=3
clin$stage[clin$stage=="stage IV"]=4



rs_merge<-distinct(rs_merge,ID,.keep_all = T) 


table(clin$cohort)
colnames(clin)
clin1=dplyr::filter(clin,cohort == "merge")
clin1=dplyr::select(clin,ID,gender,age,stage)
aaa=na.omit(rs_merge)
aaa= inner_join(rs_merge,clin1,by="ID")


rt=aaa
str(rt)
rt$stage=as.numeric(rt$stage)




#单因素独立预后分析
uniTab=data.frame()
for(i in colnames(rt[,4:ncol(rt)])){
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
uniTab=dplyr::filter(uniTab,id %in% c("RS","gender","age","stage"))
write.table(uniTab,file="test.uniCox.txt",sep="\t",row.names=F,quote=F)


############绘制森林图函数############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 10,height = 8)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
############绘制森林图函数############

bioForest(coxFile="test.uniCox.txt",forestFile=paste0("merge.test.uniForest.pdf"),forestCol="red")



