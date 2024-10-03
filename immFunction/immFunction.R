
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSVA")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSEABase")

#install.packages("ggpubr")
#install.packages("reshape2")


#???ð?
library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)
expFile="symbol.txt"             #?????????ļ?
gmtFile="immune.gmt"             #???߹??????ݼ??ļ?
riskFile="risk.txt"              #?????ļ?
socreFile="immFunScore.txt"      #???߹??ܴ??ֵ??????ļ?
setwd("C:\\Users\\lenovo\\Desktop")       #???ù???Ŀ¼

#??ȡ?????????ļ????????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
	
#??ȡ???ݼ??ļ?
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
	
#ssgsea????
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#????ssGSEA score????????
normalize=function(x){
	return((x-min(x))/(max(x)-min(x)))}
#??ssGSEA score???н???
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file=socreFile, sep="\t", quote=F, col.names=F)

#ȥ????????Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
	
#??ȡ?????ļ?
risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
	
#?ϲ?????
sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(data, risk)

#?????????ع??ܻ???????ͼ
data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
	     ylab="Score",add = "none",xlab="",palette = c("blue","red") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

#????ͼƬ?ļ?
pdf(file="immFunction.pdf", width=7, height=5.5)
print(p)
dev.off()
