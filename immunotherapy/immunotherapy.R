

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
install.packages("ggpubr")


#???ð?
library(limma)
library(reshape2)
library(ggpubr)

expFile="uniSigExp.txt"      #?????????ļ?
cluFile="ARGcluster.txt"     #???͵Ľ????ļ?
setwd("C:\\Users\\benyu\\Desktop\\ESCC")       #???ù???Ŀ¼

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

#??ȡ???͵Ľ????ļ?
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
sameSample=intersect(row.names(data), row.names(cluster))
expClu=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])

##??ȡ?????????Ļ???
#sigGene=c()
#for(i in colnames(expClu)[1:(ncol(expClu)-1)]){
#	if(sd(expClu[,i])<0.001){next}
#	if(length(levels(factor(expClu[,"ARGcluster"])))>2){
#		test=kruskal.test(expClu[,i] ~ expClu[,"ARGcluster"])
#	}else{
#		test=wilcox.test(expClu[,i] ~ expClu[,"ARGcluster"])
#	}
#	pvalue=test$p.value
#	if(pvalue<0.05){
#		sigGene=c(sigGene, i)
#	}
#}
#sigGene=c(sigGene, "ARGcluster")
#expClu=expClu[,sigGene]

#??????ת????ggplot2?????ļ?
data_melted=melt(expClu, id.vars=c("ARGcluster"))
colnames(data_melted)=c("ARGcluster", "Gene", "Expression")

#????ͼ????ɫ
bioCol=c("#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data_melted[,"ARGcluster"])))]

#????????ͼ
p=ggboxplot(data_melted, x="Gene", y="Expression", color="ARGcluster",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="ARGcluster",
	     palette = bioCol,
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=ARGcluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
	      label = "p.signif")

#????????ͼ
pdf(file="boxplot.pdf", width=10, height=4)
print(p1+theme(axis.text.x = element_text(size = 8)))
print(p1)
dev.off()




