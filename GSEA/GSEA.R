
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#???ð?
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="merge.txt"           #?????????ļ?
cluFile="ARGcluster.txt"      #???͵Ľ????ļ?
clusterType="High"               #ѡ??ͼ????չʾ?ķ???
gmtFile="c2.cp.kegg.symbols.gmt"     #???????ļ?
setwd("C:\\Users\\benyu\\Desktop")     #???ù???Ŀ¼

#??ȡ?????????ļ?,?????????ļ?????????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#??ȡ???͵Ľ????ļ?
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(cluster), colnames(data))
cluster=cluster[sameSample,,drop=F]
data=data[,sameSample,drop=F]

#????֮???Ƚ?,?õ?logFC
dataL=data[,row.names(cluster)[cluster[,"ARGcluster"]!=clusterType]]
dataH=data[,row.names(cluster)[cluster[,"ARGcluster"]==clusterType]]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)

#??ȡ???????ļ?
gmt=read.gmt(gmtFile)

#????GSEA????????
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff=1, minGSSize=15, maxGSSize=500)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)
	
#??????????ͼ??
termNum=10     #????չʾͨ·????Ŀ, չʾǰ5??????????????ͨ·
showTerm=row.names(kkTab)[1:termNum]      #??ȡչʾͨ·??????
gseaplot=gseaplot2(kk, showTerm, base_size=12, title=paste0("Enriched in Cluster ",clusterType))
pdf(file=paste0("GSEA.", clusterType, ".pdf"), width=8, height=8)
print(gseaplot)
dev.off()

