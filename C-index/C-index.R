#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("survcomp")

#install.packages("survival")
#install.packages("ggplot2")
#install.packages("ggpubr")


#???ð?
library(survival)
library(survcomp)
library(ggplot2)
library(ggpubr)

inputFile="risk.models.txt"      #?????ļ?
setwd("C:\\Users\\lenovo\\Desktop")      #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#??ģ?ͽ???ѭ??,?õ?ÿ??ģ?͵?C-indexֵ
df=data.frame()
for(i in colnames(rt)[3:ncol(rt)]){
	cindex=concordance.index(x=rt[,i], surv.time=rt$futime, surv.event=rt$fustat, method="noether")
	df=rbind(df, cbind(i,sprintf("%.03f",cindex$c.index)))
}
colnames(df)=c("signature", "cindex")
df[,"cindex"]=as.numeric(df[,"cindex"])

#????C-index????״ͼ
color=rainbow(nrow(df),alpha=0.75)
p=ggbarplot(df, x="signature", y="cindex", fill="signature",
		  xlab="", ylab="C-index", add = "none",
		  palette="lancet",
          label=T, legend="")
p=p+rotate_x_text(50)
p=p+ylim(0,round(max(df[,"cindex"])+0.15,1))
#????ͼ??
pdf(file="C-index.pdf", width=6, height=4.6)
print(p)
dev.off()
