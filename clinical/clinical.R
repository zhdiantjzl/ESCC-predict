rm(list = ls())


library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(tidyverse)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


load("RS.RDATA")
load("DATA.RDATA")
load("merge风险分组信息.RDATA")
rs_merge=rs$merge
cut1=cut1[rs_merge$ID,]
rs_merge$Risk_group=cut1$RS
rs_merge$Risk_group=ifelse(rs_merge$Risk_group == "high","High_risk","Low_risk")
table(rs_merge$Risk_group)
rs_merge=rs_merge[,c("ID","RS","Risk_group")]
clin=inner_join(rs_merge,clin,by="ID")

save(clin,file = "RS_Risk_clin.RDATA")
rt=clin
rownames(rt)=rt$ID
rt=rt[,-1]

colnames(rt)
#引用包
library(plyr)
library(ggplot2)
library(ggpubr)

trait="gender"                    #临床性状

rt$OS[rt$OS==1]="Dead"
rt$OS[rt$OS==0]="Alive"
#定义临床性状的颜色
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]

#统计高低评分组病人数目
rt1=rt[,c(trait, "Risk_group")]
rt1<-na.omit(rt1)
colnames(rt1)=c("trait", "Risk_group")
table(rt1$trait)
df=as.data.frame(table(rt1))
#计算高低评分组的百分率
df=ddply(df, .(Risk_group), transform, percent = Freq/sum(Freq) * 100)
#百分比位置
df=ddply(df, .(Risk_group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$Risk_group=factor(df$Risk_group, levels=c("Low_risk", "High_risk"))

#绘制百分率图
p=ggplot(df, aes(x = factor(Risk_group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("score")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="gender_barplot.pdf", width=4, height=5)
print(p)
dev.off()

#设置比较组
rt2=rt[,c(trait, "RS")]
rt2<-na.omit(rt2)
colnames(rt2)=c("trait", "score")
type=levels(factor(rt2[,"trait"]))
rt2$trait=factor(rt2$trait,levels =type )
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#绘制箱线图
boxplot=ggboxplot(rt2, x="trait", y="score", fill="trait",
                  xlab=trait,
                  ylab="score",
                  legend.title=trait,
                  palette=bioCol
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file="gender_boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()

colnames(rt)

