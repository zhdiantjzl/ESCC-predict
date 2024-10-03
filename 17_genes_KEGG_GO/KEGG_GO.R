
rm(list=ls())


load("DATA2.RDATA")


#基因ID转换
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(tidyverse)

genename <-c(badgene,goodgene)  #提取第一列基因名
#x - 基因组注释R包
#keys - 需要转换的基因列表
#keytype - 基因名类型
#columns - 希望返回的基因名类型数据，如返回NCBI基因ID使用 ENTREZID
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
gene_map
write.csv(as.data.frame(gene_map),"基因转换.csv",row.names =F)#导出结果至默认路径下
genelist_input<-gene_map[,2]

head(genelist_input)
genelist_input<-na.omit(genelist_input)
#GO分析
Go_result_BP <- enrichGO(genelist_input, 'org.Hs.eg.db', ont="BP", pvalueCutoff=1) #基因ID类型为ENSEMBL的ID形式，选择BP功能组，以P值0.05为界限

goplot(Go_result_BP, showCategory=5) #GO拓扑图

p1<-dotplot(Go_result_BP, showCategory=20) #气泡图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p1,filename ='GO_BP气泡图.pdf',width = 7.23,height = 8)
ggsave(p1,filename ='GO_BP气泡图.TIFF',width = 7.23,height = 8)

barplot(Go_result_BP, showCategory=20) #条形图，显示前二十个
y=as.data.frame(Go_result_BP)
y$geneID=as.character(sapply(y$geneID,function(x)paste(gene_map$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene_map$ENTREZID))],collapse="/")))
write.csv(y,"GO-BP.csv",row.names =F)#导出结果至默认路径下。
save(y,file = 'GO-BP.RDATA')

Go_result_CC <- enrichGO(genelist_input, 'org.Hs.eg.db', ont="CC", pvalueCutoff=10) #基因ID类型为ENSEMBL的ID形式，选择CC功能组，以P值0.05为界限

goplot(Go_result_CC, showCategory=5) #GO拓扑图

p1<-dotplot(Go_result_CC, showCategory=20) #气泡图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p1,filename ='GO_CC气泡图.pdf',width = 7.23,height = 8)
ggsave(p1,filename ='GO_CC气泡图.TIFF',width = 7.23,height = 8)

barplot(Go_result_CC, showCategory=20) #条形图，显示前二十个

y=as.data.frame(Go_result_CC)
y$geneID=as.character(sapply(y$geneID,function(x)paste(gene_map$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene_map$ENTREZID))],collapse="/")))
write.csv(y,"GO-CC.csv",row.names =F)#导出结果至默认路径下
save(y,file = 'GO-CC.RDATA')

Go_result_MF <- enrichGO(genelist_input, 'org.Hs.eg.db',ont="MF", pvalueCutoff=1000) #基因ID类型为ENSEMBL的ID形式，选择MF功能组，以P值0.05为界限

goplot(Go_result_MF, showCategory=5) #GO拓扑图

p1<-dotplot(Go_result_MF, showCategory=20) #气泡图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p1,filename ='GO_MF气泡图.pdf',width = 7.23,height = 8)
ggsave(p1,filename ='GO_MF气泡图.TIFF',width = 7.23,height = 8)

barplot(Go_result_MF, showCategory=20) #条形图，显示前二十个
y=as.data.frame(Go_result_MF)
y$geneID=as.character(sapply(y$geneID,function(x)paste(gene_map$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene_map$ENTREZID))],collapse="/")))

write.csv(y,"GO-MF.csv",row.names =F)#导出结果至默认路径下
save(y,file = 'GO-MF.RDATA')


## GO分析#BP MF CC一起分析
go <- enrichGO(genelist_input, OrgDb = "org.Hs.eg.db", ont="all")
library(ggplot2)
p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p <-p + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p,filename ='GO三合一气泡图.pdf',width = 7.23,height = 8)
ggsave(p,filename ='GO三合一气泡图.TIFF',width = 7.23,height = 8)
write.csv(as.data.frame(go),"GO.csv",row.names =F)#导出结果至默认路径下


#KEGG分析

KEGG_result <- enrichKEGG(genelist_input, keyType = "kegg",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500,organism = "hsa", use_internal_data=T)  #KEGG富集分析


barplot(KEGG_result, showCategory=20)#绘制条形图

p1<-dotplot(KEGG_result, showCategory=20) #气泡图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p1,filename ='KEGG气泡图.pdf',width = 7.23,height = 8)
ggsave(p1,filename ='KEGG气泡图.TIFF',width = 7.23,height = 8)


#圈图
pdf(file="KEGG_circos.pdf",width = 10,height = 7)
kkx=setReadable(KEGG_result, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kkx, showCategory = 5, circular = TRUE, colorEdge = TRUE,node_label="all")
dev.off()
x=as.data.frame(KEGG_result)
x$geneID=as.character(sapply(x$geneID,function(x)paste(gene_map$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene_map$ENTREZID))],collapse="/")))


write.csv(as.data.frame(x),"KEGG.csv",row.names =F)#导出结果至默认路径下
save(x,file = 'KEGG.RDATA')


#用P值来排序，不用padjust
#首先，我获取了富集对象x中的数据框，这是S4对象，用@符号来获取,有293行，是没有筛选过的数据
y =KEGG_result@result

#我们的横坐标有问题，是因为这里的GeneRatio是字符串，我们现在要把它变成一个数值
## 分别后去分号前面和后面的数，并变成数值
forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))
## 增加数值表示的一列GeneRatio
y$GeneRatio = forward/backward

showCategory =20
#再设定一个字体大小，大小后期可以调整
font.size =12

library(ggplot2)
library(forcats)
library(dplyr)
#复现气泡图
AAA=y %>% 
  ## 安装p值排序，选区既定数目的行
  arrange(pvalue) %>% 
  slice(1:showCategory)
## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
p1=ggplot(AAA,aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=pvalue, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(filename ='KEGG气泡图2.pdf',width = 8,height = 8)
ggsave(filename ='KEGG气泡图2.TIFF',width = 8,height = 8)

#用P值来排序，不用padjust
#首先，我获取了富集对象x中的数据框，这是S4对象，用@符号来获取,有293行，是没有筛选过的数据
y =Go_result_MF@result

#我们的横坐标有问题，是因为这里的GeneRatio是字符串，我们现在要把它变成一个数值
## 分别后去分号前面和后面的数，并变成数值
forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))
## 增加数值表示的一列GeneRatio
y$GeneRatio = forward/backward

showCategory =20
#再设定一个字体大小，大小后期可以调整
font.size =12

library(ggplot2)
library(forcats)
library(dplyr)
#复现气泡图
AAA=y %>% 
  ## 安装p值排序，选区既定数目的行
  arrange(pvalue) %>% 
  slice(1:showCategory)
## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
p1=ggplot(AAA,aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=pvalue, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(filename ='Go_result_MF2.pdf',width = 8,height = 8)
ggsave(filename ='Go_result_MF2.TIFF',width = 8,height = 8)


#用P值来排序，不用padjust
#首先，我获取了富集对象x中的数据框，这是S4对象，用@符号来获取,有293行，是没有筛选过的数据
y =Go_result_BP@result

#我们的横坐标有问题，是因为这里的GeneRatio是字符串，我们现在要把它变成一个数值
## 分别后去分号前面和后面的数，并变成数值
forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))
## 增加数值表示的一列GeneRatio
y$GeneRatio = forward/backward

showCategory =20
#再设定一个字体大小，大小后期可以调整
font.size =12

library(ggplot2)
library(forcats)
library(dplyr)
#复现气泡图
AAA=y %>% 
  ## 安装p值排序，选区既定数目的行
  arrange(pvalue) %>% 
  slice(1:showCategory)
## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
p1=ggplot(AAA,aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=pvalue, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(filename ='Go_result_BP2.pdf',width = 8,height = 8)
ggsave(filename ='Go_result_BP2.TIFF',width = 8,height = 8)


#用P值来排序，不用padjust
#首先，我获取了富集对象x中的数据框，这是S4对象，用@符号来获取,有293行，是没有筛选过的数据
y =Go_result_CC@result

#我们的横坐标有问题，是因为这里的GeneRatio是字符串，我们现在要把它变成一个数值
## 分别后去分号前面和后面的数，并变成数值
forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))
## 增加数值表示的一列GeneRatio
y$GeneRatio = forward/backward

showCategory =20
#再设定一个字体大小，大小后期可以调整
font.size =12

library(ggplot2)
library(forcats)
library(dplyr)
#复现气泡图
AAA=y %>% 
  ## 安装p值排序，选区既定数目的行
  arrange(pvalue) %>% 
  slice(1:showCategory)
## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
p1=ggplot(AAA,aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=pvalue, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(filename ='Go_result_CC2.pdf',width = 8,height = 8)
ggsave(filename ='Go_result_CC2.TIFF',width = 8,height = 8)




#GO结果画有序条形图
library(ggpubr)
library(tidyverse)
GODATA<-go@result
GODATA$"-log10(Pvalue)"<-  -log10(GODATA$pvalue)
GODATA$yyy<-  -log10(GODATA$pvalue)
colnames(GODATA)


write.csv(GODATA,"GODATA.csv",row.names =F)#导出结果至默认路径下

#分别选取BP MF CC 的前十个
aaa<-filter(GODATA,ONTOLOGY=='BP')
aaaa<-aaa[1:10,]

bbb<-filter(GODATA,ONTOLOGY=='CC')
bbbb<-bbb[1:10,]

ccc<-filter(GODATA,ONTOLOGY=='MF')
cccc<-ccc[1:10,]


drawdata<-rbind(aaaa,bbbb,cccc)


ggbarplot(drawdata, x = "Description", y = 'yyy',
          fill = "ONTOLOGY",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c('#5FB404','#01DFD7','#C238E5'),            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90           # Rotate vertically x axis texts
)
ggbarplot(drawdata, x = "Description", y = "yyy",
          fill = "ONTOLOGY",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c('#5FB404','#12B5EC','#C238E5'),            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 75, # Rotate vertically x axis texts
          ylab = '-log10(P-Value)',
          xlab = 'Pathway'
)
ggsave(filename ='GO三合一2.pdf',width = 7.23,height = 8)
ggsave(filename ='GO三合一2.TIFF',width = 7.23,height = 8)



