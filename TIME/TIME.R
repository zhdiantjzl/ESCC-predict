rm(list = ls())
#devtools::install_github("IOBR/IOBR")
library(IOBR)
names(signature_tme)[1:20]
names(signature_metabolism)[1:20]
names(signature_tumor)
names(signature_collection)[1:20]
load("DATA.RDATA")
eset_stad=outTab
eset_stad[1:5, 1:5]
#需要是标准化后的数据
range(eset_stad)

tme_deconvolution_methods

# MCPcounter
im_mcpcounter <- deconvo_tme(eset = eset_stad,
                             method = "mcpcounter"
)
## 
## >>> Running MCP-counter

# EPIC
im_epic <- deconvo_tme(eset = eset_stad,
                       method = "epic",
                       arrays = F
)
## 
## >>> Running EPIC

# xCell
im_xcell <- deconvo_tme(eset = eset_stad,
                        method = "xcell",
                        arrays = F
)

# CIBERSORT
im_cibersort <- deconvo_tme(eset = eset_stad,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)
## 
## >>> Running CIBERSORT

# IPS
im_ips <- deconvo_tme(eset = eset_stad,
                      method = "ips",
                      plot = F
)
## 
## >>> Running Immunophenoscore

# quanTIseq
im_quantiseq <- deconvo_tme(eset = eset_stad,
                            method = "quantiseq",
                            scale_mrna = T
)

# ESTIMATE
im_estimate <- deconvo_tme(eset = eset_stad,
                           method = "estimate"
)
## 

# TIMER
#需要输入肿瘤类型
im_timer <- deconvo_tme(eset = eset_stad
                        ,method = "timer"
                        ,group_list = rep("gbm",dim(eset_stad)[2])
)



dim(im_cibersort)

im_cibersort[1:4,1:4]

dim(im_xcell)

im_xcell[1:4,1:4]

library(tidyr)
# 取前12个样本做演示
#res<-cell_bar_plot(input = im_cibersort[1:50,], title = "CIBERSORT Cell Fraction")
im_estimate$ID[1:81]<-paste(substr(im_estimate$ID[1:81],1,5),gsub("-",".",substr(im_estimate$ID[1:81],6,21)),sep="")


tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
  inner_join(im_ips, by= "ID") %>% 
  inner_join(im_quantiseq, by="ID") %>% 
  inner_join(im_estimate, by= "ID") %>% 
  inner_join(im_timer, by= "ID")

tme_combine[1:4,1:4]

dim(tme_combine)

write.csv(tme_combine,file = "tme.score.csv",row.names = F)
save(tme_combine,file = "tme.score.RDATA")

#画图
rm(list = ls())
library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
#BiocManager::install("oompaBase")
library(oompaBase)
library(tidyverse)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
load("tme.score.RDATA")
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
tme_combine=as.data.frame(tme_combine)
tme_combine <- tme_combine[, colSums(is.na(tme_combine)) == 0]
rownames(tme_combine)=tme_combine$ID
tme_combine=tme_combine[,-1]
aaa=data.frame(ID = colnames(tme_combine),Type = NA)
# 使用strsplit函数分割字符串，并提取最后一个"_"后面的字符
aaa$Type <- sapply(strsplit(as.character(aaa$ID), "_"), function(x) tail(x, n = 1))
dataa=as.data.frame(t(tme_combine))
genename1=intersect(aaa$ID,rownames(dataa))
dataa=dataa[genename1,]
aaa=dplyr::filter(aaa,ID %in% genename1)

risk=clin
rownames(risk)=risk$ID
risk=risk[,-1]
sameid=intersect(colnames(dataa),rownames(risk))
dataa=dataa[,sameid]
risk=risk[sameid,]
hmdat=as.data.frame(t(dataa))
#分类信息
type=data.frame(Type =aaa$Type)
row.names(type)=aaa$ID

Typeid <- type$Type
colnames(risk)
# 创建注释
# 列注释，位于热图顶端
annCol <- data.frame(score = scale(risk$RS),
                     RiskType = risk$Risk_group,
                     OS = risk$OS,
                     Gender = risk$gender,
                     cohort = risk$cohort,
                     # 以上是risk score和risk type两种注释，可以按照这样的格式继续添加更多种类的注释信息，记得在下面的annColors里设置颜色
                     row.names = rownames(risk),
                     stringsAsFactors = F)



#基因差异分析
bbb=annCol[rownames(hmdat),]
sampleType=bbb$RiskType
sigVec=c()
for(i in colnames(hmdat)){
  test=wilcox.test(hmdat[,i] ~ sampleType)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(i, Sig))
}

colnames(hmdat)=sigVec
# 数据标准化

# 用pheatmap画图
library(pheatmap)

# 定义颜色
Type.col <- brewer.pal(n = length(unique(Typeid)),name = "Paired")


# 行注释，位于热图左侧
annRow <- data.frame(Type = factor(Typeid,levels = unique(Typeid)),
                     row.names = colnames(hmdat),
                     stringsAsFactors = F)


unique(annRow$Type)
# 为各注释信息设置颜色
annColors <- list(Type = c("MCPcounter" = Type.col[1], #行注释的颜色
                           "EPIC" = Type.col[2],
                           "xCell" = Type.col[3],
                           "CIBERSORT" = Type.col[4],
                           "IPS" = Type.col[5],
                           "quantiseq" = Type.col[6],
                           "estimate" = Type.col[7],
                           "TIMER" = Type.col[8]),
                  # 下面是列注释的颜色，可依此设置更多注释的颜色
                  "RS" = greenred(64), 
                  "RiskType" = c("High_risk" = "red","Low_risk" = "blue"))



# 样本按risk score排序
samorder <- rownames(risk[order(risk$RS),])


indata <- as.data.frame(t(hmdat))
indata[indata==0]=0.00000000000000000000000000001
#indata <- indata[,colSums(indata) > 0] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)


# pheatmap绘图
pheatmap::pheatmap(mat = as.matrix(plotdata[,samorder]), # 标准化后的数值矩阵
                   border_color = NA, # 无边框色
                   color = bluered(64), # 热图颜色为红蓝
                   cluster_rows = F, # 行不聚类
                   cluster_cols = F, # 列不聚类
                   show_rownames = T, # 显示行名
                   show_colnames = F, # 不显示列名
                   annotation_col = annCol[samorder,,drop = F], # 列注释
                   annotation_row = annRow, # 行注释
                   annotation_colors = annColors, # 注释颜色
                   gaps_col = table(annCol$RiskType)[2], # 列分割
                   gaps_row = cumsum(table(annRow$Type)), # 行分割
                   cellwidth = 0.8, # 元素宽度
                   cellheight = 10, # 元素高度
                   filename = "immune heatmap by pheatmap_score.pdf")



hmdat$ID=rownames(hmdat)
annCol$ID=rownames(annCol)
annCol=dplyr::select(annCol,ID,everything())
newdata=inner_join(annCol,hmdat,by="ID")


inputtemp<-newdata

colnames(inputtemp) <- gsub("\\**", "", colnames(inputtemp))
colnames(inputtemp) <- gsub("\\*", "", colnames(inputtemp))

gene="score"
value=colnames(inputtemp)[7:142]

for (b in 1:136) {
  test <- cor.test(inputtemp[,gene],inputtemp[,value[b]],exact=FALSE)
  paste(paste0("n = ",length(inputtemp[,gene])),
        paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
        paste0("p.value= ",round(test$p.value,4)),
        sep = ", ")
  p<-ggplot(inputtemp,aes(get(gene),get(value[b])))+
    geom_point(col="#984ea3")+
    geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
    geom_rug(col="#7fc97f")+
    theme_minimal()+
    xlab(paste("Riskscore",sep = ' '    ))+
    ylab(paste(value[b],sep = ' '    ))+
    ## 依靠函数来生成title
    labs(title = paste(
      paste0("n = ",length(inputtemp[,gene])),
      paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
      paste0("p.value= ",round(test$p.value,4)),
      sep = ", "))+
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  ggsave(p,filename = paste(gene,value[b],'相关性点图.pdf',sep=' '),width = 4.94,height = 4.72)
}




#绘制箱线图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
#bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
bioCol=bioCol[1:length(levels(factor(type[,1])))]
library(tidyverse)
library(ggpubr)


library(ggpubr)
# 自定义函数将p值转换为星号并居中
custom_label <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

genename2=str_replace_all(genename1,"-","")
colnames(inputtemp)=str_replace_all(colnames(inputtemp),"-","")


genename2=str_replace_all(genename2,"\\+","")
colnames(inputtemp)=str_replace_all(colnames(inputtemp),"\\+","")


genename2=str_replace_all(genename2,"\\(","")
colnames(inputtemp)=str_replace_all(colnames(inputtemp),"\\(","")

genename2=str_replace_all(genename2,"\\)","")
colnames(inputtemp)=str_replace_all(colnames(inputtemp),"\\)","")

inputtemp$RiskType<-factor(inputtemp$RiskType,levels = c('Low_risk','High_risk'))
my_comparisons <- list( c("Low_risk", "High_risk") )

for (i in 1:length(genename1)) {
  immunecell=genename2[i]
  # 绘制箱线图
  p <- ggboxplot(inputtemp, x = "RiskType", y = immunecell,
                 fill = "RiskType",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
    theme(legend.position='none')+
    ylab(label = paste(immunecell))+
    stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test',aes(label = ..p.signif..))
  ggsave(p,filename = paste(immunecell,'差异.pdf',sep=' '),width = 5.6,height = 4.22)
}
