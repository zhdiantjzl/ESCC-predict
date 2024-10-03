rm(list = ls())
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(survival)
library(survminer)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#save.image("DATA3.RDATA")
load("RS.RDATA")
score_t=rs$merge
#score_t$OS.time=score_t$OS.time/365
cut <- surv_cutpoint(score_t,'OS.time','OS','RS')
cut
plot(cut)
cat <- surv_categorize(cut)
fit1 <- survfit(Surv(OS.time,OS)~RS,cat)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black")) ## 自定义主题

ggsurvplot(fit1,cat,
           palette = 'jco',
           size=1.3,
           pval=T,
           legend.labs=c("High","Low"), 
           legend.title='Score',
           xlab="Time(years)",
           ylab='Overall survival',
           ggtheme = mytheme,
           break.time.by=1,
           conf.int=T,
           risk.table=TRUE,
           risk.table.title="",
           risk.table.height=.25)
dev.off()

cut1 <- surv_categorize(cut)
row.names(cut1)=score_t$ID
save(cut1,file = "merge风险分组信息.RDATA")













