#install.packages('survival')


library(survival)      #???ð?
setwd("C:\\Users\\benyu\\Desktop")      #???ù???Ŀ¼

############????ɭ??ͼ?ĺ???############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	#??ȡ?????ļ?
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#????ͼ??
	pdf(file=forestFile, width=6.6, height=4.5)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#????ɭ??ͼ???ߵ??ٴ???Ϣ
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
		
	#?????ұߵ?ɭ??ͼ
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
	axis(1)
	dev.off()
}
############ɭ??ͼ????############

############??��Ԥ??????????#############
indep=function(riskFile=null, cliFile=null, project=null){
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #??ȡ?ٴ??ļ?
	
	#???ݺϲ?
	sameSample=intersect(row.names(cli),row.names(risk))
	risk=risk[sameSample,]
	cli=cli[sameSample,]
	rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
	
	#?????ض?��Ԥ??????
	uniCoxFile=paste0(project,".uniCox.txt")
	uniCoxPdf=paste0(project,".uniCox.pdf")
	uniTab=data.frame()
	for(i in colnames(rt[,3:ncol(rt)])){
		 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
		 coxSummary = summary(cox)
		 uniTab=rbind(uniTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
	}
	write.table(uniTab,file=uniCoxFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=uniCoxFile, forestFile=uniCoxPdf, forestCol="#6699FFFF")
	
	
	#?????ض?��Ԥ??????
	multiCoxFile=paste0(project,".multiCox.txt")
	multiCoxPdf=paste0(project,".multiCox.pdf")
	uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
	rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
	multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
	multiCoxSum=summary(multiCox)
	multiTab=data.frame()
	multiTab=cbind(
	             HR=multiCoxSum$conf.int[,"exp(coef)"],
	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	multiTab=cbind(id=row.names(multiTab),multiTab)
	write.table(multiTab, file=multiCoxFile, sep="\t", row.names=F, quote=F)
	bioForest(coxFile=multiCoxFile, forestFile=multiCoxPdf, forestCol="#EE0000FF")
}
############??��Ԥ??????????#############

#???ú????????ж?��Ԥ??????
indep(riskFile="risk.TCGAall.txt", cliFile="clinical.txt", project="all")



