Args    =       commandArgs(TRUE)
inputF1 =       Args[1];#all.plink_QC3.tdt
output1 =       Args[2];#all.plink_QC3.tdt.manhattan.jpg
output2 =       Args[3];#all.plink_QC3.tdt.QQ.jpg

library(qqman)
tdt_unadj <- read.table(inputF1,header=TRUE)
tdt_unadj_clean <- tdt_unadj[which(is.na(tdt_unadj$P) == FALSE & is.na(tdt_unadj$OR) == FALSE & is.na(tdt_unadj$CHISQ) == FALSE),]
chr <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
#for (i in chr){
#        jpeg(paste(output1,i,".jpg",sep=""))
#        manhattan(subset(tdt_unadj_clean,CHR==i),main = "Manhattan Plot",ylim =c(0,12),col=c("#9a0dea"))
#        dev.off()
#}

jpeg(output1)
manhattan(tdt_unadj_clean,main = "Manhattan Plot",ylim =c(0,12))
dev.off()
jpeg(output2)
qq(tdt_unadj_clean$P,main = "Q-Q Plot",ylim = c(0,12))
dev.off()
median((tdt_unadj_clean$CHISQ)^2)/0.455 #lambda_GC
median((tdt_unadj$CHISQ)^2)/0.455 #lambda_GC
quit()

