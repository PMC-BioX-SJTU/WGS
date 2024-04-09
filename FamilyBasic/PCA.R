Args    =       commandArgs(TRUE);
inputF1 =       Args[1];#all.plink_QC4_indep.eigenvec
output1 =       Args[2];#all.plink_QC4_indep.pca.pdf
output2 =       Args[3];#all.plink_QC4_indep.remove

#setwd("/home/limo/Project_data/Fudan-BA/All-variants")
#getwd()
pc <- read.table(inputF1,header=FALSE)
pdf(output1)
plot(pc[,3],pc[,4],pch=20,xlab="PC1",ylab="PC2",main="PCA plot-PC1-PC2")
dev.off()
write.table(pc[which(abs(pc[,4]) > (mean(pc[,4]) + 3*sd(pc[,4]))),c(1,2)],output2,append=FALSE,quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
q()

