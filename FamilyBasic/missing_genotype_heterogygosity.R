Args    =       commandArgs(TRUE);
inputF1 =       Args[1];#all.plink_QC1.imiss
inputF2 =       Args[2];#all.plink_QC1.het
output1 =       Args[3];#all.plink.imiss.het.pdf
output2 =       Args[4];#all.plink.imiss.het.pdf

#setwd("/home/limo/Project_data/Fudan-BA/All-variants")
#getwd()
imiss=read.table(inputF1,header=TRUE)
imiss$logF_MISS = log10(imiss[,6])
het=read.table(inputF2,header=TRUE)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
pdf(output1)
plot(imiss$logF_MISS,het$meanHet,col="black",pch=20,xlab="Proporty of missing genotypes",ylab="Heterozygosity rate",axes=F,xlim=c(-3,0),ylim=c(0,0.5))
axis(2,at=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
abline(h=mean(het$meanHet)-(3*sd(het$meanHet)),col="RED",lty=2)
abline(h=mean(het$meanHet)+(3*sd(het$meanHet)),col="RED",lty=2)
abline(v=log10(0.05),col="RED",lty=2)
dev.off()
upplimit <- mean(het$meanHet)+(3*sd(het$meanHet))
lowlimit <- mean(het$meanHet)-(3*sd(het$meanHet))
het.remove <- het[which(het$meanHet < lowlimit | het$meanHet > upplimit) , c("FID","IID")]
imiss.remove <- imiss[which(imiss$F_MISS > 0.05),c("FID","IID")]
het.remove
imiss.remove
write.table(rbind(het.remove,imiss.remove),output2,append=FALSE,quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
q()

