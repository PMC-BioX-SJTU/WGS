Args    =       commandArgs(TRUE);
inputF1 =       Args[1];#all.plink_clean.frq
output1 =       Args[2];#all.plink_clean.maf.pdf


#setwd("/home/limo/Project_data/Fudan-BA/All-variants")
#getwd()
frq <- read.table(inputF1,header=TRUE)
pdf(output1)
hist(frq$MAF,xlab="Minor allele frequency",breaks=100,main="Histogram of MAF")
dev.off()
q()

