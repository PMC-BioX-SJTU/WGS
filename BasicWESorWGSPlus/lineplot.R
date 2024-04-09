Args    =       commandArgs(TRUE);
inputF1 =       Args[1];#input file; header, sample names;
output1 =       Args[2];#output pdf file name

a=read.table(inputF1,sep="\t",header=T)

png(output1)

sum <- (a[,13]+a[,14]+a[,15]+a[,16]+a[,17])
plot(a[,1],a[,13]/sum,type="l",col=2,xlab="read position",ylab="% of total(per read position)",lty=1,lwd=3,ylim=c(0,0.5),main="Nucletide Distribution")

for(i in 14:17){
    lines(a[,1],a[,i]/sum,type="l",col=(i-11),lwd=3)
}
legend("topright",c("A","C","G","T","N"),col=2:6,lty=1,box.col="white",lwd=3)
dev.off()
