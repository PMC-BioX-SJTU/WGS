#!/usr/bin/perl -w
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $outdir=dirname($ARGV[0]);
open DB,$ARGV[0] or die $!;
open DF,">$outdir/$ARGV[1]" or die $!;
my %depth=();
my ($i,$Initial_bases_on_targe)=(0,0);
while(<DB>){
	chomp;
	$i++;
	if($i>12){
		next if($_ eq "");
		my @tmp=split;
		$depth{$tmp[0]}=$tmp[1];
		$Initial_bases_on_target+=$tmp[1];
	}
	if($i==8){
		my @tmp=split;
		$Average_sequencing_depth_on_target=$tmp[22];
	}
}
close(DB);
my $maxCov=0;
	
foreach my $depth (sort {$a<=>$b} keys %depth){
	next if($depth==0);
	my $per=$depth{$depth}/$Initial_bases_on_target;
	$maxCov = $per if($per > $maxCov);
	print DF "$depth\t$per\t$depth{$depth}\n";
}
close(DF);
	
open CU,">$outdir/$ARGV[2]" or die $!;
print CU "Depth\tTRPercent\n";
my @depth= sort {$a<=>$b} keys %depth;
my $index=0;
foreach my $d (sort {$a<=>$b} keys %depth){
	my $tmp=0;
	$tmp=$Initial_bases_on_target-$index;
	$index+=$depth{$d};
	$rate = $tmp/$Initial_bases_on_target;
	print CU "$d\t$rate\n";
}
close(CU);
### setting for plots
my $ylim = 100*$maxCov;
my ($xbin,$ybin);
$ylim= int($ylim) + 1;
if($ylim <= 3){
	$ybin = 0.5;
}else{
	$ybin=1;
}
my $xlim=0;
if($Average_sequencing_depth_on_target<50){
	$xlim=200;
	$xbin=20;
}elsif($Average_sequencing_depth_on_target < 100){
	$xlim=400;
	$xbin=20;
}elsif($Average_sequencing_depth_on_target < 220){
	$xlim=600;
	$xbin=50;
}else{
	$xlim=1000;
	$xbin=100;
}
histPlot($outdir,$ARGV[1],$ylim,$ybin,$xlim,$xbin, $ARGV[3]);
cumuPlot($outdir,$ARGV[2],$xlim,$xbin, $ARGV[3]);




sub cumuPlot {
        my ($outdir, $dataFile, $xlim, $xbin, $sample) = @_;
        my $figFile = "$outdir/$sample\_cumuPlot.pdf";
        my $Rline=<<Rline;
        pdf(file="$figFile",w=8,h=6)
        rt <- read.table("$outdir/$dataFile",header=T)
        opar <- par()
        x <- rt\$Depth[1:($xlim+1)]
        y <- 100*rt\$TRPercent[1:($xlim+1)]
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
        xpos <- seq(0,$xlim,by=$xbin)
        ypos <- seq(0,100,by=20)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Cumulated fraction of Reads > Sequence depth(%)",side=2, line=3, at=median(ypos), cex=1.5 )
        mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
        par(opar)
        dev.off()
        png(filename="$outdir/$sample\_cumuPlot.png",width = 480, height = 360)
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=3, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
        xpos <- seq(0,$xlim,by=$xbin)
        ypos <- seq(0,100,by=20)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
        mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.5)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
        par(opar)
        dev.off()

Rline
        open (ROUT,">$figFile.R");
        print ROUT $Rline;
        close(ROUT);

        system("R CMD BATCH $figFile.R");
}


sub histPlot {
	my ($outdir, $dataFile, $ylim, $ybin, $xlim, $xbin, $sample) = @_;
	my $figFile = "$outdir/$sample\_histPlot.pdf";
	my $Rline=<<Rline; 
	pdf(file="$figFile",w=8,h=6)
	rt <- read.table("$outdir/$dataFile")
	opar <- par()
	t=sum(rt\$V2[($xlim+1):length(rt\$V2)])
	y=c(rt\$V2[1:$xlim],t)
	y <- y*100
	x <- rt\$V1[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/$sample\_histPlot.png",width = 480, height = 360)
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("R CMD BATCH $figFile.R");
}

sub Plot_insert_size {
	my ($outdir, $dataFile,$sample) = @_;
	my $figFile = "$outdir/$sample\_insert_size.pdf";
	my $Rline=<<Rline; 
	data <- scan("$dataFile")
	pdf("$figFile",width=8)
	if(range(data)[2] > 1000){
        	plot(table(data), xlab="Insert Size", ylab="Reads Number",xlim=c(0,1000))
        	data <- data[data<1000]
    		hist(data, xlim=c(0,1000))
	}else{
        	plot(table(data), xlab="Insert Size", ylab="Reads Number")
        	data <- data[data<1000]
		hist(data)
	}
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);
	system("R CMD BATCH $figFile.R");
}

