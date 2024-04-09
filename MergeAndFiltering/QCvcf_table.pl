#! /usr/bin/perl -w
##perl QualityControlBam.pl name [.verifybamid.selfSM] [.hs.raw] [.alignment.raw] [qualigyYield.raw] [datatype]>QualitControlBam.txt
#use File::Basename;
#my $outdir=dirname($ARGV[1]);

%tab=('CompOverlap',[4,5,6,7,8,9,10],'CountVariants',[19,21,26],'IndelSummary',[5,6,7,8,13,14,15,16,17,18,19,20,21,22],'MetricsCollection',[12]);
@tit=("CompOverlap","CountVariants","IndelSummary","MetricsCollection");
Process_QCBam($ARGV[0],\%tab,\@tit);

sub Process_QCBam{
	my ($file,$index,$par)=@_;
	my $i=0;
	open IN,$file||die;
	while(<IN>){
		chomp;
		if(/^#/||$_ eq ""){$i=0;next};
		my @items=split/\s+/;
		foreach my $key(@$par){
			if($items[0] eq $key){
				$i++;
				foreach my  $j(@{$$index{$key}}){
					push @{$info{$i}},$items[$j-1];
				}
			}
		}
	}
	close IN;
	foreach my $s(sort{$a<=>$b} keys %info ){
		@out=@{$info{$s}};$ss=join("\t",@out);print "$ss\n";
		#foreach my $k (0..$#out){
		#	print "$out[$k]\n";
		#}
	}	

}


