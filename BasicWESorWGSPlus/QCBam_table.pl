#! /usr/bin/perl -w
##perl QualityControlBam.pl name [.verifybamid.selfSM] [.hs.raw] [.alignment.raw] [qualigyYield.raw] [datatype]>QualitControlBam.txt
use File::Basename;
my $outdir=dirname($ARGV[1]);
@wgs=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25);
@hs=(6,7,8,10,15,16,17,19,20,23,24,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42);
@alig= (2,3,4,7,9,13,14,15,18,20,21);
@fq=(1,3,4,6,'_rate',8,'_rate');
@mix=(4,5,6,7);
Process_QCBam($ARGV[1],1,1,$ARGV[0],'CrossCont',@mix);

if($ARGV[5] eq 'WGS' ){
	Process_QCBam($ARGV[2],7,1,$ARGV[0],'Hs',@wgs);
}else{
	Process_QCBam($ARGV[2],7,1,$ARGV[0],'Hs',@hs);
}
Process_QCBam($ARGV[3],7,3,$ARGV[0],'Alig',@alig);
Process_QCBam($ARGV[4],7,1,$ARGV[0],'Q20Q30',@fq);

sub Process_QCBam{
	my ($file,$line_count,$flank,$name,$outfile,@index)=@_;
	my $i=0;
	open IN,$file||die;
	while(<IN>){
		chomp;
		@items=split/\s+/;
		$i++;
		if($i==$line_count ){
			@tit=@items;
		}elsif($i==$line_count+$flank){
			@val=@items;	
		}
	}
	close IN;
	open OUT,">$outdir/$name.$outfile.final.txt"||die;
	print OUT "Sample_ID\t$name\n";
	foreach my  $i(0..$#index){
		if($index[$i] =~ /\d+/){
			$tit[$index[$i]-1]=~s/BAIT/TARGET/g;
			print OUT  "$tit[$index[$i] - 1]\t$val[$index[$i] - 1]\n";
		}else{
			$rate=$val[$index[$i-1]]/$val[$index[2]];
			$rate=sprintf "%.5f",$rate;
			print OUT "$tit[$index[$i-1]-1]_RATE\t$rate\n";
		}
	}
	close OUT;
}

