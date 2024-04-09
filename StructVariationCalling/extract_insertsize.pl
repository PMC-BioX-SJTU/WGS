#! /usr/bin/perl -w
##perl QualityControlBam.pl name [insert_sizefile]  [bamfile]>QualitControlBam.txt
use File::Basename;
@wgs=(1);
Process_QCBam($ARGV[1],7,1,$ARGV[0],'Q20Q30',@wgs);

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
	foreach my  $i(0..$#index){
		if($index[$i] =~ /\d+/){
			$tit[$index[$i]-1]=~s/BAIT/TARGET/g;
			print   "$ARGV[2]\t$val[$index[$i] - 1]\t$ARGV[0]\n";
		}else{
			$rate=$val[$index[$i-1]]/$val[$index[2]];
			$rate=sprintf "%.5f",$rate;
			print OUT "$tit[$index[$i-1]-1]_RATE\t$rate\n";
		}
	}
	close OUT;
}

