#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	my @tt=split/\s+/;
	if(/^#/){
		print "$_\n";
	}else{
		if($tt[7]=~/hiConfDeNovo=(.*)/){
			$tt[2]=$1;
			$out=join("\t",@tt);
			print "$out\n";
		}
	
	}
}
