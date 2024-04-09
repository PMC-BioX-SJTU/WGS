#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	$line=$_;
	$line=~s/^\s+//g;
	if($line=~/^FID/){
		print "$line\n";
		next;
	}
	@tt=split(/\s+/,$line);
	if($tt[0] eq $tt[2]){
		print "$_\n";
	}


}
