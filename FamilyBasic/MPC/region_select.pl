#! /usr/bin/perl -w
open IN,$ARGV[0]||die; ## table S4
open OUT,">$ARGV[0].cand.txt"||die;
while(<IN>){
	chomp;
	if(/transcript/){
		print OUT "$_\n";
		next;
	}
	@tt=split/\t+/;
	if($tt[9] >10.38 and $tt[8]<0.6){
		print OUT "$_\n";
		$gene{$tt[1]}++;
	}
}
open IN1,$ARGV[1]||die; ## table S2
open OUT1,">$ARGV[1].cand.txt"||die;
while(<IN1>){
        chomp;
	if(/transcript/){
                print OUT1 "$_\n";
                next;
        }
        @tt=split/\t+/;
	if($tt[11]<0.6){
		next if(defined $gene{$tt[1]});
        	print OUT1 "$_\n";
	}
}
