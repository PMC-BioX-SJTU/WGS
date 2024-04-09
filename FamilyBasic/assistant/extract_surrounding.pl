#! /usr/bin/perl
$/=">";
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	next if($_ eq "");
	@ss=split(/\n/,$_,2);
	$ss[1]=~s/\n//g;
        $seq{$ss[0]}=$ss[1];
}
$/="\n";
open IN1,$ARGV[1]||die;
while(<IN1>){
        chomp;
	@tt=split/\s+/;
	$seq=substr($seq{$tt[0]},$tt[1]-300,600+$tt[2]-$tt[1]);
	$len=300+$tt[2]-$tt[1];
	$se="300:$len";
	print "$_\t$seq\t$se\n";
}
