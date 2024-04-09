#! /usr/bin/perl -w
#perl restrict2bed.pl [.bed] [.vcf] > [restricted.vcf]
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	@tt=split/\s+/;
	push @{$info{$tt[0]}},[$tt[1],$tt[2]];
}
my %out;
open IN1,$ARGV[1]||die;
while(<IN1>){
	chomp;
	if(/^#/){print "$_\n";next};
	$line=$_;
        @tt=split(/\s+/,$line);
	$key=$tt[0];	
	if(defined $info{$key}){
		@ss=@{$info{$key}};
	}else{
		print "$line\n";
		next;
	}
	$index=0;
	foreach my $i(0..$#ss){
		if($tt[1]>$ss[$i][0] && $tt[1]<$ss[$i][1] ){
			$index++;
		}
	}
	if($index==0){
		print "$line\n";
	}
}
