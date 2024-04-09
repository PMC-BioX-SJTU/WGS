#! /usr/bin/perl -w
#perl restrict2bed.pl [.bed] [.vcf] > [restricted.vcf]
my @tit;
open IN,$ARGV[0]||die;### vcf file
while(<IN>){
	chomp;
	if(/^##/){
		next
	}elsif(/^#CH/){
		@tit=split(/\t/,$_);	
		last;
	}
}

open IN1,$ARGV[1]||die;
while(<IN1>){
        chomp;
	my @tt=split(/\t/);
	if(/Gene.refGene/){
		$len=$#tt-1;
		$out=join("\t",(@tt[0..$len],@tit));
		print "$out\n";
	}else{
		$out=join("\t",(@tt[0..$len],@tt[$len+4..$#tt]));
		print "$out\n";
	}		
}
