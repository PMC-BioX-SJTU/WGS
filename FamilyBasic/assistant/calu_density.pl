#! /usr/bin/perl -w
%chr_len=('chr1'=>249250621,'chr2'=>243199373, 'chr3'=>198022430,   'chr4'=>191154276,'chr5'=>180915260,'chr6'=>171115067 , 'chr7'  =>  159138663, 'chr8' =>   146364022 , 'chr9' =>  141213431,'chr10'=>   135534747,'chr11' =>  135006516,'chr12' =>  133851895,'chr13' =>  115169878 , 'chr14' =>  107349540  , 'chr15' =>  102531392  , 'chr16' =>  90354753 , 'chr17' =>  81195210 , 'chr18' =>  78077248 ,  'chr19' =>  59128983 ,  'chr20' =>  63025520 , 'chr21' =>  48129895 ,  'chr22' =>  51304566 ,  'chrX'  =>  155270560, 'chrY'  =>  59373566 , 'chrM'  =>  16571);
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	next if(/^#/);
	@tt=split/\s+/;
	$count{$tt[0]}++;
}
foreach my $i(keys %count){
	$density=$count{$i}*1000000/$chr_len{$i};
	print "$i\t$chr_len{$i}\t$count{$i}\t$density\n"
}
