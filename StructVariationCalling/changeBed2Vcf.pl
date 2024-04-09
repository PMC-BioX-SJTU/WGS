#! /usr/bin/perl -w
$/=">";
open IN,$ARGV[0]||die;## hg19 genome
while(<IN>){
        chomp;
        next if($_ eq "");
        @ss=split(/\n/,$_,2);
        $ss[1]=~s/\n//g;
        $seq{$ss[0]}=$ss[1];
	print "$ss[0]\n";
}
$/="\n";
$i=0;
%chr=();
open IN1,"$ARGV[1]"||die;## tab file
while(<IN1>){
	chomp;
	$i++;
	@tt=split/\s+/;
	$id="DEL".$i;	
	push @{$chr{$tt[0]}},[$tt[1],$tt[2],$id];
}
@tit=();
open IN2,"$ARGV[2]"||die;## header file
while(<IN2>){
        chomp;
	push @tit,$_;
}
my $header=join("\n",@tit);
foreach my $k(keys %chr){
	open OUT,">all.merged.DEL.$k.vcf"||die;
	print OUT "$header\n";
	my @list=@{$chr{$k}};
	foreach my $m(0..$#list){
		$snp=substr($seq{$k},$list[$m][0]-1,1);
		print OUT "$k\t$list[$m][0]\t$list[$m][2]\t$snp\t<DEL>\t.\tPASS\tCHR2=$k\;END=$list[$m][1]\n";
	}
	close OUT;
}
