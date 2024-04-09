#! /usr/bin/perl -w
open IN,$ARGV[0]||die;##TDT
while(<IN>){
	chomp;
	next if(/^CHR/);
	$line=$_;
	$line=~s/^\s+//g;
	@tt=split(/\s+/,$line);
	if($tt[11]<0.0001 and ($tt[8]>1.4 or $tt[9]<0.5 )){
		$key="chr".$tt[0]."\t".$tt[2];
		$pos{$key}=$line;
		push @{$pos1{"chr".$tt[0]}},[$tt[2],$line];
	}
}
open IN1,$ARGV[1]||die;##comphet
while(<IN1>){
        chomp;
        @tt=split/\s+/;
	next if(/^Chr/);
	$key=$tt[0]."\t".$tt[1];
	if(defined $pos{$key}){
		print "$_\t0\t$pos{$key}\n";
	}
	my @chr=@{$pos1{$tt[0]}};
	@chr=sort {$a->[0] <=> $b->[0]} @chr;
	foreach my $i(0..$#chr){
		if($chr[$i][0]-200000 < $tt[1] and $tt[1] < $chr[$i][0]+200000 ){
			$len=$chr[$i][0]-$tt[1];
			print "$_\t$len\t$chr[$i][1]\n";
		}
	}
}
