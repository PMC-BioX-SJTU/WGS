#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	@tt=split/\t/;
	next if(/transcript/);
	foreach my $i(5..8){
		if($tt[$i]=~/NA/){
			$tt[$i]=-100;
		}
	}
	$missense{$tt[1]}=(10**$tt[5])*2;	
	$nonsense{$tt[1]}=(10**$tt[6])*2;
	$splice{$tt[1]}=(10**$tt[7])*2;
	$frameshift{$tt[1]}=(10**$tt[8])*2;
	$miss+=(10**$tt[5])*2;
	$fram+=(10**$tt[6])*2;
	$fram+=(10**$tt[7])*2;
	$fram+=(10**$tt[8])*2;
	$sys+=(10**$tt[4])*2;
}
print "$sys\t$miss\t$fram\n";
open IN1,$ARGV[1]||die;
while(<IN1>){
        chomp;
        @tt=split/\t/;
	next if(/^Chr/);
	if(/nonsynonymous/){	
		$nonsynony+=$missense{$tt[7]};	
		$nonsyn_trio_count{$tt[6]}++;
	}elsif(/frameshift/){
		$loF+=$frameshift{$tt[7]};
		$loF_trio_count{$tt[6]}++;
	}elsif(/stopgain/){
		$loF+=$frameshift{$tt[7]};
		$loF_trio_count{$tt[6]}++;
	}elsif(/stoploss/){
		$loF+=$frameshift{$tt[7]};
		$loF_trio_count{$tt[6]}++;
	}elsif(/splicing/){
		$loF+=$frameshift{$tt[7]};
		$loF_trio_count{$tt[6]}++;
	}
}
@count_nonsy=keys %nonsyn_trio_count;
$exp_nonsyn=$nonsynony/$#count_nonsy;
@count_loF=keys %loF_trio_count;
$exp_loF=$loF/$#count_loF;
print "$nonsynony\t$exp_nonsyn\t$loF\t$exp_loF\n";
