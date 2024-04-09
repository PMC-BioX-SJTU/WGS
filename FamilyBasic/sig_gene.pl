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
        $syn{$tt[1]}=(10**$tt[4])*2;
        $missense{$tt[1]}=(10**$tt[5])*2;
        $nonsense{$tt[1]}=(10**$tt[6])*2;
        $splice{$tt[1]}=(10**$tt[7])*2;
        $frameshift{$tt[1]}=(10**$tt[8])*2;
        $gene{$tt[1]}++;
}
open IN1,$ARGV[1]||die;
while(<IN1>){
        chomp;
	my $exp=0;
        @tt=split/\t/;
	@type=split(/,/,$tt[1]);
	foreach(@type){
		if($_=~/missense/){
			$exp+=$missense{$tt[0]};
		}elsif($_=~/frameshift/){
			$exp+=$frameshift{$tt[0]};
		}elsif($_=~/nonsense/){
			$exp+=$nonsense{$tt[0]}
		}
	}
	$exp=$exp*$ARGV[2];
	print "$tt[0]\t$tt[1]\t$exp\n";
}
