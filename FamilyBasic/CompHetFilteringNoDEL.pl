#! /usr/bin/perl -w
open IN, $ARGV[0]||die;
open OUT,">$ARGV[3].LoF.txt"||die;
open OUT1,">$ARGV[3].all.txt"||die;
while(<IN>){
	chomp;
	@tt=split/\t/;
	if(/^Chr/){
		print OUT "$_\t$ARGV[2]\n";
		print OUT1 "$_\t$ARGV[2]\n";
		$x=column('1000g2015aug_eas',@tt);
                $y=column('ExAC_EAS',@tt);
                $z=column('ExonicFunc.refGene',@tt);
		$p=column('Func.refGene',@tt);	
		next;
	}
	my ($allele,$gene)=split(/\;/,$tt[82]);
	push @{$line{$gene}},[$_,$allele];
}
close IN;

foreach my $k(keys %line){
	@items=@{$line{$k}};
	my ($mat,$mat1)=(0,0);
	my ($pat,$pat1)=(0,0);
	my @maternal=();
	my @paternal=();
	my @maternal1=();
	my @paternal1=();
	foreach my $k(0..$#items){
		@tt=split(/\t/,$items[$k][0]);
		if($tt[$x] eq "\."){$tt[$x]=0;}
		if($tt[$y] eq "\."){$tt[$y]=0;}
		if($tt[$x]<$ARGV[1] and $tt[$y]<$ARGV[1] and ($tt[$z]=~/nonsynonymous SNV|stopgain|stoploss|frameshift/ or $tt[$p]=~/splicing/)){
			if($items[$k][1] eq 'Maternal'){
				$mat1++;
                                push @maternal1,$items[$k][0]."\t".$ARGV[2];
                        }else{
                                $pat1++;
                                push @paternal1,$items[$k][0]."\t".$ARGV[2];
                        }
		}
		if($tt[$x]<$ARGV[1] and $tt[$y]<$ARGV[1] and $tt[$z]!~/nonframeshift/ and ($tt[$z]=~/stopgain|stoploss|frameshift/ or $tt[$p]=~/splicing/)){
			if($items[$k][1] eq 'Maternal'){
				$mat++;
				push @maternal,$items[$k][0]."\t".$ARGV[2];
			}else{
				$pat++;
				push @paternal,$items[$k][0]."\t".$ARGV[2];
			}
		}	
	}
	if($mat>=1 && $pat>=1){
		$out=join("\n",(@maternal,@paternal));	print OUT "$out\n";
	}
	if($mat1>=1 && $pat1>=1){
                $out1=join("\n",(@maternal1,@paternal1));  print OUT1 "$out1\n";
        }

}
sub column{
        my ($col,@tit)=@_;
        foreach my $i(0..$#tit){
                if($tit[$i] eq $col){
                        return $i
                }
        }
}

