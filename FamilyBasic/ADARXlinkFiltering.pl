#! /usr/bin/perl -w
use List::MoreUtils qw(uniq);
open IN, $ARGV[0]||die;
open OUT,">$ARGV[3].LoF.txt"||die;
open OUT1,">$ARGV[3].all.txt"||die;
open OUTHQ,">$ARGV[3].LoFHQ.txt"||die;
open OUTHQ1,">$ARGV[3].allHQ.txt"||die;
open OUTM,">$ARGV[3].Med.txt"||die;
while(<IN>){
	chomp;
	@tt=split/\t/;
	if(/^Chr/){
		print OUT "$_\t$ARGV[2]\n";
		print OUT1 "$_\t$ARGV[2]\n";
		print OUTHQ "$_\t$ARGV[2]\n";
                print OUTHQ1 "$_\t$ARGV[2]\n";
		print OUTM "$_\t$ARGV[2]\n";
		$x=column('1000g2015aug_eas',@tt);
                $y=column('ExAC_EAS',@tt);
                $z=column('ExonicFunc.refGene',@tt);
		$p=column('Func.refGene',@tt);	
		$q=column('FILTER',@tt);
		$s=column('Polyphen2_HDIV_pred',@tt);
		$t=column('Polyphen2_HVAR_pred',@tt);
		$id=column('ID',@tt);
		next;
	}
	my ($gene,$child,$allele)=split(/\;/,$tt[$id]);###need to change order later
	push @{$line{$gene}},[$_,$allele]
}
close IN;

foreach my $k(keys %line){
	@items=@{$line{$k}};
	my ($mat,$mat1,$matHQ,$matHQ1,$matM)=(0,0,0,0,0);
	my @maternal=();
	my @maternal1=();
	my @maternalHQ=();
	my @maternalHQ1=();
	my @maternalM=();
	foreach my $k(0..$#items){
		@tt=split(/\t/,$items[$k][0]);
		if($tt[$x] eq "\."){$tt[$x]=0;}
		if($tt[$y] eq "\."){$tt[$y]=0;}
		if($tt[$x]<$ARGV[1] and $tt[$y]<$ARGV[1] and $tt[$z]!~/nonframeshift/){
		### raw filtering
		if($tt[$z]=~/nonsynonymous SNV|stopgain|stoploss|frameshift/ or $items[$k][0]=~/\<DEL|\<DUP|<INV/ or ($tt[$p]!~/ncRNA/ and $tt[$p]=~/splicing/)){
                        $mat1++;
			push @maternal1,$items[$k][0]."\t".$ARGV[2];
			if($tt[$q] eq 'PASS'){
				$matHQ1++;
				push @maternalHQ1,$items[$k][0]."\t".$ARGV[2];
			}      
		}
		### medium filtering
		if($tt[$z]=~/nonsynonymous SNV|stopgain|stoploss|frameshift/ or ($items[$k][0]=~/\<DEL/ and $tt[$z]=~/frameshift/)  or($tt[$p] !~/ncRNA/ and $tt[$p]=~/splicing/)){
			next if($tt[$z]=~/nonsynonymous SNV/ and $tt[$s] =~ /B|\./ and $tt[$t] =~ /B|\./); 
			if($tt[$q] eq 'PASS'){
				$matM++;
				push @maternalM,$items[$k][0]."\t".$ARGV[2];
			}
		}
		## high quality filtering
		if($tt[$z]=~/stopgain|stoploss|frameshift/ or  ($items[$k][0]=~/\<DEL/ and $tt[$z]=~/frameshift/ )or ($tt[$p]=~/splicing/ and $tt[$p]!~/ncRNA/)){
			$mat++;
			push @maternal,$items[$k][0]."\t".$ARGV[2];
			if($tt[$q] eq 'PASS'){
				$matHQ++;
				push @maternalHQ,$items[$k][0]."\t".$ARGV[2];
			}
		}	
	}}
	if($mat>=1 ){
		$out=join("\n",uniq(@maternal));	print OUT "$out\n";
	}
	if($mat1>=1){
                $out1=join("\n",uniq(@maternal1));  print OUT1 "$out1\n";
        }
	if($matHQ>=1){
                $outHQ=join("\n",uniq(@maternalHQ));  print OUTHQ "$outHQ\n";
        }
        if($matHQ1>=1){
                $outHQ1=join("\n",uniq(@maternalHQ1));  print OUTHQ1 "$outHQ1\n";
        }
	if($matM>=1){
		$outM=join("\n",uniq(@maternalM));  print OUTM "$outM\n";
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
