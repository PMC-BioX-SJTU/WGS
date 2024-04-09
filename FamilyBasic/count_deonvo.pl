#! /usr/bin/perl -w
### perl ~/Codes/Pipelines/FamilyBasic/count_deonvo.pl all.annovar.txt 1.07 with >result.txt
### perl ~/Codes/Pipelines/FamilyBasic/count_deonvo.pl all.annovar.txt 1.07 without >result.txt
open IN,$ARGV[0]||die;
open OUT,">$ARGV[3]"||die;
while(<IN>){
	chomp;
	@tt=split/\t/;
	### header title
	if(/#CH/){
		$x=column('ID',@tt);
		$y=column('INFO',@tt);
		$x1=column('1000g2015aug_eas',@tt);
		$x2=column('1000g2015aug_all',@tt);
                $y1=column('ExAC_EAS',@tt);
		$y2=column('ExAC_ALL',@tt);
                $z=column('ExonicFunc.refGene',@tt);
                $p=column('Func.refGene',@tt);
		print OUT "$_\n";
		next
	};
	### without filtering
	if($ARGV[2] eq  "without"){
		if($tt[$x]=~/\,|\;/){
			@mm=split(/\,|\;/,$tt[$x]);
			foreach $m (@mm){
				$m=~s/.*\=//g;
				$info{$m}++;
			}
		}else{
			$info{$tt[$x]}++;
		}
	}else{
	### with filtering
		next if($tt[$y]=~/VQSLOD=(.*)\;culprit/ and $1<$ARGV[1]);
		if($tt[$x]=~/\,|\;/){
                	@mm=split(/\,|\;/,$tt[$x]);
                	foreach $m (@mm){
                        	$m=~s/.*\=//g;
                        	$info{$m}++;
                	}
        	}else{
                	$info{$tt[$x]}++;
        	}
	}	
	###count functional category
	if($tt[$x1] eq "\."){$tt[$x1]=0;}
        if($tt[$y1] eq "\."){$tt[$y1]=0;}
	if($tt[$x2] eq "\."){$tt[$x2]=0;}
        if($tt[$y2] eq "\."){$tt[$y2]=0;}
	if($tt[$x1]< 0.01 and $tt[$y1]<0.01 and $tt[$x2]< 0.01 and $tt[$y2]<0.01){
		next if($tt[$y]=~/VQSLOD=(.*)\;culprit/ and $1<$ARGV[1]);
		### Loss of Func
		if( $tt[$z]!~/nonframeshift/ and ($tt[$z]=~/stopgain|stoploss|frameshift/ or $tt[$p]=~/splicing/)){
			print OUT "$_\n";
			$lof++;
		}
		### missense
		if($tt[$z]=~/nonsynonymous SNV/){
			print OUT "$_\n";
			$mis++;
		}
		### synonymous
		if($tt[$z]=~/synonymous SNV/ and $tt[$z]!~/nonsynonymous SNV/){
			$sys++;
			print OUT "$_\n";
		}
	}
}
print "FamNo\tCount\n";
foreach $i(keys %info){
	print "$i\t$info{$i}\n";
}
print "\n\n";
print "synonymous\tmissense\tLoFunc\n";
print "$sys\t$mis\t$lof\n";
sub column{
        my ($col,@tit)=@_;
        foreach my $i(0..$#tit){
                if($tit[$i] eq $col){
                        return $i
                }
        }
}

