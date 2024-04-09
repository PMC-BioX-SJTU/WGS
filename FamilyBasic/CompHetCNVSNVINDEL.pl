#! /usr/bin/perl 
use List::MoreUtils qw(uniq);
open IN,$ARGV[0]||die; ## ped file
while(<IN>){
        chomp;
        my @tt=split/\s+/;
        if($tt[2] ne '0'){
                $child=$tt[1];
                $father=$tt[2];
                $mother=$tt[3];
                $fam=$tt[0];
                $sex=$tt[4];
        }
}
open IN1,$ARGV[1]||die;  ### snp and indel vcf file
open CH,">$ARGV[3]"||die;
while(<IN1>){
        chomp;
	if(/^##/){
		print CH "$_\n";
		next;
	};
	my @tt=split/\s+/;
        if(/^#CHROM/){
                foreach my $i(9..$#tt){
                        if($tt[$i] eq $child){
                                $c=$i;
                        }elsif($tt[$i] eq $father){
                                $f=$i;
                        }elsif($tt[$i]eq $mother){
                                $m=$i;
                        }
                } 
		my $out=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
		print CH "$out\n";
                next;
        }
        $gpc=$tt[$c];$gpc=~s/\:.*//g;$gpc=~s/\///g;$gpc=~s/\|//g; $gpc=~s/10/01/g;
        $gpf=$tt[$f];$gpf=~s/\:.*//g;$gpf=~s/\///g;$gpf=~s/\|//g;$gpf=~s/10/01/g;
        $gpm=$tt[$m];$gpm=~s/\:.*//g;$gpm=~s/\///g;$gpm=~s/\|//g;$gpm=~s/10/01/g;
	if($tt[0] !~/X|Y/ && $gpc eq '01' && $gpf eq '00' && $gpm eq '01'){
		$tt[2]="Maternal";
		my $line=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
		$ARM{$tt[0]."\t".$tt[1]}=$line; 
	}elsif($tt[0] !~/X|Y/ && $gpc eq '01' && $gpf eq '01' && $gpm eq '00'){
                $tt[2]="Paternal";
		my $line=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
		$ARF{$tt[0]."\t".$tt[1]}=$line;
	}elsif($tt[0] !~/X|Y/ && $gpc eq '11' && $gpf eq '01' && $gpm eq '00'){
		$tt[2]="Paternal";
		my $line=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
		$ARF2{$tt[0]."\t".$tt[1]}=$line;
	}elsif($tt[0] !~/X|Y/ && $gpc eq '11' && $gpf eq '00' && $gpm eq '01'){
		$tt[2]="Maternal";
		my $line=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
		$ARM2{$tt[0]."\t".$tt[1]}=$line;
	}
}
open IN3,"$ARGV[4]"||die;### Deletion only
while(<IN3>){
        chomp;
        if(/^##/){
                next;
        }
        my @tt=split/\s+/;
        if(/^#CHROM/){
                foreach my $i(9..$#tt){
                        if($tt[$i] eq $child){
                                $c=$i;
                        }elsif($tt[$i] eq $father){
                                $f=$i;
                        }elsif($tt[$i]eq $mother){
                                $m=$i;
                        }
                }
                my $out=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
                next;
        }
        $gpc=$tt[$c];$gpc=~s/\:.*//g;$gpc=~s/\///g;$gpc=~s/\|//g; $gpc=~s/10/01/g;
        $gpf=$tt[$f];$gpf=~s/\:.*//g;$gpf=~s/\///g;$gpf=~s/\|//g;$gpf=~s/10/01/g;
        $gpm=$tt[$m];$gpm=~s/\:.*//g;$gpm=~s/\///g;$gpm=~s/\|//g;$gpm=~s/10/01/g;
	if($tt[7]=~/END=(\d+)\;/){
		$end=$1;
	}
        if($tt[0] !~/X|Y/ && $gpc eq '01' && $gpf eq '00' && $gpm eq '01'){
                $tt[2]="Maternal";
                my $line=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
                $ARM1{$tt[0]}{$tt[1]."\t".$end}=$line;
        }elsif($tt[0] !~/X|Y/ && $gpc eq '01' && $gpf eq '01' && $gpm eq '00'){
                $tt[2]="Paternal";
                my $line=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
                $ARF1{$tt[0]}{$tt[1]."\t".$end}=$line;
        }
}


open IN2,$ARGV[2]||die; #refGene.txt
while(<IN2>){
	chomp;
	next if(/^#/);
	my @tt=split/\s+/;
	my ($indexF,$indexM,$indexF1,$indexM1)=(0,0,0,0);
	my @outF=();
	my @outM=();
	my @outF1=();
	my @outM1=();
	foreach my $i($tt[4]..$tt[5]){
		my $key=$tt[2]."\t".$i;
		##for snp and indel only
		if(defined $ARM{$key}){
			push @outM,$ARM{$key};
			$indexM++;	
		}
		if(defined $ARF{$key}){
			push @outF,$ARF{$key};
			$indexF++;	
		}
		##for DEL only
		if(defined $ARM2{$key}){
                        push @outM1,$ARM2{$key};
                        $indexM1++;
                }
                if(defined $ARF2{$key}){
                        push @outF1,$ARF2{$key};
                        $indexF1++;
                }
	}
	## for DEL and snp, indel
	foreach my $pos (keys %{$ARF1{$tt[2]}}){
		my ($beg,$end)=split(/\t/,$pos);
		if(($tt[4]>=$beg and $tt[4]<=$end)or($tt[5]>=$beg and $tt[5]<=$end)){
			push @outF1,$ARF1{$tt[2]}{$pos};
			$indexF1++;
		}
	}
	foreach my $pos (keys %{$ARM1{$tt[2]}}){
                my ($beg,$end)=split(/\t/,$pos);
                if(($tt[4]>=$beg and $tt[4]<=$end)or($tt[5]>=$beg and $tt[5]<=$end)){
                        push @outM1,$ARM1{$tt[2]}{$pos};
                        $indexM1++;
                }
        }
	## snp and inde output
	if($indexF>0 && $indexM>0){
		my $outF=join("\n",uniq(@outF));
		my $outM=join("\n",uniq(@outM));
		$outF=~s/Paternal/$tt[12]\;$child\;Paternal/g;
		$outM=~s/Maternal/$tt[12]\;$child\;Maternal/g;	
		print CH "$outM\n$outF\n";
	}
	## Deletion output
	if($indexF1>0 && $indexM1>0){
                my $outF1=join("\n",uniq(@outF1));
                my $outM1=join("\n",uniq(@outM1));
                $outF1=~s/Paternal/$tt[12]\;$child\;Paternal/g;
                $outM1=~s/Maternal/$tt[12]\;$child\;Maternal/g;
                print CH "$outM1\n$outF1\n";
        }

}
