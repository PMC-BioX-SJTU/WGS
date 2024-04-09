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
open IN1,$ARGV[1]||die;  ### pileup vcf file
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
	}
}

open IN2,$ARGV[2]||die; #refGene.txt
while(<IN2>){
	chomp;
	next if(/^#/);
	my @tt=split/\s+/;
	my ($indexF,$indexM)=(0,0);
	my @outF=();
	my @outM=();
	foreach my $i($tt[4]..$tt[5]){
		my $key=$tt[2]."\t".$i;
		if(defined $ARM{$key}){
			push @outM,$ARM{$key};
			$indexM++;	
		}
		if(defined $ARF{$key}){
			push @outF,$ARF{$key};
			$indexF++;	
		}
	}
	#if($indexF>0 && $indexM>0){
	#	my $outF=join("\n",uniq(@outF));
	#	my $outM=join("\n",uniq(@outM));
	#	$outF=~s/Paternal/$tt[12]\;$child\;Paternal/g;
	#	$outM=~s/Maternal/$tt[12]\;$child\;Maternal/g;	
	#	print CH "$outM\n$outF\n";
	#}
	if($indexF>0){
               my $outF=join("\n",uniq(@outF));
               $outF=~s/Paternal/$tt[12]\;$child\;Paternal/g;
               print CH "$outF\n";
	}
        if($indexM>0){
               my $outM=join("\n",uniq(@outM));
               $outM=~s/Maternal/$tt[12]\;$child\;Maternal/g;
               print CH "$outM\n";
        }

}
