#! /usr/bin/perl 
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
#open CH,">$ARGV[2]"||die;
while(<IN1>){
        chomp;
	if(/^##/){
		print "$_\n";
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
		print  "$out\n";
                next;
        }
        $gpc=$tt[$c];$gpc=~s/\:.*//g;$gpc=~s/\///g;$gpc=~s/\|//g; $gpc=~s/10/01/g;
        $gpf=$tt[$f];$gpf=~s/\:.*//g;$gpf=~s/\///g;$gpf=~s/\|//g;$gpf=~s/10/01/g;
        $gpm=$tt[$m];$gpm=~s/\:.*//g;$gpm=~s/\///g;$gpm=~s/\|//g;$gpm=~s/10/01/g;
	if($tt[0] !~/X|Y|M/ && $gpc eq '01' && $gpf eq '00' && $gpm eq '00'){
		$tt[2]=$fam;
		my $line=join("\t",(@tt[0..8],$tt[$f],$tt[$m],$tt[$c]));
		print "$line\n";
		#$ARM{$tt[0]."\t".$tt[1]}=$line; 
	}
}

