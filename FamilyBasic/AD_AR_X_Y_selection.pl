#! /usr/bin/perl -w
### this script defaultly based on the affected mutation is 1
### if u want to make the mutation 0 as affected,  add $ARGV[3]
my ($pat1,$mat1,$son1,$dau1)=("","","","");
my ($pat2,$mat2,$son2,$dau2)=("","","","");
open IN,$ARGV[0]||die; ## ped file
open IND,">$ARGV[2].inm"||die;
while(<IN>){
        chomp;
        my @tt=split/\s+/;
	if($tt[2] ne '0' && $tt[4]==1 && $tt[5] ==1){
		$son1=$tt[1];
	}elsif($tt[2] ne '0' && $tt[4]==1 && $tt[5] ==2){
		$son2=$tt[1];
	}elsif($tt[2] ne '0' && $tt[4]==2 && $tt[5] ==1){
		$dau1=$tt[1];
	}elsif($tt[2] ne '0' && $tt[4]==2 && $tt[5] ==2){
		$dau2=$tt[1];
	}elsif($tt[2]==0 && $tt[4]==1 && $tt[5] ==1){
		$pat1=$tt[1];
	}elsif($tt[2]==0 && $tt[4]==1 && $tt[5] ==2){
		$pat2=$tt[1];
	}elsif($tt[2]==0 && $tt[4]==2 && $tt[5] ==1){
		$mat1=$tt[1];
	}elsif($tt[2]==0 && $tt[4]==2 && $tt[5] ==2){
		$mat2=$tt[1];
	}
}
if($son2 ne "" && $pat1 ne "" && $mat1 ne ""){
	  print IND "AR\tXLNKR\tCompHet\n";
}elsif($dau2 ne "" && $pat1 ne "" && $mat1 ne ""){
	  print IND "AR\tCompHet\n";
}elsif($son2 ne "" && $pat1 ne "" && $mat2 ne ""){
	  print IND "AD\tAR\tXLNKD\tXLNKR\n";
}elsif($dau2 ne "" && $pat1 ne "" && $mat2 ne ""){
          print IND "AD\tAR\tXLNKD\n";
}elsif($son2 ne "" && $pat2 ne "" && $mat1 ne ""){
          print IND "AD\tAR\tXLNKR\tYLNK\n";
}elsif($dau2 ne "" && $pat2 ne "" && $mat1 ne ""){
          print IND "AD\tAR\tXLNKD\tXLNKR\n";
}
close IN;

my ($p1,$m1,$s1,$d1)=(0,0,0,0);
my ($p2,$m2,$s2,$d2)=(0,0,0,0);
open IN1,$ARGV[1]||die;  ### pileup vcf file
open AD,">$ARGV[2].primary.AD.vcf"||die;
open AR,">$ARGV[2].primary.AR.vcf"||die;
open XLNKD,">$ARGV[2].primary.XLNKD.vcf"||die;
open XLNKR,">$ARGV[2].primary.XLNKR.vcf"||die;
open YLNK,">$ARGV[2].primary.YLNK.vcf"||die;


if($ARGV[3] ==0){
	open AD2,">$ARGV[2].secondary.AD.vcf"||die;
	open AR2,">$ARGV[2].secondary.AR.vcf"||die;
	open XLNKD2,">$ARGV[2].secondary.XLNKD.vcf"||die;
	open XLNKR2,">$ARGV[2].secondary.XLNKR.vcf"||die;
	open YLNK2,">$ARGV[2].secondary.YLNK.vcf"||die;
}
while(<IN1>){
        chomp;
        $line=$_;
        if(/^##/){
		print AD "$line\n";	
		print AR "$line\n";
		print XLNKD "$line\n";
		print XLNKR "$line\n";
		print YLNK "$line\n";

		if($ARGV[3] ==0){
			print AD2 "$line\n";	
			print AR2 "$line\n";
			print XLNKD2 "$line\n";
			print XLNKR2 "$line\n";
			print YLNK2 "$line\n";
		}
		next;
	};
        my @tt=split/\s+/;
	 if(/^#CHROM/){
		foreach my $i(9..$#tt){
			if($tt[$i] eq $son1){
				$s1=$i;
			}elsif($tt[$i] eq $son2){
				$s2=$i;
			}elsif($tt[$i]eq $dau1){
				$d1=$i;
			}elsif($tt[$i]eq $dau2){
				$d2=$i;
			}elsif($tt[$i] eq $pat1){
				$p1=$i;
			}elsif($tt[$i] eq $pat2){
				$p2=$i;
			}elsif($tt[$i]eq $mat1){
				$m1=$i;
			}elsif($tt[$i]eq $mat2){
				$m2=$i;
			}
		}
        print AD "$line\n";	
		print AR "$line\n";
		print XLNKD "$line\n";
		print XLNKR "$line\n";
		print YLNK "$line\n";		
		
		if($ARGV[3] ==0){		
			print AD2 "$line\n";	
			print AR2 "$line\n";
			print XLNKD2 "$line\n";
			print XLNKR2 "$line\n";
			print YLNK2 "$line\n";
		}	
		next;
	}
	if($son2 ne "" && $pat1 ne "" && $mat1 ne ""){
		my ($gpc,$gpf,$gpm)=genotype($s2,$p1,$m1,@tt);#print "$gpc\t$gpm\t$gpf\n";
		if($tt[0] !~/X|Y/ && $gpc eq '11' && $gpf eq '01' && $gpm eq '01'){## AR
                	print AR "$line\n";
        	}
		if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '00' && $gpf eq '01' && $gpm eq '01'){## AR
                	print AR2 "$line\n";
        	}
		if($tt[0] =~/X/ && $gpc eq '11' && $gpf eq '00' && $gpm eq '01'){ ## X_link D.
			print  XLNKR "$line\n";
		}
        	if($ARGV[3] ==0 && $tt[0] =~/X/ && $gpc eq '00' && $gpf eq '11' && $gpm eq '01'){ ## X_link D.
			print  XLNKR2 "$line\n";
		}
	}elsif($dau2 ne "" && $pat1 ne "" && $mat1 ne ""){
		my ($gpc,$gpf,$gpm)=genotype($d2,$p1,$m1,@tt);#print "$gpc\t$gpm\t$gpf\n";
		if($tt[0] !~/X|Y/ && $gpc eq '11' && $gpf eq '01' && $gpm eq '01'){## AR
                	print AR "$line\n";
        	}
        	if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '00' && $gpf eq '01' && $gpm eq '01'){## AR
                	print AR2 "$line\n";
        	}
	}elsif($son2 ne "" && $pat1 ne "" && $mat2 ne ""){
		my ($gpc,$gpf,$gpm)=genotype($s2,$p1,$m2,@tt);
		if($tt[0] !~/X|Y/ && $gpc eq '01' && $gpf eq '00' && ($gpm eq '01'||$gpm eq '11')){## AD
                	print AD "$line\n";
        	}
		if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '01' && $gpf eq '11' && ($gpm eq '01'||$gpm eq '00')){## AD
                	print AD2 "$line\n";
        	}
		if($tt[0] !~/X|Y/ && $gpc eq '11' && $gpf eq '01' && $gpm eq '11'){## AR
                	print AR "$line\n";
        	}
		if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '00' && $gpf eq '01' && $gpm eq '00'){## AR
                	print AR2 "$line\n";
        	}
		if($tt[0] =~/X/ && $gpc eq '11' && $gpf eq '00' && ($gpm eq '01'||$gpm eq '11')){ ## X_link D.
			print  XLNKD "$line\n";
		}
		if($ARGV[3] ==0 && $tt[0] =~/X/ && $gpc eq '00' && $gpf eq '11' && ($gpm eq '01'||$gpm eq '00')){ ## X_link D.
			print  XLNKD2 "$line\n";
		}
		if($tt[0] =~/X/ && $gpc eq '11' && $gpf eq '00' && $gpm eq '11'){ ## X_link R.
			print  XLNKR "$line\n";
		}
		if($ARGV[3] ==0 && $tt[0] =~/X/ && $gpc eq '00' && $gpf eq '11' && $gpm eq '00'){ ## X_link R.
			print  XLNKR3 "$line\n";
		}
	}elsif($dau2 ne "" && $pat1 ne "" && $mat2 ne ""){
		my ($gpc,$gpf,$gpm)=genotype($d2,$p1,$m2,@tt);
		if($tt[0] !~/X|Y/ && $gpc eq '01' && $gpf eq '00' && ($gpm eq '01'||$gpm eq '11')){## AD
                	print AD "$line\n";
        	}
		if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '01' && $gpf eq '11' && ($gpm eq '01'||$gpm eq '00')){## AD
                	print AD2 "$line\n";
        	}
		
		if($tt[0] =~/X/ && $gpc eq '01' && $gpf eq '00' && ($gpm eq '01'||$gpm eq '11')){ ## X_link D.
			print  XLNKD "$line\n";
		}
		if($ARGV[3] ==0 && $tt[0] =~/X/ && $gpc eq '01' && $gpf eq '11' && ($gpm eq '01'||$gpm eq '00')){ ## X_link D.
			print  XLNKD2 "$line\n";
		}
		if($tt[0] !~/X|Y/ && $gpc eq '11' && $gpf eq '01' && $gpm eq '11'){## AR
                	print AR "$line\n";
        	}
		if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '00' && $gpf eq '01' && $gpm eq '00'){## AR
                	print AR2 "$line\n";
        	}
	}elsif($son2 ne "" && $pat2 ne "" && $mat1 ne ""){
		my ($gpc,$gpf,$gpm)=genotype($s2,$p2,$m1,@tt);
		if($tt[0] !~/X|Y/ && $gpc eq '01' && ($gpf eq '01' ||$gpf eq '11') && $gpm eq '00'){## AD
                	print AD "$line\n";
        	}
         	if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '01' && ($gpf eq '01' ||$gpf eq '00') && $gpm eq '11'){## AD
                	print AD2 "$line\n";
        	}
         
		if($tt[0] =~/Y/ && $gpc eq '11' && $gpf eq '11' ){ ## X_link D.
			print  YLNK "$line\n";
		}
		if($ARGV[3] ==0 && $tt[0] =~/Y/ && $gpc eq '00' && $gpf eq '00' ){ ## X_link D.
			print  YLNK2 "$line\n";
		}
		if($tt[0] !~/X|Y/ && $gpc eq '11' && $gpf eq '11' && $gpm eq '01'){## AR
                	print AR "$line\n";
        	}
		if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '00' && $gpf eq '00' && $gpm eq '01'){## AR
                	print AR2 "$line\n";
        	}
		
		if($tt[0] =~/X/ && $gpc eq '11' && $gpf eq '11' && $gpm eq '01'){ ## X_link R.
			print  XLNKR "$line\n";
		}
		if($ARGV[3] ==0 && $tt[0] =~/X/ && $gpc eq '00' && $gpf eq '00' && $gpm eq '01'){ ## X_link R.
			print  XLNKR2 "$line\n";
		}
	}elsif($dau2 ne "" && $pat2 ne "" && $mat1 ne ""){
		my ($gpc,$gpf,$gpm)=genotype($d2,$p2,$m1,@tt);
		if($tt[0] !~/X|Y/ && $gpc eq '01' && ($gpf eq '01'||$gpf eq '11') && $gpm eq '00'){## AD
                	print AD "$line\n";
        	}
		if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '01' && ($gpf eq '01'||$gpf eq '00') && $gpm eq '11'){## AD
                	print AD2 "$line\n";
        	}
		if($tt[0] !~/X|Y/ && $gpc eq '11' && $gpf eq '11' && $gpm eq '01'){## AR
                	print AR "$line\n";
        	}
		if($ARGV[3] ==0 && $tt[0] !~/X|Y/ && $gpc eq '00' && $gpf eq '00' && $gpm eq '01'){## AR
                	print AR2 "$line\n";
        	}
		if($tt[0] =~/X/ && $gpc eq '01' && $gpf eq '11' && $gpm eq '00'){ ## X_link D.
			print  XLNKD "$line\n";
		}
		if($ARGV[3] ==0 && $tt[0] =~/X/ && $gpc eq '01' && $gpf eq '00' && $gpm eq '11'){ ## X_link D.
			print  XLNKD2 "$line\n";
		}
		if($tt[0] =~/X/ && $gpc eq '11' && $gpf eq '11' && $gpm eq '01'){ ## X_link R.
			print  XLNKR "$line\n";
		}
		if($ARGV[3] ==0 && $tt[0] =~/X/ && $gpc eq '00' && $gpf eq '00' && $gpm eq '01'){ ## X_link R.
			print  XLNKR3 "$line\n";
		}
	}
}
close IN1;
sub genotype{
	my ($c,$f,$m,@tt)=@_;
	foreach my $i ($c,$f,$m){
		$tt[$i]=~s/\:.*//g;
		$tt[$i]=~s/\///g;
		$tt[$i]=~s/\|//g;
		$tt[$i]=~s/10/01/g;
	}
	return ($tt[$c],$tt[$f],$tt[$m]);
}
