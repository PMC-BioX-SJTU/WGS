#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	@tt=split/\t/;
	if(/^Chr/){
		print "$_\n";
		$x=column('1000g2015aug_eas',@tt);
		$y=column('ExAC_EAS',@tt);
		$z=column('ExonicFunc.refGene',@tt);
		$w=column('ExAC_ALL',@tt);
		$p=column('1000g2015aug_all',@tt);
		next;
	}
	if($tt[$x] eq "\."){$tt[$x]=0;}
	if($tt[$y] eq "\."){$tt[$y]=0;}
	if($tt[$w] eq "\."){$tt[$w]=0;}
	if($tt[$p] eq "\."){$tt[$p]=0;}
	if($tt[$p]< $ARGV[1] and $tt[$w]< $ARGV[1] and $tt[$x]<$ARGV[1] and $tt[$y]<$ARGV[1] and $tt[$z] ne '.'){
		print "$_\n";
	}
}close IN;

sub column{
	my ($col,@tit)=@_;
	foreach my $i(0..$#tit){
		if($tit[$i] eq $col){
			return $i
		}
	}
}
