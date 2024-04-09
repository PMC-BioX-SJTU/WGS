#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	@tt=split/\s+/;
	push @{$fam{$tt[0]}},$_;
}
foreach(keys %fam){
	$key=$_;
	@trio=@{$fam{$key}};
#	if($#trio==2){
#		open OUT,">$ARGV[1]/$key.ped"||die;
#		foreach(@trio){
#			print OUT "$_\n";
#		}
#		close OUT;
#	}
#	elsif($#trio!=2){
		#%child={};
		@parent=();
		foreach(@trio){
			@tt=split/\s+/;
			if($tt[2] ne '0'){
				$child{$tt[1]}=$_;
			}else{
				push @parent,$_;
			}
		}
		foreach (keys %child){
			$k=$_;
			$out=join("\n",(@parent,$child{$_}));
			open OUT,">$ARGV[1]/$k.ped"||die;
			print OUT "$out\n";
			close OUT;
			delete $child{$k};
		}
#	}
}
