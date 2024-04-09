open IN, $ARGV[0]||die;
while(<IN>){
	chomp;
	next if(/NA/);
	if(/^CHR/){
		print "$_\n";
	}else{
		@tt=split/\s+/;
		push @lst,[$_,$tt[12]];
	}
}
close IN;
@lst= sort {$a->[1] <=> $b->[1]} @lst;
foreach my $i(0..$#lst){
	print "$lst[$i][0]\n"; 
}
