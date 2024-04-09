#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	if(/^#/){print "$_\n";next};
	$line=$_;
        @tt=split(/\s+/,$line);
	$index=0;
	foreach my $i(9..$#tt){
		@ss=split(":",$tt[$i]);
		if(! defined $ss[2] or $ss[2] eq "\."){$ss[2]=0}
		if($ss[2]>$ARGV[1]){
			$index++;
		}
	}
	if($index> $ARGV[2]){
		print "$line\n";
	}
}
