#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	### header title
	@items=split(/\t/);
	if(/#CH/){
		$x=column('ID',@items);
		next;
	}
	@tt=split(/\;/,$items[$x]);
	my $k=$tt[0]."\t".$tt[1];
	if(!defined $fam{$k}){
		$gene{$tt[0]}++;
		$fam{$k}++;
	}
}
print "gene\tcount\n";
foreach my $key(keys %gene){
	print "$key\t$gene{$key}\n";

}
sub column{
        my ($col,@tit)=@_;
        foreach my $i(0..$#tit){
                if($tit[$i] eq $col){
                        return $i
                }
        }
}
