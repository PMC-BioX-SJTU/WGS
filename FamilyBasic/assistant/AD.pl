#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	next if(/^#/);
	@tt=split/\s+/;
	my ($index1,$index2)=(0,0);
	foreach($tt[9]){
		@ss=split/:/;
		if($ss[0]eq '0/0'){
			$index1++;
		}
	}
	foreach($tt[10],$tt[11],$tt[12]){
                @ss=split/:/;
                if($ss[0]eq '1/1'){
                        $index2++;
                }
        }
	if($index1==1 && $index2==3){
		print "$_\n";
		
	}
	$index1=0;$index2=0;
}
