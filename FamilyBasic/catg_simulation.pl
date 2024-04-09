#! /usr/bin/perl -w
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	@tt=split/\t/;
	next if(/transcript/);
	foreach my $i(5..8){
		if($tt[$i]=~/NA/){
			$tt[$i]=-100;
		}
	}
	$syn{$tt[1]}=(10**$tt[4])*2;
	$missense{$tt[1]}=(10**$tt[5])*2;	
	$nonsense{$tt[1]}=(10**$tt[6])*2;
	$splice{$tt[1]}=(10**$tt[7])*2;
	$frameshift{$tt[1]}=(10**$tt[8])*2;
	$gene{$tt[1]}++;
}

open IN1,$ARGV[1]||die;
while(<IN1>){
        chomp;
        @tt=split/\t/;
	if(/^Chr/){
		$x=column('Gene.refGene',@tt);
		next;
	}
	next if (/MUC/);
	if(! defined $gene{$tt[$x]}){
		$syn{$tt[$x]}=(10**(-100))*2;
		$missense{$tt[$x]}=(10**(-100))*2;
		$nonsense{$tt[$x]}=(10**(-100))*2;
		$splice{$tt[$x]}=(10**(-100))*2;
		$frameshift{$tt[$x]}=(10**(-100))*2;
	}
	if(/nonsynonymous/){
		$nonsyn{$tt[$x]}=$missense{$tt[$x]};
		$both{$tt[$x]}=$missense{$tt[$x]};
		$nonsyn2{$tt[$x]}++;
		$both2{$tt[$x]}++;
	}elsif(/frameshift/){
		$loF{$tt[$x]}=$frameshift{$tt[$x]};
		$both{$tt[$x]}=$frameshift{$tt[$x]};
		$loF2{$tt[$x]}++;
		$both2{$tt[$x]}++;
	}elsif(/stop/){
		$loF{$tt[$x]}=$nonsense{$tt[$x]};
		$both{$tt[$x]}=$nonsense{$tt[$x]};
		$both2{$tt[$x]}++;
		$loF2{$tt[$x]}++;
	}elsif(/splicing/){
		$loF{$tt[$x]}=$splice{$tt[$x]};
		$both{$tt[$x]}=$splice{$tt[$x]};
		$both2{$tt[$x]}++;
		$loF2{$tt[$x]}++;
	}elsif(/synonymous/){
		$synony{$tt[$x]}=$syn{$tt[$x]};
		$synony2{$tt[$x]}++;
	}
}
&observation(\%synony2,"synonymous");
&observation(\%nonsyn2,"missense");
&observation(\%loF2,"loF");
&observation(\%both2,"loF+missense");

&simulation(\%synony,"synonymous");
&simulation(\%nonsyn,"missense");
&simulation(\%loF,"loF");
&simulation(\%both,"loF+missense");

sub observation{
	my $type=shift @_;
	my $obs_type=shift @_;
	my $obs=0;
	foreach my $i (keys %$type){
		if($$type{$i} >=2){
			$obs++;
		}
	}
	print "$obs_type\t$obs\n";
}


sub simulation{
	my $catg=shift @_;
	my $type=shift @_;
	my @every_time=();
	my $time=0.01;
	my %gene2=();
	foreach my $i (1..10000) {
		my ($count,$max)=(0,0);
		for my $k (keys %$catg){
			$count=&random_value($$catg{$k},102);
			if($count>=2){
				$gene2{$k}++;
				$time++;
				$max++;
			}
		}
		push @every_time,$max;
	}
	### average expected genes with 2+ DNMs
	my @key=keys %gene2;
	my $average=($#key+1)/$time;

	### Maximum expected genes with 2+ DNMs
	@every_time=sort {$b <=> $a} @every_time;
	print "$type\t$average\t$every_time[0]\n";
}

sub random_value{
	my $pvalue=shift @_;
	my $trio_num=shift @_;
	my $count=0;
	srand();
	foreach(1..$trio_num){
		$prand=rand(1);
		if($prand<=$pvalue){
			$count++;
		}
	}
	return $count;
}

sub column{
        my ($col,@tit)=@_;
        foreach my $i(0..$#tit){
                if($tit[$i] eq $col){
                        return $i
                }
        }
}

