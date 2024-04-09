#! /usr/bin/perl -w
## perl missense_anno.pl miss_baddness.txt gamma.txt segment.txt Denovo_filtered.txt > MPC_cand.txt
$aa{'Ala'}='A';		
$aa{'Arg'}='R';	
$aa{'Asn'}='N';		
$aa{'Asp'}='D';	
$aa{'Cys'}='C';
$aa{'Gln'}='Q';
$aa{'Glu'}='E';		
$aa{'Gly'}='G';
$aa{'His'}='H';		
$aa{'Ile'}='I';	
$aa{'Leu'}='L';
$aa{'Lys'}='K';
$aa{'Met'}='M';
$aa{'Phe'}='F';
$aa{'Pro'}='P';
$aa{'Ser'}='S';
$aa{'Thr'}='T';
$aa{'Trp'}='W';
$aa{'Tyr'}='Y';
$aa{'Val'}='V';
%mis_baddness=();
%gamma=();
open IN,$ARGV[0]||die;##S10
while(<IN>){
	chomp;
	next if(/from/);
	@tt=split(/\t+/);
	$key=$aa{$tt[0]}.$aa{$tt[1]};
	$mis_baddness{$key}=$tt[9];

}
open IN1,$ARGV[1]||die;##gamma.txt.cand.txt
while(<IN1>){
        chomp;
        next if(/transcript/);
        @tt=split(/\t+/);
	push @{$gamma{'chr'.$tt[2]}},[$tt[4],$tt[5],$tt[11]];
}
open IN2,$ARGV[2]||die;##segment.txt.cand.txt
while(<IN2>){
        chomp;
        next if(/transcript/);
        @tt=split(/\s+/);
        push @{$gamma{'chr'.$tt[2]}},[$tt[4],$tt[5],$tt[8]];
}
open IN3,$ARGV[3]||die;## denovo list
while(<IN3>){
        chomp;
        my @tt=split(/\t+/);
	$line=$_;
	if(/Chr/){
		$x=column('AAChange.refGene',@tt);
		$y=column('Polyphen2_HVAR_score',@tt);
		print "$line\tMPC\n";
		next;
	}
	next if($line !~ /nonsynonymous SNV/);
	my @cand=@{$gamma{$tt[0]}};
	@cand= sort {$a->[0] <=> $b->[0]} @cand;
	foreach my $i(0..$#cand){
		if($tt[1]>=$cand[$i][0] and $tt[2]<=$cand[$i][1]){
			if($tt[$x]=~/p\.(.)\d+(.)$/){
				$key=$1.$2;
			}
			$mpc=$cand[$i][2]+$mis_baddness{$key}+$tt[$y];
			print "$line\t$mpc\n";		
		}
	}
}

sub column{
        my ($col,@tit)=@_;
        foreach my $i(0..$#tit){
                if($tit[$i] eq $col){
                        return $i
                }
        }
}

