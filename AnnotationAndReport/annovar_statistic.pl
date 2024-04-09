#!/usr/bin/perl -w
use strict;
use lib '/usr/local/share/perl5/';
use lib '/home/zhanghk/software/perl/lib/share/perl5';
use lib '/opt/perl/lib/site_perl/5.14.2';
use Statistics::R;
use Getopt::Long;
use File::Basename;

my ($annovar, $stat,  $help);

GetOptions(
	'a=s'  => \$annovar,
	's=s' => \$stat,
	'help' => \$help
);

my $usage=<<INFO;
Usage:
    perl $0 [options]
Options:
           -a   <file>    input annovar annotation result, SampleName_.*multianno.txt
           -s   <file>    output file name of statistic information
Example:
	perl annovar_statistic.pl -a cnv.txt -s cnv.stat -t CNV
Authro:
	
INFO

if($help || !$annovar || !$stat){
	die $usage;
}
my $dir = dirname $annovar;
my $sample_name = $stat;

my (%func_hash_snp,%func_hash_indel,%draw_snp,%reference_snp,%reference_indel) = ();

### setting function location and fuction in exonic region
my @func_loc=('exonic','splicing','exonic;splicing','intergenic','intronic','UTR3','UTR5','UTR5;UTR3','ncRNA_exonic','ncRNA_splicing','ncRNA_exonic;splicing','ncRNA_intronic','ncRNA_UTR5','ncRNA_UTR3','upstream','downstream','upstream;downstream','unknown');
my @exonic_func=('synonymous SNV','nonsynonymous SNV','frameshift deletion','frameshift insertion','nonframeshift deletion','nonframeshift insertion','stopgain','stoploss');


my ($total_snp,$total_indel,$hom_snp,$hom_indel,$het_snp,$het_indel,$ti,$tv) = (0,0,0,0,0,0,0,0);
my (@indel,@insertion,@deletion);

open IN,"<$annovar" or die "input annovar result open error! $!\n";
my $head = <IN>;
$head =~ /avsnp(\d+).*/;
my $dbsnp_version = "dbsnp$1";
my @ref_info=("$dbsnp_version","1000G","esp6500","ExAC","Reference_all","Novel");

while(my $line = <IN>){
	chomp $line;
	my @fields = split /\t+/,$line;
	my ($chr, $start, $end, $ref, $alt, $func, $gene, $exonic, $g1000, $dbsnp,$exac,$esp6500) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4], $fields[5], $fields[6], $fields[8], $fields[18], $fields[11],$fields[27],$fields[24]);

	if($ref ne "-" && $alt ne "-"){	
		$func_hash_snp{$func}++; ## func.refGene
		$total_snp++;
		### compared with 1000G and dbsnp
		if($g1000 ne "." || $dbsnp ne "." || $exac ne "." || $esp6500 ne "."){$reference_snp{"Reference_all"}++;}
		if($g1000 ne "."){$reference_snp{"1000G"}++;}
		if($exac ne "."){$reference_snp{"ExAC"}++;}
		if($esp6500 ne "."){$reference_snp{"esp6500"}++;}
		if($dbsnp ne "."){$reference_snp{"$dbsnp_version"}++;}
		if($dbsnp eq "." and $g1000 eq "." and $esp6500 eq "." and $exac eq "."){$reference_snp{"Novel"}++;}

		### Exonic is one type of func.refGene, Exonic function compared with refGene
		if($exonic ne "."){
			$func_hash_snp{$exonic}++;
		}
		
		if($line =~ /0\/0\:?/ or $line =~ /1\/1:?/){
			$hom_snp++;
		}else{
			$het_snp++;
		}

		## calculate ti tv ratio
		if((($ref eq "A") and ($alt eq "G")) or (($ref eq "G") and ($alt eq "A")) or (($ref eq "C") and ($alt eq "T")) or (($ref eq "T") and ($alt eq "C"))){
			$ti++;
		}else{
			$tv++;
		}
		### SNP change
		if((($ref eq "T") and ($alt eq "C")) or (($ref eq "A") and ($alt eq "G"))){
			$draw_snp{"T:A->C:G"}++;
		}
		if((($ref eq "T") and ($alt eq "A")) or (($ref eq "A") and ($alt eq "T"))){
			$draw_snp{"T:A->A:T"}++;
		}
		if((($ref eq "T") and ($alt eq "G")) or (($ref eq "A") and ($alt eq "C"))){
			$draw_snp{"T:A->G:C"}++;
		}
		if((($ref eq "C") and ($alt eq "T")) or (($ref eq "G") and ($alt eq "A"))){
			$draw_snp{"C:G->T:A"}++;
		}
		if((($ref eq "C") and ($alt eq "G")) or (($ref eq "G") and ($alt eq "C"))){
			$draw_snp{"C:G->G:C"}++;
		}
		if((($ref eq "C") and ($alt eq "A")) or (($ref eq "G") and ($alt eq "T"))){
			$draw_snp{"C:G->A:T"}++;
		}
	}
	if($ref eq "-" || $alt eq "-"){
		$func_hash_indel{$func}++;
		$total_indel++;
		### compared with 1000G and dbsnp
                if($g1000 ne "." || $dbsnp ne "." || $exac ne "." || $esp6500 ne "."){$reference_indel{"Reference_all"}++;}
                if($g1000 ne "."){$reference_indel{"1000G"}++;}
                if($exac ne "."){$reference_indel{"ExAC"}++;}
                if($esp6500 ne "."){$reference_indel{"esp6500"}++;}
                if($dbsnp ne "."){$reference_indel{"$dbsnp_version"}++;}
                if($dbsnp eq "." and $g1000 eq "." and $esp6500 eq "." and $exac eq "."){$reference_indel{"Novel"}++;}

		### Exonic is one type of func.refGene, Exonic function compared with refGene
		if($exonic ne "."){
			$func_hash_indel{$exonic}++;
		}
		### Hom and Het
		if($line =~ /0\/0\:?/ or $line =~ /1\/1:?/){
			$hom_indel++;
		}else{
			$het_indel++;
		}
		### Indel length distribution
		if($ref eq "-"){
			my $len = length($alt);
			$insertion[$len]++;
			$indel[$len]++;
		}elsif($alt eq '-'){
			my $len = length($ref);
			$deletion[$len]++;
			$indel[$len]++;
		}else{
			my $len1=length($ref);
			my $len2=length($alt);
			if($len1 < $len2){
				my $len=$len2-$len1;
				$insertion[$len]++;
				$indel[$len]++;
			}else{
				my $len=$len1-$len2;
				$deletion[$len]++;
				$indel[$len]++;
			}
		}
	}
}
close IN;

open OUT,">$dir/$stat\_statistic.xls" or die "STAT create error! $!\n";
print OUT "Type\tSNP\tINDEL\n";
print OUT "Total Variants\t$total_snp\t$total_indel\n";
my $rate1=100*$het_snp/$total_snp;
my $rate2=100*$hom_snp/$total_snp;
my $rate3=100*$het_indel/$total_indel;
my $rate4=100*$hom_indel/$total_indel;
printf OUT "%s\t%d(%4.2f%%)\t%d(%4.2f%%)\n", "Het", $het_snp, $rate1,$het_indel,$rate3;
printf OUT "%s\t%d(%4.2f%%)\t%d(%4.2f%%)\n", "Hom", $hom_snp, $rate2,$hom_indel,$rate4;


foreach (@ref_info){
	if (! defined $reference_snp{$_}){
		$reference_snp{$_}=0;
	}
	if (! defined $reference_indel{$_}){
		$reference_indel{$_}=0;
	}
	my $rate_snp = 100*$reference_snp{$_}/$total_snp;
	my $rate_indel = 100*$reference_indel{$_}/$total_indel;
	printf OUT "%s\t%d(%4.2f%%)\t%d(%4.2f%%)\n", $_, $reference_snp{$_}, $rate_snp,$reference_indel{$_}, $rate_indel;
}

foreach (@func_loc){
	if (! defined $func_hash_snp{$_}){
                $func_hash_snp{$_}=0;
        }
	if (! defined $func_hash_indel{$_}){
                $func_hash_indel{$_}=0;
        }
	my $rate_snp = 100*$func_hash_snp{$_}/$total_snp;
	my $rate_indel = 100*$func_hash_indel{$_}/$total_indel;
	printf OUT "%s\t%d(%4.2f%%)\t%d(%4.2f%%)\n", $_, $func_hash_snp{$_}, $rate_snp,$func_hash_indel{$_}, $rate_indel;
	delete $func_hash_snp{$_};
	delete $func_hash_indel{$_};

}

foreach (@exonic_func){
	if (! defined $func_hash_snp{$_}){
                $func_hash_snp{$_}=0;
        }
	if (! defined $func_hash_indel{$_}){
                $func_hash_indel{$_}=0;
        }
	my $rate_snp = 100*$func_hash_snp{$_}/$total_snp;
	my $rate_indel = 100*$func_hash_indel{$_}/$total_indel;
	printf OUT "%s\t%d(%4.2f%%)\t%d(%4.2f%%)\n", $_, $func_hash_snp{$_}, $rate_snp,$func_hash_indel{$_}, $rate_indel;
	delete $func_hash_snp{$_};
	delete $func_hash_indel{$_};
}	

foreach my $key (sort keys %func_hash_snp){
	if(defined $func_hash_indel{$key}){
		print OUT "$key\t$func_hash_snp{$key}\t$func_hash_indel{$key}\n";
		delete $func_hash_indel{$key};
	}else{
		print OUT "$key\t$func_hash_snp{$key}\t-\n";
	}
	delete $func_hash_snp{$key};
}

foreach my $key (sort keys %func_hash_indel){
	if(defined $func_hash_snp{$key}){
		print OUT "$key\t$func_hash_snp{$key}\t$func_hash_indel{$key}\n";
		delete $func_hash_snp{$key};
	}else{
		print OUT "$key\t-\t$func_hash_indel{$key}\n";
	}
	delete $func_hash_indel{$key};
}



my $rate_ti_tv = $ti/$tv;
my $rate_snp_indel=$total_snp/$total_indel;
printf OUT "Ti/Tv|SNP/INDEL Rate\t%4.2f\t%4.2f\n",$rate_ti_tv,$rate_snp_indel;

close OUT;

#my $dir = dirname $stat;
my $picture_snv = $dir."/".$sample_name."\.SNP\_Spectrum.png";
my $R = Statistics::R->new();
$R->set( 'number', [$draw_snp{"T:A->C:G"},$draw_snp{"T:A->A:T"},$draw_snp{"T:A->G:C"}, $draw_snp{"C:G->T:A"},$draw_snp{"C:G->G:C"},$draw_snp{"C:G->A:T"}] );
$R->set( 'name', ["T:A->C:G", "T:A->A:T", "T:A->G:C", "C:G->T:A", "C:G->G:C", "C:G->A:T"]);
$R->run(qq`options(bitmapType='cairo')`); 
$R->run(qq`png("$picture_snv", width=700,heigh=700)`);
$R->run(qq`bar=barplot(number,ylim=c(0,1.35*max(number)),names.arg=name,cex=1.3,xlab="Mutation type",ylab="Mutation number",main="Mutation Spectrum",col=c(rgb(186,85,211,max=255),rgb(65,105,225,max=255),rgb(105,105,105,max=255),rgb(154,205,50,max=255),rgb(0,139,139,max=255),rgb(218,165,32,max=255)),cex.lab=1.4,cex.axis=1.3,cex.main=2,space=0.9,font.lab=2)`);
$R->run(qq`legend("topright",legend=name,pch=15,cex=1.3,bty="n",col=c(rgb(186,85,211,max=255),rgb(65,105,225,max=255),rgb(105,105,105,max=255),rgb(154,205,50,max=255),rgb(0,139,139,max=255),rgb(218,165,32,max=255)))`);
$R->run(qq`text(bar, number,number,adj=c(100,100),cex=1.3,font=2,pos=3, )`);
$R->run(qq`dev.off()`);
$R->stop();




my @len;
my @indels;
my @insertions;
my @deletions;
my $flag = 1;
for(my $i = 1; $i<scalar(@indel) and $flag <=20; $i++){
#	unless($indel[$i] and $insertion[$i] and $deletion[$i]){
#		next;
#	}
	$flag++;
	push @len, $i;
	push @indels, $indel[$i]?$indel[$i]:0;
	push @insertions, $insertion[$i]?$insertion[$i]:0;
	push @deletions, $deletion[$i]?$deletion[$i]:0;
}
my $picture_indel = $dir."/".$sample_name."\.InDel_Distribution.png";
$R = Statistics::R->new();
$R->set('Len', [@len]);
$R->set('Indel', [@indels]);
$R->set('Insertion', [@insertions]);
$R->set('Deletion', [@deletions]);
$R->run(qq`options(bitmapType='cairo')`);
$R->run(qq`png("$picture_indel",width=800,heigh=500)`);
$R->run(qq`par(mar=c(5,5,5,5),mgp=c(3.5,1,0))`);
$R->run(qq`barplot(rbind(Insertion,Deletion,Indel),beside=T,names.arg=Len,col=c(rgb(255,102,204,max=255),rgb(51,204,204,max=255),rgb(0,102,204,max=255)),ylim=c(0,max(Indel)*1.2),cex.lab=1.5,font.lab=2,cex.axis=1.2,cex.main=2,las=1,xlab="InDel length(bp)",ylab="Number",main="InDel length distribution (All)")`);
$R->run(qq`legend("right",c("Insertion","Deletion","Indels"),lwd=3,lty=1,bty="n",cex=1.5,col=c(rgb(255,102,204,max=255),rgb(51,204,204,max=255),rgb(0,102,204,max=255)))`);
$R->run(qq`dev.off()`);

