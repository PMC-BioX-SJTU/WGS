#! /usr/bin/perl -w

$proj='bjxk';
$pedigree_name='bjxk_lung_cancer';
print "#creator /proj/100genomes/bin/generateGenotype.pl -pedFile $proj-hg19.ped -reference /proj/100genomes/resources/refGenomes/hg19.unmasked.2bit\n";
print "#sw_version 0.3-14-g8e4b69a (Mar 07 2011 13:28)\n";
print "#timestamp Fri Mar 11 16:53:18 2011\n";
print "#contact bobama.systemsbiology.net:\n";
print "#pedigree_name $pedigree_name\n"; 

## if there is a relationship between inviduals they are listed as "Invidual_ID:Father_ID,Mother_ID"
print "#pedigree yanfu yanguifen yanguimei yanguiying yanming yanru yanzhenping\n"; 
print "#sex yanfu:m yanguifen:f yanguimei:f yanguiying:f yanming:m yanru:m yanzhenping:m\n";
print "#pedigree_version 20170501\n";
print "#pedigree_description Beijing xiongke Lung Cancer pedigree Study\n";
print "#reference_genome hg19\n";
print "#numerical_base zero\n";
print "#allele_representation refseq_abstraction\n";
print "#genotypes unphased\n";
print "#\n";
print "#header chromosome,start_position,reference,genotype_pattern,ind1_allele1_base,ind1_allele1_score,ind1_allele2_base,ind1_allele2_score,ind2_allele1_base,ind2_allele1_score,ind2_allele2_base,ind2_allele2_score,ind3_allele1_base,ind3_allele1_score,ind3_allele2_base,ind3_allele2_score,ind4_allele1_base,ind4_allele1_score,ind4_allele2_base,ind4_allele2_score,ind5_allele1_base,ind5_allele1_score,ind5_allele2_base,ind5_allele2_score,ind6_allele1_base,ind6_allele1_score,ind6_allele2_base,ind6_allele2_score,ind7_allele1_base,ind7_allele1_score,ind7_allele2_base,ind7_allele2_score\n";
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	next if(/^#/);
	my ($output,$geno)=("","","");
	my @items=split/\s+/;
	## to remove multiple allellic positoins
	next if($items[4]=~/,/||$items[3]=~/,/);
	my $first=(split(//,$items[3]))[0];
	$output=$items[0].','.$items[1].','.$first.',';
	my $len1=length($items[3]);
	my $len2=length($items[4]);
	if($len1>1 && $len2==1){
		$index=1;
	}elsif($len1==1 && $len2>1){
		$index=2;
	}elsif($len1>1 && $len2>1){
		next
	}else{
		$index=0
	}
	my @allele=();
	foreach my $i(9..$#items){
		my @info=split(/:/,$items[$i]);
		if($info[0] eq '1/1'){
			###  SNP
			if ($index==0){
				push @allele,($items[4],'',$items[4],'');
			### Deletion
			}elsif($index==1){
				push @allele,('-','','-','');
			### Insertion
			}elsif($index==2){
				push @allele,($items[3].'+','',$items[3].'+','');
			}
			$output.='bb';
		}elsif($info[0] eq '0/1'){
			if ($index==0){
				push @allele,($items[3],'',$items[4],'');
			}elsif($index==1){
				push @allele,($first,'','-','');
			}elsif($index==2){
				push @allele,($items[3],'',$items[3].'+','');
			}
			$output.='ab';
		}elsif($info[0] eq '0/0'){
			if ($index==0 ||$index==2){
				push @allele,($items[3],'',$items[3],'');
			}elsif($index==1){
				push @allele,($first,'',$first,'');
			}
			$output.='aa';
		}else{
			push @allele,('N','','N','');
			$output.='nn';
		}
	}
	$geno=join(",",@allele);
	$output.=','.$geno;
	print "$output\n";
}

