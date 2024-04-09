#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use lib dirname(__FILE__).'/./';
use ParseConfig;
use FilterSV;
use MergeSV;

use vars "%opts";

getopts('c:o:n:g:m:', \%opts);

my $configfile = $opts{c} || die;
my $rawcallsdir = $opts{o} || die;
my $software = $opts{s} || die;
my $gaps = $opts{g} || die;
my $centel = $opts{m} || die;
my $name  = $opts{n} || die<<END;


This script is used to filter raw SV calls by score 
and by location, given one or more tab-delimited 
coordinate files with regions to avoid. The files 
must be in the locations specified by makeNewProject.sh

Usage: $0 -c configfile -n name -s software -o output4rpoject

END

#my %params = ParseConfig::getParams($configfile);
my %params=();
my @config=split(/\s+/,$configfile);

foreach my $i (@config){
        if($i=~/(.*)\=(.*)/){
                $params{$1}=$2;
        }
}


#-- Parse SV caller output --

$params{chrs} =  ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'];

my $output4step=$rawcallsdir;

my $parse = new FilterSV();
my $merge = new MergeSV();
# If filtering calls near cen/tel, gaps, other
my %gaps = ();
$params{gaps} = $opts{g} ;
$params{centel} = $opts{m};

if ($params{centel} || $params{gaps} || $params{filterOther}) {
	my @list = ();
	my $b1 = $params{centelBuffer} || '0',
	my $b2 = $params{gapsBuffer} || '0',
	my $b3 = $params{filterOtherBuffer} || '0';
	push @list, [$params{centel},$b1] if $params{centel};
	push @list, [$params{gaps},$b2] if $params{gaps};
	push @list, [$params{filterOther},$b3] if $params{filterOther};
	%gaps = $parse->parseGaps(@list);

}

my %types = ( # latest pindel output files: del,large ins,inv,tand.dup.,small ins.
	'D' => 'del',
#	'LI' => 'ins',
	'INV' => 'inv',
	'TD' => 'dup',
	'SI' => 'ins',
);
# Filtering and reformatting
foreach my $type (sort keys %types) { # latest pindel output files: del,large ins,inv,tand.dup.,small ins.
	my $type2 = uc($types{$type});
	open OUT, ">$output4step/$name.$type2.filtered.tab" || die $!;
	my @lines = $parse->parsePindelVCF("$output4step/$name.merged.$type.vcf"); 
	if (@lines) {
		foreach my $line (@lines) {
			$line = $parse->pindel2tab($line,$name,$type2);
			# filter out if near cen/tel or assembly gaps
			if ($line && %gaps) {
				$line = $parse->filterGaps($line,\%gaps,$params{gapOverlap});
			}
			print OUT $line."\n" if $line;
		}
	}
	close OUT;
}
	
# Separate small and large SVs; del.small.tab, del.large.tab
#foreach my $sv (sort keys %types) {
#	next if $sv eq 'SI' || $sv eq 'LI';
#	my $type2 = uc($types{$sv});
#	$merge->splitBySize("$output4step/$name.$type2.filtered.tab",50);
#}

#`cp $output4step/$name.SMALLINS.filtered.tab  $output4step/$name.SMALLINS.filtered.small.tab`;

# Merge redundant/overlapping calls of the same SV Type
# Creates a file called $n.filter.merged.tab
foreach my $sv (sort keys %types) {
	my $annot = uc($types{$sv});
	next if $annot eq 'SMALLINS';
	$merge->merge_spans("$output4step/$name.$annot.filtered.tab",'I','',"$name.$annot\_Pindel_$name");
	#my @alllines = `cat $output4step/$params{name}.$annot.filtered.large.merged.tab`;
	#$parse->tab2bed('array',\@alllines,"$output4step/$name.${annot}.merged.bed",$name,'PD','rawMerged');
}
