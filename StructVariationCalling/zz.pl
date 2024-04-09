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
my $output4step = $opts{o} || die;
my $gaps = $opts{g} || die;
my $centel = $opts{m} || die;
my $name  = $opts{n} || die<<END;

This script is used to filter raw SV calls by score 
and by location, given one or more tab-delimited 
coordinate files with regions to avoid. The files 
must be in the locations specified by makeNewProject.sh

Usage: $0 -c configfile -n name -s software -o ouput

END

#my $rawcallsdir = "/opt/NfsDir/UserDir/sunfl/Projects/biox_qgyw/output/SVBatchCalling";
#my $output4step="$rawcallsdir/$name/$software";
#my %params = ParseConfig::getParams($configfile);

my %params=();
my @config=split(/:/,$configfile);

foreach my $i (@config){
        if($i=~/(.*)\=(.*)/){
                $params{$1}=$2;
        }
}


$params{chrs} =  ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'];


my $parse = new FilterSV();
my $merge = new MergeSV();
# If filtering calls near cen/tel, gaps, other
my %gaps = ();
$params{centel}=$centel;
$params{gaps}=$gaps;
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

# CNVnator
# Merge redundant/overlapping calls of the same SV Type
# Creates a files called $n.merged.bed and $n.merged.tab
foreach my $sv ('DEL','DUP'){
	#my $n = lc $sv;
	# separate SVs to different files
	system "grep $sv $output4step/$name.filtered.tab > $output4step/$name.$sv.filtered.tab";
	$merge->merge_spans("$output4step/$name.$sv.filtered.tab",'I','',"$sv\_CNVnator_$name");
}
# Create BED file
#my @alllines = `cat $rawcallsdir/CNVnator/$name.dup.max.filtered.merged.tab $rawcallsdir/CNVnator/$name.del.max.filtered.merged.tab`;
#$parse->tab2bed('array',\@alllines,"$rawcallsdir/CNVnator/$name.cnvnator.merged.bed",$name,'cnvnator','rawMerged');

