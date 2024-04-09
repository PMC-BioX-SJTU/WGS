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

getopts('c:n:g:m:o:', \%opts);

my $configfile = $opts{c} || die;
my $output4step = $opts{o} || die;
my $gaps = $opts{g} || die;
my $centel = $opts{m} || die;
my $name  = $opts{n} || die<<END;

This script is used to filter raw SV calls by score 
and by location, given one or more tab-delimited 
coordinate files with regions to avoid. The files 
must be in the locations specified by makeNewProject.sh

Usage: $0 -c configfile -n name  -o output4rpoject

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
$params{centel}=$centel;
$params{gaps}=$gaps;
$params{chrs} =  ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'];


my $parse = new FilterSV();
my $merge = new MergeSV();
# If filtering calls near cen/tel, gaps, other
my %gaps = ();

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

# Breakdancer
open OUT, ">$output4step/$name.DEL.INV.INS.filtered.tab" or die;
open C, "<$output4step/$name.DEL.INV.INS.output" or die;
while (<C>) {
	my $line = $parse->cleanBD($_);
	if ($line && ($params{BDscore} || $params{BDrs})) {
		$line = $parse->filterBDscore($line,$params{BDscore},$params{BDrs});
	}
	# Reformat to tab-delimeted and BED formats
	if ($line) {
		$line = $parse->bd2tab($line,$name);
		# filter out if near cen/tel or assembly gaps
		if ($line && %gaps) {
			$line = $parse->filterGaps($line,\%gaps,$params{gapOverlap});
		}
		print OUT $line."\n" if $line;
	}
}
close C;
close OUT;

# Merge redundant/overlapping calls of the same SV Type
# Creates a files called $n.merged.bed and $n.merged.tab
foreach my $sv ('INS','DEL','INV'){
	#my $n = lc $sv;
	# separate SVs to different files
	system "grep $sv $output4step/$name.DEL.INV.INS.filtered.tab > $output4step/$name.$sv.filtered.tab";
	$merge->merge_spans("$output4step/$name.$sv.filtered.tab",'I','',"$sv\_BD_$name");
}
# Create BED file
#my @alllines = `cat $output4step/$name.ins.max.filtered.merged.tab $output4step/$name.del.max.filtered.merged.tab $output4step/$name.inv.max.filtered.merged.tab`;
#$parse->tab2bed('array',\@alllines,"$output4step/$name.BD.merged.bed",$name,'BD','rawMerged');

