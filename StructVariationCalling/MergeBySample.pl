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

getopts('n:o:', \%opts);
my $outpath4step = $opts{o} || die;
my $name  = $opts{n} || die<<END;

This script is used to filter raw SV calls by score 
and by location, given one or more tab-delimited 
coordinate files with regions to avoid. The files 
must be in the locations specified by makeNewProject.sh

Usage: $0 -c configfile

END

#my %params = ParseConfig::getParams($configfile);
#ParseConfig::checkParams('filter',\%params);

my $merge = new MergeSV();

$merge->merge_spans("$outpath4step/${name}.DEL.tab",'I','',"${name}.DEL.merged.tab") ;
$merge->merge_spans("$outpath4step/${name}.DUP.tab",'I','',"${name}.DUP.merged.tab") ;
$merge->merge_spans("$outpath4step/${name}.INV.tab",'I','',"${name}.INV.merged.tab") ;
