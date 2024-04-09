#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) iscaToGermline.pl
##  Author:
##      Robert M. Hubley <rhubley@systemsbiology.org>
##  Description:
##      Converter for IBD2 blocks from ISCA to Germline format
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2011
#*
#*  This work is licensed under the Open Source License v3.0  To view a copy
#*  of this license, visit http://www.opensource.org/licenses/OSL-3.0 or
#*  see the license.txt file contained in this distribution.
#*
#*  This software is provided ``AS IS'' and any express or implied
#*  warranties, including, but not limited to, the implied warranties of
#*  merchantability and fitness for a particular purpose, are disclaimed.
#*  In no event shall the authors or the Institute for Systems Biology
#*  liable for any direct, indirect, incidental, special, exemplary, or
#*  consequential damages (including, but not limited to, procurement of
#*  substitute goods or services; loss of use, data, or profits; or
#*  business interruption) however caused and on any theory of liability,
#*  whether in contract, strict liability, or tort (including negligence
#*  or otherwise) arising in any way out of the use of this software, even
#*  if advised of the possibility of such damage.
#*
###############################################################################
# ChangeLog:
#
###############################################################################
# TO DO
#

=head1 NAME

 iscaToGermline.pl - Convert ISCA IBD2 blocks to Germline format.

=head1 SYNOPSIS

 Usage:
 ./iscaToGermline.pl --smoothed_blocks_file <file>
                    [ --version ]

=head1 DESCRIPTION

GERMLINE Format:

    Family ID 1
    Individual ID 1
    Family ID 2
    Individual ID 2
    Chromosome
    Segment start (bp)
    Segment end (bp)
    Segment start (SNP)
    Segment end (SNP)
    Total SNPs in segment
    Genetic length of segment
    Units for genetic length (cM or MB)
    Mismatching SNPs in segment
    1 if Individual 1 is homozygous in match; 0 otherwise
    1 if Individual 2 is homozygous in match; 0 otherwise

# Example:
0    HDP3.1.1    0    HDP3.1.2    chr13    19967699    20495159    0    0    990    0.52746    MB    5
0    HDP3.2.1    0    HDP3.2.3    chr13    19610385    20726390    0    0    2788    1.116005    MB    9
0    HDP3.2.1    0    HDP3.2.6    chr13    19610385    20726390    0    0    2788    1.116005    MB    9
0    HDP3.2.4    0    HDP3.2.6    chr13    19533777    20726390    0    0    3007    1.192613    MB    9
0    HDP3.2.3    0    HDP3.2.4    chr13    19533777    20726920    0    0    3009    1.193143    MB    9
0    HDP3.2.1    0    HDP3.2.6    chr13    25720684    26564650    0    0    1585    0.843966    MB    7
0    HDP3.2.3    0    HDP3.2.6    chr13    22207406    26785830    0    0    11818    4.578424    MB    24
0    HDP3.2.6    0    HDP3.2.7    chr13    20123804    26795462    0    0    16349    6.671658    MB    33
0    HDP3.2.5    0    HDP3.2.6    chr13    20123804    26802934    0    0    16360    6.67913    MB    33
0    HDP3.2.2    0    HDP3.2.6    chr13    19410519    26808596    0    0    18751    7.398077    MB    36


=head1 OPTIONS

    -version  
        print out the version and exit

    --smoothed_blocks_file
        A file in ISCA smoothed_blocks format.

=head1 SEE ALSO

=over 4

=back

=head1 COPYRIGHT

Copyright 2011 Institute for Systems Biology

=head1 AUTHORS

  Robert M. Hubley <rhubley@systemsbiology.org> 

=cut

##
## Module dependence
##
use strict;
use FindBin;
use lib $FindBin::RealBin;
use Getopt::Long;
use Data::Dumper;

# Get version from Makefile/GIT
my $releaseVersion = "0.1.8";
my $gitVersion     = "0.1.8";
my $VERSION        = "release $releaseVersion";
$VERSION .= " ( git: $gitVersion )"
    if ( $gitVersion ne "" );

#
# DEBUG mode
#
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-version',                 # print out the version and exit
                    '-smoothed_blocks_file=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $VERSION\n";
  exec "pod2text $0";
  exit( 1 );
}

if ( $options{'version'} ) {
  print "$VERSION\n";
  exit;
}

usage()
    if ( !exists $options{'smoothed_blocks_file'} );

my $smoothedBlocksFile = $options{'smoothed_blocks_file'};

my $pedigree = "Unknown";
my $child1   = "child1";
my $child2   = "child2";

if ( -e $smoothedBlocksFile ) {
  open IN, "<$smoothedBlocksFile"
      or die "Could not open $smoothedBlocksFile: $!\n";

  while ( <IN> ) {

    if ( /^#sex\s+(.*)$/ ) {
      my @children = split( /\s+/, $1 );
      shift @children;    # father
      shift @children;    # mother
      $child1 = shift @children;
      $child1 =~ s/\:[mf]//g;
      $child2 = shift @children;
      $child2 =~ s/\:[mf]//g;
    }

    next if ( /^#/ );
    my @fields = split( /\s+/ );
    if ( $fields[ 2 ] == 0 ) {
      print "$pedigree\t$child1\t$pedigree\t$child2\t$fields[0]\t$fields[6]\t"
          . "$fields[7]\t0\t0\tNA\t"
          . ( $fields[ 5 ] / 1000000 )
          . "\tMB\t\n";
    }
  }
  close IN;
}
else {
  die "Could not find file: $smoothedBlocksFile\n";
}

1;
