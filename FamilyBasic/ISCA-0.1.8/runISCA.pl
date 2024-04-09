#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) runISCA.pl
##  Author:
##      Jared Roach      <jroach@systemsbiology.org>
##      Robert M. Hubley <rhubley@systemsbiology.org>
##  Description:
##     Run a complete Inheritance State Consistency Analysis (ISCA).
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

 runISCA.pl - Run a complete Inheritance State Consistency Analysis (ISCA)

=head1 SYNOPSIS

 Usage:
 ./runISCA.pl [ --version ]
              [ --use_precomputed_values ]
              [ --exclude_partially_called_vectors ]
              [ --exclude_fully_heterozygous_vectors ]
              [ --contact <contact_info> ]
              [ --project <proj_name> ]
              [ --pedigree_name <pedigree_name> ]
              [ --output_dir <dir> ]
              --reference_genome <hg18|hg19>
              --number_of_nonfounders #
              --inclusion_order <order> 
              --file_to_analyze <genotype_vector_file>

=head1 DESCRIPTION

=head1 OPTIONS

    --version  
        print out the version and exit

    --exclude_partially_called_vectors
	Exclude loci which are not fully called ( both haplotypes )
        in all individuals.

    --exclude_fully_heterozygous_vectors
	Exlclude loci 

    --use_precomputed_values
	Read in data structures precomputed in a previous run of the
        HMM.  This is primarily designed for use by developers.


=head1 SEE ALSO

=over 4

=back

=head1 COPYRIGHT

Copyright 2011 Institute for Systems Biology

=head1 AUTHORS

  Jared Roach      <jroach@systemsbiology.org> 
  Robert M. Hubley <rhubley@systemsbiology.org> 

=cut

##
## Module dependence
##
use strict;
use FindBin;
use lib $FindBin::RealBin;
use File::Spec;
use TimeUtils;
use Cwd;
use Getopt::Long;
use Data::Dumper;

# Get version from Makefile/GIT
my $releaseVersion = "0.1.8";
my $gitVersion     = "0.1.8";
my $VERSION        = "release $releaseVersion";
$VERSION .= " ( git: $gitVersion )"
    if ( $gitVersion ne "" );

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    "version",
                    "debug",
                    "number_of_nonfounders=i",
                    "reference_genome=s",
                    "use_precomputed_values",
                    "exclude_partially_called_vectors",
                    "exclude_fully_heterozygous_vectors",
                    "pedigree_name=s",
                    "contact=s",
                    "project=s",
                    "file_to_analyze=s",
                    "inclusion_order=s",
                    "output_dir=s"
);

my %options = ();
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $VERSION\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$VERSION\n";
  exit;
}

# Validate options
if (    !defined $options{'file_to_analyze'}
     || !-s $options{'file_to_analyze'} )
{
  print "\n\nCould not find any data files to process!\n\n";
  usage();
}
elsif ( !defined $options{'reference_genome'}
        || $options{'reference_genome'} !~ /hg18|hg19/i )
{
  print "\n\nMissing or invalid reference genome. Currently\n"
      . "supported values are \"hg18\" or \"hg19\".\n\n";
  usage();
}
elsif ( !defined $options{'number_of_nonfounders'} ) {
  print "\n\nMissing number_of_nonfounders parameter!\n\n";
  usage();
}
elsif ( !defined $options{'inclusion_order'} ) {
  print "\n\nMissing inclusion_order parameter!\n\n";
  usage();
}

# If --output_dir is not used assume the current working directory
# is the output directory.
my $outputDir = $options{'output_dir'};
if ( !defined $options{'output_dir'} ) {
  $outputDir = cwd();
}

my $contact = "Unknown";
$contact = $options{'contact'}
    if ( $options{'contact'} );

my $project = "Unknown";
$project = $options{'project'}
    if ( $options{'project'} );

my $pedigree_prefix = "Unknown";
$pedigree_prefix = $options{'pedigree_name'}
    if ( $options{'pedigree_name'} );

my $cmd;

# Start a timer for the whole workflow
TimeUtils::elapsedTime( 0 );

print "\n\nISCA Analysis - $VERSION\n\n";

##
## Run HMM
##
# Start a timer for this subtask
TimeUtils::elapsedTime( 1 );
my $msg = ""
    . timestamp( TimeUtils::FMT_YYYYMMDD_TIME )
    . ": Running Suttonian HMM\n";
print $msg;
system( "echo '$msg' > $outputDir/isca.log" );

$cmd =
      "$FindBin::RealBin/suttonian_hmm.pl "
    . "--number_of_nonfounders=$options{'number_of_nonfounders'} "
    . "--file_to_analyze=$options{'file_to_analyze'} "
    . "--inclusion_order=$options{'inclusion_order'} "
    . "--output_dir=$outputDir "
    . "--contact=$contact "
    . "--project=$project "
    . "--pedigree_name=$pedigree_prefix ";

$cmd .= "--exclude_fully_heterozygous_vectors "
    if ( defined $options{'exclude_fully_heterozygous_vectors'} );

$cmd .= "--exclude_partially_called_vectors "
    if ( defined $options{'exclude_partially_called_vectors'} );

$cmd .= "--use_precomputed_values "
    if ( defined $options{'use_precomputed_values'} );

$cmd .= " >> $outputDir/isca.log 2>&1";
print "Running: $cmd\n" if ( $options{'debug'} );
system( $cmd );
if ( $? ) { die "command $cmd failed: $!\n"; }
$msg = "  - Completed, Elapsed Time: " . TimeUtils::elapsedTime( 1 ) . "\n";
print $msg;
system( "echo '$msg' >> $outputDir/isca.log" );

##
## Post-process blocks
##
$msg = ""
    . timestamp( TimeUtils::FMT_YYYYMMDD_TIME )
    . ": Post-processing blocks\n";
print $msg;
system( "echo '$msg' >> $outputDir/isca.log" );
$cmd =
      "$FindBin::RealBin/suttonian_block_nibbler.pl "
    . "--infile_directory=$outputDir ";
$cmd .= " >> $outputDir/isca.log 2>&1";
print "Running: $cmd\n" if ( $options{'debug'} );
system( $cmd );
if ( $? ) { die "command $cmd failed: $!\n"; }
$msg = "  - Completed, Elapsed Time: " . TimeUtils::elapsedTime( 1 ) . "\n";
print $msg;
system( "echo '$msg' >> $outputDir/isca.log" );

##
## Calculate block statistics
##
$msg = ""
    . timestamp( TimeUtils::FMT_YYYYMMDD_TIME )
    . ": Calculating block statistics\n";
print $msg;
system( "echo '$msg' >> $outputDir/isca.log" );
$cmd =
      "$FindBin::RealBin/statistics_for_HMM_Suttonian_blocks.pl "
    . "--reference_genome=$options{'reference_genome'} "
    . "--infile_directory=$outputDir ";
$cmd .= " >> $outputDir/isca.log 2>&1";
print "Running: $cmd\n" if ( $options{'debug'} );
system( $cmd );
if ( $? ) { die "command $cmd failed: $!\n"; }
$msg = "  - Completed, Elapsed Time: " . TimeUtils::elapsedTime( 1 ) . "\n";
print $msg;
system( "echo '$msg' >> $outputDir/isca.log" );

#
# At this time these post-processing scripts don't work on
# pedigrees larger than quads
#
if ( $options{'number_of_nonfounders'} <= 2 ) {

  ##
  ## Smooth blocks
  ##
  $msg =
      "" . timestamp( TimeUtils::FMT_YYYYMMDD_TIME ) . ": Smoothing blocks\n";
  print $msg;
  system( "echo '$msg' >> $outputDir/isca.log" );
  $cmd =
        "$FindBin::RealBin/parse_and_smooth_HMM_Suttonian_blocks.pl "
      . "--infile_directory=$outputDir ";
  $cmd .= " >> $outputDir/isca.log 2>&1";
  print "Running: $cmd\n" if ( $options{'debug'} );
  system( $cmd );
  if ( $? ) { die "command $cmd failed: $!\n"; }
  $msg = "  - Completed, Elapsed Time: " . TimeUtils::elapsedTime( 1 ) . "\n";
  print $msg;
  system( "echo '$msg' >> $outputDir/isca.log" );

  ##
  ## Add partially called blocks
  ##
  $msg = ""
      . timestamp( TimeUtils::FMT_YYYYMMDD_TIME )
      . ": Intercalating partially called blocks\n";
  print $msg;
  system( "echo '$msg' >> $outputDir/isca.log" );
  $cmd =
        "$FindBin::RealBin/intercalate_partial_binary_blocks.pl "
      . "--infile_directory=$outputDir ";
  $cmd .= " >> $outputDir/isca.log 2>&1";
  print "Running: $cmd\n" if ( $options{'debug'} );
  system( $cmd );
  if ( $? ) { die "command $cmd failed: $!\n"; }
  $msg = "  - Completed, Elapsed Time: " . TimeUtils::elapsedTime( 1 ) . "\n";
  print $msg;
  system( "echo '$msg' >> $outputDir/isca.log" );
}

$msg =
      "\nRun Completed, Elapsed Time: "
    . TimeUtils::elapsedTime( 0 )
    . "\nLog of run saved to isca.log\n";
print $msg;
system( "echo '$msg' >> $outputDir/isca.log" );

1;
