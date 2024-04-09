#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) suttonian_hmm.pl
##  Author:
##      Jared Roach      <jroach@systemsbiology.org>
##      Robert M. Hubley <rhubley@systemsbiology.org>
##  Description:
##     hidden markov model (HMM)  for finding inheritance states
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
# created 8/31/09
# modified 8/31/09-8/30/2010 v0.02-v0.25  See archives for changes (v 0.25
#   was last version that was hard-coded for a nuclear family of four).
#   No versions v0.26-v0.29 as major change to version 0.30
# modified 10/20/2010 v0.30-v0.36  overhaul for pedigrees with an
#   arbitrary number of children (hard coded for a nuclear family, but
#   allows the number of chlidren to vary)
# modified 10/20/2010 v0.37  added Storable to precompute up-front calculations
# modified 10/22/2010 v0.38  minor changes
# modified 10/22/2010 v0.39  preload genotype vector file
# modified 10/22/2010 v0.40  minor
# modified 10/27/2010 v0.41  fixed bug in algorithm for computing Hamming
#   distance between states (state_to_state_Hamming_distances)
# modified 10/28/2010 v0.42  improved state-to-state transition probabilities
#   (load_transition_matrix)
# modified 10/28/2010 v0.43  ignore MIEs, deletions in the input file, and
#   adjacent emissions that are identical
# modified 10/29/2010 v0.44  increase pseudocount for Mendelian frequencies,
#   downweight ababababab vectors, make transition probabilities a little harder
# modified 10/29/2010 v0.44  minor: keep a record of total maternal and
#   paternal recombinations; allow taking subsets of individuals in input file
# modified 11/10/2010 v0.45  change code for dealing with linkage - keep
#   a running tab of emission counts rather than just screen adjacent
#   identical counts
# modified 11/12/2010 v0.47  removed old code reltaed to T1DGC
# modified 11/15/2010 v0.48  fixed a huge bug in v47 in the
#   "state_consistencies" subroutine -problem with not using altered
#   genotype vector for each parental pattern
# modified 11/15/2010 v0.49  minor clean up of code
# modified 11/17/2010 v0.50  minor -> improved output headers
# modified 11/17/2010 v0.52  minor
# modified 11/20/2010 v0.53  no longer exclude fully heterozygous vectors
# modified 11/21/2010 v0.54  exclude fully heterozygous vectors, but allow
#   as an option
# modified 11/21/2010 v0.55  record command line in output header
# modified 12/17/2010 v0.56  minor changes
#
###############################################################################
# TO DO
#   fix the code for correctly dealing with pseudoautosomes, X, and Y
#

=head1 NAME

 suttonian_hmm.pl - hidden markov model for finding inheritance states

=head1 SYNOPSIS

 Usage:
 ./suttonian_hmm.pl --file_to_analyze <genotype_vector_file>
                    --number_of_nonfounders # 
                    --inclusion_order <order>
                    [--reference_genome <hg18|hg19> ]
                    [ --use_precomputed_values ]
                    [ --exclude_partially_called_vectors ]
                    [ --exclude_mie_vectors ]
                    [ --exclude_fully_heterozygous_vectors ]
                    [ --contact <contact_info> ]
                    [ --project <proj_name> ]
                    [ --pedigree_name <pedigree_name> ]
                    [ --output_dir <dir> ]
                    [ --version ]

=head1 DESCRIPTION

=head1 OPTIONS

    -version  
        print out the version and exit

    --file_to_analyze <genotype_vector_file>
        A file in ISB's genotype vector format

    --number_of_nonfounders #
        The number of children in a nuclear family analysis.
        At this time this algorithm is tuned for quartets and
        this is only guaranteed to work for 2 nonfounders.

    --inclusion_order father,mother,child1,child2
        The inclusion order specifies which individuals from
        a (possibly) larger set of individuals in a genotype
        file to use.  The genotype file individuals are indexed
        from left to right as 0..n.  Likewise the inclusion_order
        parameter should be comma separated list of the genotype
        indexes in family order.  The founders may be omitted
        by specifying an "x" rather than the individual index.

        ie. Given a genotype vector file with a father, grandfather, 
        child1, child2.  The inclusion order for the mother and two
        children will be:

             --inclusion_order 0,x,2,3 

    --exclude_partially_called_vectors
	Ignore vectors for which any of the individuals is not
        fully called.  NOTE: This exclusion only applies to
        individuals included in the analysis.  Missing founders
        and individuals not specified in the inclusion order
        parameter are not considered when eliminating partially
        called vectors.

    --exclude_fully_heterozygous_vectors
        Ignore vectors for which all individuals are fully 
        heterozygous.  I.e "abababab".   Fully heterozygous vectors 
        are enriched in CNVs and compressions

    --exclude_mie_vectors
        Use only for really noisy data.

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
use Cwd;
use Getopt::Long;
use Data::Dumper;
use warnings FATAL => 'all';
use TimeUtils qw(timestamp);
use Storable qw(nstore retrieve);
use GenotypeEmissions;

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

# Save the unaldulterated command line
my $command_line_command                        = $0;
my $command_line_arguments                      = join( ' ', @ARGV );
my $invoking_command_to_record_on_output_header =
    $command_line_command . " " . $command_line_arguments;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my $number_of_nonfounders              = 0;
my $use_precomputed_values             = 0;
my $exclude_partially_called_vectors   = 0;
my $exclude_fully_heterozygous_vectors = 0;
my $pedigree_name                      = "Unknown";
my $project                            = "Unknown";
my $file_to_analyze                    = "";
my $inclusion_order                    = "";
my $print_version                      = 0;
my $ignore_MIEs                        = 0;
my $outputDir                          = "";
my $contact                            = "Unknown";
my $reference_genome                   = "Unknown";

my @getopt_args = (
   "version"                            => \$print_version,
   "number_of_nonfounders=i"            => \$number_of_nonfounders,
   "use_precomputed_values"             => \$use_precomputed_values,
   "exclude_mie_vectors"                => \$ignore_MIEs,
   "exclude_partially_called_vectors"   => \$exclude_partially_called_vectors,
   "exclude_fully_heterozygous_vectors" => \$exclude_fully_heterozygous_vectors,
   "reference_genome=s"                 => \$reference_genome,
   "pedigree_name=s"                    => \$pedigree_name,
   "project=s"                          => \$project,
   "contact=s"                          => \$contact,
   "file_to_analyze=s"                  => \$file_to_analyze,
   "inclusion_order=s"                  => \$inclusion_order,
   "output_dir=s"                       => \$outputDir
);

unless ( GetOptions( @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $VERSION\n\n";
  exec "pod2text $0";
  exit;
}

if ( $print_version ) {
  print "$VERSION\n";
  exit;
}

if ( !-s $file_to_analyze ) {
  print "\n\nCould not find any data files to process!\n\n";
  usage();
}

# If --output_dir is not used assume the current working directory
# is the output directory.
if ( $outputDir eq "" ) {
  $outputDir = cwd();
}

# Breakup the data file into directory and filename ( portable )
my ( $volOfGenotypesFile, $pathOfGenotypesFile, $genotypesFilename ) =
    File::Spec->splitpath( $file_to_analyze );
$pathOfGenotypesFile = cwd() if ( $pathOfGenotypesFile eq "" );

# Define output filename prefixes
my $dataCacheFilename = "hmm_DataCache.dat";
if ( $use_precomputed_values ) {
  if ( -s "$pathOfGenotypesFile/$dataCacheFilename" ) {
    $dataCacheFilename = "$pathOfGenotypesFile/$dataCacheFilename";
  }
  elsif ( -s "$outputDir/$dataCacheFilename" ) {
    $dataCacheFilename = "$outputDir/$dataCacheFilename";
  }
  else {
    if ( $pathOfGenotypesFile eq $outputDir ) {
      die "\n\nCould not locate data cache file $dataCacheFilename in\n"
          . "the $pathOfGenotypesFile directory.\n\n";
    }
    else {
      die "\n\nCould not locate data cache file $dataCacheFilename in\n"
          . "either the $pathOfGenotypesFile\n"
          . "or the $outputDir directories\n\n";
    }
  }
}
else {
  $dataCacheFilename = "$outputDir/$dataCacheFilename";
}

my $lastDir;
my ( $vol, $path, $filename ) = File::Spec->splitpath( $file_to_analyze );
if ( $path ) {
  my @dirs = File::Spec->splitdir( $path );
  $lastDir = pop @dirs;
}

# TODO: Grab this from the genotype file
unless ( defined( $project ) ) {
  if ( $lastDir ne "" ) {
    $project = $lastDir;
  }
  else {
    $project = "unknown";
  }
}

# TODO: Grab this from the genotype file
unless ( defined( $pedigree_name ) ) {
  if ( $lastDir ne "" ) {
    $pedigree_name = $lastDir;
  }
  else {
    $pedigree_name = "unknown";
  }
}

unless ( defined( $number_of_nonfounders ) ) {
  $number_of_nonfounders =
      2;    #synonymous with number of children for a nuclear family
}

my $skip_deletions      = 1;
my $skip_insertions     = 1;
my $linkage_suppression = 1;

if (    $inclusion_order =~ /^[\d+],x,.*/
     || $inclusion_order =~ /^x,.*/ )
{
  print STDERR "WARNING: Turning off exclude_fully_heterozygous_vectors."
      . "  This is required when one or both parents are missing!\n";
  $exclude_fully_heterozygous_vectors = 0;
}
if ( $inclusion_order =~ /^x,x,.*/ ) {
  print STDERR "WARNING: Turning off linkage supression.  "
      . "This is required when both parents are missing!\n";
  $linkage_suppression = 0;
}

#################### SETUP HMM PARAMETERS ########################

my $E           = 2.71828183;
my $pseudocount = 0.00000000001;
my $mode        = "Sutton";

# These counters will sum over all chromosomes
my $total_paternal_recombinations = 0;
my $total_maternal_recombinations = 0;

# Fully called positions (2,182,329,631 total positions)
# Hg18 has 2,855,343,769 positions that are not Ns; in this paper,
# when we refer to hg18 we mean this set of "unambiguous" positions.
# The 2,855,343,769 positions include Y chromosome positions.
# Excluding Y chromosome positions would give 2,832,359,852 positions.
my $genome_size = 2400000000;

#
# Estimated background MIE rate
#
#   This value will be used to set the emission probabilities
#   for all states on all MIE patterns.  See GenotypeEmission.pm
#   TODO: This should probably move to GenotypeEmission.pm
#   altogether.
#

# I saw 506502 MIEs across 4738970 patterns including those
# with Ns. (note that there are no patterns with n's in the
# UCSD version 2 dataset)
my $number_of_observed_MIEs                      = 513585;
my $number_of_markers                            = 4651898;
my $estimated_fraction_of_genome_in_error_blocks = 0.01;

# Smaller values here ( ie. 5 ) give too many short error states in some
# datasets.
my $definition_of_ratio_of_error_block_rate_to_background_rate = 10;

my $length_of_genome_estimated_to_be_background =
    ( 1 - $estimated_fraction_of_genome_in_error_blocks ) * $number_of_markers;

my $length_of_genome_estimated_to_be_error_prone =
    $estimated_fraction_of_genome_in_error_blocks * $number_of_markers;

my $estimated_background_MIE_rate =
    $number_of_observed_MIEs / ( $length_of_genome_estimated_to_be_background +
                   $definition_of_ratio_of_error_block_rate_to_background_rate *
                   $length_of_genome_estimated_to_be_error_prone );

# This is essentially a pseudocount; however it should be tuned for the
# actual observed rate of PCEs.
#   -- I think the rate of 0.005 is driving the HMM results to flutter
#      into very short states in specific datasets.
# This value should either be determined iteratively and empirically or
# be set as a function of the MIE observed rate.
my $pce_emission_prob = 0.03;

if ( $DEBUG ) {
  print STDERR
      "Estimated background MIE rate = $estimated_background_MIE_rate\n";
  print STDERR "PCE emission probability = $pce_emission_prob\n";
}

#
# Setup possible genotypes
#
# ie "aa" = 0, "ab" = 1 .. "nn" = 6
my @possible_genotypes = ( "aa", "ab", "an", "bb", "bn", "nn" );
my $integer = 0;
my %genotype_for_integer;
foreach my $genotype ( @possible_genotypes ) {
  $genotype_for_integer{$integer} = $genotype;
  $integer++;
}

my $number_of_possible_genotypes = scalar @possible_genotypes;

# Number of meioses in pedigree is also the length of the inheritance
# vector.
my $number_of_meioses_in_pedigree = 2 * $number_of_nonfounders;

my $number_of_states         = &getNumberOfStates( $number_of_nonfounders );
my $compsci_number_of_states =
    $number_of_states - 1;    #compsci based number of states
my $number_of_individuals_in_pedigree =
    2 + $number_of_nonfounders
    ; #this is the length of the genotype vector (still assumes 2 parents and their children)
my $compsci_number_of_phasings_of_nonfounder_genotypes =
    2**$number_of_nonfounders - 1;

my $fully_heterozygous_genotype_vector =
    "ab" x $number_of_individuals_in_pedigree;

my $completely_permissive_state_consistency_vector = "";
my $MIE_state_consistency_vector                   = "";
for ( my $n = 0 ; $n < $number_of_states ; $n++ ) {
  $completely_permissive_state_consistency_vector .= "1";
  $MIE_state_consistency_vector                   .= "0";
}

# Initialize a hash with the state indices
#the state indices are the first n whole numbers (and their binary synoynms)
my @state_indices = ( 0 .. $compsci_number_of_states );

my @list_of_founder_alleles = ( 0, 1, 2, 3 );    #hard coded for nuclear family
my @parental_pattern_indices =
    ( 1, 2, 3, 4 )
    ; #these cycle through the four different inheritance vectors that together form a state by flipping the parental bits
my @list_of_possible_phasings_of_nonfounders_in_genotype_vector =
    ( 0 .. $compsci_number_of_phasings_of_nonfounder_genotypes );
my $maximum_genotype_vector_index =
    ( $number_of_possible_genotypes**$number_of_individuals_in_pedigree ) - 1;

#
# More globals
#
my %genotype_vector_consistent_with_state;
my %genotype_content_index;
my %input_file_header_info;
my $lowestIntegerForEquivGenotypeVectorsHashRef;
my %Mendelian_prob_of_genotype_vector_emission_given_state_and_founder_genotypes;
my $emissionProbHashRef;
my $possibleEmissionsArrayRef;
my @state_names;
my %states_are_adjacent;
my %parental_transition;
my %observed_emission;
my $MIEHashRef;
my $total_genome_number_of_informative_Suttonian_locations = 0;

##
## Begin in earnest
##
print STDERR "\n\nVersion $VERSION of $0, running in mode $mode "
    . "with datasource $file_to_analyze and timestamp "
    . timestamp() . "\n";
print STDERR "Individuals\t$number_of_individuals_in_pedigree\n"
    . "Top genotype vector index\t$maximum_genotype_vector_index\n";

if ( $use_precomputed_values ) {
  my $tmpDataCache = load_previously_compiled_data( $dataCacheFilename );

  # NOTE: This is a very inefficient way to break the hash into
  #       the discrete variables. This will double the memory
  #       requirements doing it this way ( at least temporarily )
  #       TODO: Fix.
  %genotype_vector_consistent_with_state =
      %{ $tmpDataCache->{'genotype_vector_consistent_with_state'} };
  %genotype_for_integer   = %{ $tmpDataCache->{'genotype_for_integer'} };
  %genotype_content_index = %{ $tmpDataCache->{'genotype_content_index'} };
  %input_file_header_info = %{ $tmpDataCache->{'input_file_header_info'} };
  $lowestIntegerForEquivGenotypeVectorsHashRef =
      $tmpDataCache->{'lowest_integer_for_equivalent_genotype_vector'};
  %Mendelian_prob_of_genotype_vector_emission_given_state_and_founder_genotypes
      = %{
    $tmpDataCache->{
'Mendelian_prob_of_genotype_vector_emission_given_state_and_founder_genotypes'
        }
      };
  $emissionProbHashRef       = $tmpDataCache->{'emission_prob'};
  $possibleEmissionsArrayRef = $tmpDataCache->{'possible_emissions'};
  @state_names               = @{ $tmpDataCache->{'state_names'} };
  %states_are_adjacent       = %{ $tmpDataCache->{'states_are_adjacent'} };
  %parental_transition       = %{ $tmpDataCache->{'parental_transition'} };
  %observed_emission         = %{ $tmpDataCache->{'observed_emission'} };
  $MIEHashRef                = $tmpDataCache->{'MIE'};

  # Free up memory -- no need to keep a duplicate
  #undef $tmpDataCache;
}
else {

  # Create data structures
  (
    $possibleEmissionsArrayRef, $MIEHashRef,
    $lowestIntegerForEquivGenotypeVectorsHashRef,
    $emissionProbHashRef
      )
      = GenotypeEmissions::initializeTwoAlleleEmissionAlphabet( 2,
                        $pce_emission_prob, $estimated_background_MIE_rate, 0 );

  my $stateNamesRef = name_the_states( $number_of_states );

  # Recast for now
  @state_names = @$stateNamesRef;

  my ( $parentalTransitionRef, $statesAreAdjacentRef ) =
      state_to_state_Hamming_distances( $number_of_states,
                                        $number_of_nonfounders );

  # Recast for now
  %parental_transition = %$parentalTransitionRef;
  %states_are_adjacent = %$statesAreAdjacentRef;

  my ( $inputFileHdrInfoRef, $observedEmissionsRef, $informativePositions ) =
      load_genotype_vectors(
                             $file_to_analyze,
                             $skip_insertions,
                             $skip_deletions,
                             $exclude_fully_heterozygous_vectors,
                             $exclude_partially_called_vectors,
                             $maximum_genotype_vector_index,
                             $ignore_MIEs,
                             $MIEHashRef,
                             $lowestIntegerForEquivGenotypeVectorsHashRef,
                             $inclusion_order
      );

  # Recasts ( for now )
  %input_file_header_info = %$inputFileHdrInfoRef;
  %observed_emission      = %$observedEmissionsRef;
  $total_genome_number_of_informative_Suttonian_locations =
      $informativePositions;

  # Make a combined datastructure to limit the number
  # of data cache files.
  my %dataCache = (
    'MIE'                                   => $MIEHashRef,
    'genotype_vector_consistent_with_state' =>
        \%genotype_vector_consistent_with_state,
    "genotype_for_integer"                          => \%genotype_for_integer,
    "genotype_content_index"                        => \%genotype_content_index,
    "lowest_integer_for_equivalent_genotype_vector" =>
        $lowestIntegerForEquivGenotypeVectorsHashRef,
"Mendelian_prob_of_genotype_vector_emission_given_state_and_founder_genotypes"
        => \
        %Mendelian_prob_of_genotype_vector_emission_given_state_and_founder_genotypes,
    "emission_prob"          => $emissionProbHashRef,
    "possible_emissions"     => $possibleEmissionsArrayRef,
    "state_names"            => \@state_names,
    "states_are_adjacent"    => \%states_are_adjacent,
    "parental_transition"    => \%parental_transition,
    "observed_emission"      => \%observed_emission,
    "input_file_header_info" => \%input_file_header_info,
  );

  # write out data cache
  nstore \%dataCache, $dataCacheFilename;
  my $local_timestamp = timestamp();
  print STDERR "$local_timestamp: Done saving pre-computation "
      . "values to $dataCacheFilename\n";
}

my $number_of_possible_emissions = scalar @$possibleEmissionsArrayRef;
print STDERR "There are $number_of_possible_emissions characters in "
    . "the emission alphabet (plus one for tri- and quad- allelic "
    . "states, if any).\n";

my @labels_of_chrs_to_study = ( 1 .. 22 );

#my @labels_of_chrs_to_study = ( 22 );

my @chromosomes_to_study = map { "chr" . $_ } @labels_of_chrs_to_study;

my (
     $average_length_of_a_block,
     $probability_of_a_normal_recombination_per_base,
     $number_of_bad_regions_in_genome,
     $average_length_of_bad_region,
     $number_of_compression_regions_in_genome,
     $average_length_of_compression_region,
     $probability_of_transitioning_into_bad,
     $probability_of_transitioning_out_of_bad,
     $probability_of_transitioning_into_compression,
     $probability_of_transitioning_out_of_compression,
     $probability_of_staying_in_same_normal_state
    )
    = getSuttonianGenomeParameters( $genome_size,
                                    $number_of_meioses_in_pedigree );

my $transitionMatrixRef = &load_transition_matrix(
                                $number_of_states,
                                $probability_of_a_normal_recombination_per_base,
                                $probability_of_staying_in_same_normal_state
);
my %transition_matrix = %$transitionMatrixRef;

my %total_length_of_blocks_in_state = ();

#initialize a hash
for ( my $n = 0 ; $n < $number_of_states ; $n++ ) {
  $total_length_of_blocks_in_state{$n} = 0;
}

# Don't put timestamp in filename to make it easier for worflows
my $outfile                 = "$outputDir/hmm_Suttonian_intervals.txt";
my $header_for_output_files = make_header_for_output_files( $inclusion_order );
open( SYNOPSIS, ">$outfile" );
print SYNOPSIS "$header_for_output_files";
print SYNOPSIS
"#header\tchromosome\tblock_number\tmost_likely_path_state\tbinary state\tnumber_of_markers_in_block\tblock_start\tblock_end\tblock_length\trecombination meiosis\tlength of crossover interval\tcrossover start location\tpastable_UCSC_string\n";

#####################################################################
######################### M A I N   L O O P #########################
#####################################################################

my $initial_state_matrix;
my %count_of_emissions_per_state = ();

my @positions;
my @seq;
my $sequence_length;
my %most_likely_path_state;
my %block_start;
my %block_end;
my %block_state;
foreach my $chromosome_to_study ( @chromosomes_to_study ) {
  print STDERR "Working on chromosome $chromosome_to_study.\n";
  $outfile = "$outputDir/hmm_Suttonian_states_$chromosome_to_study.txt";
  open( OUT, ">$outfile" );

  my ( $positionsRef, $emissionsRef, $numberOfEmissions ) =
      get_sequence_for_this_chromosome(
                                        $chromosome_to_study,
                                        \%observed_emission,
                                        $linkage_suppression,
                                        $maximum_genotype_vector_index,
                                        $MIEHashRef
      );

  # Recast for now
  @positions       = @$positionsRef;
  @seq             = @$emissionsRef;
  $sequence_length = $numberOfEmissions;

  $total_genome_number_of_informative_Suttonian_locations += $numberOfEmissions;

  my %intermarker_distance = compute_intermarker_distances( @positions );

  # Going to calculate probabilities dynamically for now
  my ( $emissionProbMatrixRef ) =
      calculate_probability_of_each_emission_for_each_state_matrix(
                                                           $emissionProbHashRef,
                                                           \@seq );

  # Recast for now
  my %emission_probability_matrix = %$emissionProbMatrixRef;

  # Likelihood of starting in each state
  my ( $initialStateMatrixRef ) = getInitialStateMatrix( $number_of_states );

  # Recast for now
  my %initial_state_matrix = %$initialStateMatrixRef;

  # Viterbi algorithm
  # Compute the probability of the most likely path that ends at
  # each state at each position
  my ( $transitionToMaxProbRef, $probabilityOfMostLikelyPathMatrixRef ) =
      Viterbi_probability_of_most_likely_path_matrix(
        \%emission_probability_matrix,
        \%intermarker_distance, \%initial_state_matrix, \%transition_matrix, $E,
        $number_of_states, \%states_are_adjacent );

  # Dereference ( for the time being )
  my %transition_to_max_prob                 = %$transitionToMaxProbRef;
  my %probability_of_most_likely_path_matrix =
      %$probabilityOfMostLikelyPathMatrixRef;

  # Output the most likely path
  my $mostLikelyPathStateRef = Viterbi_most_likely_path(
                                       \%transition_to_max_prob,
                                       $sequence_length,
                                       \%probability_of_most_likely_path_matrix,
                                       $number_of_states
  );
  %most_likely_path_state = %$mostLikelyPathStateRef;

  print_outfile( $chromosome_to_study );
}    #end of main block for each chromosome

print STDERR
    "\nEMISSION STATITICS\nstate\temission\tcanonical genotype vector\tcount\n";

# Print out emission statistics for use in refining emission probabilities
for ( my $n = 0 ; $n < $number_of_states ; $n++ ) {
  foreach my $emission ( @{$possibleEmissionsArrayRef} ) {
    my ( $genotype_vector, $genotype_content ) = decodeEmissionIndex(
                                             $emission,
                                             $number_of_individuals_in_pedigree,
                                             $maximum_genotype_vector_index,
                                             \%genotype_for_integer
    );
    if ( defined( $count_of_emissions_per_state{$n}{$emission} ) ) {
      print STDERR "$n\t$emission\t$genotype_vector\t"
          . "$count_of_emissions_per_state{$n}{$emission}\n";
    }
    else {
      print STDERR "$n\t$emission\t$genotype_vector\t0\n";
    }
  }
}

print STDERR "\n\nState\ttotal_length_of_blocks_in_state\n";
for ( my $n = 0 ; $n < $number_of_states ; $n++ ) {
  print STDERR "$n\t$total_length_of_blocks_in_state{$n}\n";
}

print STDERR "\n\nTotal number of informative pedigree positions: "
    . "$total_genome_number_of_informative_Suttonian_locations\n";
my $maternal_paternal_ratio = "";
if ( $total_paternal_recombinations > 0 ) {
  $maternal_paternal_ratio = sprintf( "%.3f",
                                      $total_maternal_recombinations /
                                          $total_paternal_recombinations );
}
else {
  $maternal_paternal_ratio = 0;
}
my $total_recombinations =
    $total_maternal_recombinations + $total_paternal_recombinations;
print STDERR "Total maternal and paternal recombinations (and ratio): "
    . "$total_maternal_recombinations\t$total_paternal_recombinations\t"
    . "(total $total_recombinations)\t$maternal_paternal_ratio\n";

close SYNOPSIS;

exit;

#END MAIN BLOCK

#################### GLOBAL S U B R O U T I N E S #######################

## Globals:
##     $chromosome_to_study
##     %count_of_emissions_per_state
##     $chromosome
##     $emission
##     $block_number
##     @state_names
##     $sequence_length
##     $header_for_output_files
##     %block_start
##     $state
##     %most_likely_path_state
##     %probability_of_state
##     $E
##     @seq
##     @positions
##     $number_of_states
##     $probability
##     $number_of_markers_in_block
##     $m ( Used by output_block )
sub print_outfile {
  my $chromosome_to_study         = shift;
  my $compute_state_probabilities = 0;
  print OUT "$header_for_output_files";
  print OUT "#chromosome\t$chromosome_to_study\n";
  my $chromosome = $chromosome_to_study;
  print OUT "#header\tindex\tchromosome\tposition\t"
      . "emission\tstate of most likely path";
  if ( $compute_state_probabilities ) {
    my $n = 0;
    while ( $n < $number_of_states ) {
      print OUT "\tprob of being in state $n $state_names[$n]";
      $n++;
    }
  }
  print OUT "\n";

  if ( !@positions ) {
    warn "Missing positions for $chromosome_to_study!";
    close OUT;
    return;
  }

  my $block_number               = 1;
  my $number_of_markers_in_block = 0;
  $most_likely_path_state{0} = $most_likely_path_state{1};
  $block_start{$block_number} = $positions[ 0 ];
  my $m = 1;
  while ( $m <= $sequence_length ) {
    my $emission = $seq[ $m - 1 ];
    my $state    = $most_likely_path_state{$m};

    if ( defined( $count_of_emissions_per_state{$state}{$emission} ) ) {
      $count_of_emissions_per_state{$state}{$emission}++;
    }
    else {
      $count_of_emissions_per_state{$state}{$emission} = 1;
    }

    print OUT "$m\t$chromosome\t$positions[$m-1]\t$emission\t$state";

    if ( $compute_state_probabilities ) {
      my $n = 0;
      while ( $n < $number_of_states ) {
        my $probability = "Feature_not_implemented";

        #my $probability = $E**$probability_of_state{$n}{$m};
        print OUT "\t$probability";
        $n++;
      }
    }
    print OUT "\n";
    if ( $most_likely_path_state{$m} ne $most_likely_path_state{ $m - 1 } ) {
      ( $block_number, $number_of_markers_in_block ) =
          output_block( $number_of_markers_in_block,
                        $chromosome, $block_number, $m );
    }
    $number_of_markers_in_block++;
    $m++;
  }
  ( $block_number, $number_of_markers_in_block ) =
      output_block( $number_of_markers_in_block,
                    $chromosome, $block_number, $m );
  close OUT;
}    #end print_outfile

## Globals:
##     $number_of_bad_regions_in_genome
##     $linkage_suppression
##     $exclude_partially_called_vectors
##     $ratio_of_mean_to_sigma
##     $file_to_analyze
##     $average_length_of_a_block
##     $tightness_of_deletion_state
##     $project
##     %input_file_header_info
##     $tightness_of_bad_state
##     $exclude_fully_heterozygous_vectors
##     $average_length_of_bad_region
##     $genome_size
##     $skip_insertions
##     $contact
##     $pseudocount
##     $invoking_command_to_record_on_output_header
##     $probability_of_a_normal_recombination_per_base
##     $reference_genome
##     $skip_deletions
##     $mode
##     $VERSION
##     $pedigree_name
sub make_header_for_output_files {
  my $inclusion_order = shift;

  #TO DO - ADD in passthrough reporting for these genotype vector
  #        header file lines of info
  #creator
  #sw_version
  #timestamp
  #contact
  #pedigree_name
  #pedigree
  #sex
  #pedigree_version
  #pedigree_description
  #numerical_base
  #allele_representation
  #genotypes

  my $header = "";

  my $file_string = $input_file_header_info{'input_file'};
  $file_string =~ s/Hundred_Genomes_Project\///;
  $header .= "#input_file\t$file_string\n";
  $header .= "#creator\t$0 version $VERSION\n";
  $header .=
      "#note\tcommand_line\t$invoking_command_to_record_on_output_header\n";
  $header .= "#timestamp\t" . timestamp() . "\n";
  $header .= "#contact\t$contact\n";

  my @inclusionOrder = split( /,/,   $inclusion_order );
  my @pedigreeInds   = split( /\s+/, $input_file_header_info{'pedigree'} );
  my @indSex         = split( /\s+/, $input_file_header_info{'sex'} );
  $header .= "#pedigree ";
  for ( my $i = 0 ; $i <= $#inclusionOrder ; $i++ ) {
    if ( $inclusionOrder[ $i ] =~ /\d+/ ) {
      $header .= $pedigreeInds[ $inclusionOrder[ $i ] ] . " ";
    }
    else {
      $header .= "Absent ";
    }
  }
  $header .= "\n";
  $header .= "#sex ";
  for ( my $i = 0 ; $i <= $#inclusionOrder ; $i++ ) {
    if ( $inclusionOrder[ $i ] =~ /\d+/ ) {
      $header .= $indSex[ $inclusionOrder[ $i ] ] . " ";
    }
    else {
      $header .= "Absent ";
    }
  }
  $header .= "\n";

  if ( defined( $input_file_header_info{'reference_genome'} ) ) {
    $header .=
        "#reference_genome\t$input_file_header_info{'reference_genome'}\n";
  }
  else {
    $header .= "#reference_genome\t$reference_genome\n";
  }
  $header .= "#pedigree_name\t$pedigree_name\n";
  $header .= "#numerical_base\tzero\n";
  $header .= "#project\t$project\n";
  $header .= "#chromosome_nomenclature\tconventional\n";
  $header .= "#inheritance_vector_format\tbinary\n";
  $header .= "#inheritance_state_representation\tarabic_zero_based\n";
  $header .= "#parameter\tgenome_size\t$genome_size\n";
  $header .=
      "#parameter\taverage_length_of_a_block\t$average_length_of_a_block\n";
  $header .=
"#parameter\tprobability_of_a_normal_recombination_per_base\t$probability_of_a_normal_recombination_per_base\n";

#$header .=
#"#parameter\tnumber_of_bad_regions_in_genome\t$number_of_bad_regions_in_genome\n";
#$header .=
#"#parameter\taverage_length_of_bad_region\t$average_length_of_bad_region\n";
  $header .= "#parameter\tpseudocount\t$pseudocount\n";
  $header .= "#parameter\tmode\t$mode\n";
  $header .= "#parameter\tdatasource\t$file_to_analyze\n";
  $header .=
"#parameter\texclude_fully_heterozygous_vectors\t$exclude_fully_heterozygous_vectors\n";
  $header .=
"#parameter\texclude_partially_called_vectors\t$exclude_partially_called_vectors\n";
  $header .= "#parameter\tskip_deletions\t$skip_deletions\n";
  $header .= "#parameter\tskip_insertions\t$skip_insertions\n";
  $header .= "#parameter\tlinkage_suppression\t$linkage_suppression\n";

 #if ( defined( $ratio_of_mean_to_sigma ) )
 #{
 #  print OUT "#parameter\tratio_of_mean_to_sigma\t$ratio_of_mean_to_sigma\n";
 #}
 #
 #if ( defined( $tightness_of_deletion_state ) )
 #{
 #  print OUT
 #    "#parameter\ttightness_of_deletion_state\t$tightness_of_deletion_state\n";
 #}
 #if ( defined( $tightness_of_bad_state ) )
 #{
 #  print OUT "#parameter\ttightness_of_bad_state\t$tightness_of_bad_state\n";
 #}

  return $header;
}    #end make_header_for_output_files

## Globals:
##     %block_state
##     $m
##     $chromosome
##     $UCSC_browser_start
##     $block_number
##     %total_length_of_blocks_in_state
##     $distance_between_blocks
##     $UCSC_browser_end
##     %block_start
##     %most_likely_path_state
##     %block_end
##     $total_maternal_recombinations
##     $block_length
##     @positions
##     $total_paternal_recombinations
##     $number_of_markers_in_block
##     $crossover_start
##     %parental_transition
##     %states_are_adjacent
sub output_block {
  my $number_of_markers_in_block = shift;
  my $chromosome                 = shift;
  my $block_number               = shift;
  my $m                          = shift;
  my $state;
  my $binary_state;
  my $pastable_UCSC_string;

  $block_end{$block_number} = $positions[ $m - 2 ];
  my $block_length =
      $block_end{$block_number} - $block_start{$block_number} + 1;
  $block_state{$block_number} = $most_likely_path_state{ $m - 1 };
  my $recombination_meiosis;
  my $distance_between_blocks;
  my $crossover_start;
  if ( $block_number == 1 ) {
    $recombination_meiosis   = "beginning_of_chromosome";
    $distance_between_blocks = "NA";
    $crossover_start         = "NA";
  }
  else {
    $recombination_meiosis = recombination_meiosis(
                                              $block_state{ $block_number - 1 },
                                              $block_state{$block_number},
                                              \%parental_transition,
                                              \%states_are_adjacent
    );
    $distance_between_blocks =
        $block_start{$block_number} - $block_end{ $block_number - 1 };
    $crossover_start = $block_end{ $block_number - 1 };

    if ( $recombination_meiosis =~ /paternal/ ) {
      $total_paternal_recombinations++;
    }
    elsif ( $recombination_meiosis =~ /maternal/ ) {
      $total_maternal_recombinations++;
    }
  }

  my $UCSC_browser_start = $block_start{$block_number} + 1;
  my $UCSC_browser_end   = $block_end{$block_number} + 1;
  $pastable_UCSC_string = "$chromosome:$UCSC_browser_start-$UCSC_browser_end";
  $state                = $block_state{$block_number};
  $binary_state = string_binary_representation_of_inheritance_state( $state,
                                               $number_of_meioses_in_pedigree );

  print SYNOPSIS
"$chromosome\t$block_number\t$state\t$binary_state\t$number_of_markers_in_block\t$block_start{$block_number}\t$block_end{$block_number}";
  print SYNOPSIS
"\t$block_length\t$recombination_meiosis\t$distance_between_blocks\t$crossover_start\t$pastable_UCSC_string\n";

  if (
    defined( $total_length_of_blocks_in_state{ $block_state{$block_number} } ) )
  {
    $total_length_of_blocks_in_state{ $block_state{$block_number} } +=
        $block_length;
  }
  else {
    $total_length_of_blocks_in_state{ $block_state{$block_number} } =
        $block_length;
  }

  $block_number++;

  # Note that there will be a block start for a non existent
  # block at the very end of each chromosome - a coding hack, if you will
  $block_start{$block_number} = $positions[ $m - 1 ];
  $number_of_markers_in_block = 0;
  return ( $block_number, $number_of_markers_in_block );
}    #end output_block

#################### S U B R O U T I N E S #######################

##-------------------------------------------------------------------------##
## Use: $recombStr = recombination_meiosis( $recom1, $recom2,
##                                        $parentalTransitionHash,
##                                        $statesAreAdjacentHash );
##
##  Returns
##     $recombStr
##
##  Globals: None
##-------------------------------------------------------------------------##
sub recombination_meiosis {
  my $recom1                 = shift;
  my $recom2                 = shift;
  my $parentalTransitionHash = shift;
  my $statesAreAdjacentHash  = shift;

  my $hamming_distance_between_the_states =
      $statesAreAdjacentHash->{$recom1}{$recom2};
  if ( $hamming_distance_between_the_states == 1 ) {

    # Is it maternal or paternal?
    return $parentalTransitionHash->{$recom1}{$recom2};
  }
  else {
    return "$hamming_distance_between_the_states recombinations";
  }

  return "uncharacterized_transition";
}    #end recombination_meiosis

##-------------------------------------------------------------------------##
## Use:   my ( $transitionToMaxProb,
##             $probabilityOfMostLikelyPathMatrix ) =
##                 Viterbi_probability_of_most_likely_path_matrix(
##                     \%emission_probability_matrix,
##                     \%intermarker_distance,
##                     \%initial_state_matrix,
##                     \%transition_matrix,
##                     $E,
##                     $number_of_states,
##                     \%states_are_adjacent );
##
##  Returns
##      \%transitionToMaxProb
##      \%probabilityOfMostLikelyPathMatrix
##
##  Globals: None
##-------------------------------------------------------------------------##
sub Viterbi_probability_of_most_likely_path_matrix {
  my $emissionProbabilityMatrix = shift;
  my $intermarkerDistance       = shift;
  my $initialStateMatrix        = shift;
  my $transitionMatrix          = shift;
  my $E                         = shift;
  my $numberOfStates            = shift;
  my $statesAreAdjacent         = shift;

  my %transitionToMaxProb               = ();
  my %probabilityOfMostLikelyPathMatrix = ();

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $date = timestamp();
  print STDERR "[$subroutine] Begin at $date...\n";
  my $transition;
  my $n = 0;

  # Setup prob of entering from each state
  while ( $n < $numberOfStates ) {
    $probabilityOfMostLikelyPathMatrix{$n}{0} = $initialStateMatrix->{$n};
    $n++;
  }

  # Build matrix
  $n = 1;
  while ( defined $intermarkerDistance->{$n} ) {
    my $interDist = $intermarkerDistance->{$n};

    #print STDERR "[$subroutine] Working on position $n\n";
    my $m = 0;
    while ( $m < $numberOfStates ) {
      my $max_prob = undef;

      # Calculate transition probabilities
      my $p                      = 0;
      my $transitionOutProb      = 0;
      my %dynamicTransitionProbs = ();
      while ( $p < $numberOfStates ) {
        my $transition_probability_per_base = $E**$transitionMatrix->{$p}{$m};
        die if ( $transition_probability_per_base >= 1 );    #debug
        die unless ( $transition_probability_per_base >= 0 );    #debug
             # Probability of a recombination over this intermarker distance
        my $prob_of_single_event =
            1 - ( ( 1 - $transition_probability_per_base )**$interDist );
        if ( $statesAreAdjacent->{$m}{$p} == 1 ) {

          # The states are separated by exactly one crossover
          $transitionOutProb += $prob_of_single_event;
          $dynamicTransitionProbs{$m}{$p} = $prob_of_single_event;
        }
        elsif ( $statesAreAdjacent->{$m}{$p} != 0 ) {

          # The states are separated by double/multiple crossovers
          #
          # Approximated as probability of at least one event
          # I had been thinking that it is OK to use the probability of a
          # single event as the probability of at least one event.
          # This is probably OK. To do it right, I would need to sum all
          # the probabilities for an odd number of crossovers.
          # However, the probability of transitions between states 1 & 4
          # and between 2&3 is not the probability of a double crossover
          # (or more precisely, of an even number of crossovers). It is
          # half this probability. Because half of the double crossovers
          # end up back where they started.
          # So next time I run the code, I should probably divide the
          # double crossover rate by 2
          #
          # Still a bit of an approximation (see above). The farther apart
          # the states are, the higher the exponent and the lower
          # the probability
          $dynamicTransitionProbs{$m}{$p} =
              ( $prob_of_single_event**$statesAreAdjacent->{$m}{$p} ) / 2;
          $transitionOutProb += $dynamicTransitionProbs{$m}{$p};
        }
        elsif ( $m != $p ) {
          die "[$subroutine]: Expecting in-state-transition "
              . "for distance == 0\n";
        }
        if ( $m != $p && $dynamicTransitionProbs{$m}{$p} == 0 ) {
          $dynamicTransitionProbs{$m}{$p} = -10e99;
          warn "Using default transition probability! "
              . "$dynamicTransitionProbs{$m}{$p}\n";
        }
        $p++;
      }

      # Calculate the probability of remaining in our state
      $dynamicTransitionProbs{$m}{$m} = 1 - $transitionOutProb;
      if ( $dynamicTransitionProbs{$m}{$m} == 0 ) {
        $dynamicTransitionProbs{$m}{$m} = -10e99;
        warn "[$subroutine]: Using default transition probability! "
            . "$dynamicTransitionProbs{$m}{$m}\n";
      }

      # Multiply by prob of transition
      $p = 0;
      while ( $p < $numberOfStates ) {

        #print "Transition Prob ( before log ) : " .
        #       $dynamicTransitionProbs{$m}{$p}\n";
        my $transition_probability = log( $dynamicTransitionProbs{$m}{$p} );

        my $prob = $probabilityOfMostLikelyPathMatrix{$p}{ $n - 1 } +
            $transition_probability;
        if ( ( not defined( $max_prob ) ) or $prob > $max_prob ) {
          $max_prob   = $prob;
          $transition = $p;
        }
        $p++;
      }

      #print "PROB AT $n:  max_prob $max_prob ( " . exp($max_prob) .
      #      " ) state $transition\n";

      # Multiply by prob of emission at new state
      $probabilityOfMostLikelyPathMatrix{$m}{$n} =
          $max_prob + $emissionProbabilityMatrix->{$m}{$n};
      $transitionToMaxProb{$m}{$n} = $transition;
      $m++;
    }
    $n++;
  }
  $date = timestamp();
  print STDERR "[$subroutine] Done at $date.\n";
  return ( \%transitionToMaxProb, \%probabilityOfMostLikelyPathMatrix );
}    # Viterbi_probability_of_most_likely_path_matrix

##-------------------------------------------------------------------------##
## Use:   my $mostLikelyPathStateRef = Viterbi_most_likely_path(
##                   \%transition_to_max_prob,
##                   $sequence_length,
##                   \%probability_of_most_likely_path_Matrix,
##                   $number_of_states );
##
##
##  Returns
##      \%mostLikelyPathState
##
##  Globals: None
##-------------------------------------------------------------------------##
sub Viterbi_most_likely_path {
  my $transitionToMaxProb               = shift;
  my $sequenceLength                    = shift;
  my $probabilityOfMostLikelyPathMatrix = shift;
  my $numberOfStates                    = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  my $date = timestamp();
  print STDERR "[$subroutine] Begin at $date...\n";

  my $transition;
  my %mostLikelyPathState = ();

  # Figure out which state to start from
  my $max_prob = undef;
  my $p        = 0;
  while ( $p < $numberOfStates ) {
    my $prob = $probabilityOfMostLikelyPathMatrix->{$p}{$sequenceLength};
    if ( ( not defined( $max_prob ) ) or $prob > $max_prob ) {
      $max_prob   = $prob;
      $transition = $p;
    }
    $p++;
  }
  $mostLikelyPathState{$sequenceLength} = $transition;

  my $n = $sequenceLength - 1;
  while ( $n >= 1 ) {
    my $next_state = $mostLikelyPathState{ $n + 1 };

    #print STDERR "At position $n, the next state in the most " .
    #             "likely path is $next_state.\n";

    $mostLikelyPathState{$n} = $transitionToMaxProb->{$next_state}{ $n + 1 };
    $n--;
  }
  $date = timestamp();
  print STDERR "[$subroutine] Done at $date.\n";
  return ( \%mostLikelyPathState );
}    # Viterbi_most_likely_path

##-------------------------------------------------------------------------##
## Use: my ( $average_length_of_a_block,
##           $probability_of_a_normal_recombination_per_base,
##           $number_of_bad_regions_in_genome,
##           $average_length_of_bad_region,
##           $number_of_compression_regions_in_genome,
##           $average_length_of_compression_region,
##           $probability_of_transitioning_into_bad,
##           $probability_of_transitioning_out_of_bad,
##           $probability_of_transitioning_into_compression,
##           $probability_of_transitioning_out_of_compression,
##           $probability_of_staying_in_same_normal_state ) =
##           getSuttonianGenomeParameters( $genomeSize,
##                                         $numOfMeiosesInPedigree,
##                                         [$useCompressionAndBadRegions] );
##
##  Returns
##      Multiple values ( see above )
##
##  Globals: None
##-------------------------------------------------------------------------##
sub getSuttonianGenomeParameters {
  my $genomeSize                  = shift;
  my $numOfMeiosesInPedigree      = shift;
  my $useCompressionAndBadRegions = shift;

  my $average_length_of_a_block                      = 1000000000000000;
  my $probability_of_a_normal_recombination_per_base =
      1 / $average_length_of_a_block;

  my $probability_of_staying_in_same_normal_state = 0;

  my $number_of_bad_regions_in_genome                 = 0;
  my $average_length_of_bad_region                    = 0;
  my $number_of_compression_regions_in_genome         = 0;
  my $average_length_of_compression_region            = 0;
  my $probability_of_transitioning_into_bad           = 0;
  my $probability_of_transitioning_out_of_bad         = 0;
  my $probability_of_transitioning_into_compression   = 0;
  my $probability_of_transitioning_out_of_compression = 0;

  if ( defined $useCompressionAndBadRegions ) {

    # HMM seems to go into error regions way too easily and stay in
    # them too short a time, so I am dialing down the number and
    # increasing the length to see if I get better performance
    $number_of_bad_regions_in_genome         = 0.00000008;
    $average_length_of_bad_region            = 100000;
    $number_of_compression_regions_in_genome = 0.00000008;
    $average_length_of_compression_region    = 400000;

    $probability_of_transitioning_into_bad =
        $number_of_bad_regions_in_genome / $genomeSize;
    $probability_of_transitioning_out_of_bad =
        1 / $average_length_of_bad_region;

    $probability_of_transitioning_into_compression =
        $number_of_compression_regions_in_genome / $genomeSize;
    $probability_of_transitioning_out_of_compression =
        1 / $average_length_of_compression_region;

    # If using bad regions, and hard-coded for 4 person nuclear family
    $probability_of_staying_in_same_normal_state =
        1 - 2 * $probability_of_a_normal_recombination_per_base -
        $probability_of_a_normal_recombination_per_base *
        $probability_of_a_normal_recombination_per_base -
        $probability_of_transitioning_into_bad -
        $probability_of_transitioning_into_compression;
  }
  else {

    # If not using bad regions
    if ( $numOfMeiosesInPedigree == 4 ) {

      #  Two child families have a special symmetry property not
      #  present in bigger fmailies
      $probability_of_staying_in_same_normal_state =
          1 - 2 * $probability_of_a_normal_recombination_per_base
          ;  #assume the probability of a double transition is vanishingly small
    }
    else {

      # The number of meioses is the number of different states that
      # can be reached via a single recombination
      $probability_of_staying_in_same_normal_state =
          1 - $numOfMeiosesInPedigree *
          $probability_of_a_normal_recombination_per_base
          ;  #assume the probability of a double transition is vanishingly small
    }
  }

  return (
           $average_length_of_a_block,
           $probability_of_a_normal_recombination_per_base,
           $number_of_bad_regions_in_genome,
           $average_length_of_bad_region,
           $number_of_compression_regions_in_genome,
           $average_length_of_compression_region,
           $probability_of_transitioning_into_bad,
           $probability_of_transitioning_out_of_bad,
           $probability_of_transitioning_into_compression,
           $probability_of_transitioning_out_of_compression,
           $probability_of_staying_in_same_normal_state
  );
}

##-------------------------------------------------------------------------##
## Use: my $transitionMatrixRef = &load_transition_matrix( $number_of_states,
##                         $probability_of_a_normal_recombination_per_base,
##                         $probability_of_staying_in_same_normal_state );
##
##  Returns
##      \%transitionMatrix
##
##  Globals: None
##-------------------------------------------------------------------------##
sub load_transition_matrix {
  my $numberOfStates                           = shift;
  my $probabilityOfANormalRecombinationPerBase = shift;
  my $probabilityOfStayingInSameNormalState    = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $local_timestamp = timestamp();
  print STDERR "[$subroutine] Begin at $local_timestamp.\n";

  # Probabilities should be logs
  my $log_of_probability_of_staying_in_same_normal_state =
      log $probabilityOfStayingInSameNormalState;
  my $log_of_probability_of_a_normal_recombination_per_base =
      log $probabilityOfANormalRecombinationPerBase;

  # Note that we are storing probabilities in log space
  my $from_state       = 0;
  my $to_state         = 0;
  my %transitionMatrix = ();
  while ( $from_state < $numberOfStates ) {
    $to_state = $from_state;
    while ( $to_state < $numberOfStates ) {
      if ( $from_state == $to_state ) {
        $transitionMatrix{$from_state}{$to_state} =
            $log_of_probability_of_staying_in_same_normal_state;
      }
      else {
        $transitionMatrix{$from_state}{$to_state} =
            $log_of_probability_of_a_normal_recombination_per_base
            ; #note that this gets overid in the core code to deal with length between markers
        $transitionMatrix{$to_state}{$from_state} =
            $log_of_probability_of_a_normal_recombination_per_base
            ; #note that this gets overid in the core code to deal with length between markers
      }
      $to_state++;
    }
    $from_state++;
  }

  if ( $DEBUG ) {

    #print the transition matrix
    my $n = 0;
    while ( $n < $numberOfStates ) {
      my $m = 0;
      while ( $m < $numberOfStates ) {
        print "from $n to $m: $transitionMatrix{$n}{$m}\n";
        $m++;
      }
      $n++;
    }
  }

  $local_timestamp = timestamp();
  print STDERR "[$subroutine] End at $local_timestamp...\n";
  return ( \%transitionMatrix );
}    # load_transition_matrix

##-------------------------------------------------------------------------##
## Use: my ( $inputFileHdrInfoRef, $observedEmissionsRef,
##           $informativePositions ) = load_genotype_vectors(
##                  $file_to_analyze,
##                  $skip_insertions,
##                  $skip_deletions,
##                  $exclude_fully_heterozygous_vectors,
##                  $exclude_partially_called_vectors,
##                  $maximum_genotype_vector_index,
##                  $ignore_MIEs,
##                  \%MIE,
##                  \%lowest_integer_for_equivalent_genotype_vector,
##                  $inclusion_order );
##
##  Read an ISB genotype vector file and translate into HMM
##  emissions.
##
##  NOTES:
##     - If inclusion_order is provided the routine will
##       take a slice of the gentoype patterns.
##     - Filters out insertion/deletion loci, fully heterozygous loci,
##       and or partially called loci upon request.
##     - Currently only handles loci with at most four alleles ( a-d ).
##
##  Returns
##
##  Globals: None
##-------------------------------------------------------------------------##
sub load_genotype_vectors {
  my $file                            = shift;
  my $skipInsertions                  = shift;
  my $skipDeletions                   = shift;
  my $excludeFullyHeterozygousVectors = shift;
  my $excludePartiallyCalledVectors   = shift;
  my $maximumGenotypeVectorIndex      = shift;
  my $excludeMIEs                     = shift;
  my $MIEHashRef                      = shift;
  my $genotypeToEmissionHashRef       = shift;
  my $inclusionOrderStr               = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $local_timestamp = timestamp();
  print STDERR "[$subroutine] Begin at $local_timestamp.\n";

  my %inputFileHdrInfo = ();
  $inputFileHdrInfo{'input_file'} = $file;

  if ( $file =~ /.*\.bz2/ ) {
    open IN, "bunzip2 -c $file |"
        or die "[$subroutine]: Could not open $file with bunzip2!: $!\n";
  }
  elsif ( $file =~ /.*\.gz/ ) {
    open IN, "gunzip -c $file |"
        or die "[$subroutine]: Could not open $file with gunzip!: $!\n";
  }
  else {
    open IN, "<$file"
        or die "[$subroutine]: Could not open $file: $!\n";
  }

  my @inclusionOrderArray;
  if ( defined( $inclusionOrderStr ) ) {
    @inclusionOrderArray = split( ",", $inclusionOrderStr );
    print STDERR "Inclusion order:\t", "$inclusionOrderStr\n";
  }

  my $fullyCalledRE        = "";
  my $fullyHeterozygousPat = "";
  foreach my $ind ( @inclusionOrderArray ) {
    if ( $ind eq "x" ) {
      $fullyCalledRE        .= "nn";
      $fullyHeterozygousPat .= "nn";
    }
    else {
      $fullyCalledRE        .= "[ab][ab]";
      $fullyHeterozygousPat .= "ab";
    }
  }

  my $genotypeVectorStr = "";
  my @input_genotype_vector;
  my @genotype_vector;
  my ( $allele1, $allele2 );
  my $informative_positions = 0;
  my $line_of_data          = 1;
  my %observedEmissions     = ();
  while ( <IN> ) {
    s/[\n\r]+//g;
    next unless ( /\S/ );

    # Make hash of header info
    if ( /^#/ ) {

      # Skip blank header lines
      next unless ( /^\#\s*\S/ );

      # Skip keys with no values
      unless ( /^\#\s*\S+\s+\S+/ ) {
        print STDERR "$subroutine(): This header line has a key with no "
            . "value:\n$_\n";
        next;
      }

      if ( /^\#\s*(\S*)\s+(.*)$/ ) {
        $inputFileHdrInfo{$1} = $2;
      }
      else {
        print STDERR "$_ is a bad header line!\n";
        next;
      }
      next;
    }

    my @data = split( /,/, $_ );

    $genotypeVectorStr = $data[ 3 ];

    if ( defined( $inclusionOrderStr ) ) {
      my @inputGenotypeVectorArray = split( //, $genotypeVectorStr );
      my @genotypeVectorArray = ();
      while ( @inputGenotypeVectorArray ) {
        $allele1 = shift @inputGenotypeVectorArray;
        $allele2 = shift @inputGenotypeVectorArray;
        push @genotypeVectorArray, $allele1 . $allele2;
      }

      # Get genotypes for individuals if we need to
      if ( $skipDeletions || $skipInsertions ) {
        my $genotypes = "";
        foreach my $ind ( @inclusionOrderArray ) {

          # RMH: Nomenclature for missing person
          next if ( $ind eq "x" );
          $genotypes .=
              $data[ ( $ind * 4 ) + 4 ] . $data[ ( $ind * 4 ) + 4 + 2 ];
        }
        if ( $skipDeletions ) {
          next if ( $genotypes =~ /-/i );
        }
        if ( $skipInsertions ) {
          next if ( $genotypes =~ /\+/i );
        }
      }

      # RMH: Allow for uncalled missing parents
      $genotypeVectorStr = "";
      foreach my $ind ( @inclusionOrderArray ) {
        if ( $ind eq "x" ) {
          $genotypeVectorStr .= "nn";
        }
        else {
          $genotypeVectorStr .= $genotypeVectorArray[ $ind ];
        }
      }
    }
    else {

      # Assumption here is that all individuals in this
      # genotype file are going to be considered for analysis.
      if ( $skipDeletions ) {
        next if ( $_ =~ /-/i );
      }

      if ( $skipInsertions ) {
        next if ( $_ =~ /\+/i );
      }
    }

    #
    # Sanitize the genotype patterns.  Since we allow slicing of the
    # genotype file we may receive genotype patterns with gaps in them.
    #
    # ie.  bbbbbbbb => aaaaaaaa
    #      bcbcbcbc => abababab
    #      cdddcdcd => abbbabab
    # It doesn't matter to the hmm which is the reference allele
    # so we just assign the allele characters on a first come
    # first serve basis.  We also need to keep the allele characters
    # sorted on a per individual basis.
    #
    my %alleleAlphabet    = ();
    my @newAlphabet       = ( "a", "b", "c", "d" );
    my @indAlleles        = ();
    my $newGenotypeVector = "";
    foreach my $alleleChar ( split( //, $genotypeVectorStr ) ) {
      if ( $alleleChar eq "n" ) {
        push @indAlleles, "n";
      }
      elsif ( defined $alleleAlphabet{$alleleChar} ) {
        push @indAlleles, $alleleAlphabet{$alleleChar};
      }
      else {
        if ( @newAlphabet ) {
          $alleleAlphabet{$alleleChar} = shift( @newAlphabet );
          push @indAlleles, $alleleAlphabet{$alleleChar};
        }
        else {
          die "$subroutine(): Currently cannot handle genotype vector "
              . "patterns with greater than four alleles.\n";

          # Sooner or later there will be such a vector perhaps with
          # deletion/insertion characters.  Think about this.
        }
      }
      if ( ( scalar @indAlleles ) == 2 ) {
        $newGenotypeVector .= join( "", sort @indAlleles );
        undef @indAlleles;
      }
    }
    $genotypeVectorStr = $newGenotypeVector;

    # The fully heterozygous vector seems to be responsible for the
    # most false short blocks (likely to be compression or CNV blocks)
    if ( $excludeFullyHeterozygousVectors ) {
      next if ( $genotypeVectorStr eq $fullyHeterozygousPat );
    }

    if ( $excludePartiallyCalledVectors ) {
      next unless ( $genotypeVectorStr =~ /$fullyCalledRE/i );
    }

    # Only consider sites with at least two called alleles.
    next
        unless ( $genotypeVectorStr =~ /a/ and $genotypeVectorStr =~ /b/ )
        ;    #only wish to consider informative positions

    unless ( defined( $genotypeToEmissionHashRef->{$genotypeVectorStr} ) ) {

      # Do some overrides here for genotype vectors that contain
      # non-{a,b,n} characters.  I.e quad/tri-allelic sites.
      if (     $genotypeVectorStr =~ /a/
           and $genotypeVectorStr =~ /b/
           and $genotypeVectorStr =~ /c/ )
      {

        # Emission of a quad- or tri-allelic state
        $genotypeToEmissionHashRef->{$genotypeVectorStr} =
            $maximumGenotypeVectorIndex + 1;
      }
      else {
        die "$subroutine(): Odd that there is no emission assigned "
            . "to this pattern:\t$genotypeVectorStr\n";
      }
    }
    my $observedEmission = $genotypeToEmissionHashRef->{$genotypeVectorStr};

    if ( $observedEmission == $maximumGenotypeVectorIndex ) {

      # These mostly non-informative positions can overwhelm the HMM.
      # They provide little to no information on state boundaries, and
      # very little information otherwise.
      #
      # For example, abnnabnn is essentially worthless in information and
      # has been previously been converted into the equivalent of all
      # "n"s. ie.:
      # maximum_genotype_vector_index = "nnnnnnnn.....nn"
      next;
    }

    if ( $excludeMIEs ) {
      if ( $genotypeVectorStr !~ /c/ ) {
        my $genotypeVectorIndex =
            GenotypeEmissions::encodeTwoAlleleEmissionPat( $genotypeVectorStr );
        if ( $MIEHashRef->{$genotypeVectorIndex} ) {
          next;
        }
      }
    }

    #print "Genotype Passed Filters: $genotypeVectorStr " .
    #      "( $observedEmission )\n";
    my $chr      = $data[ 0 ];
    my $position = $data[ 1 ];
    $observedEmissions{$chr}{$position} = $observedEmission;
    $informative_positions++;

    $line_of_data++;
  }
  close IN;

  ## TODO Check that filters were not too restrictive and dataset
  ##      is empty.

  $local_timestamp = timestamp();
  print STDERR "[$subroutine] end at $local_timestamp\n";

  return ( \%inputFileHdrInfo, \%observedEmissions, $informative_positions );

}    #end load_genotype_vectors

##-------------------------------------------------------------------------##
## Use: $binary_state =
##         string_binary_representation_of_inheritance_state(
##                           $state, $number_of_meioses_in_pedigree );
##
##  Returns
##      $binary_state
##
##  Globals: None
##-------------------------------------------------------------------------##
sub string_binary_representation_of_inheritance_state {
  my $decimal                   = shift;
  my $numberOfMeiosesInPedigree = shift;

  # Note that this only works with (Using sprintf) Perl 5.6+  If
  # using an earlier version , need to use another recipe for
  # converting to binary
  my $binary = sprintf( "%b", $decimal );

  # Now need to pad with leading zeros
  while ( length $binary < $numberOfMeiosesInPedigree ) {
    $binary = "0" . $binary;
  }
  return $binary;
}    # string_binary_representation_of_inheritance_state

##-------------------------------------------------------------------------##
## Use: $binary_state =
##         binary_representation_of_inheritance_state(
##                           $decimal, $number_of_meioses_in_pedigree );
##
##  Returns
##      $binary_state
##
##  Globals: None
##-------------------------------------------------------------------------##
sub binary_representation_of_inheritance_state {
  my $decimal                   = shift;
  my $numberOfMeiosesInPedigree = shift;
  my $binary                    =
      sprintf( "%b", $decimal )
      ; #note that this only works with (Using sprintf) Perl 5.6+  If using an earlier version , need to use another recipe for converting to binary
        #now need to pad with leading zeros
  while ( length $binary < $numberOfMeiosesInPedigree ) {
    $binary = "0" . $binary;
  }
  my @binary = split( //, $binary );
  return @binary;
}    #end binary_representation_of_inheritance_state

##-------------------------------------------------------------------------##
## Use: my ( $positionsRef, $emissionsRef, $numberOfEmissions ) =
##          get_sequence_for_this_chromosome( $chromosome_to_study,
##                                    \%observed_emission,
##                                    $linkage_suppression,
##                                    $maximum_genotype_vector_index,
##                                    \%MIE );
##
##  Returns
##
##  Globals: None
##-------------------------------------------------------------------------##
sub get_sequence_for_this_chromosome {
  my $chr                        = shift;
  my $observedEmissionsRef       = shift;
  my $linkageSuppression         = shift;
  my $maximumGenotypeVectorIndex = shift;
  my $MIEHashRef                 = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $local_timestamp = timestamp();
  print STDERR "[$subroutine] Begin at $local_timestamp.\n";

  my @positions = sort { $a <=> $b } keys %{ $observedEmissionsRef->{$chr} };
  print STDERR "[$subroutine] Starting with "
      . scalar( @positions )
      . " positions\n";

  # Need to make the various ignoring algorithms optional on the command line

  # The order of suppressing MIEs and correcting for linakge makes
  # a difference; I am not sure which is better

  my $compsci_number_of_positions = scalar @positions - 1;

  # Suppress highly linked emissions
  my %suppress_position;
  if ( $linkageSuppression ) {

    # 20 for $flanking_length_for_linkage_suppression and 4 for
    # $flanking_count_for_suppression seem to produce bizarre results
    # so make these much less agressive (2 & 1) - pretty much just
    # smoothing out identical adjacent emissions
    # with good post-processing as now planned, linkage is better
    # handled in post-processing than by agressive code here
    # so in the long run get rid of this code block entirely
    my $flanking_length_for_linkage_suppression = 2;
    my $flanking_count_for_suppression          = 1;

    # If emission is flanked on both sides by lots of identical
    # emissions, then ignore it

    my $position_now = 0;
    foreach my $position ( @positions ) {
      my $emission = $observedEmissionsRef->{$chr}{$position};

      if ( $emission == $maximumGenotypeVectorIndex ) {

        # These mostly non-informative positions can overwhelm the
        # HMM. They provide little to no information on state
        # boundaries, and very little information otherwise.
        # NOTE: This is caught already...and shouldn't appear here
        die;
      }

      # Count identical positions in leading positions
      my $n                      = 1;
      my $leading_flanking_count = 0;
      while ( ( $position_now - $n >= 0 )
              and $n <= $flanking_length_for_linkage_suppression )
      {
        my $other_position = $positions[ $position_now - $n ];
        my $other_emission = $observedEmissionsRef->{$chr}{$other_position};
        if ( $other_emission == $emission ) { $leading_flanking_count++ }
        $n++;
      }
      $n = 1;
      my $trailing_flanking_count = 0;
      while ( ( $position_now + $n <= $compsci_number_of_positions )
              and $n <= $flanking_length_for_linkage_suppression )
      {
        my $other_position = $positions[ $position_now + $n ];
        my $other_emission = $observedEmissionsRef->{$chr}{$other_position};
        if ( $other_emission == $emission ) { $trailing_flanking_count++ }
        $n++;
      }
      if (     ( $leading_flanking_count >= $flanking_count_for_suppression )
           and ( $trailing_flanking_count >= $flanking_count_for_suppression ) )
      {

        # Do not decouple this code from the code for actually using this value
        $suppress_position{$position} = 1;
      }

      $position_now++;
    }
  }

  # Ignore MIEs and optionally Linked Sites
  my @new_positions = ();
  foreach my $position ( @positions ) {
    next
        if ( defined $suppress_position{$position} );
    my $emission = $observedEmissionsRef->{$chr}{$position};
    unless ( $MIEHashRef->{$emission} ) {
      push @new_positions, $position;
    }
  }

  my @emissions = map { $observedEmissionsRef->{$chr}{$_} } @new_positions;
  my $numberOfEmissions = scalar @emissions;

  $local_timestamp = timestamp();
  print STDERR "[$subroutine] End loading sequence of "
      . "$numberOfEmissions positions at $local_timestamp.\n";

  return ( \@new_positions, \@emissions, $numberOfEmissions );
}    #end get_sequence_for_this_chromosome

##-------------------------------------------------------------------------##
## Use: my %intermarker_distance =
##           compute_intermarker_distances( @positions );
##
##  Returns
##     %intermarker_distance
##
##  Globals: None
##-------------------------------------------------------------------------##
sub compute_intermarker_distances {
  my @positions = @_;

  my %intermarkerDistance = ( 1 => 1 );
  my $n = 2;
  while ( $n <= ( $#positions + 1 ) ) {
    $intermarkerDistance{$n} = $positions[ $n - 1 ] - $positions[ $n - 2 ];
    $n++;
  }
  return ( %intermarkerDistance );
}    # compute_intermarker_distances

##-------------------------------------------------------------------------##
## Use:   my ( $emissionProbMatrixRef ) =
##         calculate_probability_of_each_emission_for_each_state_matrix(
##                                             \%emission_prob, \@seq );
##
##  Returns
##      $emissionProbMatrixRef
##
##  Globals: None
##-------------------------------------------------------------------------##
sub calculate_probability_of_each_emission_for_each_state_matrix {
  my $emissionProbRef = shift;
  my $emissionSeqRef  = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  print STDERR "[$subroutine] Compiling emission probabilities...\n";

  my $numberOfStates    = scalar( keys( %$emissionProbRef ) );
  my $numberOfPositions = $#$emissionSeqRef;

  my %emission_probability_matrix = ();

  my $n = 0;
  while ( $n < $numberOfStates ) {
    my $m = 1;
    while ( $m <= ( $numberOfPositions + 1 ) ) {
      $emission_probability_matrix{$n}{$m} = sprintf( "%.6f",
                    log( $emissionProbRef->{$n}{ $emissionSeqRef->[ $m - 1 ] } )
      );
      $m++;
    }
    $n++;
  }

  if ( $DEBUG ) {
    print STDERR "probability of emitting  states";
    my $n = 0;
    while ( $n < $numberOfStates ) {
      print STDERR "\t#$n";
      $n++;
    }
    my $m          = 1;
    my $maxToPrint = 10;
    while ( $m <= ( $numberOfPositions + 1 ) ) {
      print STDERR "probabilities of emitting $emissionSeqRef->[$m-1] "
          . "at position $m:";
      my $n = 0;
      while ( $n < $numberOfStates ) {
        print STDERR "\t" . exp( $emission_probability_matrix{$n}{$m} );
        $n++;
      }
      $m++;
    }
    last if ( !$maxToPrint-- );
  }

  print STDERR "[$subroutine] Done.\n";
  return ( \%emission_probability_matrix );
}    # calculate_probability_of_each_emission_for_each_state_matrix

##-------------------------------------------------------------------------##
## Use: my ( $initialStateMatrixRef ) =
##          getInitialStateMatrix( $number_of_states,
##                                $base_probability );
##
##  Returns
##      $initialStateMatrixRef
##
##  Globals: None
##-------------------------------------------------------------------------##
sub getInitialStateMatrix {
  my $numberOfStates = shift;

  # Originally the code handled two bad states
  #  $initial_state_matrix{1} = 0.25;
  #  $initial_state_matrix{2} = 0.25;
  #  $initial_state_matrix{3} = 0.25;
  #  $initial_state_matrix{4} = 0.25;
  #  $initial_state_matrix{5} = $probability_of_transitioning_into_bad;
  #  $initial_state_matrix{6} = $probability_of_transitioning_into_bad;

  my $baseProbability    = log( 1 / $numberOfStates );
  my $m                  = 0;
  my %initialStateMatrix = ();
  while ( $m < $numberOfStates ) {
    $initialStateMatrix{ $m++ } = $baseProbability;
  }

  return ( \%initialStateMatrix );
}    # initial_state_matrix

##-------------------------------------------------------------------------##
## Use: my ( $genotype_vector, $genotype_content ) =
##            decodeEmissionIndex( $emission,
##                                 $number_of_individuals_in_pedigree,
##                                 $maximum_genotype_vector_index,
##                                 \%genotype_for_integer );
##
##    NOTE: This function only works for two allele genotype patterns.
##
##  Returns
##     $genotype_vector
##     $genotype_content
##
##  Globals: None
##-------------------------------------------------------------------------##
sub decodeEmissionIndex {
  my $index                         = shift;
  my $numberOfIndividualsInPedigree = shift;
  my $maximumGenotypeVectorIndex    = shift;
  my $genotypeForIntegerRef         = shift;

  if ( $index == $maximumGenotypeVectorIndex + 1 ) {

    # Very odd genotype vector with more than 2 variants at this
    # position -> probably a compression
    return ( "tri_quad", "" );
  }

  my $n = $numberOfIndividualsInPedigree - 1;

  my %genotype_content = ();
  my $genotype_vector  = "";
  my $remainder        = $index;
  for ( my $i = $n ; $i >= 0 ; $i-- ) {

    # Convert from base 10 to base 6
    my $baseSixDigit = int( $remainder / ( 6**$i ) );
    $remainder = $index % ( 6**$i );

    # Tabulate the genotype ( ie. 0 = "aa", 1 = "ab" ) content
    # of this emission.
    $genotype_content{$baseSixDigit}++;

    # Decode base-6 into genotype strings ie "0" => "aa"
    $genotype_vector =
        $genotypeForIntegerRef->{$baseSixDigit} . $genotype_vector;
  }

  my $genotype_content =
        ( $genotype_content{0} || "0" )
      . ( $genotype_content{1} || "0" )
      . ( $genotype_content{2} || "0" )
      . ( $genotype_content{3} || "0" )
      . ( $genotype_content{4} || "0" )
      . ( $genotype_content{5} || "0" );

  return ( $genotype_vector, $genotype_content );
}    # decodeEmissionIndex

##-------------------------------------------------------------------------##
## Use:  my $numberOfStates = getNumberOfStates( $numberOfNonfounders );
##
##   For a nuclear family of four there will be 4 inheritance states.
##
##   This subroutine that will need to be fleshed out for pedigrees
##   other than nuclear families.
##
##   Future Sex Chromosome Changes:
##
##    Paternal Chromosomes YX  Here Y = 0 and X = 1
##    Maternal Chromosomes XX
##
##    P/M Sex P/M Sex X-State     PAR-state   Y-State
##    -----------------------------------------------
##    0 0 m   0 0 m   Identical   Identical   Identical
##    1 0 f   0 0 m   Identical   Haplo_Mat   Null
##    0 0 m   0 1 m   Non-Ident   Haplo_Pat   Identical
##    1 0 f   0 1 m   Non-Ident   Non-Ident   Null
##    0 0 m   1 0 f   Haplo_Mat   Haplo_Mat   Null
##    1 0 f   1 0 f   Identical   Identical   Null
##    0 0 m   1 1 f   Non-Ident   Non-Ident   Null
##    1 0 f   1 1 f   Haplo_Pat   Haplo_Pat   Null
##                    -------------------------------
##                     4 states   4 states    2 states
##
##  Returns
##
##  Globals Used: None
##-------------------------------------------------------------------------##
sub getNumberOfStates {
  my $numberOfNonfounders = shift;

  return ( 4**( $numberOfNonfounders - 1 ) );
}    # number_of_states

##-------------------------------------------------------------------------##
## Use: my $hammingDistance = Hamming_distance_arrays( \@array1,
##                                                     \@array2 );
##
##  Returns
##     $hammingDistance
##
##  Globals: None
##-------------------------------------------------------------------------##
sub Hamming_distance_arrays {
  my ( $first, $second ) = @_;
  no warnings;    # silence spurious -w undef complaints
  die unless @$first == @$second;
  my $Hamming_distance = 0;
  for ( my $i = 0 ; $i < @$first ; $i++ ) {
    $Hamming_distance++ if $first->[ $i ] ne $second->[ $i ];
  }
  return $Hamming_distance;
}    # end Hamming_distance_arrays

##-------------------------------------------------------------------------##
## Use: my ( \%parental_transition, \%states_are_adjacent ) =
##                 state_to_state_Hamming_distances(
##                         $numberOfStates,
##                         $numberOfNonFounders );
##
##  Returns
##
##  Globals: None
##-------------------------------------------------------------------------##
sub state_to_state_Hamming_distances {
  my $numberOfStates      = shift;
  my $numberOfNonfounders = shift;

  my $numberOfMeiosesInPedigree = 2 * $numberOfNonfounders;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  my $local_timestamp = timestamp();
  print STDERR "[$subroutine] Begin at $local_timestamp.\n";

  # List of states next to each other (i.e., separated by one recombination)

  my $Hamming_distance;
  my $lowest_Hamming_distance;
  my @binary_representation_of_to_inheritance_state;
  my $from_state = 0;
  my $to_state;
  my $n;
  my @from_state1;
  my @from_state2;
  my @from_state3;
  my @from_state4;
  my %parental_transition = ();
  my %states_are_adjacent = ();

  # Need to test sixteen distances and take the shortest (there
  # may be a considerable speedup I haven't bothered to find)
  #   four representations of from state
  #   four representations of to state
  my @four_from_states = ();

  # OK, there is an obvious speedup (probably more) - I only
  # need to check four combinations and can hold the second fixed
  while ( $from_state < $numberOfStates ) {
    my @binary_representation_of_from_inheritance_state =
        binary_representation_of_inheritance_state( $from_state,
                                                   $numberOfMeiosesInPedigree );
    @four_from_states = ();
    push @four_from_states, \@binary_representation_of_from_inheritance_state;
    @from_state1 = @binary_representation_of_from_inheritance_state;
    push @four_from_states, \@from_state1;
    @from_state2 = @binary_representation_of_from_inheritance_state;
    $n           = 0;
    while ( $n < $numberOfMeiosesInPedigree ) {

      # Recall that $numberOfMeiosesInPedigree = the length of the
      # inheritance vector

      if ( $from_state2[ $n ] == 0 ) {

        # Zero means the paternal allele of this child is the first
        # allele of the the father
        $from_state2[ $n ] = 1;
      }
      else {

        # If not zero, must be 1
        $from_state2[ $n ] = 0;
      }

      # Skip the alleles of the other parent
      $n += 2;
    }
    push @four_from_states, \@from_state2;
    @from_state3 = @binary_representation_of_from_inheritance_state;
    $n           = 1;
    while ( $n < $numberOfMeiosesInPedigree ) {

      # Recall that $numberOfMeiosesInPedigree = the length of the
      # inheritance vector
      if ( $from_state3[ $n ] == 0 ) {

        # Zero means the paternal allele of this child is the first
        # allele of the the father
        $from_state3[ $n ] = 1;
      }
      else {    #if not zero, must be 1
        $from_state3[ $n ] = 0;
      }
      $n += 2;    #skip the alleles of the other parent
    }
    push @four_from_states, \@from_state3;
    @from_state4 = @binary_representation_of_from_inheritance_state;
    $n           = 0;
    while ( $n < $numberOfMeiosesInPedigree ) {

      # Recall that $numberOfMeiosesInPedigree = the length of the
      # inheritance vector
      if ( $from_state4[ $n ] == 0 ) {

        # Zero means the paternal allele of this child is the first
        # allele of the the father
        $from_state4[ $n ] = 1;
      }
      else {

        # If not zero, must be 1
        $from_state4[ $n ] = 0;
      }
      $n += 2;    #skip the alleles of the other parent
    }
    push @four_from_states, \@from_state4;

    $to_state = 0;
    while ( $to_state < $numberOfStates ) {
      @binary_representation_of_to_inheritance_state =
          binary_representation_of_inheritance_state( $to_state,
                                                   $numberOfMeiosesInPedigree );
      $lowest_Hamming_distance = 10e99;
      foreach my $array1 ( @four_from_states ) {
        $Hamming_distance = Hamming_distance_arrays( $array1,
                              \@binary_representation_of_to_inheritance_state );
        if ( $lowest_Hamming_distance > $Hamming_distance ) {
          $lowest_Hamming_distance = $Hamming_distance;
        }
      }

      $states_are_adjacent{$from_state}{$to_state} = $lowest_Hamming_distance;
      if ( $lowest_Hamming_distance == 1 ) {
        $parental_transition{$from_state}{$to_state} =
            parental_transition_arrays(
                              $numberOfNonfounders,
                              \@binary_representation_of_from_inheritance_state,
                              \@binary_representation_of_to_inheritance_state
            );

        #print STDERR "$from_state\t$to_state\t$lowest_Hamming_distance" .
        #             "\t$parental_transition{$from_state}{$to_state}\n";
      }
      elsif ( $lowest_Hamming_distance == 0 ) {

        #print STDERR "$from_state\t$to_state\t$lowest_Hamming_distance" .
        #             "\tfrom and to are same state\n";
      }
      else {

        #print STDERR "$from_state\t$to_state\t$lowest_Hamming_distance" .
        #             "\tnot allowed\n";
      }

      $to_state++;
    }
    $from_state++;
  }

  $local_timestamp = timestamp();
  print STDERR "[$subroutine] Done at $local_timestamp.\n";
  return ( \%parental_transition, \%states_are_adjacent );
}    #end state_to_state_Hamming_distances

##-------------------------------------------------------------------------##
## Use: my \@state_names = name_the_states( $numberOfStates );
##
##  Returns
##     \@state_names
##
##  Globals: None
##-------------------------------------------------------------------------##
sub name_the_states {
  my $numberOfStates = shift;

  # May need better names
  @state_names = ();
  my $m = 0;
  while ( $m < $number_of_states ) {
    push @state_names, $m;
    $m++;
  }

  return ( \@state_names );
}    #end name_the_states

##-------------------------------------------------------------------------##
## Use: my \%tmpDataCache = load_previously_compiled_data( $dataFile );
##
##  Returns
##     \%tmpDataCache
##
##  Globals: None
##-------------------------------------------------------------------------##
sub load_previously_compiled_data {
  my $dataCacheFile = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  my $local_timestamp = timestamp();
  print STDERR "[$subroutine] Begin at $local_timestamp.\n";

  die "load_previously_compiled_data(): Could not locate "
      . "$dataCacheFile file! Either it's missing or it's empty.\n"
      unless ( -s $dataCacheFile );
  my $age_of_file = sprintf( "%.1f", -M $dataCacheFile );
  print STDERR "\tOpening $dataCacheFile, which is "
      . "$age_of_file days old.\n";

  # Read in the data
  my $tmpDataCache = retrieve( $dataCacheFile );

  return ( $tmpDataCache );
}    #end load_previously_compiled_data

##-------------------------------------------------------------------------##
## Use: my $parentalTransitionStr = parental_transition_arrays(
##                                         $numberOfNonfounders,
##                                         \@array1, \@array2 );
##
##  Returns
##     $parentalTransitionStr
##
##  Globals: None
##-------------------------------------------------------------------------##
sub parental_transition_arrays {
  my $numberOfNonfounders = shift;
  my ( $first, $second ) = @_;

  die unless @$first == @$second;

  my $n;
  my $i;
  my $recombinant_individual;
  my $numberOfMeiosesInPedigree = 2 * $numberOfNonfounders;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  #$local_timestamp = timestamp();
  #print STDERR "[$subroutine] Begin at $local_timestamp.\n";
  #print STDERR join(",",@$first),"\tfirst array in parental naming\n";
  #print STDERR join(",",@$second),"\tsecond array in parental naming\n";
  #encoding 1 (default, with the first two bits of both arrays equal to zero)

  my $paternal = 0;
  undef $recombinant_individual;
  for ( $i = 0 ; $i < @$first ; $i += 2 ) {
    if ( $first->[ $i ] ne $second->[ $i ] ) {
      $paternal++;
      unless ( defined( $recombinant_individual ) ) {
        $recombinant_individual = $i;
      }
    }
  }
  if ( $paternal == 1 ) {
    if ( $numberOfNonfounders == 2 ) {
      return "paternal";
    }
    else {
      return "paternal " . ( $recombinant_individual / 2 );
    }
  }

  my $maternal = 0;
  undef $recombinant_individual;
  for ( $i = 1 ; $i < @$first ; $i += 2 ) {
    if ( $first->[ $i ] ne $second->[ $i ] ) {
      $maternal++;
      unless ( defined( $recombinant_individual ) ) {
        $recombinant_individual = $i;
      }
    }
  }
  if ( $maternal == 1 ) {
    if ( $numberOfNonfounders == 2 ) {
      return "maternal";
    }
    else {
      return "maternal " . ( ( $recombinant_individual - 1 ) / 2 );
    }
  }

  if ( $numberOfNonfounders == 2 ) {
    die;    #should not be possible to get here
  }

  #need to check three alternate encodings
  #encoding 2
  $n = 0;
  while ( $n < $numberOfMeiosesInPedigree )
  { #recall that $numberOfMeiosesInPedigree = the length of the inheritance vector
    if ( $$first[ $n ] == 0 )
    { #zero means the paternal allele of this child is the first allele of the the father
      $$first[ $n ] = 1;
    }
    else {    #if not zero, must be 1
      $$first[ $n ] = 0;
    }
    $n += 2;    #skip the alleles of the other parent
  }

  $paternal = 0;
  undef $recombinant_individual;
  for ( $i = 0 ; $i < @$first ; $i += 2 ) {
    if ( $first->[ $i ] ne $second->[ $i ] ) {
      $paternal++;
      unless ( defined( $recombinant_individual ) ) {
        $recombinant_individual = $i;
      }
    }
  }
  if ( $paternal == 1 ) {
    return "paternal " . ( $recombinant_individual / 2 );
  }
  $maternal = 0;
  undef $recombinant_individual;
  for ( $i = 1 ; $i < @$first ; $i += 2 ) {
    if ( $first->[ $i ] ne $second->[ $i ] ) {
      $maternal++;
      unless ( defined( $recombinant_individual ) ) {
        $recombinant_individual = $i;
      }
    }
  }
  if ( $maternal == 1 ) {
    return "maternal " . ( ( $recombinant_individual - 1 ) / 2 );
  }

  #encoding 3
  $n = 1;
  while ( $n < $numberOfMeiosesInPedigree )
  { #recall that $numberOfMeiosesInPedigree = the length of the inheritance vector
    if ( $$first[ $n ] == 0 )
    { #zero means the paternal allele of this child is the first allele of the the father
      $$first[ $n ] = 1;
    }
    else {    #if not zero, must be 1
      $$first[ $n ] = 0;
    }
    $n += 2;    #skip the alleles of the other parent
  }

  $paternal = 0;
  undef $recombinant_individual;
  for ( $i = 0 ; $i < @$first ; $i += 2 ) {
    if ( $first->[ $i ] ne $second->[ $i ] ) {
      $paternal++;
      unless ( defined( $recombinant_individual ) ) {
        $recombinant_individual = $i;
      }
    }
  }
  if ( $paternal == 1 ) {
    return "paternal " . ( $recombinant_individual / 2 );
  }

  $maternal = 0;
  undef $recombinant_individual;
  for ( $i = 1 ; $i < @$first ; $i += 2 ) {
    if ( $first->[ $i ] ne $second->[ $i ] ) {
      $maternal++;
      unless ( defined( $recombinant_individual ) ) {
        $recombinant_individual = $i;
      }
    }
  }
  if ( $maternal == 1 ) {
    return "maternal " . ( ( $recombinant_individual - 1 ) / 2 );
  }

  #encoding 4
  $n = 0;
  while ( $n < $numberOfMeiosesInPedigree )
  { #recall that $numberOfMeiosesInPedigree = the length of the inheritance vector
    if ( $$first[ $n ] == 0 )
    { #zero means the paternal allele of this child is the first allele of the the father
      $$first[ $n ] = 1;
    }
    else {    #if not zero, must be 1
      $$first[ $n ] = 0;
    }
    $n += 2;    #skip the alleles of the other parent
  }

  $paternal = 0;
  undef $recombinant_individual;
  for ( $i = 0 ; $i < @$first ; $i += 2 ) {
    if ( $first->[ $i ] ne $second->[ $i ] ) {
      $paternal++;
      unless ( defined( $recombinant_individual ) ) {
        $recombinant_individual = $i;
      }
    }
  }
  if ( $paternal == 1 ) {
    return "paternal " . ( $recombinant_individual / 2 );
  }

  $maternal = 0;
  undef $recombinant_individual;
  for ( $i = 1 ; $i < @$first ; $i += 2 ) {
    if ( $first->[ $i ] ne $second->[ $i ] ) {
      $maternal++;
      unless ( defined( $recombinant_individual ) ) {
        $recombinant_individual = $i;
      }
    }
  }
  if ( $maternal == 1 ) {
    return "maternal " . ( ( $recombinant_individual - 1 ) / 2 );
  }
  else {
    die
        ; #should not be possible not to have found which parent the recombination occured in
  }

}    #end parental_transition_arrays

1;
