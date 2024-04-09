#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) statistics_for_HMM_Suttonian_blocks.pl
##  Author:
##      Jared Roach      <jroach@systemsbiology.org>
##      Robert M. Hubley <rhubley@systemsbiology.org>
##  Description:
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
#   modified 11/12/10   v0.02 moved subroutines to modules
#   modified 11/12/10   v0.03 much more sophisticated routines for assessing
#       quality of blocks
#   modified 11/20/10   v0.05 deal with partially called genotype vectrs
#   modified 11/28/10   v0.06 downweight length so that weak centromerees don't
#       get called incorrectly (no actual changes as were gradually made in v5
#       since last version)
#   modified 11/28/10   v0.06 if the top emission has more than one heterozygous
#       founder (out of two) it is automatically weak
#   modified 11/28/10   v0.06 fixed bug that ADDED strength to blocks over a
#       certain strength
#   modified 12/16/10   v0.07 rather than using raw length of blocks, reduce
#       by the number of Ns in the block - makes the algorithm better able to
#       handle centromeres, for example
#   BOBAMA version (note that this is a bastardized version of the code and
#       so versioning is turned off - for current code version and best debugged,
#       see the non-bobama version)
#   modified 1/21/2011  updated strength algorithm; better handling of input
#       reference_genome
#
###############################################################################
# TO DO
#   this to-do items seems really neat and relevant and elegant
#     check for statistical deviation from Hardy-Weinberg equilibrium of
#     non-founder genotypes given observed founder genotypes in block (i.e.,
#     test for the presence of lots of unllinked markers in block) If they are
#     all linked, then don't trust state prediction so much.
#

=head1 NAME

 statistics_for_HMM_Suttonian_blocks.pl - 

=head1 SYNOPSIS

 Usage:
 ./statistics_for_HMM_Suttonian_blocks.pl --reference_genome <hg18|hg19>
                                          --infile_directory <DIR>

=head1 DESCRIPTION 

=head1 OPTIONS

    -version  
        print out the version and exit

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
#use strict;
use Data::Dumper;
use FindBin;
use lib $FindBin::RealBin;
use Storable qw(nstore retrieve);
use Getopt::Long;
use warnings FATAL => 'all'
    ; #like to trap all warnings when running in workflow environment becuase otherwise it is easy to miss the warning
use TimeUtils qw(timestamp);
use RSutils_for_HMM_Statistics;
use InheritanceStates
    qw(load_raw_state_for_each_position load_raw_state_for_each_position_on_a_chromosome load_block_file genotype_vector_for_index);

my $releaseVersion = "0.1.8";
my $gitVersion     = "0.1.8";
my $VERSION        = "release $releaseVersion";
my $program        = "statistics_for_HMM_Suttonian_blocks_bobama.pl";
print STDERR "Version $VERSION of $program, running with timestamp "
    . timestamp() . "\n";

my $reference_genome = "";
GetOptions( "reference_genome=s" => \$reference_genome,
            "infile_directory=s" => \$infile_directory,
    	    "gap_file=s" => \$gap_file);

if ( $reference_genome =~ /(hg18|hg19)/i ) {
  $reference_genome = lc( $1 );
}
else {
  print STDERR "I do not recognize the reference_genome type!\n";
  usage();
  die;
}

unless ( defined( $infile_directory ) ) {
  print STDERR "No infile directory defined!/\n";
  usage();
  die;
}

sub usage {
  print "$0 - $VERSION\n\n";
  exec "pod2text $0";
  exit;
}

my $statsCacheFilename = "$infile_directory/stats_DataCache.dat";

@labels_of_chrs_to_study = ( 1 .. 22 );
@chromosomes_to_study    = map { "chr" . $_ } @labels_of_chrs_to_study;

my $refGenomeResourceDirectory = "$FindBin::RealBin/refGenomes";

print STDERR
    "Directory for finding the raw HMM output files is $infile_directory.\n";

$short_filename = "nibbled_raw_blocks.txt";
$long_filename  = "$infile_directory/$short_filename";

load_previously_compiled_data( "$infile_directory/hmm_DataCache.dat" );

@possible_genotypes = ( "aa", "ab", "an", "bb", "bn", "nn" );
my $integer = 0;
foreach $genotype ( @possible_genotypes ) {
  $genotype_for_integer{$integer} = $genotype;
  $integer++;
}

$number_of_nonfounders =
    3;    #synonymous with number of children for a nuclear family
$number_of_nonfounders =
    2;    #synonymous with number of children for a nuclear family
$number_of_possible_genotypes = scalar @possible_genotypes;

$number_of_individuals_in_pedigree =
    2 + $number_of_nonfounders
    ; #this is the length of the genotype vector (still assumes 2 parents and their children)

$maximum_genotype_vector_index =
    ( $number_of_possible_genotypes**$number_of_individuals_in_pedigree ) - 1;

$compsci_number_of_states = 15;
$compsci_number_of_states = 3;
@state_indices            = ( 0 .. $compsci_number_of_states );

( $emission, $positions_of_bounding_markers ) =
    load_raw_state_for_each_position( $infile_directory );
%emission                      = %$emission;
%positions_of_bounding_markers = %$positions_of_bounding_markers;

%input_file_header_info = ();
( $block, $number_of_markers, $start, $end, $input_file_header_info ) =
    load_block_file( $long_filename );
%block             = %$block;
%number_of_markers = %$number_of_markers;
%start             = %$start;
%end               = %$end;

%input_file_header_info = %$input_file_header_info;

if ( $reference_genome eq "" ) {
  $frz_as_defined_by_input_header = $input_file_header_info{'reference_genome'};
  if ( ( $frz_as_defined_by_input_header =~ /^hg1[89]$/ ) ) {
    $reference_genome = $frz_as_defined_by_input_header;
  }
  else {
    print STDERR
"The reference genome defined by the headeer of the input file $input_file_header_info{'reference_genome'} is not in canonical format (i.e., hg18 or 19)!\n";
    die;
  }
}

$outfile = $infile_directory . "/" . "statistics_for_raw_blocks.txt";

$header_for_output_files = make_header_for_output_files();
open STATISTICS, ">$outfile" or die "Could not open $outfile for writing!\n";
print STATISTICS $header_for_output_files;
print STATISTICS
"#header\tchr\tunsmoothed_block_number\tstate\tnumber of markers\tstart\tend\tlength\n";

initialize_gap_compilation_for_correcting_lengths();

%emission_counts = ();   #probably can delete this line of code, as is redundant

## RMH: This is for jared's emission prob model
#$pce_emission_prob = 0.03;
#my $probability_cutoff = 3.1 * $pce_emission_prob;
## TODO: We have a problem with this when we use the new emission_probabilities...might
##       need to pass along the emission_prob of PCEs as a parameter
## RMH: This is for my emission prob model
#my $probability_cutoff = $emission_prob{1}{3};
my $probability_cutoff = 0.0045;
print "Setting prob cutoff to = $probability_cutoff\n";

foreach $chr ( @chromosomes_to_study ) {
  print ".";

  #print uc $chr, "\n";
  @block_numbers =
      sort { $a <=> $b }
      keys %{ $block{$chr} };
  @positions_on_this_chromosome = sort { $a <=> $b } keys %{ $emission{$chr} };

  foreach $block_number ( @block_numbers ) {
    %emission_counts = ();
    $block_type      = $block{$chr}{$block_number};

    $length = gap_corrected_length( $chr,
                                    $start{$chr}{$block_number},
                                    $end{$chr}{$block_number} );

    foreach $position ( @positions_on_this_chromosome ) {
      next
          unless (     ( $position >= $start{$chr}{$block_number} )
                   and ( $position <= $end{$chr}{$block_number} ) );
      $emission = $emission{$chr}{$position};
      if ( $emission_counts{$emission} ) {
        $emission_counts{$emission}++;
      }
      else {
        $emission_counts{$emission} = 1;
      }
    }
    print STATISTICS
"$chr\t$block_number\t$block_type\t$number_of_markers{$chr}{$block_number}\t"
        . "$start{$chr}{$block_number}\t$end{$chr}{$block_number}\t$length\n";

    @seen_emissions =
        sort { $emission_counts{$b} <=> $emission_counts{$a} }
        keys %emission_counts;

    #need to count fully called and partially called emissions
    %state_support_counts = ();

    #initialize
    # RMH: Ooops...dependency on 4 states!
    foreach $state ( 0, 1, 2, 3 ) {
      $state_support_counts{$block_number}{$state} = 0;
    }

    %state_support_counts_including_partially_called_vectors = ();

    #initialize
    # RMH: Ooops...dependency on 4 states!
    foreach $state ( 0, 1, 2, 3 ) {
      $state_support_counts_including_partially_called_vectors{$block_number}
          {$state} = 0;
    }

    @fully_called_emissions = ();

    # RMH:
    @all_emissions = ();
    foreach $emission ( @seen_emissions ) {

      # RMH: Ooops...hard coded
      if ( $emission == 1295 ) {    #debug for family of four
        print STDERR
"Statistics shows a noninformative genotype vector nnnnnnnn emission. "
            . "Not sure how it got here: $chr\t$block_number\t$block_type\t"
            . "$number_of_markers{$chr}{$block_number}\t$start{$chr}{$block_number}\t"
            . "$end{$chr}{$block_number}\t$length\n";

        die;
      }

      $genotype_vector = genotype_vector_for_index(
                                             $emission,
                                             $maximum_genotype_vector_index,
                                             $number_of_individuals_in_pedigree,
                                             \%genotype_for_integer
      );
	print "$genotype_vector\n";
      if ( $genotype_vector eq "nnnnnnnn" ) {    #debug for family of four
        print STDERR
"Statistics shows a noninformative genotype vector nnnnnnnn emission. "
            . "Not sure how it got here: $chr\t$block_number\t$block_type\t"
            . "$number_of_markers{$chr}{$block_number}\t$start{$chr}{$block_number}\t"
            . "$end{$chr}{$block_number}\t$length\n";
        die;
      }

      unless ( length $genotype_vector == 8 ) {
        print STDERR "The statistics algorithm is currently hard coded four a "
            . "nuclear family of four!\n";
        die;
      }
      @likely_states = ();

      # RMH:
      #print "Pushing $emission ( $genotype_vector )\n";
      push @all_emissions, $emission;

    #be careful of using partially called vectors for statistics, so ignore them
      if ( $genotype_vector =~ /n/ ) {
        foreach $state ( @state_indices ) {
          if ( $emission_prob{$state}{$emission} > $probability_cutoff ) {
            push @likely_states, $state;
          }
        }
        $likely_states = join( ", ", @likely_states );
        print STATISTICS
"\t$genotype_vector\t$emission\t$emission_counts{$emission}\t$likely_states\n";

        foreach $state ( @likely_states ) {
          $state_support_counts_including_partially_called_vectors{
            $block_number }{$state} += $emission_counts{$emission};
        }
      }
      else {
        push @fully_called_emissions, $emission;

        #print STDERR "$emission\n";
        #figure out the states that are likely to emit this emission
        #$emission_prob{$m}{$possible_emission} = $pce_emission_prob;
        foreach $state ( @state_indices ) {
          if ( $emission_prob{$state}{$emission} > $probability_cutoff ) {
            push @likely_states, $state;
          }
        }
        $likely_states = join( ", ", @likely_states );
        print STATISTICS
"\t$genotype_vector\t$emission\t$emission_counts{$emission}\t$likely_states\n";

        foreach $state ( @likely_states ) {
          $state_support_counts{$block_number}{$state} +=
              $emission_counts{$emission};
          $state_support_counts_including_partially_called_vectors{
            $block_number }{$state} += $emission_counts{$emission};
        }
      }

    }

    compute_statistics_for_block_quality( $block_number );
  }
}
print "\n";
dump_values_for_later_use( $statsCacheFilename );

close STATISTICS;

exit;

#END MAIN BLOCK

sub compute_statistics_for_block_quality {
  my $genotype_vector;
  my @seen_states;

  my %dp;
  $dp{0} = getDotProductOfObservedEmissionsVsStateEmissionProbabilities( 0,
                                                              \@all_emissions );
  $dp{1} = getDotProductOfObservedEmissionsVsStateEmissionProbabilities( 1,
                                                              \@all_emissions );
  $dp{2} = getDotProductOfObservedEmissionsVsStateEmissionProbabilities( 2,
                                                              \@all_emissions );
  $dp{3} = getDotProductOfObservedEmissionsVsStateEmissionProbabilities( 3,
                                                              \@all_emissions );

  my @top = sort { $dp{$b} <=> $dp{$a} } keys( %dp );
  $newStrength = $dp{ $top[ 0 ] } / ( $dp{0} + $dp{1} + $dp{2} + $dp{3} );
  print STATISTICS "New Strength = $newStrength\n";

  #make a single hash of statistics to store as a Storable
  $block_statistics{'top_state'}{$chr}{$block_number} = $top[ 0 ];
  $block_statistics{'consistent_top_emissions'}{$chr}{$block_number} = "NA";
  $block_statistics{'number_of_markers'}{$chr}{$block_number}        =
      scalar( @all_emissions );
  $block_statistics{'second_best_state'}{$chr}{$block_number} = $top[ 1 ];
  $block_statistics{'states_excluding_top_two_emissions_are_same_as_top_state'}
      {$chr}{$block_number} = "NA";
  $block_statistics{'ratio_of_support_for_top_two_states'}{$chr}
      {$block_number} = "NA";
  $block_statistics{'consistent_top_emissions'}{$chr}{$block_number} = "NA";
  $block_statistics{'length'}{$chr}{$block_number}                   = "NA";
  $block_statistics{'state_supported_after_deleting_top_emissions'}{$chr}
      {$block_number} = "NA";
  $block_statistics{'strength'}{$chr}{$block_number} = $newStrength;
  $block_statistics{'hmm_state_and_crass_state_prediction_are_the_same'}{$chr}
      {$block_number} = "NA";

  return;
}

#end compute_statistics_for_block_quality

## BLOCK STRENGTH ROUTINE MOVED TO ATTIC

sub dump_values_for_later_use {
  my $dataCacheFilename = shift;
  my %dataCache = (
              'block_statistics'              => \%block_statistics,
              'positions_of_bounding_markers' => \%positions_of_bounding_markers
  );

  nstore \%dataCache, $dataCacheFilename;
  $local_timestamp = timestamp();
  print STDERR "$local_timestamp: Done saving "
      . "values to $dataCacheFilename\n";
}

#end dump_values_for_later_use

sub load_previously_compiled_data {
  my $dataCacheFile = shift;

  my $subroutine      = "load_previously_compiled_data";
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

  # NOTE: This is a very inefficient way to break the hash into
  #       the discrete variables. This will double the memory
  #       requirements doing it this way ( at least temporarily )
  #       TODO: Fix.
  %genotype_for_integer = %{ $tmpDataCache->{'genotype_for_integer'} };
  %emission_prob        = %{ $tmpDataCache->{'emission_prob'} };

  # Free up memory -- no need to keep a duplicate
  undef $tmpDataCache;
}

#end load_previously_compiled_data

sub getDotProductOfObservedEmissionsVsStateEmissionProbabilities {
  my $state         = shift;
  my $emissions_ref = shift;

  my @obsVector = ();
  foreach my $emission ( @{$emissions_ref} ) {
    $obsVector[ $emission ] = $emission_counts{$emission};
  }

  my @theoreticalVector = ();
  foreach my $emission ( keys( %{ $emission_prob{$state} } ) ) {
    $theoreticalVector[ $emission ] = $emission_prob{$state}{$emission};
  }

  for ( my $i = 0 ; $i <= $#theoreticalVector ; $i++ ) {
    $theoreticalVector[ $i ] = 0 if ( !defined $theoreticalVector[ $i ] );
    $obsVector[ $i ]         = 0 if ( !defined $obsVector[ $i ] );
  }

  #print "Dumper obs:\n" . Dumper( \@obsVector ) . "\n";
  #print "Dumper the:\n" . Dumper( \@theoreticalVector ) . "\n";

  my $dp = vectorDotProduct( \@obsVector, \@theoreticalVector );

  return ( $dp );

}

sub vectorDotProduct {
  my $vectorA = shift;
  my $vectorB = shift;

  if ( @$vectorA != @$vectorB ) {
    warn "vectorDotProduct: Vectors are different sizes!\n";
    return;
  }

  my $total = 0;
  for ( my $i = 0 ; $i <= $#{$vectorA} ; $i++ ) {
    $total += $vectorA->[ $i ] * $vectorB->[ $i ];
  }

  return ( $total );
}

sub are_the_top_emissions_the_characteristic_emissions_for_this_state {

  #my %characteristic_emissions;
  my @test_vector;

  #BEGIN {
  @{ $characteristic_emissions{0} } = ( 1,  6, 253, 258 );    #equal weighting
  @{ $characteristic_emissions{1} } = ( 42, 1, 253 );         #42 first
  @{ $characteristic_emissions{2} } = ( 37, 6, 258 );         #37 first
  @{ $characteristic_emissions{3} } = ( 37, 42 );             #equal weighting
                                                              #}

  unless ( $number_of_individuals_in_pedigree == 4 ) {
    return "not computed";
  }
  my $state                                = shift;
  my $emissions_ref                        = shift;
  my @emissions_sorted_high_support_to_low = @$emissions_ref;

  if ( $state == 0 ) {

    if ( scalar @emissions_sorted_high_support_to_low < 4 ) {
      return 0;
    }

    @test_vector =
        sort { $a <=> $b } @emissions_sorted_high_support_to_low[ 0 .. 3 ];

    if ( scalar @test_vector < 4 ) {
      return 0;
    }

    while ( scalar @test_vector ) {
      unless ( shift @test_vector == shift @{ $characteristic_emissions{0} } ) {
        return 0;
      }
    }
    return 1;
  }

  if ( $state == 1 or $state == 2 ) {
    if ( scalar @emissions_sorted_high_support_to_low < 3 ) {
      return 0;
    }

    @test_vector = @emissions_sorted_high_support_to_low[ 0 .. 2 ];

    if ( scalar @test_vector < 3 ) {
      return 0;
    }

    unless (
            shift @test_vector == shift @{ $characteristic_emissions{$state} } )
    {
      return 0;
    }

    @test_vector = sort { $a <=> $b } @test_vector;

    while ( scalar @test_vector ) {
      unless (
            shift @test_vector == shift @{ $characteristic_emissions{$state} } )
      {
        return 0;
      }
    }
    return 1;
  }

  if ( $state == 3 ) {

    if ( scalar @emissions_sorted_high_support_to_low < 2 ) {
      return 0;
    }

    @test_vector =
        sort { $a <=> $b } @emissions_sorted_high_support_to_low[ 0 .. 1 ];

    if ( scalar @test_vector < 2 ) {
      return 0;
    }

    while ( scalar @test_vector ) {
      unless (
            shift @test_vector == shift @{ $characteristic_emissions{$state} } )
      {
        return 0;
      }
    }
    return 1;
  }

  die;

#here is an example of a good state that did not have consistent top emissions because it had a long non-reference haplotype
#3	2	1325	17557056	18980841	1423785
#abaaabaa	37	347	2, 3
#aaababab	258	213	0, 2
#abbbbbab	343	168	2, 3
#abababaa	43	159	1, 2
#aaabaaaa	6	135	0, 2
#ababbbab	331	102	1, 2
#bbaaabab	255	60	0, 1, 2, 3
#bbababab	261	44	0, 2
#aaababaa	42	28	1, 3
#bbabbbbb	765	26	0, 2
#abaaaaaa	1	22	0, 1
#ababaaaa	7	7	0
#bbabbbab	333	6	1, 3
#ababbbaa	115	3	3
#tri_quad	1296	2	0, 1, 2, 3
#abaaabab	253	2	0, 1
#abbbbbbb	775	1	0, 1
#number of markers	1325
#consistent_top_emissions 	0
#number_of_different_supported_states	4	Number worth caring about: 1
#difference_between_support_for_top_two_states	642	; Ratio: 0.488853503184713
#difference_between_counts_for_top_two_emissions	134	; Ratio: 0.613832853025937
#top state	2
#second_best_state	3
#best_state_excluding_top_emission	2
#best_state_excluding_top_two_emissions	2
#states_excluding_top_two_emissions_are_same_as_top_state	1
#state_supported_after_deleting_top_emissions	2
#strength	0.220845580640625

}

#end are_the_top_emissions_the_characteristic_emissions_for_this_state

sub make_header_for_output_files {
  my $header = "";

  #$header .= "#parameter\tmode\t$mode\n";
  $header .= "#parameter\tinput_file\t$short_filename\n";
  $header .=
"#parameter\tinput_file to suttonian_hmm.pl i.e. the input file for the input file\t$input_file_header_info{'input_file'}\n";
  $header .= "#creator\t$program version $VERSION\n";
  $header .= "#timestamp\t" . timestamp() . "\n";
  $header .= "#contact\tJared Roach jroach\@systemsbiology.org\n";
  $header .= "#reference_genome\t$input_file_header_info{'reference_genome'}\n";

#$header .= "#pedigree_name Illumina Family 1 (hard coded in header subroutine so could be wrong if was not changed appropriately)\n";
#$header .= "#pedigree_name Miller syndrome (hard coded in header subroutine so could be wrong if was not changed appropriately)\n";
  $header .= "#pedigree_name\t$input_file_header_info{'pedigree_name'}\n";

#$header .= "#numerical_base zero (hard coded in header subroutine so could be wrong if was not changed appropriately)\n";
  $header .= "#numerical_base\t$input_file_header_info{'numerical_base'}\n";

#$header .= "#project compression (hard coded in header subroutine so could be wrong if was not changed appropriately)\n";
  $header .= "#project\t$input_file_header_info{'project'}\n";

#$header .= "#chromosome_nomenclature conventional (hard coded in header subroutine so could be wrong if was not changed appropriately)\n";
  $header .=
"#chromosome_nomenclature\t$input_file_header_info{'chromosome_nomenclature'}\n";

  #$header .=  "#inheritance_vector_format binary";
  $header .= "#inheritance_state_representation\tarabic_zero_based\n";

#$header .= "#parameter\tcutoff_for_good_quality_block = $cutoff_for_good_quality_block\n";
#$header .= "#note\tMany header fields hard coded in suttonian-hmm header subroutine so could be wrong if was not changed appropriately\n";

  $header .=
"#the following header lines are automatically parsed from the header of the input file\n";
  foreach $key ( sort keys %input_file_header_info ) {

    #if (defined()) {
    $header .= "#$key\t$input_file_header_info{$key}\n";

    #}
  }
  $header .=
"#the previous header lines were automatically parsed from the header of the input file\n";

  return $header;
}

#end make_header_for_output_files

sub initialize_gap_compilation_for_correcting_lengths {

 #$frz_as_defined_by_input_header = $input_file_header_info{'reference_genome'};
 #$frz = "hg18";  #override for IF1
 #$frz = "hg19";
  my $frz = $reference_genome;
  print STDERR "Using $frz freeze for gap correction\n";

#if (($frz_as_defined_by_input_header =~ /hg/) and !($frz_as_defined_by_input_header =~ /$frz/)) {
#print STDERR "Incompatibility between reference in input header and value of refefrnece being used in initialize_gap_compilation_for_correcting_lengths!\n";
#die;
#}

  #use lib "/proj/famgen/bin/gglusman";
  $RS = new RSutils;

  #/Databases/Gustavo_gaps_in_references/hg19_gap.txt
  #/Databases/Gustavo_gaps_in_references/hg18_gap.txt
  #my $file = "$refGenomeResourceDirectory/$frz" . "_gap.txt";
 my $file =$gap_file;
  #open GAPS, "gunzip -c /proj/famgen/resources/annotation/$frz/gap.txt.gz |";
  open GAPS, $file or die "Could not open $file!\n";
  while ( <GAPS> ) {
    chomp;
    ( undef, $chrom, $start, $end ) = split /\t/;
    push @{ $gaps{$chrom} }, "$start\t$end";
  }
  close GAPS;
}

#end initialize_gap_compilation_for_correcting_lengths

sub gap_corrected_length {
  my $chrom = shift;
  my $start = shift;
  my $end   = shift;
  $end++;    #becuase data was zero based
  my $span = 0;
  foreach
      my $range ( @{ $RS->RSsubtraction( [ "$start\t$end" ], $gaps{$chrom} ) } )
  {
    my ( $rangeStart, $rangeEnd ) = split /\t/, $range;
    $span += $rangeEnd - $rangeStart;
  }

#not sure why RSsubtraction not returning defined values in some circumstances, so this is an hack override
#unless (defined($span)) {
#$span = $end-$start;
#print STDERR "RSsubtraction not working properly with $chrom $start $end.\n";
#}

  return $span;
}

#end gap_corrected_length

