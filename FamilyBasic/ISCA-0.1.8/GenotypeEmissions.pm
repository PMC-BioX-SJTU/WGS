#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) GenotypeEmissions.pm
##  Author:
##      Jared Roach      <jroach@systemsbiology.org>
##      Robert M. Hubley <rhubley@systemsbiology.org>
##  Description:
##      Data structures and routines for interpreting genotype
##      patterns and converting them into HMM emissions.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2011
#*
#*  This work is licensed under the Open Source License v3.0.  To view a copy
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
#   fix the code for correctly dealing with pseudoautosomes, X, and Y
#

=head1 NAME

 GenotypeEmissions.pm - 

=head1 SYNOPSIS

 Usage:

=head1 DESCRIPTION

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
package GenotypeEmissions;
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

# Get version from Makefile/GIT
my $releaseVersion = "0.1.8";
my $gitVersion     = "0.1.8";
my $VERSION        = "release $releaseVersion";
$VERSION .= " ( git: $gitVersion )"
    if ( $gitVersion ne "" );

# Module wide debug flag
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Use:  my ( $possible_emissions, $MIE,
##            $lowest_integer_for_equivalent_genotype_vector,
##            $emission_prob ) =
##            initializeTwoAlleleEmissionAlphabetOrig( $number_of_nonfounders,
##                                      $pce_emission_prob,
##                                      $estimated_background_MIE_rate );
##
##  Jared's original method of producing emission probabilities ( or
##  close to it ).
##
##  Supported genotypes are: aa, ab, an, bb, bn, nn
##
##  Generate list of possible emissions.  A genotype vector (after
##  applying symmetries) is encoded as an integer. The vector of all "a"s
##  is zero.  A vector of all "n"s is the maximum emmision integer.
##  These are the possible genotypes:  aa,ab,an,bb,bn,nn.
##
##  Returns
##     \@possible_emissions
##     \%MIE
##     \%lowest_integer_for_equivalent_genotype_vector
##     \%emission_prob
##
##  Globals: None
##-------------------------------------------------------------------------##
sub initializeTwoAlleleEmissionAlphabetOrig {
  my $number_of_nonfounders         = shift;
  my $pce_emission_prob             = shift;
  my $estimated_background_MIE_rate = shift;

  my $subroutine = "initializeTwoAlleleEmissionAlphabet";
  print STDERR "[$subroutine] Begin...\n" if ( $DEBUG );

  # Note currently we only allow two alleles per locus.
  #    ie: aa, ab, an, bb, bn, nn
  #        maybe in the future: a- b- n-
  # So six genotypes in all....
  my $number_of_genotypes = 6;

  #
  # For example...given a nuclear family of four...each possible genotype
  # vector has a unique base-6 integer with the lowest being 0 == aaaaaaaa
  # (base-10 integer = 0) and the highest being 5555 == nnnnnnnn (base-10
  # integer = 1295). So there are 1295 such emissions.  It makes sense
  # to collapse nearly identical emissions into equivalence classes.
  #
  # Any genotype vector with the same genotype content (same count of
  # each genotype) and the same consistencies with the 4 main states
  # is given the emission integer of the lowest such genotype vector.
  # Furthermore, any vector with any Ns (even one) that is consistent
  # with all 4 states is given the value 5555 (integer = 1295) or the
  # equivalent to "nnnnnnnn..nn" in the larger families.
  #
  my %seen;
  print STDERR "[$subroutine] Tabulating consistencies and Mendelian "
      . "probs for genotype vectors:\n"
      if ( $DEBUG );

  my $maximum_genotype_vector_index =
      ( $number_of_genotypes**( $number_of_nonfounders + 2 ) ) - 1;

  # Data structures about to be populated
  my %genotype_vector_consistent_with_state;
  my %genotype_content_index;
  my $index_for_fully_heterozygous_genotype_vector;
  my %lowest_integer_for_equivalent_genotype_vector;
  my %lowest_integer_for_character_and_content;
  my %lowest_integer_for_state_weight_profile;

  #my %number_of_heterozygote_genotypes;
  my %MIE;

  # A more general datastructure to hold the initial weights assigned to
  # states for each genotype pattern.  This will be used to develop the
  # emission probabilities.
  my %genotypeIndexToStateWeights = ();

  my $fully_heterozygous_genotype_vector =
      "ab" x ( $number_of_nonfounders + 2 );

  my $number_of_states = &getNumberOfStates( $number_of_nonfounders );
  my $completely_permissive_state_consistency_vector =
      "1" x ( $number_of_states );
  my $MIE_state_consistency_vector = "0" x ( $number_of_states );

  #
  # Calculate State Consistencies
  #    - Find lowest index for allele content and consistencies
  #
  my $genotype_vector_index = 0;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    if ( $DEBUG && ( $genotype_vector_index % 250 ) == 0 ) {
      print STDERR "\t$genotype_vector_index";
    }

    # Call this for every emission index possible in this sized
    # pedigree.
    my ( $genotype_vector, $genotype_content ) =
        &decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                       $number_of_nonfounders + 2 );

    #
    # The genotype content is a string of numerals indicating how
    # many of each genotype pattern ie "aa", "ab" .. is observed
    # in the vector.  So:  index=0 in a family of 4 will have
    # a vector of "aaaaaaaa" and a content of 400000. NOTE: This
    # will be problematic with families with greater than 9
    # individuals. Not sure who needs this.
    #
    $genotype_content_index{$genotype_vector_index} = $genotype_content;

    # Backwards calculating a global for the index
    # of the "abababab..ab" equiv.
    if ( $genotype_vector eq $fully_heterozygous_genotype_vector ) {
      $index_for_fully_heterozygous_genotype_vector = $genotype_vector_index;
    }

    # Compute consistency class for the genotype vector
    print "DEBUG: " . Dumper( getStateConsistencies( $genotype_vector ) ) . "\n"
        if ( $DEBUG );

    my ( $state_consistency_vector, $compatableStates,
         $probOfGenotypeEmissionGivenStateAndFounderGenotypes )
        = &getStateConsistencies( $genotype_vector );

    $genotype_vector_consistent_with_state{$genotype_vector_index} =
        $compatableStates;

    # 'ab' is the $genotype_content{1} value
    #  for purposes of the CNV state probabilities, we tabulate
    #  the number of heterozygote positions for each genotype vector
    #my $num_of_heterozygote_genotypes = substr( $genotype_content, 1, 1 );
    #$number_of_heterozygote_genotypes{$genotype_vector_index} =
    #    $num_of_heterozygote_genotypes;

    # Calculate the Mendelian probability of this state emitting this
    # vector, given the list of founders.  - JARED's method -
    @{ $genotypeIndexToStateWeights{$genotype_vector_index} } =
        &Mendelian_probs_of_genotype_vector_emissions_given_founder_genotypes(
                                                       $genotype_vector,
                                                       $number_of_nonfounders );

    #
    #  Create a datastructure which contains pointers to equivalent
    #  genotype indexes.  This particular method collapses genotypes
    #  into clusters which are equivalent both in genotype_content
    #  and state_consistencies.
    #
    if ( $seen{$state_consistency_vector}{$genotype_content} ) {
      $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} =
          $lowest_integer_for_character_and_content{$state_consistency_vector}
          {$genotype_content};
    }
    else {
      $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} =
          $genotype_vector_index;
      $lowest_integer_for_character_and_content{$state_consistency_vector}
          {$genotype_content} = $genotype_vector_index;
      $seen{$state_consistency_vector}{$genotype_content} = 1;
    }
    if (
         ( $genotype_vector =~ /n/ )
         and ( $state_consistency_vector eq
               $completely_permissive_state_consistency_vector )
        )
    {
      $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} =
          $maximum_genotype_vector_index;
    }

    # Define a hash of all indexes which are MIEs
    if ( $state_consistency_vector eq $MIE_state_consistency_vector ) {
      $MIE{$genotype_vector_index} = 1;
    }
    else {
      $MIE{$genotype_vector_index} = 0;
    }

    $genotype_vector_index++;
  }    # while ( $genotype_vector_index <= $max...

  if ( $DEBUG ) {
    print STDERR "\t$genotype_vector_index";
    print STDERR "\n";
    my $local_timestamp = timestamp();
    print STDERR "[$subroutine] Done with state compatibilities and "
        . "Mendelian probs of states at $local_timestamp.\n";
  }

  ## PRINT OUT STATE WEIGHTS
  if ( $DEBUG ) {
    print "\nGENOTYPE STATE WEIGHTS\n";
    for ( my $i = 0 ; $i <= $maximum_genotype_vector_index ; $i++ ) {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $i, $number_of_nonfounders + 2 );
      my $possibleEmissionIdx =
          $lowest_integer_for_equivalent_genotype_vector{$genotype_vector};
      my $cv = join( "",
                     @{ $genotype_vector_consistent_with_state{$i} }
                         { 0 .. ( $number_of_states - 1 ) } );
      my ( $equivVector, $foo ) =
          decodeTwoAlleleEmissionIndex( $possibleEmissionIdx,
                                        $number_of_nonfounders + 2 );
      print "$i\t$genotype_vector\t$cv\t$possibleEmissionIdx($equivVector)\t"
          . join( "\t", @{ $genotypeIndexToStateWeights{$i} } ) . "\n";
    }
  }

  #
  # Develop emission probabilties
  #

  # This sets the pseudocount for the probability of emission of all states
  # V43:
  #   $emission_pseudocount = $pce_emission_prob;
  # V44: Increased to this in v44 to allow states to not flip into short
  #      states that happen to have a high SCE content.
  #   $emission_pseudocount = 2 * $pce_emission_prob;
  # V44-Revised: Increased to this in a revised v44 to allow states to not flip
  #              into short states that happen to have a high SCE content
  #              Consider completing downweighting vectors that are specific for
  #              a single state as they may get overweighted
  my $emission_pseudocount = 3 * $pce_emission_prob;

  ## Initialize all possible emission indexes to the emission_pseudocount
  $genotype_vector_index = 0;
  my %possible_emission;
  my %emission_prob;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );
    my $possibleEmissionIdx =
        $lowest_integer_for_equivalent_genotype_vector{$genotype_vector};
    $possible_emission{$possibleEmissionIdx} = 1;
    my $m = 0;
    while ( $m < $number_of_states ) {
      $emission_prob{$m}{$possibleEmissionIdx} = $emission_pseudocount;
      $m++;
    }
    $genotype_vector_index++;
  }

  my @possible_emissions = sort { $a <=> $b } keys %possible_emission;

  # Now either generate theoretical emission probs for each possible character
  # in the emission alphabet, or read in empirical probs for very rough theory,
  # a state has a high prob of emitting all vectors consistent with it, and a
  # low prob of emitting all vectors inconsistent.
  #   Low prob is $probability_of_emitting_a_Suttonian_error
  #   High prob depends on parental genotypes -> {ab,aa} much more frequent
  #   than {ab,ab} -> about 15-20% of all emissions from a state will come
  #   from parent patterns {ab,ab}
  # Also divide the prob of an emission by the number of states it is consistent
  # with (e.g., abababab ~5% of emissions from identical because it also emits
  # from non-identical, whereas ababaaaa is ~10% of identical emissions becuase
  # it only emits from there)
  # Prob of generating N depends on sequencing technology and filters - tricky
  # to estimate theoretically

  # For each state
  # list all compatible emissions and weight theme and assign normalized probs

  # Code that follows -> use the emission probs from the genotype vectors
  # conditioned on the founder genotypes
  # Foreach genotype vector class, use the emission prob per state of whichever
  # genotype vector in that class is highest -> not necessarily the same
  # genotype vector within a class for each state
  $genotype_vector_index = 0;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );
    my $emission =
        $lowest_integer_for_equivalent_genotype_vector{$genotype_vector};
    my $state = 0;
    while ( $state < $number_of_states ) {
      my $emission_prob =
          ${ $genotypeIndexToStateWeights{$genotype_vector_index} }[ $state ];

  # If this prob is greater than the greatest prob seen so far for this class of
  # genotype vectors, then use this prob for this class
      if ( $emission_prob > $emission_prob{$state}{$emission} ) {
        $emission_prob{$state}{$emission} = $emission_prob;
      }
      $state++;
    }
    $genotype_vector_index++;
  }

  # Now overrule the emission probs for all MIEs
  foreach my $emission ( @possible_emissions ) {
    if ( $MIE{$emission} ) {
      my $state = 0;
      while ( $state < $number_of_states ) {
        $emission_prob{$state}{$emission} = $estimated_background_MIE_rate;
        $state++;
      }
    }
  }

  if ( $DEBUG ) {
    for ( my $genotype_vector_index = 0 ;
          $genotype_vector_index <= $maximum_genotype_vector_index ;
          $genotype_vector_index++ )
    {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                        $number_of_nonfounders + 2 );

      print "$genotype_vector_index\t$genotype_vector\t";

      # Only consider emissions
      if ( $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} ==
           $genotype_vector_index )
      {
        for ( my $state = 0 ; $state < $number_of_states ; $state++ ) {
          print "" . $emission_prob{$state}{$genotype_vector_index} . "\t";
        }
      }
      print "\n";
    }
  }

  # Treat tri- and quad- allelic as MIEs
  my $state = 0;
  while ( $state < $number_of_states ) {
    $emission_prob{$state}{ $maximum_genotype_vector_index + 1 } =
        $estimated_background_MIE_rate;
    $state++;
  }
  $lowest_integer_for_equivalent_genotype_vector{'tri_quad'} =
      $maximum_genotype_vector_index + 1;

  #print "FINAL EMISSION PROB: " . Dumper( \%emission_prob ) . "\n";

  if ( $DEBUG ) {
    my $local_timestamp = timestamp();
    print STDERR "[$subroutine] Done at $local_timestamp.\n";
  }

  if ( $DEBUG ) {
    for ( my $genotype_vector_index = 0 ;
          $genotype_vector_index <= $maximum_genotype_vector_index + 1 ;
          $genotype_vector_index++ )
    {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                        $number_of_nonfounders + 2 );

      print "$genotype_vector_index\t$genotype_vector\t"
          . "$lowest_integer_for_equivalent_genotype_vector{$genotype_vector}\t";

      # Only consider emissions
      if ( $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} ==
           $genotype_vector_index )
      {
        for ( my $state = 0 ; $state < $number_of_states ; $state++ ) {
          print "" . $emission_prob{$state}{$genotype_vector_index} . "\t";
        }
      }
      print "\n";
    }
  }

  return ( \@possible_emissions, \%MIE,
           \%lowest_integer_for_equivalent_genotype_vector,
           \%emission_prob );

}    #end initializeTwoAlelleEmissionAlphabetOrig

##-------------------------------------------------------------------------##
## Use:  my ( $possible_emissions, $MIE,
##            $lowest_integer_for_equivalent_genotype_vector,
##            $emission_prob ) =
##            initializeTwoAlleleEmissionAlphabet( $number_of_nonfounders,
##                                      $pce_emission_prob,
##                                      $estimated_background_MIE_rate );
##
##  Supported genotypes are: aa, ab, an, bb, bn, nn
##
##  Generate list of possible emissions.  A genotype vector (after
##  applying symmetries) is encoded as an integer. The vector of all "a"s
##  is zero.  A vector of all "n"s is the maximum emmision integer.
##  These are the possible genotypes:  aa,ab,an,bb,bn,nn.
##
##  NOTES:
##    Currently doesn't use $estimated_background_MIE_rate and instead
##    simply uses the psedocounts for these patterns.
##
##  Returns
##     \@possible_emissions
##     \%MIE
##     \%lowest_integer_for_equivalent_genotype_vector
##     \%emission_prob
##
##  Globals: None
##-------------------------------------------------------------------------##
sub initializeTwoAlleleEmissionAlphabet {
  my $number_of_nonfounders              = shift;
  my $pce_emission_prob                  = shift;
  my $estimated_background_MIE_rate      = shift;
  my $useOriginalWeightCalculationMethod = shift;

  my $subroutine = "initializeTwoAlleleEmissionAlphabet";
  print STDERR "[$subroutine] Begin...\n" if ( $DEBUG );

  # Note currently we only allow two alleles per locus.
  #    ie: aa, ab, an, bb, bn, nn
  #        maybe in the future: a- b- n-
  # So six genotypes in all....
  my $number_of_genotypes = 6;

  # This sets the pseudocount for the probability of emission of all states
  # V43:
  #   $emission_pseudocount = $pce_emission_prob;
  # V44: Increased to this in v44 to allow states to not flip into short
  #      states that happen to have a high SCE content.
  #   $emission_pseudocount = 2 * $pce_emission_prob;
  # V44-Revised: Increased to this in a revised v44 to allow states to not flip
  #              into short states that happen to have a high SCE content
  #              Consider completing downweighting vectors that are specific for
  #              a single state as they may get overweighted
  #my $emission_pseudocount = 3 * $pce_emission_prob;
  #my $emission_pseudocount =  $pce_emission_prob;

  ## As the emission_pseudocount gets small the number of state transitions goes up
  ## and many small states are produced
  ##
  ## This is what I settled on for family of 4 fully_called & exclude_fully_het alphabet
  my $emission_pseudocount = 4 * $pce_emission_prob;

  #
  # For example...given a nuclear family of four...each possible genotype
  # vector has a unique base-6 integer with the lowest being 0 == aaaaaaaa
  # (base-10 integer = 0) and the highest being 5555 == nnnnnnnn (base-10
  # integer = 1295). So there are 1295 such emissions.  It makes sense
  # to collapse nearly identical emissions into equivalence classes.
  #
  # Any genotype vector with the same genotype content (same count of
  # each genotype) and the same consistencies with the 4 main states
  # is given the emission integer of the lowest such genotype vector.
  # Furthermore, any vector with any Ns (even one) that is consistent
  # with all 4 states is given the value 5555 (integer = 1295) or the
  # equivalent to "nnnnnnnn..nn" in the larger families.
  #
  my %seen;
  print STDERR "[$subroutine] Tabulating consistencies and Mendelian "
      . "probs for genotype vectors:\n"
      if ( $DEBUG );

  my $maximum_genotype_vector_index =
      ( $number_of_genotypes**( $number_of_nonfounders + 2 ) ) - 1;

  # Data structures about to be populated
  my %genotype_vector_consistent_with_state;
  my %genotype_content_index;
  my $index_for_fully_heterozygous_genotype_vector;
  my %lowest_integer_for_equivalent_genotype_vector;
  my %lowest_integer_for_character_and_content;
  my %lowest_integer_for_state_weight_profile;
  my %MIE;

  # A more general datastructure to hold the initial weights assigned to
  # states for each genotype pattern.  This will be used to develop the
  # emission probabilities.
  my %genotypeIndexToStateWeights = ();

  my $fully_heterozygous_genotype_vector =
      "ab" x ( $number_of_nonfounders + 2 );

  my $number_of_states = &getNumberOfStates( $number_of_nonfounders );
  my $completely_permissive_state_consistency_vector =
      "1" x ( $number_of_states );
  my $MIE_state_consistency_vector = "0" x ( $number_of_states );

  #
  # Calculate State Consistencies
  #    - Find lowest index for allele content and consistencies
  #
  my $genotype_vector_index = 0;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    if ( $DEBUG && ( $genotype_vector_index % 250 ) == 0 ) {
      print STDERR "\t$genotype_vector_index";
    }

    # Call this for every emission index possible in this sized
    # pedigree.
    my ( $genotype_vector, $genotype_content ) =
        &decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                       $number_of_nonfounders + 2 );

    #
    # The genotype content is a string of numerals indicating how
    # many of each genotype pattern ie "aa", "ab" .. is observed
    # in the vector.  So:  index=0 in a family of 4 will have
    # a vector of "aaaaaaaa" and a content of 400000. NOTE: This
    # will be problematic with families with greater than 9
    # individuals. Not sure who needs this.
    #
    $genotype_content_index{$genotype_vector_index} = $genotype_content;

    # Backwards calculating a global for the index
    # of the "abababab..ab" equiv.
    if ( $genotype_vector eq $fully_heterozygous_genotype_vector ) {
      $index_for_fully_heterozygous_genotype_vector = $genotype_vector_index;
    }

    # Compute consistency class for the genotype vector
    print "DEBUG: " . Dumper( getStateConsistencies( $genotype_vector ) ) . "\n"
        if ( $DEBUG );

    my ( $state_consistency_vector, $compatableStates,
         $probOfGenotypeEmissionGivenStateAndFounderGenotypes )
        = &getStateConsistencies( $genotype_vector );

    $genotype_vector_consistent_with_state{$genotype_vector_index} =
        $compatableStates;

    @{ $genotypeIndexToStateWeights{$genotype_vector_index} } =
        &getStateWeightsForTwoAlleleGenotypePattern( $genotype_vector,
                                      \%genotype_vector_consistent_with_state );

    # Handle MIEs
    if ( $state_consistency_vector eq $MIE_state_consistency_vector ) {
      $MIE{$genotype_vector_index} = 1;
    }
    else {
      $MIE{$genotype_vector_index} = 0;
    }

    $genotype_vector_index++;
  }    # while ( $genotype_vector_index <= $max...

  ## PRINT OUT STATE WEIGHTS
  if ( $DEBUG ) {
    print "\nGENOTYPE STATE WEIGHTS ( ALL )\n";
    for ( my $i = 0 ; $i <= $maximum_genotype_vector_index ; $i++ ) {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $i, $number_of_nonfounders + 2 );
      my $cv = join( "",
                     @{ $genotype_vector_consistent_with_state{$i} }
                         { 0 .. ( $number_of_states - 1 ) } );
      print "$i\t$genotype_vector\t$cv\t"
          . join( "\t", @{ $genotypeIndexToStateWeights{$i} } ) . "\n";
    }
  }

  ## Develop final emission alphabet
  ##   Cluster emissions by their state weights.
  $genotype_vector_index = 0;
  my %possible_emission;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );

    #
    #  Create a datastructure which contains pointers to equivalent
    #  genotype indexes.  This method collapses genotypes
    #  into clusters which share the same state weights.
    #
    my $weightStr = join(
                          ",",
                          @{
                            $genotypeIndexToStateWeights{$genotype_vector_index}
                              }
    );
    print "$genotype_vector_index $weightStr\n" if ( $DEBUG );
    if ( exists $lowest_integer_for_state_weight_profile{$weightStr} ) {
      $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} =
          $lowest_integer_for_state_weight_profile{$weightStr};
    }
    else {
      $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} =
          $genotype_vector_index;
      $lowest_integer_for_state_weight_profile{$weightStr} =
          $genotype_vector_index;
      $possible_emission{$genotype_vector_index} = 1;
    }
    $genotype_vector_index++;
  }
  my @possible_emissions = sort { $a <=> $b } keys %possible_emission;

  ## PRINT OUT STATE WEIGHTS
  if ( $DEBUG ) {
    print "\nGENOTYPE STATE WEIGHTS ( ALL POSSIBLE )\n";
    for ( my $i = 0 ; $i <= $maximum_genotype_vector_index ; $i++ ) {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $i, $number_of_nonfounders + 2 );
      my $possibleEmissionIdx =
          $lowest_integer_for_equivalent_genotype_vector{$genotype_vector};
      my $cv = join( "",
                     @{ $genotype_vector_consistent_with_state{$i} }
                         { 0 .. ( $number_of_states - 1 ) } );
      my ( $equivVector, $foo ) =
          decodeTwoAlleleEmissionIndex( $possibleEmissionIdx,
                                        $number_of_nonfounders + 2 );
      print "$i\t$genotype_vector\t$cv\t$possibleEmissionIdx($equivVector)\t"
          . join( "\t", @{ $genotypeIndexToStateWeights{$i} } ) . "\n";
    }
  }

  if ( $DEBUG ) {
    print STDERR "\t$genotype_vector_index";
    print STDERR "\n";
    my $local_timestamp = timestamp();
    print STDERR "[$subroutine] Done with state compatibilities and "
        . "Mendelian probs of states at $local_timestamp.\n";
  }

  ##
  ## Develop emission probabilties
  ##

  # Total weights+pseudocounts for each state restricted
  # to possible emissions.
  $genotype_vector_index = 0;
  my %emission_prob     = ();
  my @totalStateWeights = ();
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );
    my $possibleEmissionIdx =
        $lowest_integer_for_equivalent_genotype_vector{$genotype_vector};

    # Canonical Emission
    if ( $possibleEmissionIdx == $genotype_vector_index ) {
      my $state = 0;
      while ( $state < $number_of_states ) {
        $totalStateWeights[ $state ] +=
            $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] +
            $emission_pseudocount;
        $state++;
      }
    }
    $genotype_vector_index++;
  }

  if ( $DEBUG ) {
    print "totalStateWeight[ 0 ] = $totalStateWeights[0]\n";
    print "totalStateWeight[ 1 ] = $totalStateWeights[1]\n";
    print "totalStateWeight[ 2 ] = $totalStateWeights[2]\n";
    print "totalStateWeight[ 3 ] = $totalStateWeights[3]\n";
    print "Total Possible Emissions = " . scalar( @possible_emissions ) . "\n";
  }

  # Now either generate theoretical emission probs for each possible character
  # in the emission alphabet, or read in empirical probs for very rough theory,
  # a state has a high prob of emitting all vectors consistent with it, and a
  # low prob of emitting all vectors inconsistent.
  #   Low prob is $probability_of_emitting_a_Suttonian_error
  #   High prob depends on parental genotypes -> {ab,aa} much more frequent
  #   than {ab,ab} -> about 15-20% of all emissions from a state will come
  #   from parent patterns {ab,ab}
  # Also divide the prob of an emission by the number of states it is consistent
  # with (e.g., abababab ~5% of emissions from identical because it also emits
  # from non-identical, whereas ababaaaa is ~10% of identical emissions becuase
  # it only emits from there)
  # Prob of generating N depends on sequencing technology and filters - tricky
  # to estimate theoretically

  # For each state
  # list all compatible emissions and weight theme and assign normalized probs

  # Code that follows -> use the emission probs from the genotype vectors
  # conditioned on the founder genotypes
  # Foreach genotype vector class, use the emission prob per state of whichever
  # genotype vector in that class is highest -> not necessarily the same
  # genotype vector within a class for each state

  # Normalize across state
  $genotype_vector_index = 0;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {

    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );

    # Only consider emissions
    if ( $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} ==
         $genotype_vector_index )
    {
      my $state = 0;
      while ( $state < $number_of_states ) {
        $emission_prob{$state}{$genotype_vector_index} =
            ( $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] +
              $emission_pseudocount ) / $totalStateWeights[ $state ];
        $state++;
      }
    }
    $genotype_vector_index++;
  }

  if ( $DEBUG ) {
    for ( my $genotype_vector_index = 0 ;
          $genotype_vector_index <= $maximum_genotype_vector_index ;
          $genotype_vector_index++ )
    {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                        $number_of_nonfounders + 2 );

      print "$genotype_vector_index\t$genotype_vector\t";

      # Only consider emissions
      if ( $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} ==
           $genotype_vector_index )
      {
        for ( my $state = 0 ; $state < $number_of_states ; $state++ ) {
          print "" . $emission_prob{$state}{$genotype_vector_index} . "\t";
        }
      }
      print "\n";
    }
  }

  # Treat tri- and quad- allelic as MIEs
  my $state = 0;
  while ( $state < $number_of_states ) {
    $emission_prob{$state}{ $maximum_genotype_vector_index + 1 } =
        $estimated_background_MIE_rate;
    $state++;
  }
  $lowest_integer_for_equivalent_genotype_vector{'tri_quad'} =
      $maximum_genotype_vector_index + 1;

  if ( $DEBUG ) {
    my $local_timestamp = timestamp();
    print STDERR "[$subroutine] Done at $local_timestamp.\n";
  }

  if ( $DEBUG ) {
    for ( my $genotype_vector_index = 0 ;
          $genotype_vector_index <= $maximum_genotype_vector_index + 1 ;
          $genotype_vector_index++ )
    {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                        $number_of_nonfounders + 2 );

      if ( $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} !=
           $genotype_vector_index )
      {
        print "$genotype_vector_index\t$genotype_vector\t"
            . "$lowest_integer_for_equivalent_genotype_vector{$genotype_vector}\n";
        next;
      }

      print "$genotype_vector_index\t$genotype_vector\t"
          . "$lowest_integer_for_equivalent_genotype_vector{$genotype_vector}\t";

      # Only consider emissions
      if ( $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} ==
           $genotype_vector_index )
      {
        for ( my $state = 0 ; $state < $number_of_states ; $state++ ) {
          print "" . $emission_prob{$state}{$genotype_vector_index} . "\t";
        }
      }
      print "\n";
    }
  }

  return ( \@possible_emissions, \%MIE,
           \%lowest_integer_for_equivalent_genotype_vector,
           \%emission_prob );

}    #end initializeTwoAlelleEmissionAlphabet

##-------------------------------------------------------------------------##
## Use:  my ( $possible_emissions, $MIE,
##            $lowest_integer_for_equivalent_genotype_vector,
##            $emission_prob ) =
##            initializeTwoAlleleEmissionAlphabetExp( $number_of_nonfounders,
##                                      $pce_emission_prob,
##                                      $estimated_background_MIE_rate );
##
##  ****Experimental version****
##
##  Supported genotypes are: aa, ab, an, bb, bn, nn
##
##  Generate list of possible emissions.  A genotype vector (after
##  applying symmetries) is encoded as an integer. The vector of all "a"s
##  is zero.  A vector of all "n"s is the maximum emmision integer.
##  These are the possible genotypes:  aa,ab,an,bb,bn,nn.
##
##  Returns
##     \@possible_emissions
##     \%MIE
##     \%lowest_integer_for_equivalent_genotype_vector
##     \%emission_prob
##
##  Globals: None
##-------------------------------------------------------------------------##
sub initializeTwoAlleleEmissionAlphabetExp {
  my $number_of_nonfounders              = shift;
  my $pce_emission_prob                  = shift;
  my $estimated_background_MIE_rate      = shift;
  my $useOriginalWeightCalculationMethod = shift;

  my $subroutine = "initializeTwoAlleleEmissionAlphabet";
  print STDERR "[$subroutine] Begin...\n" if ( $DEBUG );

  # Note currently we only allow two alleles per locus.
  #    ie: aa, ab, an, bb, bn, nn
  #        maybe in the future: a- b- n-
  # So six genotypes in all....
  my $number_of_genotypes = 6;

  # This sets the pseudocount for the probability of emission of all states
  # V43:
  #   $emission_pseudocount = $pce_emission_prob;
  # V44: Increased to this in v44 to allow states to not flip into short
  #      states that happen to have a high SCE content.
  #   $emission_pseudocount = 2 * $pce_emission_prob;
  # V44-Revised: Increased to this in a revised v44 to allow states to not flip
  #              into short states that happen to have a high SCE content
  #              Consider completing downweighting vectors that are specific for
  #              a single state as they may get overweighted
  #my $emission_pseudocount = 3 * $pce_emission_prob;
  #my $emission_pseudocount =  $pce_emission_prob;

  ## As the emission_pseudocount gets small the number of state transitions goes up
  ## and many small states are produced
  ##
  ## This is what I settled on for family of 4 fully_called & exclude_fully_het alphabet
  #my $emission_pseudocount = 4 * $pce_emission_prob;
  my $emission_pseudocount = 0.001;

  #
  # For example...given a nuclear family of four...each possible genotype
  # vector has a unique base-6 integer with the lowest being 0 == aaaaaaaa
  # (base-10 integer = 0) and the highest being 5555 == nnnnnnnn (base-10
  # integer = 1295). So there are 1295 such emissions.  It makes sense
  # to collapse nearly identical emissions into equivalence classes.
  #
  # Any genotype vector with the same genotype content (same count of
  # each genotype) and the same consistencies with the 4 main states
  # is given the emission integer of the lowest such genotype vector.
  # Furthermore, any vector with any Ns (even one) that is consistent
  # with all 4 states is given the value 5555 (integer = 1295) or the
  # equivalent to "nnnnnnnn..nn" in the larger families.
  #
  my %seen;
  print STDERR "[$subroutine] Tabulating consistencies and Mendelian "
      . "probs for genotype vectors:\n"
      if ( $DEBUG );

  my $maximum_genotype_vector_index =
      ( $number_of_genotypes**( $number_of_nonfounders + 2 ) ) - 1;

  # Data structures about to be populated
  my %genotype_vector_consistent_with_state;
  my %genotype_content_index;
  my $index_for_fully_heterozygous_genotype_vector;
  my %lowest_integer_for_equivalent_genotype_vector;
  my %lowest_integer_for_character_and_content;
  my %lowest_integer_for_state_weight_profile;
  my %MIE;

  # A more general datastructure to hold the initial weights assigned to
  # states for each genotype pattern.  This will be used to develop the
  # emission probabilities.
  my %genotypeIndexToStateWeights = ();

  my $fully_heterozygous_genotype_vector =
      "ab" x ( $number_of_nonfounders + 2 );

  my $number_of_states = &getNumberOfStates( $number_of_nonfounders );
  my $completely_permissive_state_consistency_vector =
      "1" x ( $number_of_states );
  my $MIE_state_consistency_vector = "0" x ( $number_of_states );

  #
  # Calculate State Consistencies
  #    - Find lowest index for allele content and consistencies
  #
  my $genotype_vector_index = 0;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    if ( $DEBUG && ( $genotype_vector_index % 250 ) == 0 ) {
      print STDERR "\t$genotype_vector_index";
    }

    # Call this for every emission index possible in this sized
    # pedigree.
    my ( $genotype_vector, $genotype_content ) =
        &decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                       $number_of_nonfounders + 2 );

    #
    # The genotype content is a string of numerals indicating how
    # many of each genotype pattern ie "aa", "ab" .. is observed
    # in the vector.  So:  index=0 in a family of 4 will have
    # a vector of "aaaaaaaa" and a content of 400000. NOTE: This
    # will be problematic with families with greater than 9
    # individuals. Not sure who needs this.
    #
    $genotype_content_index{$genotype_vector_index} = $genotype_content;

    # Backwards calculating a global for the index
    # of the "abababab..ab" equiv.
    if ( $genotype_vector eq $fully_heterozygous_genotype_vector ) {
      $index_for_fully_heterozygous_genotype_vector = $genotype_vector_index;
    }

    # Compute consistency class for the genotype vector
    print "DEBUG: " . Dumper( getStateConsistencies( $genotype_vector ) ) . "\n"
        if ( $DEBUG );

    my ( $state_consistency_vector, $compatableStates,
         $probOfGenotypeEmissionGivenStateAndFounderGenotypes )
        = &getStateConsistencies( $genotype_vector );

    $genotype_vector_consistent_with_state{$genotype_vector_index} =
        $compatableStates;

    @{ $genotypeIndexToStateWeights{$genotype_vector_index} } =
        &getStateWeightsForTwoAlleleGenotypePattern( $genotype_vector,
                                      \%genotype_vector_consistent_with_state );

    # Handle MIEs
    if ( $state_consistency_vector eq $MIE_state_consistency_vector ) {
      $MIE{$genotype_vector_index} = 1;

      # Override the weights for MIEs
      #my $state = 0;
      #while ( $state < $number_of_states )
      #{
      #  $genotypeIndexToStateWeights{$genotype_vector_index}->[$state] =
      #       $estimated_background_MIE_rate - $emission_pseudocount;
      #  $state++;
      #}

    }
    else {
      $MIE{$genotype_vector_index} = 0;
    }

    $genotype_vector_index++;
  }    # while ( $genotype_vector_index <= $max...

  ## PRINT OUT STATE WEIGHTS
  if ( $DEBUG ) {
    print "\nGENOTYPE STATE WEIGHTS ( ALL )\n";
    for ( my $i = 0 ; $i <= $maximum_genotype_vector_index ; $i++ ) {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $i, $number_of_nonfounders + 2 );
      my $cv = join( "",
                     @{ $genotype_vector_consistent_with_state{$i} }
                         { 0 .. ( $number_of_states - 1 ) } );
      print "$i\t$genotype_vector\t$cv\t"
          . join( "\t", @{ $genotypeIndexToStateWeights{$i} } ) . "\n";
    }
  }

  $genotype_vector_index = 0;
  my @cntOfEmission              = ();
  my $totalCntOfEmission         = 0;
  my $totalFilteredCntOfEmission = 0;
  open IN, "</09TB_2/famgen/git/isca/pOfEmission.txt"
      or die "Could not open pOfEmission.txt";
  while ( <IN> ) {

    if ( /([abn]+)\s+(\d+)\s+([\d\.e-]+)/ ) {
      my $genotype_vector = $1;
      my $count           = $2;
      my $index           = encodeTwoAlleleEmissionPat( $genotype_vector );
      if (    $genotype_vector =~ /^[ab][ab][ab][ab][ab][ab][ab][ab]$/
           && $genotype_vector !~ /^abababab$/ )
      {
        $cntOfEmission[ $index ] = $count;
        $totalFilteredCntOfEmission += $count;
      }
      $totalCntOfEmission += $count;
    }
  }
  close IN;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {

    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );

    if (    $genotype_vector =~ /^[ab][ab][ab][ab][ab][ab][ab][ab]$/
         && $genotype_vector !~ /^abababab$/ )
    {
      my $probOfEmission =
          ( $cntOfEmission[ $genotype_vector_index ] /
            $totalFilteredCntOfEmission );
      print "** $genotype_vector => $probOfEmission [ "
          . "$cntOfEmission[$genotype_vector_index] "
          . "/ $totalFilteredCntOfEmission ] "
          if ( $DEBUG );
      foreach my $state ( 0, 1, 2, 3 ) {

        # Add pseudocount & renormalize
        $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] =
            ( $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] +
              $emission_pseudocount ) / ( 1 + ( $emission_pseudocount * 4 ) );
        print
            "$genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ]("
            if ( $DEBUG );

        # Calculate inverse of conditional probablities
        $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] =
            $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] *
            ( $probOfEmission / 0.25 );
        print
"$genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ])\t"
            if ( $DEBUG );
      }
      print "\n";
    }
    else {
      foreach my $state ( 0, 1, 2, 3 ) {
        $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] = 0;
      }
    }
    $genotype_vector_index++;
  }
  ## END TESTING

  ## Develop final alphabet - The possible emissions are dependent on
  ## the filters used by the program.  This section applies some typical
  ## filters to the data and then determines the alphabet by clustering
  ## similar weight patterns together.
  $genotype_vector_index = 0;
  my %possible_emission;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );

    #if ( 1 )
    #if ( $genotype_vector =~ /^[ab][ab]nn[ab][ab][ab][ab]$/ &&
    #$genotype_vector !~ /^abnnabab$/ )
    #if ( $genotype_vector =~ /^nnnn[ab][ab][ab][ab]$/ )
    #if ( $genotype_vector =~ /^[ab][ab][ab][ab][ab][ab][ab][ab]$/ )
    if (    $genotype_vector =~ /^[ab][ab][ab][ab][ab][ab][ab][ab]$/
         && $genotype_vector !~ /^abababab$/ )
    {

      #
      #  Create a datastructure which contains pointers to equivalent
      #  genotype indexes.  This method collapses genotypes
      #  into clusters which share the same state weights.
      #
      my $weightStr = join(
                            ",",
                            @{
                              $genotypeIndexToStateWeights{
                                $genotype_vector_index}
                                }
      );

      #print "$genotype_vector_index $weightStr\n";
      if ( exists $lowest_integer_for_state_weight_profile{$weightStr} ) {
        $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} =
            $lowest_integer_for_state_weight_profile{$weightStr};
      }
      else {
        $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} =
            $genotype_vector_index;
        $lowest_integer_for_state_weight_profile{$weightStr} =
            $genotype_vector_index;
        $possible_emission{$genotype_vector_index} = 1;
      }
    }
    else {
      $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} = 1296;
    }
    $genotype_vector_index++;
  }
  my @possible_emissions = sort { $a <=> $b } keys %possible_emission;

  ## PRINT OUT STATE WEIGHTS
  if ( $DEBUG ) {
    print "\nGENOTYPE STATE WEIGHTS ( ALL POSSIBLE )\n";
    for ( my $i = 0 ; $i <= $maximum_genotype_vector_index ; $i++ ) {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $i, $number_of_nonfounders + 2 );
      my $possibleEmissionIdx =
          $lowest_integer_for_equivalent_genotype_vector{$genotype_vector};
      my $cv = join( "",
                     @{ $genotype_vector_consistent_with_state{$i} }
                         { 0 .. ( $number_of_states - 1 ) } );
      my ( $equivVector, $foo ) =
          decodeTwoAlleleEmissionIndex( $possibleEmissionIdx,
                                        $number_of_nonfounders + 2 );
      print "$i\t$genotype_vector\t$cv\t$possibleEmissionIdx($equivVector)\t"
          . join( "\t", @{ $genotypeIndexToStateWeights{$i} } ) . "\n";
    }
  }

  if ( $DEBUG ) {
    print STDERR "\t$genotype_vector_index";
    print STDERR "\n";
    my $local_timestamp = timestamp();
    print STDERR "[$subroutine] Done with state compatibilities and "
        . "Mendelian probs of states at $local_timestamp.\n";
  }

  #
  # Develop emission probabilties
  #

  ## Initialize all possible emission indexes to the emission_pseudocount
  $genotype_vector_index = 0;
  my %emission_prob;
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );
    my $possibleEmissionIdx =
        $lowest_integer_for_equivalent_genotype_vector{$genotype_vector};
    my $m = 0;
    while ( $m < $number_of_states ) {
      $emission_prob{$m}{$possibleEmissionIdx} =
          $genotypeIndexToStateWeights{$possibleEmissionIdx}->[ $m ];
      $m++;
    }
    $genotype_vector_index++;
  }

  # Now either generate theoretical emission probs for each possible character
  # in the emission alphabet, or read in empirical probs for very rough theory,
  # a state has a high prob of emitting all vectors consistent with it, and a
  # low prob of emitting all vectors inconsistent.
  #   Low prob is $probability_of_emitting_a_Suttonian_error
  #   High prob depends on parental genotypes -> {ab,aa} much more frequent
  #   than {ab,ab} -> about 15-20% of all emissions from a state will come
  #   from parent patterns {ab,ab}
  # Also divide the prob of an emission by the number of states it is consistent
  # with (e.g., abababab ~5% of emissions from identical because it also emits
  # from non-identical, whereas ababaaaa is ~10% of identical emissions becuase
  # it only emits from there)
  # Prob of generating N depends on sequencing technology and filters - tricky
  # to estimate theoretically

  # For each state
  # list all compatible emissions and weight theme and assign normalized probs

  # Code that follows -> use the emission probs from the genotype vectors
  # conditioned on the founder genotypes
  # Foreach genotype vector class, use the emission prob per state of whichever
  # genotype vector in that class is highest -> not necessarily the same
  # genotype vector within a class for each state

  $genotype_vector_index = 0;
  my @totalStateWeights = ();
  while ( $genotype_vector_index <= $maximum_genotype_vector_index ) {
    my ( $genotype_vector, $genotype_content ) =
        decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                      $number_of_nonfounders + 2 );

    # Only consider emissions
    if ( $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} ==
         $genotype_vector_index )
    {
      my $state = 0;
      while ( $state < $number_of_states ) {
        if (
          $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] > 0 )
        {
          $totalStateWeights[ $state ] +=
              $genotypeIndexToStateWeights{$genotype_vector_index}->[ $state ] +
              $emission_pseudocount;
        }
        else {
          $totalStateWeights[ $state ] += $emission_pseudocount;
        }
        $state++;
      }
    }
    $genotype_vector_index++;
  }

  if ( $DEBUG ) {
    print "totalStateWeight[ 0 ] = $totalStateWeights[0]\n";
    print "totalStateWeight[ 1 ] = $totalStateWeights[1]\n";
    print "totalStateWeight[ 2 ] = $totalStateWeights[2]\n";
    print "totalStateWeight[ 3 ] = $totalStateWeights[3]\n";
    print "Total Possible Emissions = " . scalar( @possible_emissions ) . "\n";
  }

  # Treat tri- and quad- allelic as MIEs
  my $state = 0;
  while ( $state < $number_of_states ) {
    $emission_prob{$state}{ $maximum_genotype_vector_index + 1 } =
        $estimated_background_MIE_rate;
    $state++;
  }
  $lowest_integer_for_equivalent_genotype_vector{'tri_quad'} =
      $maximum_genotype_vector_index + 1;

  #print "FINAL EMISSION PROB: " . Dumper( \%emission_prob ) . "\n";

  if ( $DEBUG ) {
    my $local_timestamp = timestamp();
    print STDERR "[$subroutine] Done at $local_timestamp.\n";
  }

  if ( $DEBUG ) {
    print "Total Possible Emissions = " . scalar( @possible_emissions ) . "\n";
    for ( my $genotype_vector_index = 0 ;
          $genotype_vector_index <= $maximum_genotype_vector_index + 1 ;
          $genotype_vector_index++ )
    {
      my ( $genotype_vector, $genotype_content ) =
          decodeTwoAlleleEmissionIndex( $genotype_vector_index,
                                        $number_of_nonfounders + 2 );

      next
          if (
             $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} !=
             $genotype_vector_index );

      print
"$genotype_vector_index\t$genotype_vector\t$lowest_integer_for_equivalent_genotype_vector{$genotype_vector}\t";

      # Only consider emissions
      if ( $lowest_integer_for_equivalent_genotype_vector{$genotype_vector} ==
           $genotype_vector_index )
      {
        for ( my $state = 0 ; $state < $number_of_states ; $state++ ) {
          print "" . $emission_prob{$state}{$genotype_vector_index} . "\t";
        }
      }
      print "\n";
    }
  }

  return ( \@possible_emissions, \%MIE,
           \%lowest_integer_for_equivalent_genotype_vector,
           \%emission_prob );

}    #end initializeTwoAlelleEmissionAlphabetExp

##-------------------------------------------------------------------------##
## Use: my %vectorToState = getStateConsistencies( $genotype_vector,
##                                                 $minPloidy,
##                                                 $nonfounderGenderVector );
##
##   Given a genotype vector ie. "abbbabbb" determine which inheritance
##   states it is compatible with.  Return a hash of the states with
##   compatible states set to 1.
##
##      $genotype_vector: The genotype pattern for this family
##      $minPloidy:  Optional. Either "1" or "2"
##      $nonfounderGenderVector: Required if $minPloidy is set to "1".
##                               A bitstring indicating the gender of
##                               the non-founders where "0" = male, and
##                               "1" = female.
##
##  Returns
##      %compatibleStates{state} = 1;
##     The state consistency vector given the $genotype_vector.
##     This is a bit array where each bit corresponds to a
##     particular inheritance state.  A one indicates the pattern
##     is compatable with the state and a zero if not.  For a
##     family of four and example might be: 0110.  Indicating
##     that the genotype pattern ( vector ) is compatible with
##     states 1 & 2 but not 0 and 3.
##
##  Globals: None
##-------------------------------------------------------------------------##
sub getStateConsistencies {
  my $pattern                = shift;
  my $minPloidy              = shift;
  my $nonfounderGenderVector = shift;

  my %seen;
  my $subroutine = "getStateConsistencies";

  #
  # Check parameters
  #
  if ( defined $minPloidy ) {
    if ( $minPloidy == 1 ) {
      unless ( defined $nonfounderGenderVector
        && $nonfounderGenderVector =~ /^[01]+$/
        && length( $nonfounderGenderVector ) == ( length( $pattern ) - 4 ) / 2 )
      {
        die "$subroutine(): Error! Value for nonfounderGenderVector is "
            . "strange: $nonfounderGenderVector\n";
      }
    }
    elsif ( $minPloidy != 2 ) {
      die "$subroutine(): Error! Value for minPloidy is incorrect: "
          . "$minPloidy\n";
    }
  }
  else {

    # Default
    $minPloidy = 2;
  }

  print STDERR "[$subroutine] Begin analyzing $pattern...\n" if ( $DEBUG );

  # Create useful bitmasks if we are considering
  # the X chromosome in this analysis
  my $highBitMask               = 0;
  my $lowBitMask                = 0;
  my $maleChildrenSingleBitMask = 0;
  my $femaleChildrenTwoBitMask  = 0;
  if ( defined $nonfounderGenderVector
       && $nonfounderGenderVector ne "" )
  {
    foreach my $child ( split( //, $nonfounderGenderVector ) ) {

      # A mask containing "1" if the individual is a male.
      $maleChildrenSingleBitMask = $maleChildrenSingleBitMask << 1;
      $maleChildrenSingleBitMask |= 1 if ( $child eq "0" );

      # A mask where each individual's two bit indicator has the
      # high order bit set.
      $highBitMask = $highBitMask << 2;
      $highBitMask |= 2;

      # A mask where each individual's two bit indicator has the
      # low order bit set.
      $lowBitMask = $lowBitMask << 2;
      $lowBitMask |= 1;

      # A two bit ( per individual ) mask where "1" indicates the
      # individual is a female.
      $femaleChildrenTwoBitMask = $femaleChildrenTwoBitMask << 2;
      $femaleChildrenTwoBitMask |= 2 if ( $child eq "1" );
    }
  }

  my @genotype_vector       = split( "", $pattern );
  my $number_of_nonfounders = ( length( $pattern ) / 2 ) - 2;
  my $number_of_states      = &getNumberOfStates( $number_of_nonfounders );

  # By default, set all the consistencies to zero (we will flip them to 1
  # if we determine them to be consistent)
  my %genotype_vector_consistent_with_state;
  for ( my $i = 0 ; $i < $number_of_states ; $i++ ) {
    $genotype_vector_consistent_with_state{$i} =
        0;    # 1 is consistent; 0 is inconsistent
  }

  # For each inheritance state, we need to examine all assignments of
  # phase to all non-founders to see if the resulting phased genotypes
  # match one of the inheritance vectors in that state
  #  ie In a family of 4 there are 2 non-founders.  Each non-founder
  #  can inherit the first allele from it's father ( 0 = first allele )
  #  or the second allele ( 1 = second allele ) and likewise for her
  #  mother.  A two bit pattern can represent this information ie. "00"
  #  for first allele father, and first allele mother.  To represent
  #  all possible ways 2 children can inherit their alleles
  #  one just needs to enumerate all patterns of 4 bits ( 16 ).
  #
  #       00 00 = identical = 0
  #       01 01 = identical
  #       10 10 = identical
  #       11 11 = identical
  #
  #       00 01 = haplo-pat = 1
  #       01 00 = haplo-pat
  #       10 11 = haplo-pat
  #       11 10 = haplo-pat
  #
  #       00 10 = haplo-mat = 2
  #       10 00 = haplo-mat
  #       11 01 = haplo-mat
  #
  #       00 11 = non-ident = 3
  #       01 10 = non-ident
  #       10 01 = non-ident
  #       11 00 = non-ident
  #

  # Computationally, the founder swaps and the non-founder swaps are
  # equivalent so the total number of patterns investigated is 2**(number
  # of indiviudals in the pedigree)

  # Is there at least one parental pattern for which the genotype
  # vector matches an inheritance vector
  my %seen_parental_patterns = ();
  my %seen_phasing;

  #
  # Explore all orderings of parent alleles ( nuclear family = 4 ):
  #  ie. given parental alleles abcd:
  #   1 =  ab  cd
  #   2 =  ab  dc
  #   3 =  ba  cd
  #   4 =  ba  dc
  #
  my @parentalPatternIndices = ( 1 .. 4 );

  #
  #  ChrX is special.
  #
  #  First some conventions.  Genotype patterns for NonPAR chrX in males will
  #  always include an "n" to denote the non-existant allele.  Genotype
  #  patterns should be fed to us with sorted alleles so "n" will always
  #  appear last if it exists.  Ie. genotype of "b" on the single male X
  #  chromosome should be handed to us as "bn".
  #
  #  If we are in the non Pseudo Autosomal Regions ( nonPAR ) the
  #  minPloidy will be haploid ( $minPloidy == 1 )
  #  and the Father's allele patterns are placed in a fixed ordering
  #  of paternal-maternal.  I.e if the parent's genotype pattern is
  #  "anab" then we only allow the following two parental pattern swaps:
  #      na ab
  #      na ba
  #  Where "n" in the Father is the X allele he got from his Father and
  #  "a" is the allele he got from his Mother.  Obviously he didn't get
  #  a chrX allele from his father so this is assumed to be a non-inherited
  #  allele position.
  if ( defined $minPloidy && $minPloidy == 1 ) {
    @parentalPatternIndices = ( 3 .. 4 );
  }

  foreach my $parental_pattern_index ( @parentalPatternIndices ) {

    # Initialize the list_of_equivalent alleles from the parental
    # allele reshuffle.
    # parental_alleles = "abaa"
    #   $equivalentAlleleList[0] = "a"
    #   $equivalentAlleleList[1] = "b"
    #   ...
    my @equivalentAlleleList =
        &getParentalPattern( $parental_pattern_index, \@genotype_vector );

    # Basically the string version of the above array
    my $parental_pattern = join( "", @equivalentAlleleList );

    print "[$subroutine] Parental Pattern = $parental_pattern\n"
        if ( $DEBUG );

    # If already seen this pattern, skip it
    next if ( $seen_parental_patterns{$parental_pattern} );
    $seen_parental_patterns{$parental_pattern} = 1;

    %seen_phasing = ();

    # Is there a phasing of the children's genotypes such that the genotype
    # vector matches the inheritance vector implied by the choice of
    # parental allele arrangement?
    #
    # Because many phasings are equivalent if you instead swap the
    # parent alleles ( something we can do because the phase of each
    # locus isn't known ).
    #
    #    ie.
    #              00 00  Intrinsic phase
    #              00 01  Intrinsic phase
    #              00 10  Intrinsic phase
    #              00 11  Intrinsic phase
    #              01 00  same as 00 01 but swap mother's alleles
    #              01 01  same as 00 00 but swap mother's alleles
    #              01 10  same as 00 11 but swap mother's alleles
    #              01 11  same as 00 10 but swap mother's alleles
    #              10 00  same as 00 10 but swap father's alleles
    #              ...
    #
    my $numberOfPhasings = 2**$number_of_nonfounders;
    for ( my $phasing = 0 ; $phasing < $numberOfPhasings ; $phasing++ ) {

      # chrX - non PAR
      #   If we are considering chrX we need to handle things a little
      #   differently.  The convention for display of alleles in a haploid
      #   region is [abn]n.  Where the first allele is the single call
      #   followed by a constant "n" for the non-existant second allele.
      #
      #   The inheritance vector convention this code is following requires
      #   the ordering in haploid regions to be maternal/paternal.  For
      #   example the "X" inherited by the father from his mother would be
      #   the second allele and the "X" inherted by the father from his
      #   father ( ie. non-existent ) would be the first allele.  So
      #   if the father's pattern at a given position is "an" we need
      #   to flip it always to "na".
      next
          unless (    ( not defined $minPloidy )
                   || $minPloidy == 2
                   || ( $phasing & $maleChildrenSingleBitMask ) ==
                   $maleChildrenSingleBitMask );

      # Convert phase to binary and store in an array
      my @genotype_vector_phasing_indicators =
          &binary_representation_of_genotype_phasing( $phasing,
                                                      $number_of_nonfounders );

      # Using the binary pattern create a phased genotype vector...only
      # phasing the children...the parents are flipped elsewhere
      # NOTE: This is hardcoded for two non-founders
      my $n                      = 0;
      my @phased_genotype_vector = @genotype_vector;
      while ( $n < $number_of_nonfounders ) {

        # if "1" then flip alleles
        if ( $genotype_vector_phasing_indicators[ $n ] ) {
          @phased_genotype_vector[ 2 * $n + 4, 2 * $n + 5 ] =
              @phased_genotype_vector[ 2 * $n + 5, 2 * $n + 4 ];
        }
        $n++;
      }

      # Don't reprocess the same patterns
      my $phased_genotype = join( "", @phased_genotype_vector );

      next if ( $seen_phasing{$phased_genotype} );
      $seen_phasing{$phased_genotype} = 1;

      # Inheritance States
      #   This bit pattern defines how alleles are inherited by the
      #   children.  For a nuclear family of four, allowing for parent
      #   allele swapping and child allele swapping we have 4 distinct
      #   states:
      #
      #      00 00
      #      00 01
      #      00 10
      #      00 11
      #
      #    Where each children's two bit pattern represents inheritance
      #    of parternal alleles and maternal alleles respectively.  Also
      #    0 indicates they received the parents first allele and 1
      #    indicates they received the second.
      #
      #    ChrX - non PAR
      #
      #    This region is handeled differently.  As this is a haploid
      #    region in males we have a reduction of states available.
      #    ie.
      #       m/m  00 00 = 0
      #            00 01 = 1
      #       m/f  00 10 = 2
      #            00 11 = 3
      #       f/m  10 00 = 8
      #            10 01 = 9
      #       f/f  10 10 = 10
      #            10 11 = 11
      foreach my $inheritance_state ( 0 .. ( $number_of_states - 1 ) ) {
        my @appendedEquivAlleles = @equivalentAlleleList;

        # Don't need to calculate if was already consistent in another
        # phasing or parental pattern
        next
            if ( $genotype_vector_consistent_with_state{$inheritance_state} );

        if ( defined $minPloidy && $minPloidy == 1 ) {
          if ( $nonfounderGenderVector =~ /^1.*/ ) {
            $inheritance_state |= 1 << ( ( $number_of_nonfounders * 2 ) - 1 );
          }

          next
              if ( ( $inheritance_state & $highBitMask ) !=
                   $femaleChildrenTwoBitMask );
        }

        # Note that this only works with (Using sprintf) Perl 5.6+
        # If using an earlier version , need to use another recipe for
        # converting to binary
        my ( @binary_representation_of_inheritance_state ) =
            &binary_representation_of_inheritance_state( $inheritance_state,
                                                       $number_of_nonfounders );

        # There needs to be allelic equality between all the equivalent
        # alleles in the inheritance vector with the parental alleles and
        # with each other so get a list of all the descendant alleles from
        # each founder allele (including the founder allele itself), and
        # then compute equality on the whole list initialize with the
        # founder allele

        if ( $DEBUG ) {
          print STDERR "Parental Vector:\t$parental_pattern\n";
          print STDERR "Genotype vector:\t", join( ",", @genotype_vector ),
              "\n";
          print STDERR "Phased Genotype vector:\t",
              join( ",", @phased_genotype_vector ), "\n";
          print STDERR "Inheritance vector( $inheritance_state ):\t",
              join( ",", @binary_representation_of_inheritance_state ), "\n";
        }

        # Now add in the descendant alleles
        # Append the paternal allele of each child
        # NOTE: The $phased_genotype_vector still contains the
        #       original first for alleles of the parents. Skip these.
        $n = 0;
        while ( $n < ( 2 * $number_of_nonfounders ) ) {
          if ( $binary_representation_of_inheritance_state[ $n ] == 0 ) {

            # Zero means the paternal allele of this child is the first
            # allele of the the father
            $appendedEquivAlleles[ 0 ] .= $phased_genotype_vector[ $n + 4 ];
          }
          else {
            $appendedEquivAlleles[ 1 ] .= $phased_genotype_vector[ $n + 4 ];
          }

          # skip the alleles of the other parent
          $n += 2;
        }

        # Append the maternal allele of each child
        $n = 1;
        while ( $n < ( $number_of_nonfounders * 2 ) ) {
          if ( $binary_representation_of_inheritance_state[ $n ] == 0 ) {
            $appendedEquivAlleles[ 2 ] .= $phased_genotype_vector[ $n + 4 ];
          }
          else {
            $appendedEquivAlleles[ 3 ] .= $phased_genotype_vector[ $n + 4 ];
          }
          $n += 2;
        }

        if ( $DEBUG ) {
          print STDERR "Equiv allele list: "
              . join( ", ", @appendedEquivAlleles ) . "\n";
        }

        # Now we have 4 strings one for each parental allele.
        # The first character of the string *is* the parental allele
        # and is followed by the phased children's matched up alleles.
        # If each string contains a single pattern character ( or Ns...
        # here encoded as "0"s ) type then the phasing is compatible
        # with this state.
        my $equivalent = 1;
        foreach my $founder_allele ( 0 .. 3 ) {
          my $alleles = $appendedEquivAlleles[ $founder_allele ];

          # First remove all N's
          $alleles =~ s/n//g;

          # Make sure the remaining characters are homogenous
          unless (    $alleles eq ""
                   || $alleles =~ /^(.)(\1)*$/ )
          {
            $equivalent = 0;
          }
        }
        if ( $DEBUG ) {
          print STDERR "Equivalent status: $equivalent\n";
        }

        # If they are all equivalent, then the genotype vector
        # is compatible with the state
        if ( $equivalent ) {
          $genotype_vector_consistent_with_state{$inheritance_state} = 1;
        }
      }    # inheritance state
    }    # phasings
  }    # parental pattern

  #
  # Create a string vector of state consistencies
  #
  my $state_consistency_vector = "";
  for ( my $n = 0 ; $n < $number_of_states ; $n++ ) {
    if ( defined $genotype_vector_consistent_with_state{$n}
         && $genotype_vector_consistent_with_state{$n} == 1 )
    {
      $state_consistency_vector .= "1";
    }
    else {
      $state_consistency_vector .= "0";
    }
  }

  return ( $state_consistency_vector, \%genotype_vector_consistent_with_state );
}    # end getStateConsistencies

##-------------------------------------------------------------------------##
## Use:  my getStateWeightsForTwoAlleleGenotypePattern( $genotypePattern,
##                             $genotypeIndexToStateConsistenciesHashRef );
##
## This routine uses a method similar to Gustavo's Kwanza state weighting
## method to calculate state weights for possibly degenerate genotype
## patterns.
##
##-------------------------------------------------------------------------##
sub getStateWeightsForTwoAlleleGenotypePattern {
  my $genotypePattern                          = shift;
  my $genotypeIndexToStateConsistenciesHashRef = shift;

  # Flag to indicate that one pattern completion matches: aa aa aa (aa)n
  # which is reference homozygous ( an extremely common uninformative
  # pattern ).
  my $nonInformativeCompletion = 0;

  # Number of completions which have at least one state consistency
  my $numConsistentCompletions = 0;

  my @alpha          = ( "a", "b" );
  my @consis         = ();
  my %seen           = ();
  my $numberOfStates = 0;

  # Number of degenerate bases in the pattern ( ie. Ns )
  my ( $numNs ) = ( $genotypePattern =~ tr/n/n/ );

  if ( $numNs > 0 ) {

    # Loop over all the possible completions of a given degenerate
    # genotype pattern.
    for ( my $i = 0 ; $i < 2**$numNs ; $i++ ) {

      # Fill in the N's given the bitpattern representation
      # of $i
      my $nIdx             = 0;
      my $completeGenotype = "";
      for ( my $j = 0 ; $j < length( $genotypePattern ) ; $j++ ) {
        my $allele = substr( $genotypePattern, $j, 1 );
        if ( $allele eq "n" ) {
          if ( $i & ( 1 << $nIdx ) ) {
            $completeGenotype .= $alpha[ 1 ];
          }
          else {
            $completeGenotype .= $alpha[ 0 ];
          }
          $nIdx++;
        }
        else {
          $completeGenotype .= $allele;
        }
      }

      # Normalize the allele order in the pattern. I.e "a" should come before
      # "b" in a genotype for an individual
      for ( my $j = 0 ; $j < length( $genotypePattern ) ; $j += 2 ) {
        substr( $completeGenotype, $j, 2 ) = join( "",
                          sort { $a cmp $b }
                              split( //, substr( $completeGenotype, $j, 2 ) ) );
      }

      # Frequency adjustment
      #   The non-informative sites "aa aa aa (aa)n" are very frequent
      #   in the dataset.  Using these patterns or ones that can be
      #   completed to them would bias the weights quite a bit.
      $nonInformativeCompletion = 1
          if ( $completeGenotype =~ /^[a]+$/ || $completeGenotype =~ /^[b]+$/ );

      #
      # Alternative which includes less common non-informative completions
      # aabbabab or bbaaabab..
      #$nonInformativeCompletion = 1 if ( $completeGenotype =~ /^[a]+$/ ||
      #                           $completeGenotype =~ /^[b]+$/ ||
      #                           $completeGenotype =~ /^(aabb|bbaa)(ab)+/ );

      # Could have been seen as the simplistic completion routine above
      # distinguishes between completing "nn" as "ab" and "ba".  The
      # sorting routine then corrects the second ordering and creates
      # a duplicate.
      next if ( $seen{$completeGenotype} );
      $seen{$completeGenotype} = 1;

      my $genotype_vector_index =
          &encodeTwoAlleleEmissionPat( $completeGenotype );
      my $maxState = keys(
           %{
             $genotypeIndexToStateConsistenciesHashRef->{$genotype_vector_index}
               }
      ) - 1;
      $numberOfStates = $maxState + 1;
      my $state_consistency_vector = join(
                                     "",
                                     @{
                                       $genotypeIndexToStateConsistenciesHashRef
                                           ->{$genotype_vector_index}
                                         }{ 0 .. $maxState }
      );

      # Only consider completion if it's consistent with at least one state
      if ( $state_consistency_vector =~ /[1]/ ) {
        push @consis, $state_consistency_vector;
        $numConsistentCompletions++;
      }
    }
  }
  else {

    #
    # Not a degenerate genotype pattern
    #

    # Frequency adjustment
    #   The non-informative sites "aa aa aa (aa)n" are very frequent
    #   in the dataset.  Using these patterns or ones that can be
    #   completed to them would bias the weights quite a bit.
    $nonInformativeCompletion = 1
        if ( $genotypePattern =~ /^[a]+$/ || $genotypePattern =~ /^[b]+$/ );

    #
    # Alternative which includes less common non-informative completions
    # aabbabab or bbaaabab..
    #$nonInformativeCompletion = 1 if ( $completeGenotype =~ /^[a]+$/ ||
    #                           $completeGenotype =~ /^[b]+$/ ||
    #                           $completeGenotype =~ /^(aabb|bbaa)(ab)+/ );

    my $genotype_vector_index = &encodeTwoAlleleEmissionPat( $genotypePattern );
    my $maxState = keys(
           %{
             $genotypeIndexToStateConsistenciesHashRef->{$genotype_vector_index}
               }
    ) - 1;
    $numberOfStates = $maxState + 1;
    my $state_consistency_vector = join(
           "",
           @{
             $genotypeIndexToStateConsistenciesHashRef->{$genotype_vector_index}
               }{ 0 .. $maxState }
    );

    # Only consider pattern if it's consistent with at least one state
    if ( $state_consistency_vector =~ /[1]/ ) {
      push @consis, $state_consistency_vector;
      $numConsistentCompletions++;
    }
  }

  my @stateWeights = ();
  if ( $nonInformativeCompletion ) {
    my $equalWeight = ( 1 / $numberOfStates );
    for ( my $i = 0 ; $i < $numberOfStates ; $i++ ) {
      push @stateWeights, $equalWeight;
    }
    return @stateWeights;
  }
  @stateWeights = split( //, "0" x $numberOfStates );

  # ab nn aa aa
  #   ab aa aa aa ( ident/hap-pat )  1/2 occurrances ( assuming eq distributed )
  #   ab ab aa aa ( ident ) 1/2 occurrances
  #
  #    ident = (1/2 / 2) + 1/2 = 3/4
  #    hap-pat = 1/2 / 2 = 1/4
  #
  foreach my $consist ( @consis ) {
    my ( $numOnes ) = ( $consist =~ tr/1/1/ );
    for ( my $i = 0 ; $i < length( $consist ) ; $i++ ) {
      $stateWeights[ $i ] += ( 1 / ( $numConsistentCompletions * $numOnes ) )
          if ( substr( $consist, $i, 1 ) eq "1" );
    }
  }

  return ( @stateWeights );

}    # end getStateWeightsForTwoAlleleGenotypePattern

##-------------------------------------------------------------------------##
## Use:  my
##
##  Given a genotpe_vector...say "aaabaaaa"
##
##  How many of the possible ways I can arrange the parental alleles
##  will be consistent with state N?
##
##   aa 00 ab aa
##
##     0000   a0a0
##     0001
##
##-------------------------------------------------------------------------##
sub Mendelian_probs_of_genotype_vector_emissions_given_founder_genotypes {
  my $genotype_vector       = shift;
  my $number_of_nonfounders = shift;

  my $number_of_states = getNumberOfStates( $number_of_nonfounders );
  my $number_of_meioses_in_pedigree = ( 2 * $number_of_nonfounders );

  # Will do this foreach Mendelian state given the state and the founder
  # genotypes, probability of observing the non-founder genotypes
  # but note that I have already collapsed the genotype vectors into
  # sets of genotype vectors so we have to either average probs across
  # the class, or use the probs of the modal genotype, or some other
  # approach. It may be simplest to just use the modal genotype or
  # the "class prototype". There is a bit of an issue when there are
  # n's -> by assuming the most permissive, I may be biasing the
  # probs a bit it might be best to assume that all n's are "A"s for
  # purpose of prob calling or to mix the probs of letting n = a and
  # letting n =b based on expected background MAF allele frequency
  # by letting n's get called as the value producing the highest
  # posterior probabilities, I may be letting the rate of partial calls
  # in a region to influence whether it is called as a particular state

  $genotype_vector =~ tr/n/0/;
  my $probability;
  my @genotype_vector = split( "", $genotype_vector );
  my $n               = 0;

  my @state_consistency_counts   = ( 0, 0, 0, 0 ); # should be at least 4 states
  my @state_inconsistency_counts = ( 0, 0, 0, 0 );
  my @probs                      = ( 0, 0, 0, 0 );
  for ( $n = 5 ; $n <= $number_of_states ; $n++ ) {
    push @state_consistency_counts,   0;
    push @state_inconsistency_counts, 0;
    push @probs,                      0;
  }

  # Calculate all Mendelian possibilities given the state
  # then prob is the fraction of the possibilities for which this
  # pattern is compatible (may be slightly off for some cases of partial calls)
  # Since we are holding the state fixed, each possibility is derived from
  # combinatorial switches of the founder alleles

  # Fixed ordering for the 4 alleles of both parents.
  # These cycle through the four different inheritance vectors that
  # together form a state by flipping the parental bits.
  my @parental_pattern_indices = ( 1, 2, 3, 4 );

  my @theoretical_nonfounders_genotype_vector;
  my @nonfounders_genotype_vector = @genotype_vector;
  shift @nonfounders_genotype_vector
      ;    #hardcoded to assume the four founders alleles are at the front
  shift @nonfounders_genotype_vector
      ;    #hardcoded to assume the four founders alleles are at the front
  shift @nonfounders_genotype_vector
      ;    #hardcoded to assume the four founders alleles are at the front
  shift @nonfounders_genotype_vector
      ;    #hardcoded to assume the four founders alleles are at the front
   #should be sorted, but maybe not (perhaps I have not thought through ways in which it might not be sorted or be sorted differently)
  $n = 0;

  while ( $n < $number_of_meioses_in_pedigree ) {

    #sort the alleles of each individual
    @nonfounders_genotype_vector[ $n, $n + 1 ] =
        sort @nonfounders_genotype_vector[ $n, $n + 1 ];
    $n += 2;
  }

  foreach my $parental_pattern_index ( @parental_pattern_indices ) {
    my @reshuffledParentalPattern =
        &getParentalPattern( $parental_pattern_index, \@genotype_vector );

    foreach my $inheritance_state ( 0 .. ( $number_of_states - 1 ) ) {

      # Make the predicted childrens' genotype vector, given the
      # inheritance state
      my @binary_representation_of_inheritance_state =
          binary_representation_of_inheritance_state( $inheritance_state,
                                                      $number_of_nonfounders );

      $n                                       = 0;
      @theoretical_nonfounders_genotype_vector = ();

      # Recall that $number_of_meioses_in_pedigree = the length of
      # the inheritance vector
      while ( $n < $number_of_meioses_in_pedigree ) {

        # If this child received the allele in question, push it on the stack
        if ( $binary_representation_of_inheritance_state[ $n ] == 0 ) {

          # Zero means the paternal allele of this child is the first
          # allele of the the father
          push @theoretical_nonfounders_genotype_vector,
              $reshuffledParentalPattern[ 0 ];
        }
        else {    #if not zero, must be 1
          push @theoretical_nonfounders_genotype_vector,
              $reshuffledParentalPattern[ 1 ];
        }
        $n++;

        # Alleles of the other parent (mother, as currently hard-coded)
        if ( $binary_representation_of_inheritance_state[ $n ] == 0 ) {

          # Zero means the paternal allele of this child is the first
          # allele of the the mother
          push @theoretical_nonfounders_genotype_vector,
              $reshuffledParentalPattern[ 2 ];
        }
        else {    #if not zero, must be 1
          push @theoretical_nonfounders_genotype_vector,
              $reshuffledParentalPattern[ 3 ];
        }
        $n++;
      }

      # Now put the resulting children's genotype vector into canonical
      # unphased format (i.e., each genotype alphabetized) to see if it
      # matches the vector we are testing
      $n = 0;

      while ( $n < $number_of_meioses_in_pedigree ) {

        # Sort the alleles of each individual
        @theoretical_nonfounders_genotype_vector[ $n, $n + 1 ] =
            sort @theoretical_nonfounders_genotype_vector[ $n, $n + 1 ];

        $n += 2;
      }

      # Compare the theoretical and observed
      my $are_equal =
          compare_arrays_permitting_nocalls( \@nonfounders_genotype_vector,
                                    \@theoretical_nonfounders_genotype_vector );

      if ( $are_equal ) { $state_consistency_counts[ $inheritance_state ]++ }
      else { $state_inconsistency_counts[ $inheritance_state ]++ }
    }
  }

  #  Probabilities for each state are the number of consistent over
  #  (inconsistent+consistent) results
  $n =
      scalar @state_consistency_counts
      ; #this should be the number of states - I am not sure why I wrote it this way
  my $m = 0;
  while ( $m < $n ) {
    $probs[ $m ] = $state_consistency_counts[ $m ] /
        ( $state_consistency_counts[ $m ] + $state_inconsistency_counts[ $m ] );
    $m++;
  }
  return @probs;
}    #end Mendelian_probs_of_genotype_vector_emissions_given_founder_genotypes

##-------------------------------------------------------------------------##
## Use:  my
## Here nocalls are coded as "0"
##-------------------------------------------------------------------------##
sub compare_arrays_permitting_nocalls {
  my ( $first, $second ) = @_;
  return 0 unless @$first == @$second;
  for ( my $i = 0 ; $i < @$first ; $i++ ) {

    # If one of the alleles being compared is a no-call,
    # consider the arrays equal at that position
    next if (    $first->[ $i ] eq '0'
              || $second->[ $i ] eq '0' );

    return 0 if ( $first->[ $i ] ne $second->[ $i ] );
  }
  return 1;
}    # end compare_arrays_permitting_nocalls

##-------------------------------------------------------------------------##
## Use:  my
##-------------------------------------------------------------------------##
sub encodeTwoAlleleEmissionPat {
  my $genotypePattern = shift;

  my %alleles = (
                  "aa" => 0,
                  "ab" => 1,
                  "an" => 2,
                  "bb" => 3,
                  "bn" => 4,
                  "nn" => 5
  );

  my @g_array = $genotypePattern =~ /../sg;
  my $result  = 0;
  my $idx     = 0;
  foreach my $pat ( @g_array ) {
    $result += $alleles{$pat} * ( 6**$idx++ );
  }
  return ( $result );
}    # encodeTwoAlleleEmissionPat()

##-------------------------------------------------------------------------##
## Use:  my
### WARNING IN BOTH MODULES ( different name )
##-------------------------------------------------------------------------##
sub decodeTwoAlleleEmissionIndex {
  my $index           = shift;
  my $num_individuals = shift;

  my %alleles = (
                  0 => "aa",
                  1 => "ab",
                  2 => "an",
                  3 => "bb",
                  4 => "bn",
                  5 => "nn"
  );

  my $remainder             = $index;
  my $genotype_content      = "";
  my $genotype_vector       = "";
  my %genotype_content_hash = ();
  while ( $num_individuals-- ) {

    # Convert from base 10 to base 6
    my $baseSixDigit = $remainder % 6;
    $remainder = int( $remainder / 6 );

    # Build allele string ie. 0 => "aa"
    $genotype_vector .= $alleles{$baseSixDigit};

    # Tabulate the genotype content ( ie. how many "aa"s are
    # there in the genotype_vector ).
    $genotype_content_hash{$baseSixDigit}++;
  }

  for ( my $i = 0 ; $i < 6 ; $i++ ) {
    $genotype_content .= ( $genotype_content_hash{$i} || "0" );
  }

  if ( $remainder ) {
    return ( "tri_quad" );
  }
  else {
    return ( $genotype_vector, $genotype_content );
  }
}    # decodeTwoAlleleEmissionIndex()

##-------------------------------------------------------------------------##
## Use:  my
### CURRENTLY IN BOTH MODULES!
##-------------------------------------------------------------------------##
sub binary_representation_of_inheritance_state {
  my $decimal               = shift;
  my $number_of_nonfounders = shift;

  die "binary_representation_of_inheritance_state(): Error "
      . "missing number_of_nonfounders!"
      if ( !defined $number_of_nonfounders );

  # Perl 5.6+ dependency
  my $binary = sprintf( "%b", $decimal );

  # Now need to pad with leading zeros
  while ( length $binary < ( 2 * $number_of_nonfounders ) ) {
    $binary = "0" . $binary;
  }
  my @bArray = split( //, $binary );
  return @bArray;
}    #end binary_representation_of_inheritance_state

##-------------------------------------------------------------------------##
## Use:  my
##-------------------------------------------------------------------------##
sub binary_representation_of_genotype_phasing {
  my $decimal               = shift;
  my $number_of_nonfounders = shift;

  # Perl 5.6+ dependency
  my $binary = sprintf( "%b", $decimal );

  # Now need to pad with leading zeros
  while ( length $binary < $number_of_nonfounders ) {
    $binary = "0" . $binary;
  }
  my @bArray = split( //, $binary );
  return @bArray;
}    #end binary_representation_of_genotype_phasing

##-------------------------------------------------------------------------##
## Use: my @parentalPatternArray = getParentalPattern(
##                                        $patternIdx,
##                                        \@genotypeVectorArray );
##
##   Given a genotype vector array ie. ( "a", "b", "a", "a", ... )
##   return an array representing a paticular shuffle of the parental
##   alleles ( first 4 alleles in the @genotypeVectorArray ).  The
##   patternIdx's represent the following shuffle patterns:
##
##              1: 0, 1, 2, 3
##              2: 0, 1, 3, 2
##              3: 1, 0, 2, 3
##              4: 1, 0, 3, 2
##
##   Returns
##
##   The shuffled parental genotype pattern as a four position array.
##
##   Globals: None
##-------------------------------------------------------------------------##
sub getParentalPattern {
  my $patternIdx          = shift;
  my $genotypeVectorArray = shift;

  # Given an genotype vector as an array and that the first
  # 4 positions are labled "a", "b", "c", and "d", return:
  # 1: abcd
  if ( $patternIdx eq 1 ) {
    return (
             (
               $genotypeVectorArray->[ 0 ],
               $genotypeVectorArray->[ 1 ],
               $genotypeVectorArray->[ 2 ],
               $genotypeVectorArray->[ 3 ]
             )
    );
  }

  # 2: abdc
  elsif ( $patternIdx eq 2 ) {
    return (
             (
               $genotypeVectorArray->[ 0 ],
               $genotypeVectorArray->[ 1 ],
               $genotypeVectorArray->[ 3 ],
               $genotypeVectorArray->[ 2 ]
             )
    );
  }

  # 3: bacd
  elsif ( $patternIdx eq 3 ) {
    return (
             (
               $genotypeVectorArray->[ 1 ],
               $genotypeVectorArray->[ 0 ],
               $genotypeVectorArray->[ 2 ],
               $genotypeVectorArray->[ 3 ]
             )
    );
  }

  # 4: badc
  elsif ( $patternIdx eq 4 ) {
    return (
             (
               $genotypeVectorArray->[ 1 ],
               $genotypeVectorArray->[ 0 ],
               $genotypeVectorArray->[ 3 ],
               $genotypeVectorArray->[ 2 ]
             )
    );
  }

  # Undefined
  else {
    die "getParentalPattern( $patternIdx ): Error! "
        . "patternIdx out of bounds!\n";
  }

}    # end getParentalPattern

##-------------------------------------------------------------------------##
## Use:  my $numberOfStates = getNumberOfStates( $numberOfNonfounders );
##
##   For a nuclear family of four there will be 4 inheritance states.
##
##   This subroutine that will need to be fleshed out for pedigrees
##   other than nuclear families.
##
##   X&Y Chromosomes:
##    PAR Regions:
##      Paternal Allele Order: No special ordering...same as autosomes
##      Maternal Allele Order: No special ordering...same as autosomes
##
##    nonPAR Region Chromsome X:
##      Paternal Allele Order: -X  Ordering:
##                                    0 = no allele transmission due
##                                        to haploid state.
##                                    1 = allele transmission.
##      Maternal Allele Order: XX - No special ordering...same as autosomes.
##
##    nonPAR Region Chromosome Y:  Not relevant to ISCA, inherited directly
##
##
##
##    P/M Sex P/M Sex X-State     PAR-state   Y-State
##                                Auto-state  ( unused )
##    -----------------------------------------------
##    0 0 m   0 0 m   Identical   Identical   Identical
##    0 0 m   0 1 m   Non-Ident   Haplo_Pat   Identical
##    0 0 m   1 1 f   Non-Ident   Non-Ident   Null
##    0 0 m   1 0 f   Haplo_Mat   Haplo_Mat   Null
##    1 0 f   0 1 m   Non-Ident   Non-Ident   Null
##    1 0 f   0 0 m   Identical   Haplo_Mat   Null
##    1 0 f   1 0 f   Identical   Identical   Null
##    1 0 f   1 1 f   Haplo_Pat   Haplo_Pat   Null
##                    -------------------------------
##                     4 states   4 states    2 states
##
##  Returns
##
##  Globals Used: None
##
##-------------------------------------------------------------------------##
###WARNING USED IN BOTH MODULES!
sub getNumberOfStates {
  my $numberOfNonfounders = shift;

  return ( 4**( $numberOfNonfounders - 1 ) );
}    #end number_of_states

##-------------------------------------------------------------------------##
## Use:
## TODO: Restate this!
## See table in getNumberOfStates().  We can flip the maternal
## alleles for Autosomes and on X so we can fix the maternal
## indicator for the first child ( 00/10 ).  If it's autosomes we can fix
## both maternal/paternal for the first child (00).   The others can
## freely flip.  So depending on the gender of the first child we have
## 10(f) or 00 (m) in on the X...and max 11 for each additional child.
##-------------------------------------------------------------------------##
sub getMaxStateValue {
  my $numberOfNonfounders    = shift;
  my $minPloidy              = shift;
  my $nonfounderGenderVector = shift;

  my $maxStateValue = 0;

  if ( defined $minPloidy && $minPloidy == 1 ) {
    if ( defined $nonfounderGenderVector ) {
      my $indPos = 0;
      foreach my $gender ( split( //, $nonfounderGenderVector ) ) {
        $maxStateValue = $maxStateValue << 2;
        if ( $gender == 1 ) {
          if ( $indPos == 0 ) {
            $maxStateValue |= 2;
          }
          else {
            $maxStateValue |= 3;
          }
        }
        else {
          $maxStateValue |= 1;
        }
        $indPos++;
      }
    }
    else {
      $maxStateValue = 2**( $numberOfNonfounders * 2 ) - 1;
      $maxStateValue -= 2**( ( $numberOfNonfounders * 2 ) - 2 );
    }
  }
  else {
    $maxStateValue = 4**( $numberOfNonfounders - 1 ) - 1;
  }

  return ( $maxStateValue );
}

1;
