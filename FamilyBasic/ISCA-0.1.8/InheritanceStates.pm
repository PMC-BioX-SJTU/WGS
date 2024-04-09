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
#*  This work is licensed under the Open Source License v2.1.  To view a copy
#*  of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
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
#
###############################################################################
# TO DO
#

=head1 NAME

InheritanceStates

=head1 SYNOPSIS

use InheritanceStates;

Usage:

  InheritanceStates::Package_Method();

=head1 DESCRIPTION

A set of parsers for ISCA text output files.

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2011 Institute for Systems Biology

=head1 AUTHORS

  Jared Roach      <jroach@systemsbiology.org>
  Robert M. Hubley <rhubley@systemsbiology.org>

=head1 PACKAGE METHODS

=cut

##
## Module dependence
##
use strict;
use warnings;

package InheritanceStates;
use Data::Dumper;
use File::Spec;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
require Exporter;
@ISA = qw(Exporter);

@EXPORT =
    qw(load_raw_state_for_each_position load_raw_state_for_each_position_on_a_chromosome load_block_file genotype_vector_for_index);

##-------------------------------------------------------------------------##

=head2 load_raw_state_for_each_position()

  Use: my ( $emission, $positions_of_bounding_markers ) = 
       InheritanceStates::load_raw_state_for_each_position( $infile_directory );

  $emissions = Hash reference
  $positions_of_boiunding_markers = Hash reference

=cut

##-------------------------------------------------------------------------##
sub load_raw_state_for_each_position {
  my $infile_directory = shift;

  my $nextfile;
  my %seen_chromosomes;
  my $chr = "";
  my %emission;
  my %positions_of_bounding_markers;

  # Process directory looking for state files.
  # e.g. hmm_Suttonian_states_chr2.txt
  while ( $nextfile = <$infile_directory/hmm_Suttonian_states_*.txt> ) {

    #next if file is a link
    next if ( -l $nextfile );

    #next if file is a directory
    next if ( -d $nextfile );

    my ( $volume, $path, $shortfilename ) = File::Spec->splitpath( $nextfile );
    print STDERR "Working on $shortfilename\n";

    if ( $shortfilename =~ /_(chr[^_]*)\./ ) {
      $chr = $1;
      if ( $seen_chromosomes{$chr} ) {
        die "This chromosome $chr appears more than once in the input "
            . "directory for example in $shortfilename. Please fix the "
            . "error or rename the file to trump this error correction "
            . "check.\n";
      }
      $seen_chromosomes{$chr} = 1;
    }

    my ( $emissionRef, $lowestPosition, $highestPosition ) =
        load_raw_state_for_each_position_on_a_chromosome( $nextfile );

    # hacking the data structure a little bit to encode the first and
    # last marker positions on each chromosome
    $positions_of_bounding_markers{$chr}{'low'}  = $lowestPosition;
    $positions_of_bounding_markers{$chr}{'high'} = $highestPosition;

    $emission{$chr} = {%$emissionRef};

  }

  return ( \%emission, \%positions_of_bounding_markers );
}    #end load_raw_state_for_each_position

##-------------------------------------------------------------------------##

=head2 load_raw_state_for_each_position_on_a_chromosome()

  Use: my ( $emission, $positions_of_bounding_markers ) = 
       InheritanceStates::load_raw_state_for_each_position_on_a_chromosome( 
                                                           $file );

  $emissions = Hash reference
  $positions_of_boiunding_markers = Hash reference

=cut

##-------------------------------------------------------------------------##
sub load_raw_state_for_each_position_on_a_chromosome {
  my $file = shift;

  my $chr;
  my $position;
  my %emission = ();
  my $emissionVal;
  my @data;
  my $lowest_position  = 1e99;
  my $highest_position = 0;
  open IN, $file;

  while ( <IN> ) {
    chomp;
    next if /^#/;        #
    next if /^index/;    #
    next unless ( /\S/ );

    @data     = split( /\t/, $_ );
    $chr      = $data[ 1 ];
    $position = $data[ 2 ];
    if ( $position < $lowest_position )  { $lowest_position  = $position }
    if ( $position > $highest_position ) { $highest_position = $position }
    $emissionVal = $data[ 3 ];

#$block = $data[4];  #may use this later, but for now it should corresponf to the unsmoothed state in the block file
    $emission{$position} = $emissionVal;

#print STDERR "load_raw_state_for_each_position_on_a_chromosome $emissionVal $chr $position\n";
  }
  close IN;

  return ( \%emission, $lowest_position, $highest_position );

}    #end load_raw_state_for_each_position_on_a_chromosome

##-------------------------------------------------------------------------##

=head2 load_block_file()

  Use: my ( $block, $number_of_markers, $start, $end, 
            $input_file_header_info ) = 
           InheritanceStates::load_block_file( $file )

  All return values are scalar references to hashes.

=cut

##-------------------------------------------------------------------------##
sub load_block_file {
  my $long_filename = shift;

  my $index_for_duplicate_header_keys = 1;
  my @inline;
  my @data;
  my $chr;
  my $block_number;
  my ( %block, %number_of_markers, %start, %end, %seen_chr,
       %input_file_header_info );
  %input_file_header_info = ();

  unless ( -e $long_filename ) {
    print STDERR
"No file $long_filename, which was sent as input to subroutine load_block_file in InheritanceStates.pm!\n";
    die;
  }
  my $age_of_file = sprintf( "%.1f", -M $long_filename );
  print STDERR "Opening $long_filename, which is $age_of_file days old.\n";
  open VAR, $long_filename
      or die "[load_block_file]: Could not open $long_filename!\n";
  my $key;

  while ( <VAR> ) {
    chomp;

    # Keep header lines
    if ( /^#/ ) {

      # Make hash of header info
      if ( $_ =~ /^\#(\S*)\s+(.*)$/ ) {
        $key = $1;
        if ( $input_file_header_info{$key} ) {
          $key = $key . $index_for_duplicate_header_keys;
          $index_for_duplicate_header_keys++;
        }
        $input_file_header_info{$key} = $2;
      }
      else {
        print STDERR "$_ is a bad header line!\n";
        die;
      }
      next;
    }

  # Field order: chromosome	block_number	most_likely_path_state
  #              number_of_markers_in_block	block_start	block_end
  #              block_length	recombination meiosis	length of crossover interval
  #              crossover start location
  # NOT SURE WHY CHOMPING TWICE????
    chomp;

    @data = split( /\t/, $_ );

    $chr          = $data[ 0 ];
    $block_number = $data[ 1 ];

    $block{$chr}{$block_number} = $data[ 2 ];

    #print STDERR "Loading block $block_number $chr\n";

    #binary inheritance state = $data[3];

    $number_of_markers{$chr}{$block_number} = $data[ 4 ];
    $start{$chr}{$block_number}             = $data[ 5 ];
    $end{$chr}{$block_number}               = $data[ 6 ];

    #$length{$chr}{$block_number} = $data[6];
    #print STDERR "The genotype at $chr $pos is $genotype.\n";
    #$seen_chr{$chr} = 1;
    #$seen_position{$chr}{$pos} = 1;
    #die;
  }
  close VAR;

  #$/ = "\n"; #switch back in case mac-specific returns were turned on;

  return ( \%block, \%number_of_markers, \%start, \%end,
           \%input_file_header_info );
}    #end load_block_file

##-------------------------------------------------------------------------##

=head2 genotype_vector_for_index()

  Use: my $genotype_vector = InheritanceStates::genotype_vector_for_index( 
                              $index,
                              $maximum_genotype_vector_index,
                              $number_of_individuals_in_pedigree,
                              $genotype_for_integer );


=cut

##-------------------------------------------------------------------------##
sub genotype_vector_for_index {

  my $index                             = shift;
  my $maximum_genotype_vector_index     = shift;
  my $number_of_individuals_in_pedigree = shift;
  my $genotype_for_integerRef           = shift;

  my %genotype_for_integer = %$genotype_for_integerRef;
  my $genotype_integer;
  my %genotype_content = ();

  if ( $index == $maximum_genotype_vector_index + 1 ) {
    return
        "tri_quad"
        ; #very odd genotype vector with more than 2 variants at this position -> probably a compression
  }

  my @genotype_integer_vector = ();
  my $n                       = $number_of_individuals_in_pedigree - 1;
  $genotype_integer_vector[ $n ] = int( $index / ( 6**$n ) );
  my $remainder = $index % ( 6**$n );
  $n--;
  while ( $n >= 0 ) {
    $genotype_integer_vector[ $n ] = int( $remainder / ( 6**$n ) );
    $remainder = $index % ( 6**$n );
    $n--;
  }

  my $genotype_vector = "";
  while ( @genotype_integer_vector ) {
    $genotype_integer = shift @genotype_integer_vector;
    if ( defined( $genotype_for_integer{$genotype_integer} ) ) {
      $genotype_vector .= $genotype_for_integer{$genotype_integer};
    }
    else {
      print STDERR "No genotype for integer $genotype_integer!\n";
      $genotype_vector .= "00";
    }

  }
  return $genotype_vector;
}    #end genotype_vector_for_index

#### Return Value of Package Module
1;
