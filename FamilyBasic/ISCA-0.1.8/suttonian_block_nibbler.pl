#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) suttonian_block_nibbler.pl
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
#   modified 1/18/2011 v.03; crash recovery during editing of basic code
#
###############################################################################
# TO DO
#    Currently this program depends on the emission probability table.
#      - Formalize how HMM data is stored/retreived.
#

=head1 NAME

 suttonian_block_nibbler.pl - Remove ambiguous edges of HMM state output

=head1 SYNOPSIS

 Usage:
 ./suttonian_block_nibbler.pl -infile_directory <dir_where_hmm_data_is_stored>

=head1 DESCRIPTION

  This program chews away at the end of HMM-determined blocks becuase 
  the HMM will call a state even if it is not 100% certain a particular 
  position is in that state.

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
use strict;
use FindBin;
use lib $FindBin::RealBin;
use Storable qw(nstore retrieve);
use Data::Dumper;
use Getopt::Long;
use warnings FATAL => 'all';
use TimeUtils qw(timestamp);
use InheritanceStates
    qw(load_raw_state_for_each_position load_raw_state_for_each_position_on_a_chromosome load_block_file genotype_vector_for_index);

# Get version from Makefile/GIT
my $releaseVersion = "0.1.8";
my $gitVersion     = "0.1.8";
my $VERSION        = "release $releaseVersion";
$VERSION .= " ( git: $gitVersion )"
    if ( $gitVersion ne "" );

my $program = "suttonian_block_nibbler.pl";
print STDERR "Version $VERSION of $program, running with timestamp "
    . timestamp() . "\n";

my $infile_directory;
GetOptions( "infile_directory=s" => \$infile_directory );

if ( !defined $infile_directory ) {
  usage();
}

if (    !-s "$infile_directory/hmm_Suttonian_intervals.txt"
     || !-s "$infile_directory/hmm_DataCache.dat" )
{
  print STDERR
      "\n\nCannot locate input files ( ie. hmm_Suttonian_intervals.txt ..)\n\n";
  usage();
}
my $dataCacheFilename = "$infile_directory/hmm_DataCache.dat";

sub usage {
  print "$0 - $VERSION\n\n";
  exec "pod2text $0";
  exit;
}

my @labels_of_chrs_to_study = ( 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22  );

my @chromosomes_to_study = map { "chr" . $_ } @labels_of_chrs_to_study;

# TODO MAKE THIS A PARAMETER OR READ FROM HMM OUTPUT! ( and not limited to 4 )
# RMH: This is for Jared's emission model
#my @pce_emission_prob = ( 0.003, 0.003, 0.003, 0.003 );
# RMH: This is for my emission model
my @pce_emission_prob = (
                          0.00419894709779826, 0.00503551912849316,
                          0.00503551912849316, 0.00512820831385635
);

my $number_of_nonfounders =
    2;    #synonymous with number of children for a nuclear family
my @possible_genotypes = ( "aa", "ab", "an", "bb", "bn", "nn" );

my $number_of_meioses_in_pedigree =
    2 * $number_of_nonfounders;    #this is the length of the inheritance vector

print STDERR
    "Directory for finding the raw HMM output files is $infile_directory.\n";
my ( $emission, $positions_of_bounding_markers ) =
    load_raw_state_for_each_position( $infile_directory );
my %emission = %$emission;

my $short_filename = "hmm_Suttonian_intervals.txt";
my $long_filename  = "$infile_directory/$short_filename";
my %emission_prob  = ();

load_previously_compiled_data( $dataCacheFilename );

my %input_file_header_info = ();
my ( $block, $number_of_markers, $start, $end, $input_file_header_info ) =
    load_block_file( $long_filename );
my %block             = %$block;
my %number_of_markers = %$number_of_markers;
my %start             = %$start;
my %end               = %$end;

%input_file_header_info = %$input_file_header_info;

my $dir_root                = $infile_directory;
my $outfile                 = "nibbled_raw_blocks.txt";
my $header_for_output_files = make_header_for_output_files();
open( BLOCKS, ">$dir_root/$outfile" );
print BLOCKS $header_for_output_files;

#print BLOCKS "#header\tchr\tunsmoothed original block_number\tstate\tnumber of markers\tstart\tend\tlength\n";

my $number_of_nibbled_positions         = 0;
my $inconsistent_nibbled_positions      = 0;
my $dually_consistent_nibbled_positions = 0;
my $total_markers_in_considered_blocks  = 0;
my $total_completely_lost_blocks        = 0;

#print "p = " . Dumper( \@positions_on_this_chromosome ) . "\n";
foreach my $chr ( @chromosomes_to_study ) {

  #print uc $chr, "\n";
  my @block_numbers =
      sort { $a <=> $b }
      keys %{ $block{$chr} };    #$block{$chr}{$block_number} = $data[2];
  my $first_block_number = $block_numbers[ 0 ];
  my $last_block_number  = $block_numbers[ -1 ];
print "\n\n\naaaaaaaaa\n".join(":",@block_numbers)."zzzzzzzzz\n\n\n";
  my @positions_on_this_chromosome =
      sort { $a <=> $b } keys %{ $emission{$chr} };
	#print join(":",@positions_on_this_chromosome)."zzzzzzzzz\n\n\n";
  if ( $first_block_number == $last_block_number ) { #no chewing to be done here
    my $state        = $block{$chr}{1};
    my $binary_state =
        string_binary_representation_of_inheritance_state( $state );
    my @positions_in_this_block = ();

    #make a list of positions in this block
    foreach my $position ( @positions_on_this_chromosome ) {
      next
          unless (     ( $position >= $start{$chr}{1} )
                   and ( $position <= $end{$chr}{1} ) );
      push @positions_in_this_block, $position;
    }
    $start = $positions_in_this_block[ 0 ];
    $end   = $positions_in_this_block[ -1 ];
    my $length                          = $end - $start + 1;
    my $number_of_markers_in_this_block = scalar @positions_in_this_block;
    $total_markers_in_considered_blocks += $number_of_markers_in_this_block;
    print BLOCKS
"$chr\t1\t$state\t$binary_state\t$number_of_markers_in_this_block\t$start\t$end\t$length\n";
    next;
  }

  my $completely_lost_blocks = 0;
  foreach my $block_number ( @block_numbers ) {

    #%emission_counts = ();
    #$block_type = $block{$chr}{$block_number};
    my $state                   = $block{$chr}{$block_number};
    my @positions_in_this_block = ();

    #make a list of positions in this block
    foreach my $position ( @positions_on_this_chromosome ) {
      next
          unless (     ( $position >= $start{$chr}{$block_number} )
                   and ( $position <= $end{$chr}{$block_number} ) );

      #$emission = $emission{$chr}{$position};
      push @positions_in_this_block, $position;
    }

    my $chewed_positions_in_front = 0;
    unless ( $block_number == $first_block_number ) {
      my $previous_block_number   = $block_number - 1;
      my $state_of_previous_block = $block{$chr}{$previous_block_number};

      # Chew beginning of block
      # If emission is inconsistent with this block, chew block
      # If emission is compatible with both this block and the
      # previous block, chew block.
      foreach my $position ( @positions_in_this_block ) {
        my $emission                        = $emission{$chr}{$position};
        my $emission_prob                   = $emission_prob{$state}{$emission};
        my $emission_prob_of_previous_block =
            $emission_prob{$state_of_previous_block}{$emission};
        my $chewed = 0;
        if ( $emission_prob <= $pce_emission_prob[ $state ] )
        {    #emission is inconsistent with this block
          $chewed = 1;
          $inconsistent_nibbled_positions++;
        }
        elsif ( $emission_prob_of_previous_block >= $emission_prob )
        {    #emission is compatible with both this block and the previous block
          $chewed = 1;
          $dually_consistent_nibbled_positions++;
        }

        if ( $chewed ) {
          $chewed_positions_in_front++;
          $number_of_nibbled_positions++;

          #print STDERR "Chewing position $position with emission " .
          #             "$emission and prob $emission_prob with " .
          #             "prob of previous block " .
          #             "$emission_prob_of_previous_block\n";
        }
        else {
          last;    #stop chewing if reach an anchor
        }
      }
    }

    my $chewed_positions_in_back = 0;
    unless ( $block_number == $last_block_number ) {
      my $next_block_number   = $block_number + 1;
      my $state_of_next_block = $block{$chr}{$next_block_number};

      # Chew end of block
      # if emission is inconsistent with this block, chew block
      # if emission is compatible with both this block and the
      # next block, chew block.
      foreach my $position ( reverse @positions_in_this_block ) {
        my $emission                    = $emission{$chr}{$position};
        my $emission_prob               = $emission_prob{$state}{$emission};
        my $emission_prob_of_next_block =
            $emission_prob{$state_of_next_block}{$emission};
        my $chewed = 0;
        if ( $emission_prob <= $pce_emission_prob[ $state ] )
        {    #emission is inconsistent with this block
          $chewed = 1;
          $inconsistent_nibbled_positions++;
        }
        elsif ( $emission_prob_of_next_block >= $emission_prob )
        {    #emission is compatible with both this block and the previous block
          $chewed = 1;
          $dually_consistent_nibbled_positions++;
        }
        if ( $chewed ) {

          #print STDERR "Chewing position $position with emission " .
          #             "$emission and prob $emission_prob with " .
          #             "prob of next block $emission_prob_of_next_block\n";
          $chewed_positions_in_back++;
          $number_of_nibbled_positions++;
        }
        else {
          last;    #stop chewing if reach an anchor
        }
      }

    }

    #trim now

    my $number_of_markers_in_this_block = scalar @positions_in_this_block;
    my $total_markers_in_considered_blocks += $number_of_markers_in_this_block;
    my $new_end =
        $number_of_markers_in_this_block - $chewed_positions_in_back - 1;
    print STDERR "FOR BLOCK $chr:$block_number there were "
        . "$number_of_markers_in_this_block markers.\n"
        . "      Chewing $chewed_positions_in_back from "
        . "back and $chewed_positions_in_front from "
        . "beginning with new end $new_end\n";

    my $n = 0;
    while ( $n < $chewed_positions_in_front ) {
      shift @positions_in_this_block;
      $n++;
    }
    $n = 0;
    while ( $n < $chewed_positions_in_back ) {
      pop @positions_in_this_block;
      $n++;
    }

    $number_of_markers = scalar @positions_in_this_block;

    if ( $number_of_markers <= 0 ) {
      $completely_lost_blocks++;
    }
    else {
      my $start        = $positions_in_this_block[ 0 ];
      my $end          = $positions_in_this_block[ -1 ];
      my $length       = $end - $start + 1;
      my $binary_state =
          string_binary_representation_of_inheritance_state( $state );
      my $output_block_number = $block_number - $completely_lost_blocks;
      print BLOCKS
"$chr\t$output_block_number\t$state\t$binary_state\t$number_of_markers\t$start\t$end\t$length\n";
    }
  }

  $total_completely_lost_blocks += $completely_lost_blocks;

}
close BLOCKS;

print STDERR "number_of_nibbled_positions\t$number_of_nibbled_positions\n";
print STDERR
    "inconsistent_nibbled_positions\t$inconsistent_nibbled_positions\n";
print STDERR
"dually_consistent_nibbled_positions\t$dually_consistent_nibbled_positions\n";
print STDERR
    "total_markers_in_considered_blocks\t$total_markers_in_considered_blocks\n";
print STDERR "total_completely_lost_blocks\t$total_completely_lost_blocks\n";

exit;

#END MAIN BLOCK

####################### S U B R O U T I N E S #######################

sub string_binary_representation_of_inheritance_state {
  my $decimal = shift;
  my $binary  = sprintf( "%b", $decimal );

  #now need to pad with leading zeros
  while ( length $binary < $number_of_meioses_in_pedigree ) {
    $binary = "0" . $binary;
  }
  return $binary;
}    #end binary_representation_of_inheritance_state

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
  %emission_prob = %{ $tmpDataCache->{'emission_prob'} };

  # Free up memory -- no need to keep a duplicate nor all the extra stuff
  undef $tmpDataCache;
}    #end load_previously_compiled_data

sub make_header_for_output_files {
  my $header = "";
  $header .= "#parameter\tinput_file\t$short_filename\n";
  $header .=
"#parameter\tinput_file to suttonian_hmm.pl i.e. the input file for the input file\t$input_file_header_info{'input_file'}\n";
  $header .= "#creator\t$program version $VERSION\n";
  $header .= "#timestamp\t" . timestamp() . "\n";
  $header .= "#contact\tJared Roach jroach\@systemsbiology.org\n";
  $header .= "#reference_genome\t$input_file_header_info{'reference_genome'}\n";
  $header .= "#pedigree_name\t$input_file_header_info{'pedigree_name'}\n";
  $header .= "#numerical_base\t$input_file_header_info{'numerical_base'}\n";
  $header .= "#project\t$input_file_header_info{'project'}\n";
  $header .=
"#chromosome_nomenclature\t$input_file_header_info{'chromosome_nomenclature'}\n";
  $header .= "#inheritance_state_representation\tarabic_zero_based\n";
  $header .=
"#the following header lines are automatically parsed from the header of the input file\n";

  foreach my $key ( sort keys %input_file_header_info ) {
    $header .= "#$key\t$input_file_header_info{$key}\n";
  }
  $header .=
"#the previous header lines were automatically parsed from the header of the input file\n";

  return $header;
}    #end make_header_for_output_files

