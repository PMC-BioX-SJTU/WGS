#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) intercalate_partial_binary_blocks.pl
##  Author:
##      Jared Roach      <jroach@systemsbiology.org>
##      Robert M. Hubley <rhubley@systemsbiology.org>
##  Description:
##      This program puts in partially called inheritance states in between
##      fully called states.
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
#   modified 2/3/2011 v0.02; specify end positions in zero based rather than
#      one-based coordinates; this imporves compatibilty with downstream scripts
#   modified 2/3/2011 v0.03; added ambiguity character to first two bits of
#      vector rather than defaulting to 0; fixing trailing period that
#      rarely shows up
#
###############################################################################
# TO DO
#

=head1 NAME

 intercalate_partial_binary_blocks.pl - Add partially called inheritance states

=head1 SYNOPSIS

 Usage:
 ./intercalate_partial_binary_blocks.pl  --infile_directory <DIR>

=head1 DESCRIPTION

   This program puts in partially called inheritance states in between 
   fully called states.

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
use FindBin;
use lib $FindBin::RealBin;
use Getopt::Long;
use warnings FATAL => 'all'
    ; #like to trap all warnings when running in workflow environment becuase otherwise it is easy to miss the warning
use TimeUtils qw(timestamp);
use InheritanceStates
    qw(load_raw_state_for_each_position load_raw_state_for_each_position_on_a_chromosome load_block_file genotype_vector_for_index);

my $releaseVersion = "0.1.8";
my $gitVersion     = "0.1.8";
my $VERSION        = "release $releaseVersion";
my $program        = "intercalate_partial_binary_blocks.pl";
print STDERR "Version $VERSION of $program, running with timestamp "
    . timestamp() . "\n";

my $infile_directory;
GetOptions( "infile_directory=s" => \$infile_directory );
unless ( defined( $infile_directory ) ) {
  print STDERR "\n\nNo infile directory defined!\n\n";
  usage();
}

sub usage {
  print "$0 - $VERSION\n\n";
  exec "pod2text $0";
  exit;
}

my @labels_of_chrs_to_study = ( 1 .. 22 );

#@labels_of_chrs_to_study = (21);  #debug
my @chromosomes_to_study = map { "chr" . $_ } @labels_of_chrs_to_study;

my $number_of_nonfounders =
    2;    #synonymous with number of children for a nuclear family
my $number_of_meioses_in_pedigree =
    2 * $number_of_nonfounders;    #this is the length of the inheritance vector

my $short_filename = "smoothed_blocks.txt";
my $long_filename  = "$infile_directory/$short_filename";

my %input_file_header_info = ();
my ( $block, $number_of_markers, $start, $end,, $input_file_header_info_ref ) =
    load_block_file( $long_filename );
my %block             = %$block;
my %number_of_markers = %$number_of_markers;
my %start             = %$start;
my %end               = %$end;
%input_file_header_info = %$input_file_header_info_ref;

my $dir_root = $infile_directory;
my $outfile  = "smoothed_blocks_with_intercalated_partial_blocks.txt";
my $header_for_output_files = make_header_for_output_files();
open( BLOCKS, ">$dir_root/$outfile" );
print BLOCKS $header_for_output_files;

#print BLOCKS "#header\tchr\tunsmoothed original block_number\tstate\tnumber of markers\tstart\tend\tlength\n";
print BLOCKS
"#header\tchromosome\tblock_number\tsmoothed nibbled HMM state\tbinary state\tnumber of markers\tblock_start\tblock_end\tblock_length\n";

my $state;
foreach my $chr ( @chromosomes_to_study ) {
  print "Working on " . uc( $chr ) . "\n";
  my @block_numbers =
      sort { $a <=> $b }
      keys %{ $block{$chr} };    #$block{$chr}{$block_number} = $data[2];
  my $first_block_number = $block_numbers[ 0 ];
  my $last_block_number  = $block_numbers[ -1 ];
  ## RMH: Fix...we should report singleton blocks anyway
  if ( @block_numbers == 0 || $first_block_number == $last_block_number ) {
    my $bs =
        string_binary_representation_of_inheritance_state( $block{$chr}{1} );
    my $length = $end{$chr}{1} - $start{$chr}{1};
    print BLOCKS
        "$chr\t1\t$state\t$bs\tND\t$start{$chr}{1}\t$end{$chr}{1}\t$length\n";

    # No intercalating to be done here.
    next;
  }
  my $new_block_numbering = 1;

  my %old_block;
  my %new_state;
  my %new_start;
  my $new_start;
  my $new_end;
  my %new_end;
  foreach my $block_number ( @block_numbers ) {

    $state = $block{$chr}{$block_number};
    $start = $start{$chr}{$block_number};

    #$end = $end{$chr}{$block_number};
    $end =
        $end{$chr}{$block_number} + 1
        ; #the input file is not using zero-based ends, so we have to add one to get to the zero-based output required of this program

    $old_block{$chr}{$new_block_numbering} = $block_number;
    $new_state{$chr}{$new_block_numbering} = $state;
    $new_start{$chr}{$new_block_numbering} = $start;

    $new_end{$chr}{$new_block_numbering} = $end;

    last if ( $block_number == $last_block_number );

    #$previous_block_number = $block_number-1;
    my $next_block_number = $block_number + 1;

    my $next_intercalated_block = $new_block_numbering + 1;

    #$new_start = $end+1;
    $new_start = $end;

    #$new_end = $start{$chr}{$next_block_number}-1;
    $new_end = $start{$chr}{$next_block_number};

    if ( $new_end < $new_start )
    {    #intercalation length is zero (highly unlikely)
          #do nothing
      $new_block_numbering++;
    }
    else {

      #$new_start{$chr}{$next_intercalated_block} = $end+1;
      $new_start{$chr}{$next_intercalated_block} = $end;

      $old_block{$chr}{$next_intercalated_block} = $block_number + 0.5;

 #$new_end{$chr}{$next_intercalated_block} = $start{$chr}{$next_block_number}-1;
      $new_end{$chr}{$next_intercalated_block} =
          $start{$chr}{$next_block_number};

      $new_block_numbering++;
      $new_block_numbering++;
    }
  }
  @new_block_numbers = sort { $a <=> $b } keys %{ $new_start{$chr} };
  foreach my $block_number ( @new_block_numbers ) {
    if ( defined( $new_state{$chr}{$block_number} ) ) {
      $state = $new_state{$chr}{$block_number};

      $binary_state =
          string_binary_representation_of_inheritance_state( $state );

#print STDERR "$chr $block_number $state $binary_state\t$old_block{$chr}{$block_number}\n";

    }
    else {
      $state                  = "partial";
      $previous_defined_state = $block_number - 1;
      $next_defined_state     = $block_number + 1;

#print STDERR "Interpolating state between block $chr $previous_defined_state and $next_defined_state which have states $new_state{$chr}{$previous_defined_state},$new_state{$chr}{$next_defined_state}\tpending\t$old_block{$chr}{$block_number}\n";

      $binary_state =
          intercalated_binary_state( $new_state{$chr}{$previous_defined_state},
                                     $new_state{$chr}{$next_defined_state} );

#print STDERR "Interpolating state between block $chr $previous_defined_state and $next_defined_state which have states $new_state{$chr}{$previous_defined_state},$new_state{$chr}{$next_defined_state}\t$binary_state\t$old_block{$chr}{$block_number}\n";

    }

    $length = $new_end{$chr}{$block_number} - $new_start{$chr}{$block_number};
    print BLOCKS
"$chr\t$block_number\t$state\t$binary_state\tND\t$new_start{$chr}{$block_number}\t$new_end{$chr}{$block_number}\t$length\n";
  }
}
close BLOCKS;

exit;

#END MAIN BLOCK

sub intercalated_binary_state {
  my $previous_state = shift;
  die unless defined( $previous_state );
  my $next_state = shift;
  die unless defined( $next_state );
  ( $previous_state, $next_state ) =
      sort { $a <=> $b } ( $previous_state, $next_state );
  if ( $previous_state == 0 ) {
    if ( $next_state == 0 ) {
      return "0000";
    }
    elsif ( $next_state == 1 ) {

      #return "000."
      return "0.0.";
    }
    elsif ( $next_state == 2 ) {

      #return "00.0"
      return ".0.0";
    }
    elsif ( $next_state == 3 ) {

      #return "00.."
      return "....";
    }
    else {
      die;
    }
  }
  elsif ( $previous_state == 1 ) {
    if ( $next_state == 1 ) {
      return "0001";
    }
    elsif ( $next_state == 2 ) {

      #return "00.."
      return "....";
    }
    elsif ( $next_state == 3 ) {

      #return "00.1"
      return ".0.1";
    }
    else {
      die;
    }
  }
  elsif ( $previous_state == 2 ) {
    if ( $next_state == 2 ) {

      #return "0010."  #trailing period a bug prior to version v.03
      return "0010";
    }
    elsif ( $next_state == 3 ) {

      #return "001."
      return "0.1.";
    }
    else {
      die;
    }
  }
  elsif ( $previous_state == 3 ) {
    if ( $next_state == 3 ) {
      return "0011";
    }
    else {
      die;
    }
  }
  else {
    die;
  }
  return "NoResult";
}

#end intercalated_binary_state

sub string_binary_representation_of_inheritance_state {
  my $decimal = shift;
  my $binary  =
      sprintf( "%b", $decimal )
      ; #note that this only works with (Using sprintf) Perl 5.6+  If using an earlier version , need to use another recipe for converting to binary
        #now need to pad with leading zeros
  while ( length $binary < $number_of_meioses_in_pedigree ) {
    $binary = "0" . $binary;
  }
  return $binary;
}

#end binary_representation_of_inheritance_state

sub make_header_for_output_files {
  my $header = "";

#$header .= "#parameter\tmode\t$mode\n";
#$header .= "#parameter\tinput_file\t$short_filename\n";
#$header .= "#parameter\tinput_file to suttonian_hmm.pl i.e. the input file for the input file\t$input_file_header_info{'input_file'}\n";
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

