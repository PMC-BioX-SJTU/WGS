#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) parse_and_smooth_HMM_Suttonian_blocks.pl
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
#   modified 9/28/09   v0.02  added chrX mode
#   modified 10/15/10   v0.03  major changes to make aware of block statistics
#       and re-do general algorithm
#   modified 10/16/10   v0.04  also get rid of strong,w,{w*},w,same strong blocks
#   modified 10/16/10   v0.05  improve headers; add in output file for suspicious blocks
#   modified 10/17/10   v0.07  allow manual override of strength for particular blocks
#
#
################################################################################
# TO DO
#
#

=head1 NAME

 parse_and_smooth_HMM_Suttonian_blocks.pl -

=head1 SYNOPSIS

 Usage:
 ./parse_and_smooth_HMM_Suttonian_blocks [ --infile_directory <DIR> ]

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
#like to trap all warnings when running in workflow environment becuase otherwise it is easy to miss the warning
use warnings FATAL => 'all';
use FindBin;
use lib $FindBin::RealBin;
use Cwd;
use Data::Dumper;
use Getopt::Long;
use InheritanceStates
    qw(load_raw_state_for_each_position load_raw_state_for_each_position_on_a_chromosome genotype_vector_for_index);
use TimeUtils qw(timestamp);
use Storable qw(nstore retrieve);

my $infile_directory = cwd();
GetOptions( "infile_directory=s" => \$infile_directory );

my $releaseVersion = "0.1.8";
my $gitVersion     = "0.1.8";
my $VERSION        = "release $releaseVersion";
my $program        = "parse_and_smooth_HMM_Suttonian_blocks_bobama.pl";
print "Version $VERSION of $program, running with timestamp "
    . timestamp() . "\n";

#my $mode = "chrX";
my $mode = "autosomes";

my $short_filename = "$infile_directory/nibbled_raw_blocks.txt";
my $long_filename  = $short_filename;

my $hmmDataCacheFilename  = "$infile_directory/hmm_DataCache.dat";
my $statDataCacheFilename = "$infile_directory/stats_DataCache.dat";

if ( !-s $short_filename ) {
  print "\n\nMissing $short_filename!\n\n";
  usage();
}

sub usage {
  print "$0 - $VERSION\n\n";
  exec "pod2text $0";
  exit;
}

load_previously_compiled_data( $hmmDataCacheFilename, $statDataCacheFilename );
load_file( $long_filename );

#$number_of_markers_to_be_short = 500;
#$length_to_be_short = 500000;
#$absolute_minimum_CNV_or_error_block_length = 50000;
#$absolute_minimum_CNV_or_error_number_of_markers = 200;

$number_of_nonfounders =
    2;    #synonymous with number of children for a nuclear family
$number_of_meioses_in_pedigree =
    2 * $number_of_nonfounders;    #this is the length of the inheritance vector

# Jared's original
#$cutoff_for_good_quality_block = 0.4;

# New cutoff
# Both Parents
if ( $input_file_header_info{'note'} =~ /inclusion_order\s*=\s*([\d,xX]+)\s+/ )
{
  my $incOrder = $1;
  if ( $incOrder =~ /^[xX],[xX],.*/ ) {

    # Both parents missing
    $cutoff_for_good_quality_block = 0.26;
    print "Zero founders cutoff = $cutoff_for_good_quality_block\n";
  }
  elsif ( $incOrder =~ /^[xX],[\dxX],.*/ || $incOrder =~ /^[\dxX],[xX],.*/ ) {

    # One parent missing
    $cutoff_for_good_quality_block = 0.30;
    print "One founder cutoff = $cutoff_for_good_quality_block\n";
  }
  else {

    # Both parents
    $cutoff_for_good_quality_block = 0.34;
    print "Two founders cutoff = $cutoff_for_good_quality_block\n";
  }
}

$genome_length                     = 0;
$total_length_of_suspicious_blocks = 0;
$total_length_of_uncalled_blocks   = 0;

$outfile                 = "suspicious_blocks.txt";
$full_path_outfile       = "$infile_directory/$outfile";
$header_for_output_files = make_header_for_output_files();
open( SUSPICIOUS, ">$full_path_outfile" );
print SUSPICIOUS $header_for_output_files;
print SUSPICIOUS
"#header chr\tunsmoothed original block_number\tstate\tbinary state\tnumber_of_markers\tstart\tend\tlength\n";

@labels_of_chrs_to_study = ( 1 .. 22 );
@chromosomes_to_study    = map { "chr" . $_ } @labels_of_chrs_to_study;

foreach $chr ( @chromosomes_to_study ) {

  #for debugging
  #next unless ($chr eq "chr8");

  print uc $chr, "\n";
  @block_numbers =
      sort { $a <=> $b }
      keys %{ $block{$chr} };    #$block{$chr}{$block_number} = $data[2];
  $number_of_blocks =
      scalar( @block_numbers )
      ;    #should be same as last block number - why twice? - bad programming
  @blocks_to_eliminate =
      ()
      ; #these are weak blocks that are bordered by two strong blocks in the same state
  @blocks_to_consider_uncalled =
      ()
      ;  #these are weak blocks that are bordered by two different strong blocks
  print "There are $number_of_blocks blocks.\n";

  #if no strong block on entire chr report and error
  $n                 = 1;
  $strong_block_seen = 0;
  while ( $n <= $number_of_blocks ) {
    $block_strength = $block_statistics{'strength'}{$chr}{$n};
    if ( $block_strength > $cutoff_for_good_quality_block ) {
      $strong_block_seen = 1;
      last;
    }
    $n++;
  }

  if ( $strong_block_seen ) {

    # More agressive version of the above
    $n                    = 1;
    $first_block_strength = $block_statistics{'strength'}{$chr}{$n};
    print "First Block strength for $chr block $n is $first_block_strength\n";
    die unless ( defined( $first_block_strength ) );
    while ( $first_block_strength < $cutoff_for_good_quality_block ) {
      push @blocks_to_eliminate, $n;
      print "Eliminating\t$n\t$block{$chr}{$n} [weak,weak,etc.\n";
      $n++;

      last unless ( defined( $block_statistics{'strength'}{$chr}{$n} ) );

      $first_block_strength = $block_statistics{'strength'}{$chr}{$n};
      print "First block strength for $chr block $n is $first_block_strength\n";
    }
    $strong_block_seen = 0
        if ( $first_block_strength < $cutoff_for_good_quality_block );

    if ( $strong_block_seen ) {
      $n                   = $number_of_blocks;
      $last_block_strength = $block_statistics{'strength'}{$chr}{$n};
      print "Last Block strength for $chr block $n is $last_block_strength\n";
      die unless ( defined( $last_block_strength ) );
      while ( $last_block_strength < $cutoff_for_good_quality_block ) {
        push @blocks_to_eliminate, $n;
        print "Eliminating\t$n\t$block{$chr}{$n} weak,weak,etc.]\n";
        $n--;

        last unless ( defined( $block_statistics{'strength'}{$chr}{$n} ) );

        $last_block_strength = $block_statistics{'strength'}{$chr}{$n};
        print "Last block strength for $chr block $n is $last_block_strength\n";
      }

      #get rid of strong,w,{w*},w,same strong blocks
      $n                      = 1;
      $flag_for_pruning       = 0;
      @blocks_to_add_to_prune = ();
      while ( $n <= $number_of_blocks ) {
        $block_strength = $block_statistics{'strength'}{$chr}{$n};
        if ( ( $block_strength > $cutoff_for_good_quality_block )
             and $flag_for_pruning )
        {

          #prune the currently accumulated list
          #print STDERR "Pruning list added when saw block $n.\n";
          $current_block_type = $block{$chr}{$n};
          if ( $previous_block_type == $current_block_type )
          {    #flanking strong block states are the same
            @blocks_to_eliminate =
                ( @blocks_to_eliminate, @blocks_to_add_to_prune );print "bbb\t$n\t$current_block_type\n";
            foreach $q ( @blocks_to_add_to_prune ) {
              print
"Eliminating\t$q\t$block{$chr}{$q} strong,weak*,same strong\n";
            }
          }
          else {
            $Hamming_distance_between_states =
                $states_are_adjacent{$previous_block_type}{$current_block_type};print "aaaaaaa\t$n\t$Hamming_distance_between_states\n";
            if ( $Hamming_distance_between_states == 1 )
            { #go ahead and delete the intervening states because the two strong states can directly transition
               #problem here is that we do not really know how to assign the markers in the deleted states because they might belong to either of the flnaking states
               #best solution is to send the data back to the HMM (or equivalent border-detection-algorithm), permitting only two states and allowing only one transition
               #for now, we are going to declare the entire region "uncalled" for state
              @blocks_to_consider_uncalled =
                  ( @blocks_to_consider_uncalled, @blocks_to_add_to_prune );
              foreach $q ( @blocks_to_add_to_prune ) {
                print
"Considering uncalled\t$q\t$block{$chr}{$q} strong,weak*,strong one recombination away from first strong\n";
              }
            }
          }
          $previous_block_type    = $block{$chr}{$n};
          @blocks_to_add_to_prune = ();
          $n++;
          next;
        }
        elsif ( ( $block_strength > $cutoff_for_good_quality_block )
                and not $flag_for_pruning )
        {
          $flag_for_pruning       = 1;
          $previous_block_type    = $block{$chr}{$n};
          @blocks_to_add_to_prune = ();
          $n++;
          next;
        }
        else {    #weak block
          if ( $flag_for_pruning ) {

            #print STDERR "Pushing $n to prune.\n";
            push @blocks_to_add_to_prune, $n;
          }
          else {

            #do nothing
          }
          $n++;
          next;
        }
      }
    }
  }
  else {
    $n = 1;
    while ( $n <= $number_of_blocks ) {
      push @blocks_to_eliminate, $n++;
    }
  }

  #report any weak blocks remaining
  $number_of_blocks = scalar @block_numbers;
  $n                = 1;
  while ( $n <= $number_of_blocks ) {
    if ( grep { $_ eq $n } @blocks_to_eliminate ) {
      $n++;
      next;
    }
    if ( grep { $_ eq $n } @blocks_to_consider_uncalled ) {
      $n++;
      next;
    }
    $block_strength = $block_statistics{'strength'}{$chr}{$n};

#print "Block Strength = $block_strength - Cutoff = $cutoff_for_good_quality_block\n";
    if ( $block_strength < $cutoff_for_good_quality_block ) {
      print "This weak block remains $chr $n\n";
    }
    $n++;
  }

  eliminate_flanking_and_intercalated_blocks();
  set_start_of_first_block_to_1_and_end_of_last_block_to_end_of_chromosome();
}

dump_values_for_later_use();

#print STDERR "Still need to implement strong,weak,different strong & strongs compatible!!\n\n";

if ( $total_paternal_recombinations ) {
  $maternal_paternal_ratio = sprintf( "%.3f",
                                      $total_maternal_recombinations /
                                          $total_paternal_recombinations );
}
else {
  $maternal_paternal_ratio = "UNDEF";
}
$total_recombinations =
    $total_maternal_recombinations + $total_paternal_recombinations;

unless ( $genome_length ) { die }
$fraction_of_suspicious =
    sprintf( "%.4f", $total_length_of_suspicious_blocks / $genome_length );
print
"Total length of suspicious blocks:\t$total_length_of_suspicious_blocks\tout of\t$genome_length\tgenome length (\t$fraction_of_suspicious\t)\n";
print SUSPICIOUS
"#note\tTotal length of suspicious blocks:\t$total_length_of_suspicious_blocks\tout of\t$genome_length\tgenome length (\t$fraction_of_suspicious\t)\n";

$fraction_of_uncalled =
    sprintf( "%.4f", $total_length_of_uncalled_blocks / $genome_length );
print
"Total length of uncalled blocks:\t$total_length_of_uncalled_blocks\tout of\t$genome_length\tgenome length (\t$fraction_of_uncalled\t) [note that this length is included in the total of all suspicious blocks]\n";
print SUSPICIOUS
"#note\tTotal length of uncalled blocks:\t$total_length_of_uncalled_blocks\tout of\t$genome_length\tgenome length (\t$fraction_of_uncalled\t) [note that this length is included in the total of all suspicious blocks]\n";

print
"Total maternal and paternal recombinations (and ratio):\t$total_maternal_recombinations\t$total_paternal_recombinations\t(total\t$total_recombinations\t)\t$maternal_paternal_ratio\n";
print OUT
"#note\tTotal maternal and paternal recombinations (and ratio):\t$total_maternal_recombinations\t$total_paternal_recombinations\t(total\t$total_recombinations\t)\t$maternal_paternal_ratio\n";

close SUSPICIOUS;
close OUT;

exit;

#END MAIN BLOCK

sub by_integer {
  $a =~ /(\d+)/;
  $a1 = $1;
  $b =~ /(\d+)/;
  $b1 = $1;
  $a1 <=> $b1;
}

sub output_suspicious_block_to_suspicious_file {
  my $block_number                          = shift;
  my $relation_to_neighboring_strong_blocks = shift;
  $length       = $end{$chr}{$block_number} - $start{$chr}{$block_number};
  $binary_state =
      string_binary_representation_of_inheritance_state(
                                                  $block{$chr}{$block_number} );
  print SUSPICIOUS
"$chr\t$block_number\t$block{$chr}{$block_number}\t$binary_state\t$start{$chr}{$block_number}\t$end{$chr}{$block_number}\t$length\t$relation_to_neighboring_strong_blocks\n";
  $total_length_of_suspicious_blocks += $length;
  if ( $relation_to_neighboring_strong_blocks eq "transitional" ) {
    $total_length_of_uncalled_blocks += $length;
  }
}

#end output_suspicious_block_to_suspicious_file

sub set_start_of_first_block_to_1_and_end_of_last_block_to_end_of_chromosome {

#or maybe better practice just to report blocks extending from first observed marker to last observed marker
  @block_numbers     = sort { $a <=> $b } keys %{ $block{$chr} };
  $last_block_number = $block_numbers[ -1 ];

  #$start{$chr}{1} = $first_observed_marker;
  #$end{$chr}{$last_block_number} = $last_observed_marker;
  $start{$chr}{1} = $positions_of_bounding_markers{$chr}{'low'};
  $end{$chr}{$last_block_number} = $positions_of_bounding_markers{$chr}{'high'};
  $genome_length += $positions_of_bounding_markers{$chr}{'high'} -
      $positions_of_bounding_markers{$chr}{'low'};
}

#end set_start_of_first_block_to_1_and_end_of_last_block_to_end_of_chromosome

sub eliminate_flanking_and_intercalated_blocks {
  @blocks_to_eliminate = sort { $a <=> $b } @blocks_to_eliminate;
  print "Blocks to eliminate\t", join( ", ", @blocks_to_eliminate ), "\n";
  print "Blocks to consider uncalled\t",
      join( ", ", @blocks_to_consider_uncalled ), "\n";
  %new_block = ();
  %new_start = ();
  %new_end   = ();

  #%new_number_of_markers = ();
  #traverse the list of blocks, eliminating those indicated
  $n         = 1;
  $new_index = 1;
  while ( $n <= $number_of_blocks ) {
    if ( grep { $_ eq $n } @blocks_to_eliminate ) {
      output_suspicious_block_to_suspicious_file( $n, "intercalated" );
    }
    elsif ( grep { $_ eq $n } @blocks_to_consider_uncalled ) {
      output_suspicious_block_to_suspicious_file( $n, "transitional" );
    }
    else {

      #print STDERR "Outputting new block $new_index from old index $n\n";
      $new_block{$new_index} = $block{$chr}{$n};    #block_type
      $new_start{$new_index} = $start{$chr}{$n};
      $new_end{$new_index}   = $end{$chr}{$n};

#no longer computing quality statistics directly here (although might re-consider if important to do this iteratively
#$new_number_of_markers{$new_index} = $number_of_markers{$chr}{$n};  #note that I am not going to count eliminated markers as markers supporting the new block
      $new_index++;
    }
    $n++;
  }

  #now merge any adjacent identical blocks
  @block_numbers    = sort { $a <=> $b } keys %new_block;
  $number_of_blocks = scalar( @block_numbers );

#print "New blocks\t",join(", ",@block_numbers)," - number_of_blocks = $number_of_blocks\n";

  $n           = 1;
  $final_index = 1;
  %final_block = ();
  %final_start = ();
  %final_end   = ();

  #%final_number_of_markers = ();
  while ( $n < $number_of_blocks ) {
    $block_now  = $new_block{$n};          #block_type
    $next_block = $new_block{ $n + 1 };    #block_type
                                           #if they are equal, merge
                                           #else output the block
    if ( $block_now == $next_block ) {
      $new_start{ $n + 1 } = $new_start{$n};

#$new_number_of_markers{$n+1} = $new_number_of_markers{$n} + $new_number_of_markers{$n+1};
      print
"Merging new block $n into the next block because they are both $block_now type\n";
    }
    else {

      #print STDERR "Final block $final_index\n";
      $final_block{$final_index} = $block_now;
      $final_start{$final_index} = $new_start{$n};    #block_type
      $final_end{$final_index}   = $new_end{$n};      #block_type
        #$final_number_of_markers{$final_index} = $new_number_of_markers{$n};  #block_type
      $final_index++;
    }
    $n++;
  }

  $final_block{$final_index} = $new_block{$n};
  $final_start{$final_index} = $new_start{$n};    #block_type
  $final_end{$final_index}   = $new_end{$n};      #block_type

#$final_number_of_markers{$final_index} = $new_number_of_markers{$n};  #block_type
#print STDERR "Final block $final_index\n";

  #replace the block hashes with the final_block_hashes
  %{ $block{$chr} } = %final_block;
  %{ $start{$chr} } = %final_start;
  %{ $end{$chr} }   = %final_end;

  #print "Dumper block: " . Dumper( \%block ) . "\n";

  #%{$number_of_markers{$chr}} = %final_number_of_markers;
}

#end eliminate_flanking_and_intercalated_blocks

sub make_header_for_output_files {
  my $header = "";
  $header .= "#parameter\tmode\t$mode\n";
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
  $header .= "#pedigree\t$input_file_header_info{'pedigree'}\n";
  $header .= "#sex\t$input_file_header_info{'sex'}\n";

  #$header .=  "#inheritance_vector_format binary";
  $header .= "#inheritance_state_representation\tarabic_zero_based\n";
  $header .=
"#parameter\tcutoff_for_good_quality_block = $cutoff_for_good_quality_block\n";

#$header .= "#note\tMany header fields hard coded in suttonian-hmm header subroutine so could be wrong if was not changed appropriately\n";
  return $header;
}

sub dump_values_for_later_use {
  $total_paternal_recombinations = 0;
  $total_maternal_recombinations = 0;

  if ( $mode eq "chrX" ) {
    $outfile = "$infile_directory/smoothed_Xblocks_" . timestamp() . ".txt";
  }
  else {
    $outfile = "$infile_directory/smoothed_blocks_" . timestamp() . ".txt";
    $outfile = "$infile_directory/smoothed_blocks.txt";
  }

  open( OUT, ">$outfile" );
  print OUT $header_for_output_files;
  print OUT
"#header\tchromosome\tblock_number\tmost_likely_path_state\tbinary state\tnumber of markers\tblock_start\tblock_end\tblock_length\trecombination meiosis\tlength of crossover interval\tcrossover start location\n";

#print OUT "#number_of_markers_to_be_short\t$number_of_markers_to_be_short\n";
#print OUT "#length_to_be_short\t$length_to_be_short\n";
#print OUT "#absolute_minimum_CNV_or_error_block_length\t$absolute_minimum_CNV_or_error_block_length\n";
#print OUT "#absolute_minimum_CNV_or_error_number_of_markers\t$absolute_minimum_CNV_or_error_number_of_markers\n";
#print OUT "#block_pruning_order\t",join(", ",@block_pruning_order),"\n";

  foreach $chr ( @chromosomes_to_study ) {

#print OUT "#chromosome	block_number\tmost_likely_path_state\tblock_start\tblock_end\tblock_length\trecombination meiosis\tlength of crossover interval\tcrossover start location\n";
#print OUT "#chromosome\t",uc $chr,"\n";
#print OUT "#chromosome\t",$chr,"\n";
    @block_numbers =
        sort { $a <=> $b }
        keys %{ $block{$chr} };    #$block{$chr}{$block_number} = $data[2];

#print "chr = $chr; bn=" . join(",", @block_numbers ) . " ; " . Dumper( $block{$chr} ) . "\n";
    ## RMH: Ignore blocks that don't have defined states
    next
        if ( @block_numbers == 1
             && !defined $block{$chr}{ $block_numbers[ 0 ] } );

    foreach $block_number ( @block_numbers ) {
      $length = $end{$chr}{$block_number} - $start{$chr}{$block_number};
      if ( $block_number > 1 ) {
        $length_of_crossover_interval =
            $start{$chr}{$block_number} - $end{$chr}{ $block_number - 1 };
        $crossover_start_location = $end{$chr}{ $block_number - 1 };
        if ( $chr eq "chrXnonPAR" ) {
          $recombination_meiosis =
              recombination_Xmeiosis( $block{$chr}{ $block_number - 1 },
                                      $block{$chr}{$block_number} );
        }
        else {
          $recombination_meiosis =
              recombination_meiosis( $block{$chr}{ $block_number - 1 },
                                     $block{$chr}{$block_number} );
          if ( $recombination_meiosis =~ /paternal/ ) {
            $total_paternal_recombinations++;
          }
          elsif ( $recombination_meiosis =~ /maternal/ ) {
            $total_maternal_recombinations++;
          }
        }
      }
      else {
        $length_of_crossover_interval = "NA";
        $crossover_start_location     = "NA";
        $recombination_meiosis        = "NA";
      }

      $binary_state =
          string_binary_representation_of_inheritance_state(
                                                  $block{$chr}{$block_number} );
      $number_of_markers = "ND";

      # New code RMH: TODO REMOVE
      if ( 1 ) {
        print OUT
"$chr\t$block_number\t$block{$chr}{$block_number}\t$binary_state\t$number_of_markers\t$start{$chr}{$block_number}\t$end{$chr}{$block_number}\t$length\t";
        print OUT
"$recombination_meiosis\t$length_of_crossover_interval\t$crossover_start_location\n";
      }
      else {

        # Legacy
        print OUT
"$chr\t$block_number\t$block{$chr}{$block_number}\t$start{$chr}{$block_number}\t$end{$chr}{$block_number}\t$length\t";
        print OUT
"$recombination_meiosis\t$length_of_crossover_interval\t$crossover_start_location\n";
      }
    }
  }

  print "Done storing the the data (with timestamp " . timestamp() . ")!\n";
}

#end dump_values_for_later_use

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

sub recombination_meiosis {
  my $recom1 = shift;
  my $recom2 = shift;

  #($recom1,$recom2) = sort ($recom1,$recom2);
  my $hamming_distance_between_the_states =
      $states_are_adjacent{$recom1}{$recom2};
  if ( $hamming_distance_between_the_states == 1 ) {

    #is it maternal or paternal?
    #return "paternal"
    #return "maternal"
    return $parental_transition{$recom1}{$recom2};
  }
  else {
    return "$hamming_distance_between_the_states recombinations";
  }

  #note that currently will never get past this point
  return "uncharacterized_transition";
  if ( ( $recom1 == 5 ) or ( $recom2 == 5 ) ) {
    return "error_transition";
  }
  if ( ( $recom1 == 6 ) or ( $recom2 == 6 ) ) {
    return "compression_transition";
  }
}

#end recombination_meiosis

sub recombination_Xmeiosis {
  my $recom1 = shift;
  my $recom2 = shift;
  ( $recom1, $recom2 ) = sort ( $recom1, $recom2 );
  return "waiting for implementation";
}

#end recombination_meiosis

sub load_file {
  my @inline;
  my $long_filename = shift;
  die unless ( -e $long_filename );
  my $age_of_file = sprintf( "%.1f", -M $long_filename );
  print "Opening $short_filename, which is $age_of_file days old.\n";
  open VAR, $long_filename;
  while ( <VAR> ) {
    chomp;
    if ( /^#/ ) {

      #make hash of header info
      if ( $_ =~ /^\#(\S*)\s+(.*)$/ ) {
        $input_file_header_info{$1} = $2;

        #print "$1\t$2\n";
      }
      else {
        print "$_ is a bad header line!\n";
        die;
      }
      next;
    }

#next if (/^chromosome/);  #skip headerline
#print $_;
#chromosome	block_number	most_likely_path_state	number_of_markers_in_block	block_start	block_end	block_length	recombination meiosis	length of crossover interval	crossover start location

    @data = split( /\t/, $_ );
    unless ( $data[ 0 ] =~ /\S/ ) { shift @inline }
    ; #not sure why getting empty first field; this is a hack to solve porblem without understanding it; something to do with mac return issue
    unless ( $data[ -1 ] =~ /\S/ ) { pop @inline }
    ; #not sure why getting empty first field; this is a hack to solve porblem without understanding it; something to do with mac return issue
    $chr          = $data[ 0 ];
    $block_number = $data[ 1 ];

    #$chr{$block_number} = $chr;
    $block{$chr}{$block_number} = $data[ 2 ];

    #binary state = $data[3]

    #$number_of_markers{$chr}{$block_number} = $data[4];
    $start{$chr}{$block_number} = $data[ 5 ];
    $end{$chr}{$block_number}   = $data[ 6 ];

    #$length{$chr}{$block_number} = $data[6];
    #print STDERR "The genotype at $chr $pos is $genotype.\n";
    #$seen_chr{$chr} = 1;
    #$seen_position{$chr}{$pos} = 1;
    #die;
  }
  close VAR;
}

#end load_file

sub load_previously_compiled_data {
  my $hmmDataCache  = shift;
  my $statDataCache = shift;

  my %dataCache = (
              'block_statistics'              => \%block_statistics,
              'positions_of_bounding_markers' => \%positions_of_bounding_markers
  );

  my $subroutine      = "load_previously_compiled_data";
  my $local_timestamp = timestamp();
  print "[$subroutine] Begin at $local_timestamp.\n";

  die "load_previously_compiled_data(): Could not locate "
      . "$statDataCache file! Either it's missing or it's empty.\n"
      unless ( -s $statDataCache );
  my $age_of_file = sprintf( "%.1f", -M $statDataCache );
  print "\tOpening $statDataCache, which is " . "$age_of_file days old.\n";

  # Read in the data
  my $tmpDataCache = retrieve( $statDataCache );

  # NOTE: This is a very inefficient way to break the hash into
  #       the discrete variables. This will double the memory
  #       requirements doing it this way ( at least temporarily )
  #       TODO: Fix.
  %block_statistics              = %{ $tmpDataCache->{'block_statistics'} };
  %positions_of_bounding_markers =
      %{ $tmpDataCache->{'positions_of_bounding_markers'} };

  undef $tmpDataCache;

  die "load_previously_compiled_data(): Could not locate "
      . "$hmmDataCache file! Either it's missing or it's empty.\n"
      unless ( -s $hmmDataCache );
  $age_of_file = sprintf( "%.1f", -M $hmmDataCache );
  print "\tOpening $hmmDataCache, which is " . "$age_of_file days old.\n";

  # Read in the data
  $tmpDataCache = retrieve( $hmmDataCache );

  # NOTE: This is a very inefficient way to break the hash into
  #       the discrete variables. This will double the memory
  #       requirements doing it this way ( at least temporarily )
  #       TODO: Fix.
  %states_are_adjacent = %{ $tmpDataCache->{'states_are_adjacent'} };
  %parental_transition = %{ $tmpDataCache->{'parental_transition'} };

  # Free up memory -- no need to keep a duplicate
  undef $tmpDataCache;
}

#end load_previously_compiled_data

