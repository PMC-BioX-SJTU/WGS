#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) RSutils_for_HMM_Statistics.pm
##  Author:
##      Gustavo Glusman  <gglusman@systemsbiology.org>
##      Robert M. Hubley   <rhubley@systemsbiology.org>
##  Description:
##      A set of generic functions on range sets.
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
################################################################################
#
## To Do:
#
#

=head1 NAME

RSutils_for_HMM_Statistics - A set of generic functions on range sets

=head1 SYNOPSIS

use RSutils_for_HMM_Statistics;

Usage:


  my $obj = new RSutils_for_HMM_Statistics();

  $obj->Object_Method();

=head1 DESCRIPTION

 This modules is a slight modification to the work of Gustavo Glusman 
 ( "RSUtils.pm" ). 

 All methods use zero-based start, one-based end.

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2011 Institute for Systems Biology

=head1 AUTHORS

  Gustavo Glusman  <gglusman@systemsbiology.org>
  Jared Roach      <jroach@systemsbiology.org>

=head1 OBJECT METHODS

=cut

##
## Module dependence
##
package RSutils;
use strict;

sub new {
  my $package = shift;

  my $obj = {};
  bless $obj, $package;
  return $obj;
}

##-------------------------------------------------------------------------##

=head2 RSsort()

  Use: my $sortedSet = $obj->RSsort( $RSref );

=cut

##-------------------------------------------------------------------------##
sub RSsort {

  #sorts a range set by starting position and then by ending position
  #overlapping ranges are not collapsed, labels are not modified
  my ( $self, $RSref ) = @_;
  my ( @sorted, $start, $end, @starts, @ends, $i );

  foreach my $item ( @{$RSref} ) {
    ( $start, $end ) = split /\t/, $item;
    push @starts, $start;
    push @ends,   $end;
  }
  foreach $i (
    sort {
             $starts[ $a ] <=> $starts[ $b ]
          || $ends[ $a ] <=> $ends[ $b ]
    } ( 0 .. $#{$RSref} )
      )
  {
    push @sorted, $$RSref[ $i ];
  }

  return \@sorted;
}

##-------------------------------------------------------------------------##

=head2 RSunion()

  Use: my $unionSet = $obj->RSunion( $RSreflist );

=cut

##-------------------------------------------------------------------------##
sub RSunion {

  #creates the union of range sets
  #can be used on just one range set to collapse overlapping ranges
  #the resulting set is returned sorted
  #repeated labels are condensed; labels are sorted per range
  my ( $self, @RSreflist ) = @_;
  my (
       $RSref, $item,      $start,   $end,    @rest,
       %info,  $laststart, $lastend, @sorted, @info
  );
#print "bbbb\t".join("\t",@RSreflist)."\n";
  foreach $RSref ( @RSreflist ) {
    print "bbbb\t".join(",",@$RSref)."\n";
	foreach $item ( @{$RSref} ) {
      ( $start, $end, @rest ) = split /\t/, $item;
      push @{ $info{$start}{$end} }, @rest;
    }
  }

  $laststart = "x";
  foreach $start ( sort { $a <=> $b } keys %info ) {
    foreach $end ( sort { $a <=> $b } keys %{ $info{$start} } ) {
      if ( $laststart eq "x" ) {

        #start first combined range
        $laststart = $start;
        $lastend   = $end;
        @info      = @{ $info{$start}{$end} };
      }
      elsif ( $start > $lastend ) {

        #finish current combined range
        my %labels;
        foreach ( @info ) {
          $labels{$_}++;
        }
        push @sorted, join( "\t", $laststart, $lastend, sort keys %labels );

        #...and start a new one
        $laststart = $start;
        $lastend   = $end;
        @info      = @{ $info{$start}{$end} };
      }
      else {

        #extend current combined range
        $lastend = $end if $end > $lastend;
        push @info, @{ $info{$start}{$end} };
      }
    }
  }
  if ( $laststart ne "x" ) {

    #finish last combined range
    my %labels;
    foreach ( @info ) {
      $labels{$_}++;
    }
    push @sorted, join( "\t", $laststart, $lastend, sort keys %labels );
  }
  return \@sorted;
}

##-------------------------------------------------------------------------##

=head2 RSintersection()

  Use: my $intersecSet = $obj->RSintersection( $RSreflist );

=cut

##-------------------------------------------------------------------------##
sub RSintersection {

  #calculates the intersection of any number of range sets
  #if given just one set, behaves like RSunion
  #currently, it returns just the ranges, losing all extra information
  my ( $self, $joint, @RSreflist ) = @_;
  my (
       $RSref, @newjoint, $rangen, $lastrangen, $item,
       $start, $end,      $rstart, $rend,       $sorted
  );

  $joint = $self->RSunion( $joint );
  foreach $RSref ( @RSreflist ) {
    @newjoint   = ();
    $sorted     = $self->RSunion( $RSref );
    $lastrangen = 0;
    foreach $item ( @{$sorted} ) {
      ( $start, $end ) = split /\t/, $item;
      for ( $rangen = $lastrangen ; $rangen <= $#{$joint} ; $rangen++ ) {
        ( $rstart, $rend ) = split /\t/, $$joint[ $rangen ];
        last if $rstart >= $end;
        $lastrangen = $rangen;
        next if $start > $rend;
        push @newjoint,
            join( "\t",
                  $start < $rstart ? $rstart : $start,
                  $end < $rend     ? $end    : $rend );
      }
    }
    $joint = \@newjoint;
  }
  return $joint;
}

##-------------------------------------------------------------------------##

=head2 RSintersection2()

  Use: my $intersecSet = $obj->RSintersection2( $joint, $RSreflist );

=cut

##-------------------------------------------------------------------------##
sub RSintersection2 {

  #calculates the intersection of any number of range sets
  #assumes all range sets have already been passed by RSunion
  #currently, it returns just the ranges, losing all extra information
  my ( $self, $joint, @RSreflist ) = @_;
  my (
       $RSref, @newjoint, $rangen, $lastrangen, $item,
       $start, $end,      $rstart, $rend,       $sorted
  );

  #$joint = $self->RSunion($joint);
  foreach $sorted ( @RSreflist ) {
    @newjoint   = ();
    $lastrangen = 0;
    foreach $item ( @{$sorted} ) {
      ( $start, $end ) = split /\t/, $item;
      for ( $rangen = $lastrangen ; $rangen <= $#{$joint} ; $rangen++ ) {
        ( $rstart, $rend ) = split /\t/, $$joint[ $rangen ];
        last if $rstart >= $end;
        $lastrangen = $rangen;
        next if $start > $rend;
        push @newjoint,
            join( "\t",
                  $start < $rstart ? $rstart : $start,
                  $end < $rend     ? $end    : $rend );
      }
    }
    $joint = \@newjoint;
  }
  return $joint;
}

##-------------------------------------------------------------------------##

=head2 RSsubtraction()

  Use: my $remainingSet = $obj->RSsubtraction( $joint, $RSreflist );

=cut

##-------------------------------------------------------------------------##
sub RSsubtraction {

#calculates what's left of a range set when subtracting from it any number of range sets
#if given just one set, behaves like RSunion
  ##modified ranges retain labels - some labels may have become irrelevant
  my ( $self, $joint, @RSreflist ) = @_;
  my (
       $RSref,  $rangen, $lastrangen,  $item,   $start, $end,
       $rstart, $rend,   @replacement, $sorted, @labels
  );
	print "ccc\t".join("\t",@$joint)."\n";
  $joint = $self->RSunion( $joint );print "aaaaaaa\t".join("\t",@$joint)."\n";
  foreach $RSref ( @RSreflist ) {
    $sorted     = $self->RSunion( $RSref );
    $lastrangen = 0;
    foreach $item ( @{$sorted} ) {
      ( $start, $end ) = split /\t/, $item;
      for ( $rangen = $lastrangen ; $rangen <= $#{$joint} ; $rangen++ ) {
        ( $rstart, $rend, @labels ) = split /\t/, $$joint[ $rangen ];
        last if $rstart >= $end;
        $lastrangen = $rangen;
        next if $start >= $rend;
        @replacement = ();
        if ( $rstart < $start ) {
          push @replacement,
              join( "\t", $rstart, $rend < $start ? $rend : $start, @labels );
        }
        if ( $rend > $end ) {
          push @replacement, join( "\t", $end, $rend, @labels );
        }
        splice @$joint, $rangen, 1, @replacement;
        $rangen--;
      }
    }
  }
  return $joint;
}

##-------------------------------------------------------------------------##

=head2 RSsuppression()

  Use: my $set = $obj->RSsuppresion( $joint, $RSreflist );

=cut

##-------------------------------------------------------------------------##
sub RSsuppression {

  #calculates what's left of a range set when removing from it any range
  #that overlaps ranges in any number of range sets
  #if given just one set, behaves like RSunion
  my ( $self, $joint, @RSreflist ) = @_;
  my (
       $RSref, $rangen, $lastrangen, $item, $start,
       $end,   $rstart, $rend,       $sorted
  );

  $joint = $self->RSunion( $joint );
  foreach $RSref ( @RSreflist ) {
    $sorted     = $self->RSunion( $RSref );
    $lastrangen = 0;
    foreach $item ( @{$sorted} ) {
      ( $start, $end ) = split /\t/, $item;
      for ( $rangen = $lastrangen ; $rangen <= $#{$joint} ; $rangen++ ) {
        ( $rstart, $rend ) = split /\t/, $$joint[ $rangen ];
        last if $rstart >= $end;
        $lastrangen = $rangen;
        next if $start >= $rend;
        splice @$joint, $rangen, 1;
      }
    }
  }
  return $joint;
}

##-------------------------------------------------------------------------##

=head2 RSsegmentation()

  Use: my $set = $obj->RSsegmentation( $RSreflist );

=cut

##-------------------------------------------------------------------------##
sub RSsegmentation {

#segments a redundant range set (or group of range sets)
#each resulting range in the set is undivided
#it is also labeled with the identifiers of the original ranges covering it
#if an id is mentioned more than once from the same starting point, the longer version is used
  my ( $self, @RSreflist ) = @_;
  my ( %edge, %end, @ids, @result );

  foreach my $RSref ( @RSreflist ) {

    #print STDERR "\tRSref\t$RSref\n";
    foreach my $item ( @{$RSref} ) {

      #print STDERR "\titem\t$item\n";
      my ( $start, $end, @rids ) = split /\t/, $item;
      $edge{$start}++;
      $edge{$end}++;
      foreach my $id ( @rids ) {
        if ( defined $end{$start}{$id} ) {
          $end{$start}{$id} = $end if $end > $end{$start}{$id};
        }
        else {
          $end{$start}{$id} = $end;
        }
      }
    }
  }
  my @edges = sort { $a <=> $b } keys %edge;

  foreach my $i ( 0 .. $#edges - 1 ) {

    my $edge = $edges[ $i ];

    #print STDERR "$i\t$edge\n";
    while ( my ( $id, $end ) = each %{ $end{$edge} } ) {
      for ( my $j = $i ; $j < $#edges ; $j++ ) {
        last if $edges[ $j ] >= $end;
        $ids[ $j ]{$id}++;
      }
    }
    my @rids = sort keys %{ $ids[ $i ] };
    if ( @rids ) {
      push @result, join( "\t", $edge, $edges[ $i + 1 ], @rids );
    }
  }

  #print STDERR join("\n",@result),"\n";

  return \@result;
}

##-------------------------------------------------------------------------##

=head2 SRSsort()

  Use: my $set = $obj->SRSsort( $RSreflist );

=cut

##-------------------------------------------------------------------------##
sub SRSsort {

#sorts a spliced range set by starting position of the first segment and then by ending position of first segment
#overlapping ranges are not collapsed, labels are not modified
  my ( $self, $RSref ) = @_;
  my ( @sorted, $start, $end, @starts, @ends, $i );

  foreach my $item ( @{$RSref} ) {
    ( $start, $end ) = split /\t/, $item;
    ( $start ) = split /,/, $start;
    ( $end )   = split /,/, $end;
    push @starts, $start;
    push @ends,   $end;
  }
  foreach $i (
    sort {
             $starts[ $a ] <=> $starts[ $b ]
          || $ends[ $a ] <=> $ends[ $b ]
    } ( 0 .. $#{$RSref} )
      )
  {
    push @sorted, $$RSref[ $i ];
  }

  return \@sorted;
}

##-------------------------------------------------------------------------##

=head2 SRScluster()

  Use: my $result = $obj->SRScluster( $RSreflist );

=cut

##-------------------------------------------------------------------------##
sub SRScluster {

#clusters spliced ranges in one or more sets, by range overlap (e.g. cluster transcripts by exon overlap)
  my ( $self, @RSreflist ) = @_;
  my ( %info );
  foreach my $RSref ( @RSreflist ) {
    foreach my $item ( @{$RSref} ) {
      my ( $starts, $ends, $id ) = split /\t/, $item, 3;
      my @starts = split /,/, $starts;
      my @ends   = split /,/, $ends;
      my $start = $starts[ 0 ];
      my $end   = $ends[ $#ends ];

      if ( defined $info{$starts}{$ends} ) {

#we already saw such an exact spliced range definition, simply append id info to the old one
        $info{$starts}{$ends}{'id'} .= ",$id";
        next;
      }

      #register the range
      $info{$starts}{$ends}{'start'} = $start;
      $info{$starts}{$ends}{'end'}   = $end;
      $info{$starts}{$ends}{'id'}    = $id;

      #compare to previous ranges...
      foreach my $pstarts ( keys %info ) {
        foreach my $pends ( keys %{ $info{$pstarts} } ) {
          next if $info{$pstarts}{$pends}{'id'} eq $id;
          next if $start + 1 > $info{$pstarts}{$pends}{'end'};
          next if $end < $info{$pstarts}{$pends}{'start'} + 1;

          #there is genomic overlap, test exons...
          my @pstarts = split /,/, $pstarts;
          my @pends   = split /,/, $pends;
          foreach my $i ( 0 .. $#starts ) {
            foreach my $j ( 0 .. $#pstarts ) {
              last if $pstarts[ $j ] + 1 > $ends[ $i ];
              next if $pends[ $j ] < $starts[ $i ] + 1;

              #there is exon overlap!
              my $me             = join( "\t", $starts, $ends );
              my $slavestarts    = $pstarts;
              my $slaveends      = $pends;
              my $previousMaster = $info{$slavestarts}{$slaveends}{'master'};
              last if $previousMaster eq $me;
              $info{$slavestarts}{$slaveends}{'master'} = $me;
              while ( $previousMaster && ( $previousMaster ne $me ) ) {
                ( $slavestarts, $slaveends ) = split /\t/, $previousMaster;
                $previousMaster = $info{$slavestarts}{$slaveends}{'master'};
                $info{$slavestarts}{$slaveends}{'master'} = $me;
              }
            }
          }
        }
      }
    }
  }

  my %clusters;
  foreach my $starts ( keys %info ) {
    foreach my $ends ( keys %{ $info{$starts} } ) {

      #identify ultimate master
      my $master   = $info{$starts}{$ends}{'master'};
      my $umaster  = join( "\t", $starts, $ends );
      my @remaster = ( $umaster );
      while ( $master ) {
        push @remaster, $master;
        my ( $mstarts, $mends ) = split /\t/, $master;
        if ( $master = $info{$mstarts}{$mends}{'master'} ) {
          last if $master eq $umaster;
          $umaster = $master;
        }
      }

      #remaster to ultimate master and register in cluster
      foreach my $id ( @remaster ) {
        my ( $idstarts, $idends ) = split /\t/, $id;
        $info{$idstarts}{$idends}{'master'} = $umaster;
        $clusters{$umaster}{$id} = 1;
      }
    }
  }

  #merge exons
  my @result;
  foreach my $master ( keys %clusters ) {
    my ( @ranges, @ids );
    foreach my $id ( keys %{ $clusters{$master} } ) {
      my ( $starts, $ends ) = split /\t/, $id;
      my @starts = split /,/, $starts;
      my @ends   = split /,/, $ends;
      foreach my $i ( 0 .. $#starts ) {
        push @ranges, join( "\t", $starts[ $i ], $ends[ $i ] );
      }
      push @ids, $info{$starts}{$ends}{'id'};
    }
    my $joint = $self->RSunion( \@ranges );
    my ( @starts, @ends );
    foreach my $range ( @{$joint} ) {
      my ( $start, $end ) = split /\t/, $range;
      push @starts, $start;
      push @ends,   $end;
    }
    push @result,
        join( "\t",
              join( ",", @starts ),
              join( ",", @ends ),
              join( ",", @ids ) );
  }
  return \@result;
}

1;
