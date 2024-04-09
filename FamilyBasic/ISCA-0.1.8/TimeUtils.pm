#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) TimeUtils.pm
##  Author:
##      Jared Roach        <jroach@systemsbiology.org>
##      Robert M. Hubley   <rhubley@systemsbiology.org>
##  Description:
##      A set of generic functions for creating timers, and displaying
##      timestamps.
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

TimeUtils

=head1 SYNOPSIS

use TimeUtils;

Usage:

  TimeUtils::timestamp();

=head1 DESCRIPTION

A set of generic functions for creating timers, and displaying
timestamps.      

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
package TimeUtils;
use strict;
use Data::Dumper;
use Carp;
require Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

@ISA    = qw(Exporter);
@EXPORT = qw(timestamp elapsedTime);

use constant FMT_YYYYMMDD_TIME => 1;

##
## Package Variables
##
my %timeHistory = ();

##-------------------------------------------------------------------------##

=head2 timestamp()

  Use: my $timeStr = TimeUtils::timestamp( [$report_back] );

  Report the time stamp as a scalar $month.$day.$hours.$min.$sec or
  if $report_back is set to something ( other than "text" ) return
  the string formated as:

         01, 01, 2011, 03:38:20

  Or if $report_back is set to "text" print the time to STDERR
  formated like:

         The current date is 01-01-2011  03:38:20

=cut

##-------------------------------------------------------------------------##
sub timestamp {
  my $report_back = shift;
  my ( $sec, $min, $hours, $day, $month, $year ) =
      ( localtime )[ 0, 1, 2, 3, 4, 5 ];
  $month = sprintf( "%02d", $month + 1 );
  $day   = sprintf( "%02d", $day );
  $hours = sprintf( "%02d", $hours );
  $min   = sprintf( "%02d", $min );
  $sec   = sprintf( "%02d", $sec );
  if ( $report_back ) {

    if ( $report_back eq "text" ) {
      return "$month, $day, " . ( $year + 1900 ) . ", $hours:$min:$sec";
    }
    elsif ( $report_back == FMT_YYYYMMDD_TIME ) {
      return "" . ( $year + 1900 ) . "$month$day $hours:$min:$sec";
    }
    else {
      printf STDERR "The current date is %02d-%02d-%04d\t%02d:%02d:%02d\n",
          $month, $day, $year + 1900, $hours, $min, $sec;
    }
  }
  return $month . $day . $hours . $min . $sec;
}

##-------------------------------------------------------------------------##
## Use: my $string = elapsedTime( $index );
##
##   Returns
##
##      Great little utility for measuring the elapsed
##      time between one or more sections of perl code.
##
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 elapsedTime()

  Use: my $timeStr = TimeUtils::elapsedTime( $index );

  Keep a set of time indices.  The first time an index is used the
  time is stored in the index and the function returns 0.  Subsequent
  calls with the index will return the difference between the previously
  recorded time and the current time.  The difference is returned as a
  string formatted as "hours:minutes:seconds".

=cut

##-------------------------------------------------------------------------##
sub elapsedTime {
  my ( $timeHistoryIdx ) = @_;
  if ( defined $timeHistory{$timeHistoryIdx} ) {
    my $diffTime = time - $timeHistory{$timeHistoryIdx};
    $timeHistory{$timeHistoryIdx} = time;
    my $Min = int( $diffTime / 60 );
    $diffTime -= $Min * 60;
    my $Hours = int( $Min / 60 );
    $Min -= $Hours * 60;
    my $Sec = $diffTime;
    return "$Hours:$Min:$Sec";
  }
  else {
    $timeHistory{$timeHistoryIdx} = time;
    return 0;
  }
}

#### Return Value of Package Module
1;
