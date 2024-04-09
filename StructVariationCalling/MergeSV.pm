# Copyright (c) 2010, 2011, 2012 Genome Research Ltd.
#
# Author: Kim Wong <kw10@sanger.ac.uk>
# 
# This file is part of SVMerge.
# 
# SVMerge is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

package MergeSV;

use warnings;
use strict;
use Carp;
use Set::IntSpan;

=head1 NAME 

MergeSV

=head1 AUTHOR

kw10@sanger.ac.uk

=head1 DESCRIPTION

Methods for comparing, manipulating, and merging coordinate sets

=head1 METHODS

=head2 new

	Arguments   : none
	Description : creates a new MergeSV object

=cut

sub new {
	my ($class,%arg) = @_;
	my $self={};
	bless($self,$class || ref($class));
	return $self;
}

=head2 coordCheck

	Arg[1]      :  string with coordinates [1000-2000 or "1000 2000"]
	Example     :  coordCheck($line)
	Description :  checks coordinate formatting and start < end
	Returns     :  nothing if passes check

=cut
sub coordCheck {
	my ($self,$coord) = @_;
	my ($s,$e) = ($1,$2) if $coord =~ /^(\d+)-(\d+)$/ || $coord =~ /^(\d+)\s+(\d+)$/;
	if (!$e) {
		croak "Coordinates not formatted properly $coord\n";
	}
	elsif ($e<$s) {
		croak "End coordinate is less than the Start coordinate: start $s, end $e\n";
	}
}

=head2 merge_spans

	Arg[1]      :  tab-delimited coordinate file
	Arg[2]      :  use intersect (I) or union (U) of overlapping coordinates
	Arg[3]      :  buffer added to coordinates before overlap check
	Arg[4]      :  annotation to use in column 4
	Example     :  merge_spans("coords.tab",I,1000,"DEL_merge_NA18506")
	Description :  finds overlapping coordinates from a single file
	Returns     :  prints out merged set of tab-delimited SV calls; returns nothing

=cut
sub merge_spans {
	my ($self,$file,$type,$buffer,$append,$addsize) = @_;
	my $sets;
	my $sizes;
	# Must sort by coordinates first!
	open F, "sort -k1,1 -k2,2g -k3,3g $file |" or die;
	my $count=0;
	while (<F>) {
		$count++;
		print STDERR "working on $count\n" if ($count%10000==0) ;
		my ($chrom,$s,$e,$annot)=split /\s+/, $_;
		$self->coordCheck("$s-$e");
		next if !$chrom || !$s || !$e;
		my $size = $annot =~ /\S+_\S+_\S+_(\d+)$/ ? $1 : ''; # format: TYPE_METHOD_STRAIN_SIZE
		$size = 0 if !$size;
		if (!$sets->{$chrom}) {
			push @{$sets->{$chrom}}, "$s-$e";
			push  @{$sizes->{$chrom}}, $size;
			next;
		}
		my $found=0;
		# Now compare current range to previous spans
		my $j='';
		foreach ($j=$#{$sets->{$chrom}}; $j>=-1; $j--) {
			if ($type eq 'I') {
				last if ($sets->{$chrom}->[$j] =~ /-(\d+)/ && $1 < $s-1);
			}
			# if a buffer is provided for the overlap calculation:
			my $s_tmp = $s;
			if ($type eq 'U') {
				$s_tmp -=$buffer if $buffer;
			}
			my $set1 = $sets->{$chrom}->[$j];
			my ($min1,$max1) = ($1,$2) if $set1 =~ /^(\d+)-(\d+)/;
			my $int=1 if ($min1<=$s_tmp && $max1>=$e) ||
					($min1<=$s_tmp && $max1>=$s_tmp) ||
					($min1<=$e && $max1>=$e) ;
			# Replace the current span if there is overlap
			if ( $int ) {
				if ($type eq 'I') {
					my $min= $min1 > $s_tmp ? $min1 : $s_tmp;
					my $max= $max1 > $e ? $e : $max1;
					if (($max1-$min1<100 && $e-$s_tmp>=100 && $max-$min<75)||
						($max1-$min1>=100 && $e-$s_tmp<100 && $max-$min<75)) {
						# outer coordinates if low overlap and one large, one small
						my $min=$min1 > $s_tmp ? $s_tmp : $min1;
						my $max=$max1 > $e ? $max1 : $e;
						$sets->{$chrom}->[$j] = "$min-$max";
						$sizes->{$chrom}->[$j] .= ";$size" if $size > 0;
						$found++;
						last;

					}
					elsif ($max-$min<75 && $s-$min1>50 && $e-$max1>50) {
						# outer coordinates if low overlap and small coord. span
						if ($e-$s+1<200 && $max-$min+1<200) {
							my $min=$min1 > $s_tmp ? $s_tmp : $min1;
							my $max=$max1 > $e ? $max1 : $e;
							$sets->{$chrom}->[$j] = "$min-$max";
							$sizes->{$chrom}->[$j] .= ";$size" if $size > 0;
							$found++;
							last;
						}
						# keep 2 separated calls if one or both are >=200bp spans
						#  adjust coords (to make no overlap)
						else {
							my $min=$min1 < $s_tmp ? $min1 : $s_tmp;
							my $max=$min1 < $s_tmp ? $s_tmp : $min1;
							$s=$max1 < $e ? $max1 : $e;
							$e=$max1 < $e ? $e : $max1;
							if ($max==$s) {
								$max-=1;
								$s+=1;
							}
							$sets->{$chrom}->[$j] = "$min-$max";
							push @{$sets->{$chrom}}, "$s-$e";
							push  @{$sizes->{$chrom}}, $size ;
							$found++;
							last;
						}
					}
					else {
						$sets->{$chrom}->[$j] = "$min-$max"; #get the minimal coords
						$sizes->{$chrom}->[$j] .= ";$size" if $size > 0;
						$found++;
						last;
					}
				}
				elsif ($type eq 'U') {
					my $end = $e > $max1 ? $e : $max1;
					$sets->{$chrom}->[$j] = "$min1-$end" if $type eq 'U'; #get the maximum coords
					$sizes->{$chrom}->[$j] .= ";$size" if $size > 0;
					$found++;
					last;
				}
			}
		}
		if (!$found) {
			push @{$sets->{$chrom}}, "$s-$e";
			push  @{$sizes->{$chrom}}, $size ;
		}
			
	}
	close F;
	# Note: if two spans don't intersect but are adjacent, spans will merge them here
	# eg:				$sets->{$chrom}->[$_] = $int if $type eq 'I'; #get the minimal coords

	# 1       8102001 8103000 1
	# 1       8103001 8104000 1
	# will merge to 1 8102001 8104000
	if ($file =~ /(\S+).tab$/) {
		$file = $1;
	}
	open M, ">$file.merged.tab" or croak "Can't open $file.merged.tab for writing\n";
	foreach my $c (sort keys %{$sets}) {
			my $list = \@{$sets->{$c}};
		    my $printline = $self->print_spans($c,$list,$append,$addsize,$sizes) if $list;
			print M $printline;
	}
	close M; 
	return;
}


=head2 diff_spans

	Arg[1]      :  tab-delimited coordinate file
	Arg[2]      :  tab-delimited coordinate reference file
	Arg[3]      :  buffer added to coordinates before overlap check
	Arg[4]      :  annotation to use in column 4
	Example     :  diff_spans("coords.tab","coords.tab",1000,"DEL_merge_NA18506")
	Description :  finds non-overlapping coordinates unique to file in Arg[1]
	Returns     :  prints out a unique set of tab-delimited SV calls; returns nothing

=cut
sub diff_spans {
	my ($self,$file,$file2,$buffer,$append,$addsize) = @_;
	open F2, "<$file2" or die "Can't open $file2\n";
	my %refspans;
	while (<F2>) {
		my ($c,$s,$e) = $self->coord_parse($_);
		$self->coordCheck("$s-$e");
		push @{$refspans{$c}}, "$s-$e";
	}
	close F2;
	my %refobj=();
	for (keys %refspans) {
		$refobj{$_} = new Set::IntSpan @{$refspans{$_}};
	}
	open F1, "<$file" or die "Can't open $file\n";
	my %coords;
	my %annot;
	while (<F1>) {
		my ($c,$s,$e,$annot) = $self->coord_parse($_);
		$self->coordCheck("$s-$e");
		my $s2 = $s;
		$s2-=1 unless $s2==1;
		my $e2 = $e+1;
		my $set = $refobj{$c};
		if ($set) {
			my $int = intersect $set "$s2-$e2";
			next if $int;
			push @{$coords{$c}}, "$s-$e";
		}
		else {
			push @{$coords{$c}}, "$s-$e";
		}
		$annot{$c}{"$s"}[0]=$annot if $annot;
		$annot{$c}{"$s"}[1]=$e if $annot;

	}
	close F1;
	$file = `basename $file`;
	if ($file =~ /(\S+).tab$/) {
		$file = $1;
	}
	open M, ">$file.unique.tab" or croak "Can't open $file.unique.tab for writing\n";
	foreach my $chr (sort keys %coords ){
		my $set = \@{$coords{$chr}};
		my $printline = $self->print_spans($chr,$set,$append,$addsize,'',%annot);
		print M $printline;
	}
	close M;
	#return;

}


=head2 splitBySize

	Arg[1]      :  tab-delimited coordinate file
	Arg[2]      :  cutoff size (bp)
	Example     :  splitBySize("file.tab",100)
	Description :  divides coordinate spans into *.small.tab and *.large.tab files
	Returns     :  prints small and large tab-delimited SV calls; returns nothing

=cut
sub splitBySize { 
	my ($self,$file,$cutoff) = @_;
	open F, "<$file" or croak "Can't open file $file for reading\n";
	my $name = $file;
	$name =~ s/\.tab$//;
	open S, ">$name.small.tab" or croak "Can't open $name.small.tab for writing.\n";
	open L, ">$name.large.tab" or croak "Can't open $name.large.tab for writing.\n";
	while (<F>) {
		if (/^\S+\s+(\d+)\s+(\d+)/) {
			if ($2-$1 < $cutoff) {
				print S $_;
			}
			else {
				print L $_;
			}
		}
		else {
			croak "$file not in BED format:\n$_\n";
		}
	}
	close S;
	close L;
	close F;
	return;

}


sub print_spans {
	my ($self,$c,$spans,$append,$addsize,$calcsizes,%annot) = @_;
	#my ($c,$set,$append,$addsize,$calcsizes,%annot) = @_;
	#my @spans = $set->spans();
	
	my @spans = @$spans;
	my $print = '';

	for (0..$#spans) {
		my $calcsize = $calcsizes->{$c}->[$_] if $calcsizes;
		my ($crd,$end) = $spans[$_]=~ /-/ ? split "-", $spans[$_] : ($spans[$_]-1,$spans[$_]);
		die "$crd,$end,$spans[$_]\n" if !$crd||!$end;
		$print .= "$c\t$crd\t$end";
		my $tmp= "$c\t$crd\t$end";
		my $size = $end - $crd + 1;
		if ($calcsize && $calcsize =~ /;/) {
			my @val=split ";", $calcsize;
			my $ave = 0;
			for (@val) {
				$ave+=$_ if $_;
			}
			$ave = int(0.5+($ave/(scalar @val)));
			$calcsize = $ave;
		}
		my $found = '';
		if (%annot) {
			print STDERR "$c,$crd,$annot{$c}{$crd}\n" if 
				!$c || !$crd || !$annot{$c}{$crd};
			if ($annot{$c}{$crd}[1] == $end) {
				$found = 1;
				$print .= "\t";
				$print .= $annot{$c}{$crd}[0];
			}
		}
		if (!$found && $append) {
			$print .= "\t$append";
			$print .= "_$size" if $addsize;
			$print .= "_$calcsize" if $calcsize;
		}
		$print .= "\n";
	}
	return $print;

}

sub coord_parse {
	my ($self,$line) =@_;
	if ($line =~ /^[chr]*(\S+):(\d+)-(\d+)/) {
		my ($s,$e) = $self->cleanCoord($2,$3);
		return ($1,$s,$e);
	}
	elsif ($line =~ /^[chr]*(\S+)\s+(\d\S*)\s+(\d\S*)\s*(\S+)*/) {
		my ($s,$e) = $self->cleanCoord($2,$3);
		return ($1,$s,$e,$4) if $4;
		return ($1,$s,$e);
	}
	else {
		croak "Unrecognized format: $line\n";
	}
}

sub cleanCoord {
	my ($self,@coords) = @_;
	for (@coords) {
		s/,//g;
		s/chr//gi;
	}
	return @coords;
}


1;
