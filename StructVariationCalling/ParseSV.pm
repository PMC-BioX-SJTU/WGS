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

package ParseSV;

use strict;
use warnings;
use Carp;
use DBI;
use Data::Dumper; 
use File::Basename;
use lib dirname(__FILE__).'/../lib/';
use Set::IntSpan;
use ExonerateParse;

sub new {
	my ($class,%arg) = @_;
	my $self={};
	bless($self,$class || ref($class));
	return $self;
}

sub findSplit {
	my ($self,$hits) =  @_;
	# first check for contigs with more than one hit
	#
	my %seen=();
	for (@$hits) {
		chomp;
		my $line = $_;
		push @{$seen{$1}}, $line if $line =~ /^(\S+).+/;

	}
	my @res;
	for (keys %seen) {
		if (scalar @{$seen{$_}} > 1) {
			my $res = ($self->parseSplit($seen{$_}));
			push @res, @$res if $res;
		}
	}
	return \@res;

}

sub matchType {
	my ($self,$lines,$type,$svcoords,$svsize) = @_;
	my $size = 0; # actual SV size (best hit)
	my $size2 = 0; # SV size (second best hit)
	my $diff=0; # size diff between best hit and predicted
	my $found = 0; # second best hit found that doesn't match SV type
	my $ovr = 0; # best overlap size
	my ($invbuffer,$glbuffer,$delinsbuffer) = (500,1000,100); # if no Match, then search this region
	my $bestmatch; 
	my @secondmatch;
	my $matchtype;
	my $gapfound;
	my $gapovrlpsize;
	my ($c1,$c2) = split '-', $svcoords;
	$c1 =~ s/\S+://;
	foreach my $line (@$lines) {
		my @cols = split /\s+/, $line; 
		my $sizediff = abs($cols[3]-$svsize) unless $svsize eq 'mergeINS';  # expected size vs. actual size
		
		# Overlap check
		my $crds = ($c1-1)."-".($c2+1);
		my $set = new Set::IntSpan "$crds"; # predicted coords
		die "$line\n" if $cols[1]=~/-/ || $cols[2]=~/-/;
		if ($cols[1]>$cols[2]) {
			print $line;
			die;
		}
		my $ovrlp = intersect $set "$cols[1]-$cols[2]";
		my $ovrlpsize = $2-$1 if $ovrlp && $ovrlp =~ /(\S+)-(\S+)/;
		$ovrlpsize = 1 if !$ovrlpsize && $ovrlp;
		if ($cols[4] =~ /$type/ && $cols[4] !~ /GAP/) { # match based on SV type first; 
		# DELINS will match DEL OR INSl largeINS will match INS; INV will match INVCOMPLEX
			if ($type eq 'INS' && (!$size || ($size && $cols[3]>$size) || ($matchtype && $matchtype eq 'Complex' ))){ # For Ins, match based on size
				#my $buffer = $cols[4] =~ /internalm/ ? 100 : 50;
				my $buffer = 100;
				my ($c1,$c2) = ($cols[1]-$buffer,$cols[2]+$buffer);
				my $ovrlp2 = intersect $set "$c1-$c2";

				if ( $ovrlp2 ) {
					if ($bestmatch && ( $cols[4] !~ /internal/ && $cols[4] ne 'INS')) {# complex match 
						push @secondmatch, $line;
						$size = $cols[3];
					} 
					else {
						push @secondmatch, $bestmatch if $bestmatch;
						$bestmatch = $line;
						$size = $cols[3];
						$ovrlp=$ovrlp2;
						$matchtype = ($cols[4] eq $type || $cols[4] =~ /internal/) ? 'Match' : 'Complex';
					}
				}
				else {
					push @secondmatch, $line;
					$size = $cols[3];
				}
			}
			elsif ($type eq 'INV') { # for inversions, 'match' must be closer in size -lots of false inversions
				if (($cols[3]<($c2-$c1) && $cols[3]/($c2-$c1)>=0.3) ||  
					($cols[3]>($c2-$c1) && ($c2-$c1)/$cols[3]>=0.3)) {
					if ($size==0 || $diff==0 || ($diff && $sizediff < $diff) ||$cols[4] =~ /INVCOMPLEX/ ) {
						my @coords = split /\s+/, $line;
						my $index=0;
						my @svcoord = ($c1,$c2);
						for ($coords[1],$coords[2]) {
							my $coord=$_;
							my $range1=($coord-1000)."-".($coord+1000);
							my $range2 = ($svcoord[$index]-1000)."-".($svcoord[$index]+1000);
							my $ovr1 = $self->overlapCheck("$range1","$range2");
							if ($ovr1 && (!$bestmatch || ($bestmatch && $bestmatch !~ /INVCOM|INVDEL|INVINS/))) {
								push @secondmatch, $bestmatch if $bestmatch;
								$bestmatch = $line;
								$size=$cols[3];
								$diff=$sizediff;
								$matchtype = ($cols[4] eq $type || $cols[4] =~ /internal/) ? 'Match' : 'Complex';
								last;
							}
							$index++;
						}
					}
					else {
						push @secondmatch, $line;
					}
				}
			}
			# For the rest, match based on smallest difference from expected size, or best overlap
			elsif ($type ne 'INS' && ($size==0 || $diff==0 || ($diff && $sizediff < $diff) || ($size<50 && $cols[3]>=50)
					|| ($ovrlp && $ovrlpsize>$ovr))
			) {
				# Compare to existing bestmatch first
				if ($bestmatch && $ovrlp && $ovrlpsize < $ovr) {
					push @secondmatch, $line;
					$size2 = $cols[3];
				}
				# If no gap overlap AND no previous best, or if current overlap > previous overlap
				elsif ($ovr==0 || ($ovrlp && $ovrlpsize>=$ovr && $sizediff < $diff)) {
					push @secondmatch, $bestmatch if $bestmatch;
					$bestmatch = $line;
					$matchtype = ($cols[4] eq $type || $cols[4] =~ /internal/) ? 'Match' : 'Complex';
					$size = $cols[3];
					$diff=$sizediff;
					$ovr=$ovrlpsize if $ovrlpsize;
				}
				# if no gap found yet, AND no previous overlap found OR current overlap>prev.overlap
				else {
					push @secondmatch, $line;
					$size2 = $cols[3];
				}
				# If gap overlap, and current overlap > gap overlap
				if  ($ovrlpsize && $gapfound && $gapfound<$ovrlpsize) {
					$gapfound='';
				}

			}
		}
		elsif ($type eq 'INV' ){ # if no match to inversion, check for other SV nearby
			next if $cols[4] =~ /GAP/ && $cols[3]>=5000; # ignore GAP when only the INV ends are assembled
			if ((($cols[1]>$c1-$invbuffer && $cols[1]<$c1+$invbuffer) || ($cols[1]>$c2-$invbuffer && $cols[1]<$c2+$invbuffer) ||
				($cols[2]>$c1-$invbuffer && $cols[2]<$c1+$invbuffer) || ($cols[2]>$c2-$invbuffer && $cols[2]<$c2+$invbuffer))
			) {
			push @secondmatch, $line;
			$size2 = $cols[3];
			$found++;
		    }
		}
		elsif ($type eq 'GAIN' || $type eq 'LOSS' ){ # check for SVs within 1kb
			if ( (($cols[1]>$c1-$glbuffer && $cols[1]<$c1+$glbuffer) || 
				($cols[1]>$c2-$glbuffer && $cols[1]<$c2+$glbuffer) ||
				($cols[2]>$c1-$glbuffer && $cols[2]<$c1+$glbuffer) || 
				($cols[2]>$c2-$glbuffer && $cols[2]<$c2+$glbuffer))
			) {
			push @secondmatch, $line;
			$size2 = $cols[3];
			$found++;
		    }
		}
		#  if no type match for a del or ins
		elsif ( $found==0)  {
			if ($cols[4] =~ /GAP/) {
				my ($c1,$c2) = ($cols[1]-50,$cols[2]+50);
				my $ovrlp2 = intersect $set "$c1-$c2";
				my $ovrlpsize = $2-$1 if $ovrlp2 && $ovrlp2 =~ /(\d+)-(\d+)/;
				$ovrlpsize = 1 if $ovrlp2 && $ovrlp2 !~ /(\d+)-(\d+)/; 
				if ($ovrlp2 && $ovr<$ovrlpsize) {
					push @secondmatch, $line;
					$size2 = $cols[3];
					$gapfound = (($gapfound && $gapfound < $ovrlpsize) || !$gapfound) ? $ovrlpsize : $gapfound;
				}
			}
			elsif ( ($cols[1]>$c1-$delinsbuffer && $cols[1]<$c1+$delinsbuffer) || 
				($cols[1]>$c2-$delinsbuffer && $cols[1]<$c2+$delinsbuffer) ||
				($cols[2]>$c1-$delinsbuffer && $cols[2]<$c1+$delinsbuffer) || 
				($cols[2]>$c2-$delinsbuffer && $cols[2]<$c2+$delinsbuffer)
			) {
				# check for larger SV with good overlap if the only match is <50bp)
				if ($size && $size<100 && $cols[3]>=100 && $cols[3]<(($c2-$c1)*2)) {
					push @secondmatch, $bestmatch if $bestmatch;
					$bestmatch = '';
					$size = 0;
					$matchtype = '';
					$diff = 0;
				}
				push @secondmatch, $line;
				$size2 = $cols[3];

			}
		}
	}
	return ($bestmatch,\@secondmatch,$matchtype,$gapfound);
}


sub overlapCheck {
	my ($self,$span1,$span2) = @_;
	my @s1 = split /-/, $span1;
	my @s2 = split /-/, $span2;
	if (($s1[0]>=$s2[0] && $s1[0]<=$s2[1]) ||
		($s1[1]>=$s2[0] && $s1[1]<=$s2[1]) ||
	 	($s2[0]>=$s1[0] && $s2[0]<=$s1[1]) ||
		($s2[1]>=$s1[0] && $s2[1]<=$s1[1]) 
		) 
	{
		return 1;
	}
	else {
		return 0;
	}
}

sub parseNeighbours { # check different contigs only, not split contigs
	my ($self,$hits,$buffer,@predictedCoords) =  @_;
	my @starts = ();
	my %hits;
	my @res;
	my ($lowest,$highest)=(); # check missing alignments at ends
	foreach my $line (@$hits) {
		if ($line =~ /^\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\S)/) {
			my $start = $3 eq '+' ? $1 : $2;
			if (!$hits{$start}) {
				push @starts, $start;
				$hits{$start} = $line;
			}
			else {
				my $size1 = $3 eq '+' ? $2-$1 : $1-$2;
				my $size2 = $hits{$start}=~/^\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\S)/ && $3 eq '+' ? $2-$1 : $1-$2;
				if ($size2<$size1) {
					$hits{$start+1} = $line;
				}
				else {
					$hits{$start+1}=$hits{$start};
					push @starts, $start+1;
					$hits{$start}=$line;
				}
			}
		}
	}
	my @sorted = sort {$a<=>$b} @starts;
	my ($q,$qstart,$qend,$tdir,$tstart,$tend,$qsize,$ch) = ();
	my ($lqdiff,$rqdiff) = ''; # diff at left side, diff at right side  
	foreach my $i (0..$#sorted) {
		my @cols = split " ", $hits{$sorted[$i]};
		# reverses target coords if on the minus strand
		($cols[5],$cols[6]) = ($cols[6],$cols[5]) if $cols[6]<$cols[5];
		# check for missing alignment at beginning
		if (!$q) {
			if ($cols[5]>$predictedCoords[0]+500) {
				my $gapsize = $cols[5]-$predictedCoords[0];
				my $c= $predictedCoords[0]+$buffer < $cols[5] ? $predictedCoords[0]+$buffer : $predictedCoords[0];
				#my $c=$predictedCoords[0];
				$c-=1;
				my $print = "$cols[4] $c $cols[5] $gapsize GAPEND $cols[0]";
				push @res, $print;
			}
			($q,$qstart,$qend,$tdir,$tstart,$tend,$qsize,$ch) = ($cols[0],$cols[1],$cols[2],$cols[7],$cols[5],$cols[6],$cols[10],$cols[4]);
			if ($tdir eq '+') {
				$rqdiff = $qsize - $qend;
				$lqdiff = $qstart;
			}
			else {
				$lqdiff = $qsize - $qend;
				$rqdiff = $qstart;
			}
			#$lqdiff = $qstart;
			next;
		}
		next if $cols[6]<$tend; # skip internal alignments
		if ($i==$#sorted) {
			if ($cols[6]<$predictedCoords[1]-500) {
				my $gapsize = $predictedCoords[1]-$cols[6];
				#my $c = $predictedCoords[1]-$buffer;
				my $c= $predictedCoords[1]-$buffer > $cols[6] ? $predictedCoords[1]-$buffer : $predictedCoords[1];
				my $c1 = $cols[6]-1;
				my $print = "$cols[4] $c1 $c $gapsize GAPEND $cols[0]";
				push @res, $print;
			}
		}
		my $current_lqdiff = $cols[7] eq '+' ? $cols[1] : $cols[10]-$cols[2];
		my $current_rqdiff = $cols[7] eq '+' ? $cols[10]-$cols[2] : $cols[1];
		if ($q && $cols[0] eq $q) { # ignore if two hits from the same contig are adjacent
			($q,$qstart,$qend,$tdir,$tstart,$tend,$qsize,$ch) = ($cols[0],$cols[1],$cols[2],$cols[7],$cols[5],$cols[6],$cols[10],$cols[4]);
			($lqdiff,$rqdiff)=($current_lqdiff,$current_rqdiff); # shift contigs over by one
			next;
		}
		elsif ($rqdiff>10 || $current_lqdiff>10) {
			my $size = '';
			$size = $rqdiff+$current_lqdiff;
			
			my ($begin,$end,$contigs,$delgap) = ($tend,$cols[5]+1,"$q;$cols[0]","DEL");# defaults
			my $delsize=$end-$begin+1;
			if ($delsize>100) {
				if ($rqdiff>10 && $current_lqdiff<=10) {
					$contigs = $q;
					$delgap="GAP";
				}
				elsif ($rqdiff<=10 && $current_lqdiff>10) {
					$contigs = $cols[0];
					$delgap="GAP";
				}
				$begin-=1 if $begin==$end;
				$begin-=1;
				my $print = "$cols[4] $begin $end $delsize $delgap"."INS $contigs"; #delins or gapins
				push @res, $print;
			}
			else {
				if ($rqdiff>10 && $current_lqdiff>10) {
					$begin = $tend;
					$end = $cols[5]+1;
					# account for overlapping contigs
					($begin,$end) = ($end-1,$begin+1) if $end < $begin;
					$contigs = "$q;$cols[0]";
				}
				elsif ($rqdiff>10) { # prev contig only
					($begin,$end) = ($tend,$tend+1);
					$contigs = $q;
				}
				else { #current contig only
					($begin,$end) = ($cols[5],$cols[5]+1);
					$contigs = $cols[0];
				}
				$begin-=1 if $begin==$end;
				$begin-=1;
				my $print = "$cols[4] $begin $end $size INS $contigs";
				push @res, $print;
			}	
		}
		else {
			my $size = $cols[5]+1-$tend+1;
			if ($size>100) {
				$tend-=1;
				my $print = "$cols[4] $tend ".($cols[5]+1)." $size GAP $q;$cols[0]";
				push @res, $print;
			}
			elsif ($size>0) {
				$tend-=1;
				my $print = "$cols[4] $tend ".($cols[5]+1)." $size smGAP $q;$cols[0]";
				push @res, $print;
			}

		}
		($q,$qstart,$qend,$tdir,$tstart,$tend,$qsize,$ch) = ($cols[0],$cols[1],$cols[2],$cols[7],$cols[5],$cols[6],$cols[10],$cols[4]);
		($lqdiff,$rqdiff)=($current_lqdiff,$current_rqdiff); # shift contigs over by one
	}
	if ($rqdiff && $rqdiff > 0) {
		$tend-=2;
		my $print = "$ch $tend ".($tend+1)." $rqdiff INS $q\n";
		push @res, $print;
	}
	return \@res;
}

sub parseSplit {
	my ($self,$hits) = @_;
    # sort by qstart
	my @starts = ();
	my %hits;
	my @res;
	my $tothit=0;
	my $small='';
	my $exonParse = new ExonerateParse;
	foreach my $line (@$hits) {
		if ($line =~ /^\S+\s+(\d+)\s+(\d+)/) {
			$tothit+=($2-$1);
			$small=1 if $2-$1<50;
			push @starts, $1;
			$hits{$1} = $line;
		}
	}
	my @sorted = sort {$a<=>$b} @starts;
	my ($qstart,$qend,$tdir,$tstart,$tend,$qsize) = ();
	my $i = 0;
	for (0..$#sorted) {
		my $which=$_;
		my @cols = split " ", $hits{$sorted[$which]};
		if (!$tdir) {
			($qstart,$qend,$tdir,$tstart,$tend,$qsize) = ($cols[1],$cols[2],$cols[7],$cols[5],$cols[6],$cols[10]);
			$i++;
			next;
		}
		my $qdiff = $qend-$cols[1];
		# Adjust coords if there is overlap of alignments (repeats at breakpoints)
		#
		# inversion or "inverted duplication" example: NOD chr1:1_195956456_INV;1:195956456-195957256
		next if $cols[2]-$cols[1]<50; #spurious hits
		if ($tdir ne $cols[7] && $cols[1]>100 && abs($qdiff)<100 ) { # overlap is less than 50bp
			# inverted duplications (insertions)
			#current is -, previous is +
			if ($cols[7] eq '-' && $tstart==$cols[6]) { # dupilcated region is from the left
				my ($istart,$iend) = ();
				my ($c1start,$c1end) = $tstart < $tend ? ($tstart,$tend) : ($tend,$tstart);
				my ($c2start,$c2end) = $cols[5] < $cols[6] ? ($cols[5],$cols[6]) : ($cols[6],$cols[5]);
				my $ovrlp = $exonParse->spanCheck("$c1start-$c1end","$c2start-$c2end"); # check for overlapping contigs
				if ($ovrlp > 1) {
					my @sort = sort {$a<=>$b}($c1start,$c1end,$c2start,$c2end);
					($istart,$iend) = ($sort[0],$sort[-1]);
				}
				else {
					$istart = $tstart;
					$iend = $tend > $cols[5] ? $tend : $cols[5];
				}
				($istart,$iend) =  $exonParse->breakpointCheck($qdiff,$istart,$iend) if $qdiff>0 && !$ovrlp;
				my $size = $iend-$istart;
				if ($size==0) {
					$istart-=2;
					$size=1;
				}
				my $print = "$cols[4] ".($istart)." $iend $size INVDUPL $cols[0]";
				push @res, $print;
			}
			#current is +, previous is -
			elsif ($cols[7] eq '+' && $tstart==$cols[6]) {# duplicated region is from the right 
				my ($istart,$iend) = ();
				my ($c1start,$c1end) = $tstart < $tend ? ($tstart,$tend) : ($tend,$tstart);
				my ($c2start,$c2end) = $cols[5] < $cols[6] ? ($cols[5],$cols[6]) : ($cols[6],$cols[5]);
				my $ovrlp = $exonParse->spanCheck("$c1start-$c1end","$c2start-$c2end"); # check for overlapping contigs
				if ($ovrlp > 1) {
					my @sort = sort {$a<=>$b}($c1start,$c1end,$c2start,$c2end);
					($istart,$iend) = ($sort[0],$sort[-1]);
				}
				else {
					$istart = $tstart <= $cols[6] ? $tend : $cols[5];
					$iend = $tend < $cols[5] ? $tend : $cols[5];
				}
				($istart,$iend) =  $exonParse->breakpointCheck($qdiff,$istart,$iend) if $qdiff>0 && !$ovrlp;
				my $size = $iend-$istart;
				if ($size==0) {
					$istart-=2;
					$size=1;
				}
				my $print = "$cols[4] ".($istart)." $iend $size INVDUPR $cols[0]";
				push @res, $print;
			}
			# inversions
			else {
				my $type="INV";
				# 4 differnet situations for inversions; must check for deletion if
				# end of contig is elsewhere or 
				# beginning of contig is elsewhere
				# Remember hits are sorted by qstarts!
				my ($istart,$iend) = ();
				my ($c1start,$c1end) = $tstart < $tend ? ($tstart,$tend) : ($tend,$tstart);
				my ($c2start,$c2end) = $cols[5] < $cols[6] ? ($cols[5],$cols[6]) : ($cols[6],$cols[5]);
				my $ovrlp = $exonParse->spanCheck("$c1start-$c1end","$c2start-$c2end"); # check for overlapping contigs
				if ($ovrlp > 1) {
					my @sort = sort {$a<=>$b}($c1start,$c1end,$c2start,$c2end);
					($istart,$iend) = ($sort[0],$sort[-1]);
				}
				else {
					($istart,$iend) = $tstart <= $cols[6] ? ($tend,$cols[5]) : ($cols[5],$tend);
				}
				# check for insertion next to inversion

				if ( ($cols[7] eq '+' && (($cols[5]<=$tstart && $qsize-$cols[2]>50)||($cols[6]>$tend && $qstart>50))) ||
						($cols[7] eq '-' && (($cols[6]<=$tstart && $qstart>50)||($cols[5]>$tend && $qsize-$cols[2]>50)))
				) {	
					$type.='INS' if ($i==$#sorted && $qsize-$cols[2]>100)||($i==1 && $qstart>100); #make sure the previous hit is not a middle hit
				}
				# check for deletion inside inversion
				if ($i<$#sorted && (($cols[7] eq '+' && $tend-$cols[6]>100)  ||($cols[7] eq '-' && $cols[6]-$tend>100 ))) {
					#$type.='DEL' if ($i<$#sorted); # if current is a middle hit	
					my @cols2 = split " ", $hits{$sorted[$i+1]};
					if ($cols[7] eq '+' && $cols2[7] eq '+' && $cols2[5]-$cols[6]>100 && abs($cols[2]-$cols2[1])<50) {
						$type.='DEL';
					}
					elsif ($cols[7] eq '+' && $cols2[7] eq '-' && $cols[6]-$cols2[5]>100 && abs($cols[2]-$cols2[1])<50) {
						$type.='DEL';
					}
				}
				if ( $i>=2 && (($cols[7] eq "-" && $tstart-$cols[5]>100)||($cols[7] eq "+" && $cols[5]-$tstart>100))) {
					my @cols3 = split " ", $hits{$sorted[$i-2]};
					if ($cols3[7] ne $tdir) {
						$type.='DEL';
					}
				}
				my $first = $tstart<$cols[6] ? 1 :2 ;
				($istart,$iend) =  $exonParse->breakpointCheck($qdiff,$istart,$iend,$first) if $qdiff>0 && !$ovrlp;
				my $size = $iend-$istart;
				if ($size==0) {
					$istart-=2;
					$size=1;
				}
				$istart-=1 if $istart==$iend;
				my $print = "$cols[4] ".($istart)." $iend $size $type $cols[0]";
				push @res, $print;
			}

		}
		elsif ($tdir ne $cols[7] && abs($qdiff)>=100 ) { # overlap is more than 100
			my $type="INVCOMPLEX";
			my ($istart,$iend) = ();
			my ($c1start,$c1end) = $tstart < $tend ? ($tstart,$tend) : ($tend,$tstart);
			my ($c2start,$c2end) = $cols[5] < $cols[6] ? ($cols[5],$cols[6]) : ($cols[6],$cols[5]);
			($istart,$iend) = $tstart <= $cols[6] ? ($tend,$cols[5]) : ($cols[5],$tend);
			($istart,$iend) = ($iend,$istart) if $istart>$iend; # for more complicated cases
			my $size = $iend-$istart;
			if ($size==0) {
				$istart-=2;
				$size=1;
			}
			my $print = "$cols[4] ".($istart)." $iend $size $type $cols[0]";
			push @res, $print;
			


		}
		# insertion or deletion
		if ($tdir eq $cols[7]) {
			# when alignments are  --b--> --a-->  or  <--a--   <--b-- # 
			if (($tdir eq '+' && $cols[6]<$tend && $cols[5]<$tend) || ($tdir eq '-' && $cols[5]>$tstart && $cols[6]>$tstart)) {
				next;
			}
			elsif ($cols[1]<$qend && $cols[5]==$tend) {
				my $type = "INS";
				my ($istart,$iend) = ($tend-1,$tend);
				my $size = $qend-$cols[1]+1;
				$istart-=1 if $istart==$iend;
				$istart-=1;
				my $print = "$cols[4] $istart $iend $size $type $cols[0]";
				push @res, $print;
			}
			elsif ($cols[1]-$qend>50) { #items sorted by qstart
				my $type = abs($cols[5]-$tend)<50 ? "INS" : "GAP";
				my ($istart,$iend) = $cols[5]<=$tend ? ($cols[5],$tend) : ($tend,$cols[5]);
				my $first = $tstart<$cols[6] ? 1 :2 ;
				($istart,$iend) =  $exonParse->breakpointCheck($qdiff,$istart,$iend,$first) if $qdiff>0;
				my $size = $type eq 'GAP' ? $iend-$istart : abs($qdiff) ;
				if ($size==0) {
					$istart-=2;
					$size=1;
				}
				$istart-=1 if $istart==$iend;
				$istart-=1;
				my $print = "$cols[4] $istart $iend $size $type $cols[0]";
				push @res, $print;
			}
			else {	
				next if ($cols[2]-$cols[1]<25 && $which==$#sorted) || ($cols[2]-$cols[1]<50 && $which<$#sorted && $which>0); # exclude spurious matches
				my $type = abs($tend-$cols[5])>50 ? "DEL" : "GAP";
				my ($istart,$iend) = $cols[5]<=$tend ? ($cols[5],$tend) : ($tend,$cols[5]);
				my $first = $tstart<$cols[6] ? 1 :2 ;
				($istart,$iend) =  $exonParse->breakpointCheck($qdiff,$istart,$iend,$first) if $qdiff>0;
				if (($qend-$cols[1])/($qend-$qstart)<0.5) {
					my $size = $iend-$istart;
					# Cases where the target coords are the same, and the contig overlap is <50
					if ($size==0) {
						$istart-=2;
						$size=1;
					}
					$istart-=2 if $istart==$iend;
					my $print = "$cols[4] ".($istart)." $iend $size $type $cols[0]";
					push @res, $print;
				}
			}

		}
		($qstart,$qend,$tdir,$tstart,$tend) = ($cols[1],$cols[2],$cols[7],$cols[5],$cols[6]);
		$i++;
	}
	@res=$self->checkContigs(@res) if @res;
	return \@res;

}

sub checkContigs {
	my ($self,@results) = @_;
	my @res=();
	my %seen;
	my $chr = $1 if $results[0] =~/^(\S+)/;
	for (@results) {
		if ($_ !~ /INV/) {
			push @res, $_;
		}
		else {
			my $line = $_;
			my @cols = split /\s+/, $line;
			$seen{"$cols[4] $cols[5]"}{"$cols[1]-$cols[2]"}=$line;
		}
	}
	foreach my $key (keys %seen) {
		my @coords = keys %{$seen{$key}};
		if (@coords!=2) {
			for (@coords) {
				push @res, $seen{$key}{$_};
			}
		}
		else {
			my $set = new Set::IntSpan "$coords[0]";
			$set = union $set "$coords[1]";
			if ($set =~ /^(\d+)-(\d+)$/) {	
				my $size=$2-$1;
				$set=~s/-/ /;
				push @res, "$chr $set $size $key";
			}
			else {
				push @res, $seen{$key}{$coords[0]};
				push @res, $seen{$key}{$coords[1]};
			}
		}
	}
	return @res;
}

sub parseIndels {
	my ($self,$hits) = @_;
	my @res;
	my $exonParse = new ExonerateParse;
	for (@$hits) {
		next if /no hit/ || /^\s*$/;
		my @cols = $exonParse->exonerateHitParse($_);
		($cols[5],$cols[6]) = ($cols[6],$cols[5]) if $cols[5]>$cols[6];
		$cols[11] = $exonParse->reverseCIGAR if $cols[7] eq '-';
		my $gcoord=$cols[5];
		my %indels=();
		my ($merged,$change) = $exonParse->mergeCIGAR($cols[11]);
		while ($merged =~ m/(\S+)\s+(\d+)/g) {
			my $type = 'internal';
			$type .= "m" if $change;
			if ($1 eq 'M') {
				$gcoord+=$2;
			}
			elsif ($1 eq 'I') {
				push @{$indels{$cols[0]}}, [$gcoord-1,$gcoord+1,$2,"INS$type"];
			}
			elsif ($1 eq 'D') {
				push @{$indels{$cols[0]}}, [$gcoord,$gcoord+$2,$2,"DEL$type"];
				$gcoord+=$2;
			}
		}
		for (keys %indels) {
			my $contig = $_;
			for (0..$#{$indels{$contig}}) {
				my $print = "$cols[4] ";
				$print .= join(" ", @{$indels{$contig}[$_]});
				$print .= " $contig";
				push @res, $print;
			}
		}
	}
	return \@res;

}


sub filterMatches {
	my ($self,$matches) = @_;
	my %coords=();
	my $exonParse = new ExonerateParse;
	for (@$matches) {
		next if $_ !~ /^\S+\s+\d+\s+\d+\s+\d+\s+/;
		my $line = $_;
		if ($line =~ /INV/ || $line eq 'DEL' || $line eq 'DELINS') {
			$coords{$1}{"$2-$3"}=$line if $line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+/;
		}
	}
	my @pass;
	# this part filters out false INS from parseNeighbor
	if (keys %coords) {
		foreach my $match (@$matches) {
			my ($c,$s,$e) = ($1,$2,$3) if $match =~ /^(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+/;
			next if !$s || !$e;
			if ($coords{$c}{"$s-$e"} && $coords{$c}{"$s-$e"} eq $match ) {
				push @pass, $match;
				next;
			}
			elsif ($match =~ /internal/) {
				push @pass, $match;
				next;
			}
			my $over=0;
			foreach my $span (keys %{$coords{$c}}) {
				my ($a,$b) = split "-", $span;
				foreach ($a,$b) {
					my $span1 = ($_-50)."-".($_+50);
					my $int = $exonParse->spanCheck($span1,"$s-$s");
					my $int2 = $exonParse->spanCheck($span1,"$e-$e");
					$over++ if $int>0 || $int2>0;
				}
			}
			if (!$over) {
				push @pass, $match;
			}
		}
		return \@pass;	
	}
	else { return $matches }
}

1;
