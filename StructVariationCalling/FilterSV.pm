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

package FilterSV;


use strict;
use warnings;
use Carp;
use Set::IntSpan;

=head1 NAME 

FilterSV

=head1 AUTHOR

kw10@sanger.ac.uk

=head1 DESCRIPTION

Methods to filter output from various SV callers, and re-format to tab-delimited and BED formats.

=head1 METHODS

=head2 new

	Arguments   :  none
	Description :  create a new FilterSV object

=cut

sub new {
	my ($class,%arg) = @_;
	my $self={};
	bless($self,$class || ref($class));
	return $self;
}

=head2 coordCheck

	Arg[1]      :  coordinate set [1000-2000 or "1000 2000"]
	Description :  checks for proper format and start < end

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
# parse tab/space delimited files with coordinates to exclude
# when filtering SVs

=head2 parseGaps

	Arg[1]      :  array of arrays with coordinate files and buffer size
	Example     :  parseGaps(@listOfFiles)
	Description :  pushes all coordinates from the files into a single hash
	Returns     :  a hash of hashes with chromosome keys

=cut
sub parseGaps {
	my ($self,@list) = @_;
	my %gaps = ();
	my $set = 0;
	foreach (0..$#list) {
		$set++;
		open F, "<$list[$_][0]" or croak "Can't open $list[$_][0] for reading\n";
		my $buffer = $list[$_][1];
		while (<F>) {
			if (/^(\S+)\s+(\S+)\s+(\S+)/) {
				my ($c,$s,$e) = ($1,$2,$3);
				$s-=$buffer;
				$e+=$buffer;
				if (!$gaps{$set}{$c}) {
					$gaps{$set}{$c} = new Set::IntSpan "$s-$e";
				}
				else {
					$gaps{$set}{$c}->U("$s-$e");
				}
			}
		}
	}
	return %gaps;
}

# Take tab-delimited line and check for overlap with gaps

=head2 filterGaps 

	Arg[1]      :  tab or space delimited coordinates
	Arg[2]      :  reference to hash from parseGaps
	Arg[3]      :  minimum fraction of query with overlap
	Example     :  filterGaps($line,\$gaps,0.5)
	Description :  compares coordinates to coordinates in a hash
	               and returns the line only if there is less
	               overlap than specified in Arg[3]
	Returns:    :  returns original string in Arg[1]

=cut

sub filterGaps {
	my ($self,$line,$gaps,$perc) = @_;
	my ($c,$s,$e) = split /\s+/, $line;
	$self->coordCheck("$s-$e");
	for (keys %$gaps) {
		my $set = $gaps->{$_}->{$c};
		next if !$set;
		foreach my $sp (spans $set) {
			my $newset="$sp->[0]-$sp->[1]";
			$newset = new Set::IntSpan $set;
			my $check = intersect $newset "$s-$e";
			if ($check) {
				my $size1 = $e-$s+1;
				my $size2 = $check =~ /(\d+)-(\d+)/ ? $2-$1 : '1';
				#if (!$size2 || !$size1 || !$perc) {
				#	print STDERR "$size1, $size2, $perc, $line\n";
				#}
				if ($size2/$size1>$perc) { #overlap>25% of prediction
					#push @hits, $_;
					return '';
				}
			}
		}
	}
	# if it reaches here than no overlap was found
	return $line;


}

# initial filtering of BD calls with 'unknown' or extreme # or read supports
#
sub cleanBD {
	my ($self,$line) = @_;
	return if $line =~ /^#/;
	chomp $line;
	my @col = split "\t", $line;
	my $check=();
	# Filter out bad calls first and "unknown"
	if ($col[6] eq "UN" || $col[6] eq 'ITX') {
		return;
	}
	elsif ( (($col[6] eq "DEL" || $col[6] eq "INS") && $col[1] >= $col[4]) ||
		($col[6] eq "INV"  && ($col[7]<0 || $col[1] > $col[4] || $col[4]-$col[1]<200 ))) { 
		return;
	}
	# Insertion calls <100 usually artifacts
	elsif ($col[6] eq "INS" && abs($col[7])<100) {
	#elsif ($col[6] eq "INS" && abs($col[7])<250) {
		return;
	}
	elsif ($col[6] eq "DEL" && $col[4]-$col[1]<100) {
		return;
	}
	else {
		for ($col[2],$col[5]) {
			my $supports=$_;
			if ($supports=~/(\d+)\+(\d+)-/ && ($1>=1000 || $2>=1000)) {
				return;
			}
		}
	}
	return $line;
}

=head2 filterBDscore

	Arg[1]      :  output line from BreakDancerMax
	Arg[2]      :  minimum score for filtering
	Arg[3]      :  minimum read supports
	Example     :  filterBDscore($line,30,2)
	Description :  filters BreakDancerMax output by score and read supports
	Returns     :  original line in Arg[1] if it passes filtering

=cut
sub filterBDscore {
	my ($self,$line,$score,$rs) = @_;
	my @cols = split "\t", $line;
	if ($score && $rs && $cols[8] >= $score && $cols[9] >= $rs) {
		return $line;
	}
	elsif ($score && !$rs && $cols[8] >= $score) { 
		return $line;
	}
	elsif (!$score && $rs && $cols[9] >= $rs ) {
		return $line;
	}
	else {
		return '';
	}
}
=head2 filterCNVnator

        Arg[1]      :  output line from BreakDancerMax
        Arg[2]      :  minimum Pvalue
        Arg[3]      :  minimum normalized read depth
	Arg[4]      :  read depth complex index  
        Example     :  filterBDscore($line,30,2)
        Description :  filters CNVnatr output by score and P-value
        Returns     :  original line in Arg[1] if it passes filtering

=cut
sub filterCNVnator {
        my ($self,$line,$P,$q0,$RDcomplex) = @_;
        my @cols = split "\t", $line;
        if ($cols[4]<$P  && $cols[8]<0.5 and $cols[0] eq 'deletion' and $cols[3]*(1+$cols[8]) < $RDcomplex) {
	        return $line;
	}elsif($cols[4]<$P  && $cols[8]<0.5 and $cols[0] eq 'duplication'){
		return $line;
        }else {
                return '';
        }
}
=head2 filterBDcn

	Arg[1]      :  output line from BreakDancerMax
	Arg[2]      :  maximum copy number for filtering deletions
	Example     :  filterBDcn($line,2) # column 12 <= 2
	Description :  filters BreakDancerMax output by copy number (column 12)
	Returns     :  original line in Arg[1] if it passes filtering

=cut

sub filterBDcn {
	my ($self,$line,$copynum) = @_;
	my @cols = split "\t", $line;
	if ($copynum && ($cols[11] eq 'NA' || $cols[11] <= $copynum)) {
		return $line;
	}
	else {
		return '';
	}
}


sub getRGB {
	my ($self,$type) = @_;
	my %SVcolors = (
		'DEL' => '255,0,0', # red
		'INS' => '0,100,0', # d.green
		'INV' => '0,0,255', # blue
		'GAIN' => '34,139,34', # forest green
		'LOSS' => '255,0,0', # red
	);
	return $SVcolors{$type} if $SVcolors{$type};
	return '0,100,0' if $type=~/INS|/;
	return '0,0,0'; 
}

=head2 tab2bed

	Arg[1]      :  'array' or 'file'
	Arg[2]      :  reference to array, or file name
	Arg[3]      :  output file name
	Arg[4]      :  sample name
	Arg[5]      :  method [BD,PD,CND,SEC,RDX]
	Arg[6]      :  annotation
	Example     :  tab2bed('array',\@lines,"SVcalls.bed",NA18506,BD,final)
	Description :  converts a tab delimited set of SVcalls to BED format
	Returns     :  prints out BED lines; returns nothing

=cut
sub tab2bed {
	my ($self,$what,$input,$output,$name,$method,$desc) = @_;
	my @lines = ();
	if ($what eq 'file') {
		@lines = `cat $input`;
	}
	else {
		@lines = @$input;
	}
	open OUT, ">$output" or croak "Can't open $output for writing\n";
	print OUT qq(track name="$name $desc $method" description="$name $desc $method" useScore=0 itemRgb="On" visibility=2\n);
	for (@lines) {
		next if /^\s*$/;
		my $line = $_;
		chomp $line;
		my @cols = split "\t", $line;
		$cols[0] = "chr$cols[0]" if $cols[0] !~ /chr/;
		if (@cols<4) {
			croak "File must be tab-delimited with at least 4 columns\n";
		}	
		my $type = $1 if $cols[3] =~ /^([^_]+)/;
		$cols[1]-=1 unless $cols[1]==0;
		my $rangesize = $cols[2]-$cols[1];
		my $rgb = $self->getRGB($type);
		my $coord = "$cols[0]\t$cols[1]\t$cols[2]";
		my $size = $1 if $cols[3] =~ /_(\d+)$/;
		$size = 'large' if !$size;
		my $name = "$type\_$method\_$name\_$size";
		my $strand = "+"; # all + by default

		my $bed = join "\t", ($coord,$name,0,$strand,$cols[1],$cols[2],$rgb,1,$rangesize,0);
		print OUT $bed."\n";
	
	}
	close OUT;
	return;


}

=head2 bd2tab

	Arg[1]      :  BreakDancerMax line
	Arg[2]      :  ID [eg: sample name]
	Example     :  bd2tab($line,NA18506)
	Description :  Convert BreakDancer output to tab-delimited output with
	               specific IDs in column 4.
	Returns     :  tab-delimited line

=cut 
sub bd2tab {
	my ($self, $line, $id) = @_;
	my @cols = split "\t", $line;
	# chr, start, stop, SV_BD_name_size, score (out of 1000)
 	my $score = $cols[8]*10 > 1000 ? 1000 : $cols[8]*10; 
 	my $size = $cols[6] eq 'INS' ? $cols[7] : ($cols[4]-$cols[1]+1);
	my $tab = "$cols[0]\t$cols[1]\t$cols[4]\t$cols[6]\_BD_$id\_". $size;
	return $tab;
	}


=head2 cnvnator2tab

	Arg[1]      :  BreakDancerMax line
	Arg[2]      :  ID [eg: sample name]
	Example     :  bd2tab($line,NA18506)
	Description :  Convert BreakDancer output to tab-delimited output with
	               specific IDs in column 4.
	Returns     :  tab-delimited line

=cut 
sub cnvnator2tab {
	my ($self, $line, $id) = @_;
	my @cols = split "\t", $line;
	# chr, start, stop, SV_BD_name_size, score (out of 1000)
	#my $score = $cols[8]*10 > 1000 ? 1000 : $cols[8]*10; 
	my ($chr,$start,$end)=split(/:|\-/,$cols[1]);
	my $size =  ($end-$start+1);
	my $type="";
	if($cols[0] eq 'deletion'){
		$type="DEL";
	}else{
		$type="DUP"
	}
	my $tab = "$chr\t$start\t$end\t$type\_cnvnator_$id\_". $size;
	return $tab;
}
# Format BD line to BED
sub bd2bed {
	my ($self, $line, $id) = @_;
	my @cols = split "\t", $line;
	# chr, start, stop, SV_BD_name_size, score (out of 1000)
	$cols[6] = uc $cols[6] if $cols[6];
	my $rgb = $self->getRGB($cols[6]);
	my $coord = "$cols[0]\t$cols[1]\t$cols[4]";
	my $size = $cols[4]-$cols[1]+1;
	my $name = "$cols[6]\_BD_$id\_$size";
	my $score = $cols[8]*10 > 1000 ? 1000 : $cols[8]*10;
	my $strand = "+"; # all + by default

	$cols[1]-=1 unless $cols[1]==0;
	my $bed = join "\t", ($coord,$name,$score,$strand,$cols[1],$cols[4],$rgb,1,$size,0);
	return $bed;
}

=head2 parsePindel

	Arg[1]      :  SV type [del,ins,delins]
	Arg[2]      :  file name
	Example     :  parsePindel(del,del.txt)
	Description :  reformats Pindel output
	Returns     :  array with parsed Pindel lines

=cut
sub parsePindel {
	my ($self, $type, $file) = @_;
	if ($type eq 'del' || $type eq 'ins') {
		my @out =  `grep -h -i support $file | perl -ne 's/\\s+/ /g; print \$_."\n"' | cut -f 1-5,7,8,12,13,19,21 -d " "  | sort -k 5,8 -k11,11r | sort -u -k5,8`;
		return @out;
	}
	elsif ($type eq 'delins'){
		my @out = `grep -h -i support $file | perl -ne 's/\\s+/ /g; print \$_."\n"' | cut -f  1-3,6,7,9-12,18,20 -d " "  | sort -k 5,8 -k11,11r | sort -u -k5,8`;
		return @out;
	}

}

=head2 parsePindelv2

	Arg[1]      :  file name
	Example     :  parsePindelv2(del.txt)
	Description :  reformats Pindel v2 output
	Returns     :  array with parsed Pindel lines

=cut

sub parsePindelv2 {
	my ($self, $file) = @_;
	my @out = ();
	my @lines = `grep -h -i support $file`;
	for (@lines) {
		if (/^(\S+\s+\S+\s+\S+)\s+.*(ChrID\s+\S+)\s+.*BP\s+(\d+\s+\d+)\s+.*(Supports\s+\d+)\s+.*S1\s+(\d+)\s+.*SUM_MS\s+(\S+)/) {
			push @out, "$1 $2 $3 $4 $5 $6";
		}
	}
	return @out;

}


=head2 parsePindelv3

	Arg[1]      :  file name
	Example     :  parsePindelv3(pindeloutput)
	Description :  reformats Pindel v0.2.3 output
	Returns     :  array with parsed Pindel lines

=cut

sub parsePindelv3 {
	my ($self, $file) = @_;
	my @out = ();
	my @lines = `grep -h -i support $file`;
	@lines = `grep -h ChrID  $file | grep -w LI` if !@lines;
	for (@lines) {
		# v0.2.4q
		# Num LI ChrID chr pos1 + numreads pos - numreads 
		if (/^(\S+\s+LI)\s+(ChrID\s+\S+)\s+(\S+)\s+\+\s+(\S+)\s+(\S+)\s+-\s+(\S+)/) {
			my $range = $3>$5? "$5 $3" : "$3 $5";
			my $supports = $4+$6;
			push @out, "$1 NA $2 $range Supports $supports NA";
		}
		# v0.2.3
		elsif (/^(\S+\s+LI)\s+(ChrID\s+\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
			my $range = $3>$5? "$5 $3" : "$3 $5";
			my $supports = $4+$6;
			push @out, "$1 NA $2 $range Supports $supports NA";
		}
		elsif (/^(\S+\s+\S+\s+\S+)\s+.*(ChrID\s+\S+)\s+.*BP\s+(\d+\s+\d+)\s+.*(Supports\s+\d+)\s+.*S1\s+(\d+)\s+.*SUM_MS\s+(\S+)/) {
			push @out, "$1 $2 $3 $4 $5 $6";
		}
	}
	return @out;

}

=head2 parsePindelVCF

	Arg[1]      :  file name
	Example     :  parsePindelvcf(pindeloutput)
	Description :  reformats Pindel  V0.2.5b9 output
	Returns     :  array with parsed Pindel lines

=cut

sub parsePindelVCF {
	my ($self, $file) = @_;
	my @out = ();
	open IN, $file||die;
	while(<IN>){
		chomp;
		next if(/^#/);
		my @items = split(/\t/);
		my ($end,$support,$len)=(0,0,0);
		next if($items[7]=~/SVLEN=-(\d+)/ and $1<50);
		next if($items[7]=~/SVLEN=(\d+)/ and $1<50);
		if($items[7]=~/END=(\d+)/){$end=$1};
		foreach my $i (9..$#items){
			if($items[$i]=~/.*:(\d+),(\d+)/){
				$support+=$2;
			}
		}
		if($support >10){
			if($file=~/SI/){
				$len=length $items[4]
			}else{
				$len=$end-$items[1];	
			}
			push @out,"$items[0] $items[1] $end $len";
		}	
	}
	return @out;

}

=head2 parseDellyVCF

	Arg[1]      :  file name
	Example     :  parsePindelvcf(pindeloutput)
	Description :  reformats Pindel  V0.2.5b9 output
	Returns     :  array with parsed Pindel lines

=cut

sub parseDellyVCF {
	my ($self, $file) = @_;
	my @out = ();
	open IN, $file||die;
	while(<IN>){
		chomp;
		next if(/^#/);
		my @items = split(/\t/);
		my ($end,$support,$len)=(0,0,0);
		if($items[7]=~/END=(\d+)/){
			$end=$1;
			$len=$1-$items[1];
		};
		next if($len<100);
		if($items[7]=~/PE=(\d+).*SR=(\d+)/){
			$support=$1+$2;
		};
		next if($support <10);
		push @out,"$items[0] $items[1] $end $len";
	}	
	return @out;

}


=head2 filterPindelScore

	Arg[1]      :  minimum Pindel score
	Arg[2]      :  SV type (del,ins,delins)
	Arg[3]      :  array of reformatted Pindel lines (parsePindel)
	Example     :  filterPindelScore(30,del,@lines)
	Description :  reformats Pindel output
	Returns     :  array with filtered Pindel lines

=cut

sub filterPindelScore {
	my ($self,$score,$type,@lines) = @_;
	my @pass = ();
	foreach my $line (@lines) {
		my @cols = split /\s+/, $line;
		push @pass, $line if ($type eq 'ins' || $type eq 'del') && $cols[10]>=$score;
		push @pass, $line if ($type eq 'delins') && $cols[10]/$cols[8]>=$score; #sum_ms/#supports
	}
	return @pass;

}

=head2 filterPindelScorev2

	Arg[1]      :  minimum Pindel score
	Arg[2]      :  array of reformatted Pindel lines (parsePindelv2)
	Example     :  filterPindelScorev2(30,@lines)
	Description :  reformats Pindel output
	Returns     :  array with filtered Pindel lines

=cut

sub filterPindelScorev2 {
	my ($self,$score,@lines) = @_;
	my @pass = ();
	foreach my $line (@lines) {
		my @cols = split /\s+/, $line;
		push @pass, $line if $cols[10]>=$score;
	}
	return @pass;

}

=head2 filterPindelScorev3

	Arg[1]      :  minimum Pindel score
	Arg[2]      :  array of reformatted Pindel lines (parsePindelv3)
	Example     :  filterPindelScorev3(30,@lines)
	Description :  reformats Pindel v0.2.3 output
	Returns     :  array with filtered Pindel lines

=cut

sub filterPindelScorev3 {
	my ($self,$score,$supp,@lines) = @_;
	my @pass = ();
	foreach my $line (@lines) {
		my @cols = ($1,$2) if $line=~/Supports\s+(\S+)\s+(\S+)/;
		if ($score && $supp) {
			push @pass, $line if ($cols[1] eq 'NA' || $score<=$cols[1]) && $supp<=$cols[0];
		}
		elsif ($score) {
			push @pass, $line if $cols[1] eq 'NA' || $score<=$cols[1];
		}
		elsif ($supp) {
			push @pass, $line if $supp<=$cols[0];
		}
	}
	return @pass;

}

=head2 pindel2tab

	Arg[1]      :  Pindel line
	Arg[2]      :  ID [eg: sample name]
	Arg[3]      :  SV type [eg: DEL,INS,INV]
	Example     :  pindel2tab($line,NA18506,DEL)
	Description :  Convert Pindel output to tab-delimited output with
	               specific IDs in column 4.
	Returns     :  tab-delimited line

=cut 

sub pindel2tab {
	my ($self,$line,$id,$type) = @_;
	$type = uc $type if $type;
	my @cols = split /\s+/, $line;
	# chr, start, stop, SV_Pindel_name_size
	my $tab = "$cols[0]\t$cols[1]\t$cols[2]\t$type\_Pindel\_$id\_$cols[3]";
	return $tab;
}

=head2 delly2tab

	Arg[1]      :  Pindel line
	Arg[2]      :  ID [eg: sample name]
	Arg[3]      :  SV type [eg: DEL,INS,INV]
	Example     :  delly2tab($line,NA18506,DEL)
	Description :  Convert Pindel output to tab-delimited output with
	               specific IDs in column 4.
	Returns     :  tab-delimited line

=cut 

sub delly2tab {
	my ($self,$line,$id,$type) = @_;
	$type = uc $type if $type;
	my @cols = split /\s+/, $line;
	# chr, start, stop, SV_Pindel_name_size
	my $tab = "$cols[0]\t$cols[1]\t$cols[2]\t$type\_Delly\_$id\_$cols[3]";
	return $tab;
}

# Reformat CND calls to tab-delimited

=head2 cnd2tab

	Arg[1]      :  cnD output line
	Arg[2]      :  ID [eg: sample name)
	Example     :  cnd2tab($line,NA18506)
	Description :  reformats cnD output to tab-delimited output with a
	               unique ID for each call
	Returns     :  tab-delimited string

=cut
sub cnd2tab {
	my ($self,$line,$id) = @_;
	my @cols = split /\s+/, $line;
	# chr, start, stop, SV_Pindel_name_size, score (1000 for gain 500 for loss)
	my $state = $cols[3] eq 'G' ? 'GAIN' : 'LOSS';
	my $tab = "$cols[0]\t$cols[1]\t$cols[2]\t$state\_CND\_$id\_" . ($cols[2]-$cols[1]+1);
	return $tab;


}


=head2 sec2tab

	Arg[1]      :  SECluster output line
	Arg[2]      :  ID [eg: sample name)
	Arg[3]      :  SV type [eg: INS]
	Example     :  sec2tab($line,NA18506,INS)
	Description :  reformats SECluster output to tab-delimited output with a
	               unique ID for each call
	Returns     :  tab-delimited string

=cut

sub sec2tab {
	my ($self,$line,$id,$type) = @_;
	$type = uc $type if $type;
	my @cols = split /\s/, $line;
	# chr, start, stop, SV_Pindel_name_size, score = 1000
	my $chr = $1 if $line =~ /^(\S+):/;
	my ($s,$e) = split "-", $cols[2];
	my $tab = "$chr\t$s\t$e\t$type\_SEC\_$id\_na";
	return $tab;
}

=head2 retroParse 

	Arg[1]      :  RetroSeq VCF calls file
	Arg[2]      :  Output file name
	Example     :  retroParse(retro.calls.vcf,retro.vcf.parsed)
	Description :  Reformats RetroSeq VCF output to tab-delimited output with a
	               unique ID for each call; merges calls with same coordinates
	Returns     :  1 upon success; tab-delimited file is created

=cut



sub retroParse {
	my ($self,$chrlist,$in,$out,$name) = @_;
	my %alllines=();
	open I, "<$in" or croak "Can't open $out\n";
	my %chrlist;
	for(@{$$chrlist}) {
		$chrlist{$_}=1;
	}
	while (<I>) {
		chomp;
		next if /^#/;
		my @cols = split /\s/, $_;
		# vcf format
		my ($chr,$s) = ($cols[0],$cols[1]);
		next if !$chrlist{$chr};
		my $type = $cols[7]=~/MEINFO=([^,]+),/ ? $1 : 'INS';
		$type=~s/_/-/g;
		$type = "|$type" if $alllines{$chr}{$s};
		$alllines{$chr}{$s}.=$type;
	}
	close I;
    open O, ">$out" or croak "Can't open $out\n";
	foreach my $c (sort keys %alllines) {
		foreach my $s (sort {$a<=>$b} keys %{$alllines{$c}}) {
			my $e=$s+1;
			#print O2 "$c\t$s\t$e\tINS|$alllines{$c}{$s}\_RETRO\_na\n";
			print O "$c\t$s\t$e\tINS|$alllines{$c}{$s}\_merge_$name\_na\n";
		}
	}
	close O;
	return 1 if -f $out;
	return '';
}

# not used

sub retro2tab {
	my ($self,$line,$id) = @_;
	my @cols = split /\s/, $line;
	# vcf format
	my ($chr,$s) = ($cols[0],$cols[1]);
	my $e = $s+1;
	my $type = $cols[8]=~/MEINFO=([^,]+),/ ? $1 : 'INS';
	my $tab = "$chr\t$s\t$e\t$type\_RETRO\_$id\_na";
	return $tab;
}

=head2 filterRDXscore 

	Arg[1]      :  output line from RDXplorer
	Arg[2]      :  Z-score cutoff
	Example     :  filterRDXscore($line,5)
	Description :  filter RDXplorer output by Z-score
	Returns     :  line in Arg[1] if it passes score filter

=cut

sub filterRDXscore {

	my ($self,$line,$score) = @_;
	my @cols = split "\t", $line;
	if (abs($cols[4])>=$score) {
		return $line;
	}
	else {
		return '';
	}

}

=head2 rdx2tab

	Arg[1]      :  RDXplorer output line
	Arg[2]      :  ID [eg: sample name)
	Example     :  rdx2tab($line,NA18506
	Description :  reformats RDXplorer output to tab-delimited output with a
	               unique ID for each call
	Returns     :  tab-delimited string

=cut
sub rdx2tab {
	my ($self,$line,$id) = @_;
	my @cols = split /\s/, $line;
	# chr, start, stop, SV_RDX_name_size, score
	my $sv = $cols[4]<0 ? 'LOSS' : 'GAIN';
	my $tab = "$cols[0]\t$cols[1]\t$cols[2]\t$sv\_RDX\_$id\_" . ($cols[2]-$cols[1]+1);
	return $tab;

}

1;
