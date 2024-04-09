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

package ParseConfig;

use strict;
use Carp;

# Read in and parse SV config file


=head1 NAME

ParseConfig

=head1 AUTHOR

kw10@sanger.ac.uk

=head1 DESCRIPTION

Used to parse configuration files for the SVMerge pipeline and check for required parameters.

=head1 METHODS

=head2 getParams

	Arg[1]      :  configuration file
	Example     :  getParams("svmerge.config")
	Description :  parse parameters for all SVMerge parameters
	Returns:    :  A hash with parameters as keys

=head2 checkParams

	Arg[1]      :  parameters to check [svcallers,filter,assembly,parse]
	Arg[2]      :  reference to parameter hash from getParams
	Example     :  checkParams(svcallers,\%params)
	Description :  checks for required parameters for the SVMerge pipeline
	Returns     :  true if check is successful

=cut

sub getParams {
	my $file = shift @_;
	my %params = ();
	open F, "<$file" or die "Can't open configuration file $file\n";
	while (<F>) {
		next if /^#/;
		if ( /([^=]+)=(.+)/) {
			my ($par,$val)=($1,$2);
			$val=~s/\s+$//;
			$params{$par} = $val;#print "zzz\t$1\t$2\n";
		}
	}
	close F;
	#print "$params{chrrefdir}\tzzz\n";

	# Create a list of chromosomes to be analysed
	#
	my @chrList = ();
	if ($params{chrList}) {
		@chrList = `cat $params{chrList}`;
	}
	if ($params{chrRange}) {
		my @list = split ",", $params{chrRange};
		for (@list) {
			if (/-/) {
				my @range = split "-", $_;
				for ($range[0]..$range[1]) {
					push @chrList, $_;
				}
			}
			else {
				push @chrList, $_;
			}
		}

	}
	$params{chrs} = @chrList;
	if ($params{callerlist}) {
		my @callers = split /\s+/, $params{callerlist};
		for (@callers) {
			$params{$_} = 1;
		}
	}
	return %params;
}


# Check for required parameters missing

sub checkParams {
	my ($which,$params) = @_;
	my %params=%$params;
	# general parameters
	for ('name','project','version','name','svdir','exedir','projdir','callerlist','defaultQueue') {
		croak "Missing parameter $_\n" if !$params{$_};
	}
	if ($params{name} =~ /_/) {
		$params->{name} =~ s/_/-/g;
		print STDERR "name parameter contains '_'; changing to $params->{name}\n";
	}
	if ($params{chrrefdir} && $params{reffile}) {
		croak "Specify EITHER 'chrrefdir' OR 'reffile'\n";
	}
	elsif (!$params{chrrefdir} && !$params{reffile}) {
		croak "Missing reference fasta(s): specify EITHER 'chrrefdir' OR 'reffile'\n";
	}
	if (!$params{chrRange} && !$params{chrOther}) {
		croak "Must specify at least one of 'chrRange' and/or 'chrOther'\n";
	}
	if ($params{bam} && $params{bamdir}) {
		croak "Specify either a single BAM file (bam and bai) or a BAM directory with chromosome BAMs (bamdir)\n";
	}
	elsif ($params{bam} && !$params{bai}) {
		croak "Specify the the BAM index file (bai)\n";
	}
	# SV calling parameters
	if ($which eq 'svcallers') {
		if ($params{breakdancer}) {
			for ('bam2conf','bdexe','bam','bai','samtools','BDmapq','BDsd') {
				croak "Missing Breakdancer parameter $params{$_}\n" if !$params{$_};
			}
		}
		if ($params{pindel}) {
			for ('pinexe','PDinsert','samtools') {
				croak "Missing Breakdancer parameter $params{$_}\n" if !$params{$_};
			}
		}
		if ($params{cnD}) {
			for ('cnddir','CNDnohet','CNDsnprate','CNDrep','samtools') {
				croak "Missing Breakdancer parameter $params{$_}\n" if !$params{$_};
			}
		}
		if ($params{sec}) {
			for ('SECminCluster','SECmin','SECmax','samtools') {
				croak "Missing Breakdancer parameter $params{$_}\n" if !$params{$_};
			}
		}
		if ($params{rdx}) {
			croak "Missing RDXplorer paramater RDXout\n" if !$params{RDXout};
		}
	}
	# filtering parameters
	elsif ($which eq 'filter') {
		for ('gaps','centel','filterOther'){
			my $buffer = $_ . "Buffer";
			croak "Missing parameter $buffer for $_" if $params{$_} && !$params{$buffer};
		}
	}
	# config/assembly parameters
	elsif ($which eq 'assembly') {
		for ('samtools','submatrix','exonerateExe') {
			croak "Missing parameter $_\n" if !$params{$_};
		}
		# running velvet and exonerate
		if ($params{velvet}) {
			for ('velveth','velvetg','hashlen','ins_len','exp_cov','cov_cutoff') {
				croak "Missing parameter $_\n" if !$params{$_};
			}
		}
		if ($params{abyss}){
			for ('abyss-pe','kmer') {
				croak "Missing ABySS parameter $_\n" if !$params{$_};
			}
			if ($params{paired} && !$params{npairs}) {
				croak "Missing ABySS npairs parameter for paired assembly\n";
			}		
		}
#		# not running Velvet - require contig file and outdir for /exonerate/
#		elsif (!$params{exonerate}) {
#			croak "No paramaters for velvet - running alignments only. No exonerate parameters found\n";
#		}
#		elsif ($params{exonerate}) {
#			for ('contigs','exonOutdir') {
#				croak "Missing parameter $_\n" if !$param{$_};
#			}
			
		#}
	}
	elsif ($which eq 'parse') {
		if ($params{bamcheck}) {
			for ('meanCov','offset','parseSplitLines') {
				croak "Missing parameter $_\n" if !$params{$_};
			}
			if (!$params{bam} && !$params{bamdir}) {
				croak "Must specify BAM file of BAM directory\n";
			}
			elsif ($params{bam} && !$params{bai}) {
				croak "Must specify BAM index file\n";
			}
		}
		if (!$params{offset} && $params{subseq}) {
			croak "Missing parameter 'offset' (subseq=1)\n";
		}
	}

	return 1;

}
1;
