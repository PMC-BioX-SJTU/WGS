#! /usr/bin/perl -w
#argv[1]: path to fq data
#argv[2]:_1.clean.fastq.gz
#argv[3]:_2.clean.fastq.gz
$j=1;
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	$i=$_;
	$name=$i;
	$name=~s/.*\///g;
	push @name,'name'.$j.'='.$name;
	$R1=`ls $ARGV[1]/*$ARGV[2]`;
	$R1=~s/\n/ /g;
	$R2=$R1;
	$R2=~s/$ARGV[2]/$ARGV[3]/g;
	$R2=~s/\n/ /g;
	$out1= "$R1:$R2";print "$out1\tsssss\n";
	$out1=~s/ :/:/g;
	push @fq, 'fq'.$j.'='."$out1";
	push @dir,'dir'.$j.'='.$ARGV[1].'/';
	$j++;

}
$out1=join("\n",@name);
$out2=join("\n",@dir);
$out3=join("\n",@fq);
print "[Sample.name]\n$out1\n[Sample.dir]\n$out2\n[Sample.data]\n$out3\n"

