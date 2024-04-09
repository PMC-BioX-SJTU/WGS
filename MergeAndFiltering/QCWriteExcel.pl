#! /usr/bin/perl -w
use Spreadsheet::WriteExcel;
@input=($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]);
@index=('Q20Q30','Alignment','CrossContamination','Coverage');
WriteExcel(\@input,\@index,$ARGV[4]);
sub WriteExcel{
	###=============================get option
	my ($file,$sheet,$output)=@_;
	###===============================formation
	my $workbook = Spreadsheet::WriteExcel->new("$output");
	my @worksheet;
	foreach(@$sheet){
		my $key=$_;
		my $wsh = $workbook->add_worksheet("$key");
		push @worksheet,$wsh;
	}
	my $title = $workbook->add_format();
	$title->set_bold();
	$title->set_color('red');
	$title->set_align('left');
	###==============================write excel line
	my @files=@$file;
	foreach my $f(0..$#files){
		my $row = 0;
		open IN,$files[$f]||die;
		while(<IN>){
			my @info=split(/\t/);
			my @info1=($info[0]);
			for my $k (0..$#info){if($k%2==1 ){push @info1,$info[$k]}}
			if($row==0){
				###==============================write title line
				for(my $i=0;$i<=$#info1;$i++){
					$worksheet[$f]->write($i,$row,$info1[$i],$title);
				}
			}else{
				###===============================write body region
				for (my $j=0;$j<=$#info1;$j++){
					$worksheet[$f]->write($j,$row,$info1[$j]);
				}
			}
			$row++;
		}close IN;
	}
}


