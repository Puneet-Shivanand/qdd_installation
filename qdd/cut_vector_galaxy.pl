#!/usr/bin/perl -w

# name of the program: cut_vector.pl

# INPUT: 	xxx.tag (fasta file with vectors/adapters)
# 		yyy.fas fasta file to be cleaned and sorted out.

# OUPUT: yyy_m8.tbl results of blast againts te adapter file
#		yyy_tagY.fas 1 fasta file with all of the cleaned sequences that had the given tag
#		yyy_NOTAG.fas 1 fasta file with all of the cleaned sequences that did not have tag


  use warnings;
  use strict;
  use subprogramQDD;

my $syst = 'win';
my $fas_file = 'D:/QDD2.3/data_qual/dr_test.fas'; # file including path
my $outfolder = 'D:/QDD2.3/dataout/';
my $blast_path = 'C:/Program Files/NCBI/blast-2.2.25+/bin/';
my $blast_version = 'blast+';
my $Evector = '1E-3';
my $adapt_file = 'D:/QDD2.3/data/adapter.fas'; # file with path
my $outfile = 'D:/QDD2.3/dataout/dr_test.wov';
my $out_length = 'D:/QDD2.3/dataout/dr_test_wov_length.txt';
my $limit_l = 80;
my $num_threads =1;

my $qual_file = 'D:\QDD2.3\data_qual\dr_test.qual';
my $out_qual = 'D:/QDD2.3/dataout/dr_test_out.qual';
my $qual = 1;

if ( scalar@ARGV == 14 )
{
   	$syst = $ARGV[0];
   	$fas_file = $ARGV[1]; # input file
   	$outfolder = $ARGV[2]; # otpout fasta file wo vector/adapter
   	$blast_path = $ARGV[3];
   	$blast_version = $ARGV[4];
   	$Evector = $ARGV[5];
   	$adapt_file = $ARGV[6];
    $outfile = $ARGV[7];
	$out_length = $ARGV[8]; # output file with number of bases cut from each sequence
	$limit_l = $ARGV[9];
	$num_threads = $ARGV[10];
	$qual = $ARGV[11];
	$qual_file = $ARGV[12];
	$out_qual = $ARGV[13];
#print "cut_vector: Correct number of paramters\n";
}

if ($qual)
{
	open(Q, ">$out_qual") or die "Cannot open output quality file ($out_qual)\n";
}

my %bl_beg = ();
my %bl_end = ();
my %seq = ();
my %adapt = ();

my $fas_root = $fas_file;
$fas_root =~  s/.*[\\\/]//;
my $out_blast = $outfolder.$fas_root.'.blast';
my $out_noadapt = $outfolder.$fas_root.'.noadapt';


# blast sequence file againts the adaptor file; filter is switched off, 
#to be able to clean even low complexity sequences
call_blast_adapter_no_filter_1($fas_file, $out_blast, $adapt_file, $Evector, $blast_path, $blast_version, $syst, $num_threads);


%seq = read_hach_from_fasta_code_wo_space($fas_file);
my %seq_l = length_sequences(\%seq);

open_file_in($out_blast, 'BL');
while(my $line = <BL>)
{
	my $code = '';
	my @line = split(' ', $line);
	if ($line[6] < $seq_l{$line[0]} - $line[7])
	{
		if (exists $bl_beg{$line[0]})
		{
			if ($line[7]> $bl_beg{$line[0]})
			{
				$bl_beg{$line[0]} = $line[7];
			}
		}
		else
		{
			$bl_beg{$line[0]} = $line[7];
		}
	}
	else
	{
		if (exists $bl_end{$line[0]})
		{
			if ($line[6] < $bl_end{$line[0]})
			{
				$bl_end{$line[0]} = $line[6];
			}
		}
		else
		{
			$bl_end{$line[0]} = $line[6];
		}
	}

}
close BL;

foreach my $end (keys %bl_end)
{
	my $seq_t = substr($seq{$end}, 0, $bl_end{$end}-1);
	$seq{$end} = $seq_t;
}
foreach my $beg (keys %bl_beg)
{
	my $l = length $seq{$beg};
	if ($l - $bl_beg{$beg}> 0)
	{
		my $seq_t = substr($seq{$beg}, $bl_beg{$beg}, $l - $bl_beg{$beg});
		$seq{$beg} = $seq_t;
	}
	else
	{
		$seq{$beg} = '';
	}
}


open_file_out($outfile, 'FAS');
open_file_out($out_length, 'OUT');
print OUT "code\torig_length\tbases_cut_beg\tbases_cut_end\tlength_wov\n";

my %qual = ();
if ($qual)
{
	%qual = read_hach_from_qual_code_wo_space($qual_file);
}
foreach my $seq (sort keys %seq)
{
	my $base_cut_beg = 0;
	my $base_cut_end = 0;
	if (exists $bl_end{$seq})
	{
		$base_cut_end = $seq_l{$seq} - $bl_end{$seq} +1;
	}
	if (exists $bl_beg{$seq})
	{ 
		$base_cut_beg = $bl_beg{$seq};
		if (length $seq{$seq} >= $limit_l)
		{
			print FAS ">$seq\n";
			print FAS cut_up_fasta_line($seq{$seq}, '100');
			if($qual)
			{
				if($base_cut_end >0)
				{
					$qual{$seq} =~ s/([0-9]+\s+){$base_cut_end}$//;
				}
				$qual{$seq} =~ s/^([0-9]+\s+){$base_cut_beg}//;
				print Q ">$seq\n$qual{$seq}\n";
			}
		}
	}
	my $l = length $seq{$seq};
	print OUT "$seq\t$seq_l{$seq}\t$base_cut_beg\t$base_cut_end\t$l\n";
}
close OUT;
close FAS; 
if($qual)
{close Q;}


exit;

