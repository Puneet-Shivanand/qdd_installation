#!/usr/bin/perl -w

# name of the program: sort_tag.pl
# sort sequences in a fasta file according to a list of tags

# reads all, sequences from a fasta file and identify which of the given tag is present, removes tag (and adaptor)
# writes one fasta file per tag with all sequences that had the given tag.
# if blast hit is shorter then adaptor, removes the adaptor length not the length of the hit

# INPUT: 	xxx.fas (fasta file with tags)
# 		yyy.fas fasta file to be cleaned and sorted out.

# OUPUT: yyy_m8.tbl results of blast againts te adapter file
#		yyy_tagY.fas 1 fasta file with all of the cleaned sequences that had the given tag
#		yyy_NOTAG.fas 1 fasta file with all of the cleaned sequences that did not have tag


  use warnings;
  use strict;
  use subprogramQDD;

my $syst = 'win';
my $fas_file = 'D:\QDD2.3\Ea\Ea_tag.fas'; # file with path
my $outfolder = 'D:\QDD2.3\Ea_out/';
my $Evector = '1E-1';
my $word_size = '7';
my $blast_path = 'C:\Program Files\NCBI\blast-2.2.25+\bin/';
my $blast_version = 'blast+';
my $tag_file = 'D:\QDD2.3\data\tag.fas'; # file wo path
my $tempfolder = 'temp/';
my $num_threads = 1;
my $qual_file = 'D:\QDD2.3\Ea\Ea_tag.qual';
my $qual = 1;

  
if ( scalar@ARGV ==  12)
{
   	$syst = $ARGV[0];
   	$fas_file = $ARGV[1];
   	$outfolder = $ARGV[2];
   	$blast_path = $ARGV[3];
   	$blast_version = $ARGV[4];
   	$Evector = $ARGV[5];
   	$word_size = $ARGV[6];
   	$tag_file = $ARGV[7];
	$tempfolder = $ARGV[8];
	$num_threads =  $ARGV[9];
	$qual_file = $ARGV[10];
	$qual = $ARGV[11];
}

my %adapt = ();
my %seq = ();
my %sort = ();
my %qual = ();

my $fas_file_root = get_filename_root($fas_file);
my $outfile = $outfolder.$fas_file_root.'_NOTAG.fas';
my $outqual = $outfolder.$fas_file_root.'_NOTAG.qual';
my $out_blast = $tempfolder.$fas_file_root.'_m6.tbl';

 
# blast sequence file againts the adaptor file; filter is switched off, *
#to be able to clean even low complexity sequences
call_blast_adapter_no_filter_word_size_1($fas_file, $out_blast, $tag_file, $Evector, $blast_path, $blast_version, $word_size, $syst, $num_threads);
delete_tag_adapt_files($tag_file, $syst);

%seq = read_hach_from_fasta_code_wo_space($fas_file);
if($qual)
{
	%qual = read_hach_from_qual_code_wo_space($qual_file);
}
my %seq_l = length_sequences(\%seq);
%adapt = read_hach_from_fasta_code_wo_space($tag_file);
my %adapt_l = length_sequences(\%adapt);

open_file_in($out_blast, 'BL');
while(my $line = <BL>)
{
	my @line = split(' ', $line);
# hit sur brin +, hit at the begining of the sequence (+2 for to allow for 1 seq error (indel))
	if ($line[8] < $line[9] and $line[7] < $adapt_l{$line[1]} +2)
	{
		$sort{$line[1]} .= $line[0].';';
		if ($line[7]> $adapt_l{$line[1]})
		{
			$seq{$line[0]} = substr($seq{$line[0]}, $line[7], $seq_l{$line[0]} - $line[7]);
			if($qual)
			{
			$qual{$line[0]} =~ s/^([0-9]+\s+){$line[7]}//;
			}
		}
		else
		{
			$seq{$line[0]} = substr($seq{$line[0]}, $adapt_l{$line[1]}, $seq_l{$line[0]} - $adapt_l{$line[1]});
			if($qual)
			{
			$qual{$line[0]} =~ s/^([0-9]+\s+){$adapt_l{$line[1]}}//;
			}
		}
	}
}

close (BL);

#print_hachage(\%sort);

foreach my $adapt (keys %sort)
{
	my @list = split(';', $sort{$adapt});
	my $out = $outfolder.'/'.$fas_file_root.'_'.$adapt.'.fas';
	my $out_qual = $outfolder.'/'.$fas_file_root.'_'.$adapt.'.qual';
	open_file_out($out, 'OUT');
	if($qual)
	{
	open_file_out($out_qual, 'OQ');
	}
	my $i = 0;
	foreach my $code (@list)
	{
		if ($seq{$code} eq '')
		{
			print "In $code an adapter is already removed\n";
		}
		else
		{
			print OUT ">$code\n";
			print OUT cut_up_fasta_line($seq{$code}, '100');
			if($qual)
			{
				print OQ ">$code\n$qual{$code}\n";
			}
			$seq{$code} = '';
			++$i;
		}
	}
	close OUT;
	if($qual)
	{
		close OQ;
	}
	print "Number of sequences with $adapt tag: $i\n";
}
open_file_out($outfile, 'OUT');
if($qual)
{
	open_file_out($outqual, 'OQ');
}

# prints a file with the non_cleaned sequences
my $i = 0;
foreach my $code (keys %seq)
{
	unless ($seq{$code} eq '')
	{
		print OUT ">$code\n";
		print OUT cut_up_fasta_line($seq{$code}, '100');
		if($qual)
		{
			print OQ ">$code\n$qual{$code}\n";	
		}
		++$i;
	}
}
	print "Number of sequences without tag: $i\n";

close (OUT);
if($qual)
{
	close OQ;
}

exit;

