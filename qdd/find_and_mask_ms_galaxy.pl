#!/usr/bin/perl -w

# Nom du programme:find_and_mask_ms_galaxy.pl (modified from mask_ms.pl and extract_perfect_microsat_with_hach.pl)

#INPUT: xxx.fas (file with all sequences in it (adaptor already cut))

# Finds all microsat of a sequence 
# $minx is the minimum number of repeats for all motifs, where X is the lenght of the motif (e.g. $min2 is the minimum mumber of dinucleotide repeats)

# Put sequences into fasta file 
# works with multiple files with multiple sequences in it
# ONE FILE PER SPECIES

#OUTPUT:
#xxx_mask.fas # if $keep_non_ms == 1 writes all sequences into fasta file, If 0, writes only sequences with ms in it in fasta file (ms is masked by lowercase))

# xxx.ms (info on seq)
# [0] SEQ_CODE
# [1] NUMBER OF MS
# [2] SEQ_LENGTH
# [3] MS MOTIF
# [4] MS FIRTS POS
# [5] MS LAST POS
# [6] NUMBER OF REP
# [3-6] REPEATED FOR ALL MS



  use warnings;
  use strict;
  use subprogramQDD;


my $min1 = 0;
my $min2 = 0;
my $min3 = 0;
my $min4 = 0;
my $min5 = 0;
my $min6 = 0;
my $del_prev_mask = 0;
my $keep_non_ms = 0;
my $input_file = '';
my $output_fasta = '';
my $output_ms = '';
my $write_ms_file = '';
my $write_fas_file = '';

 if ( scalar@ARGV == 13 )
 {
  	$min1 = $ARGV[0];
	$min2 = $ARGV[1];
	$min3 = $ARGV[2];
	$min4 = $ARGV[3];
	$min5 = $ARGV[4];
	$min6 = $ARGV[5];
	$del_prev_mask = $ARGV[6];
	$keep_non_ms = $ARGV[7];
	$input_file = $ARGV[8];
	$output_fasta = $ARGV[9];
	$output_ms = $ARGV[10];
	$write_ms_file = $ARGV[11];
        $write_fas_file = $ARGV[12];
 }

else
{
	$min1 = 100000;
	$min2 = 5;
	$min3 = 5;
	$min4 = 5;
	$min5 = 5;
	$min6 = 5;
# if $del_prev_mask = 0; keeps previous lower case masking, else deletes previous masking
	$del_prev_mask = 1;
# if 1 writes all sequences info fasta file and ms file, If 0, writes only sequences with ms in it in fasta file (ms is masked by lowercase)) and info on ms composition to ms file
	$keep_non_ms = 1;
	$input_file = '/home/pri/galaxy-dist/tools/qdd/datain/sample.fas';
	$output_ms = '/home/pri/galaxy-dist/tools/qdd/dataout/sample.fas.ms_info.txt';
	$output_fasta = '/home/pri/galaxy-dist/tools/qdd/dataout/sample_ms.fas';
	$write_ms_file = 1; # if 1 ms_file is produced
  	$write_fas_file = 1; # if 1 masked fasta file is produced
 	print "Direct parameteres are used in find_and_mask_ms!!\n";
}


  if ($write_ms_file ==1)
  {
   open_file_out($output_ms, 'SEQ');
   print SEQ "seq_code;number_of_ms;seq_length;ms_motif;ms_first_pos;ms_last_pos;number_of_rep\n";
  }

  if ($write_fas_file ==1)
  {
	open_file_out($output_fasta, 'FAS');
  }


  open_file_in($input_file, 'DATA');
  while (my $seq = one_seq($input_file))
  {
	if ($seq eq '>')
	{
		next;
	}
	$seq =~ /.+\n/;
	my $code = $&;
	$seq =~ s/$code//;
	$seq =~ s/>//;
	$seq =~ s/\s//g;
	chomp $seq;
	$code =~ s/\s.*//;
	chomp $code;

  	my $seq_length = length $seq;
	my %ms_motif = ();
	my %ms_pos = ();

	if ($del_prev_mask == 1)
	{
		$seq = uc $seq;
	}
#	print "$seq\n";
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
	find_microsat_gb_motif_def_6($seq, \%ms_motif, \%ms_pos, $min1, $min2, $min3, $min4, $min5, $min6);
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
   	my @ms = sort {$a <=> $b} keys %ms_pos;
    	my $ms_numb = scalar @ms;
	if ($write_fas_file == 1)
	{
		if (scalar@ms > 0)
		{
			my $fragment_l = change_ms_into_lower_case($seq, \%ms_pos);
			print FAS ">$code\n";
			print FAS cut_up_fasta_line($fragment_l, '100');
		}
		elsif ($keep_non_ms == 1)
		{
			print FAS ">$code\n";
			print FAS cut_up_fasta_line($seq, '100');
		}
	}

#print info into seq file
	if ($write_ms_file == 1 and ($ms_numb > 0 or $keep_non_ms))
	{
		print SEQ "$code;$ms_numb;$seq_length;";
   	 	foreach my $ms (@ms)
    		{
			my $ms_length = $ms_pos{$ms} - $ms +1;
			my $rep_length = length $ms_motif{$ms};
			my $rep_numb = $ms_length/$rep_length;
			print SEQ "$ms_motif{$ms};$ms;$ms_pos{$ms};$rep_numb;";
    		}
    		print SEQ "\n";   
	}

# while
  }
  close DATA;
  if ($write_ms_file ==1)
  {
  close SEQ;
  }
    if ($write_fas_file ==1)
  {
  close FAS;
  }

exit;


############################################################


sub one_seq
{
#    read the data concernant one sequence

    my ($filename) = @_; 
    my ($res_seq) = ''; 

    my ($sep_entree) = $/;
    $/ = ">";
    $res_seq = <DATA>;
    $/ = $sep_entree;
    return $res_seq;
}


##########################################################
