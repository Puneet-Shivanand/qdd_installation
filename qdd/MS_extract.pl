#!/usr/bin/perl -w

# Nom du programme:MS_extract.pl
# emese.meglecz@imbe.fr

#INPUT: Fasta file(s)
#	 	-Accept several sequences in the same file
#		-Works with multiple files placed into a folder (if $batch ==1), or just a single input file (default; if $batch ==0)
#		-Can treat very long sequences (chr) 
#			if sequence length is greater than $sl_limit, the sequence is cut up into $sl_limit length fragments =>
#			(only infuences the output if there is MS in the junction between fragments. It is a very minor bias, since the default of $sl_limit is 100 000 000, and can be set higher if you have enough memory)
#		-Lines in the input fasta files must be shorter then $sl_limit


# Finds all pure microsats in sequences 
# $minx is the minimum number of repeats for all motifs, where X is the lenght of the motif (e.g. $min2 is the minimum number of repeats for dinucleotide motifs)
# pooles microsatellites that are closer then 2*$flank
# extract sequences with a single or pooled MS in the middle and $flank bases on both sides 
# if there is a more than 100 Ns in a fragment with MS, the the fragment is cut up and Ns are deleted.
# 	if the MS is at the beining or and at the end of the sequence and the 5' or 3' flanking region is shorter than $flank, the other flanking region is set longer, to have 2*$flank = 5'flank + 3'$flank


# if $del_prev_mask = 1 => deletes all previous soft masking (whole sequence is put into upper case)
# $folder (-in) if $batch ==0 (default) => name of the input file, 
#				if $batch ==1 name of input folder with files to be analysed (batch submission); in this case only files with $file_motif in the name are analysed
# $outfolder (-dir) output folder; must exists before starting the run

#OUTPUT:
#-csv file: sequence ID; fisrt position of the extracted fragment in which the ms is found; MS motif; first position of MS; last position of MS; number of repeats
#-fasta file with fragments with $flank flanking regions in both sides; if ms are closer than 2*flank bases, than they are kept in the same fragment;
# 	if the MS is at the beining or and of the end of the sequence and the 5' or 3' fkanking region is shorter than $flank, the other flanking region is set longer, to have 2*$flank = 5'flank + 3'$flank
#	MS masked in LC letters (soft masking); 

# log.txt file with the input parameters


  use warnings;
  use strict;
  use Data::Dumper;
  
# defaut settings 
my %par = define_defaults();
###########################################
sub define_defaults
{
	my %par = (-min1 => 10000, 
				-min2 => 5, 
				-min3 => 5, 
				-min4 => 5, 
				-min5 => 5, 
				-min6 => 5, 
				-in => 'D:\Tribolium\prg\datain\test.fas',
				-dir => 'dataout',
				-dmask => 1,
				-batch => 0,
				-out_fas => '', # outfile name Use only if not batch submisson
				-out_csv => '',  # outfile name Use only if not batch submisson
				-length_limit => 30,
				-flank_length => 200
				);
	return %par;
}
########################################
if (scalar @ARGV)
{
	read_params(\@ARGV, \%par);
}
else
{
	print_usage(\%par);
}

#Parameteres that can be defined from the commande line : 
my $batch = $par{-batch}; # if 0 -in is a the name of a single input file; if 1, -in is the name of a folder that contains all input files to be treated
my $folder = $par{-in}; #-in; input filename OR input folder name if several file are to be run, all of them should be in the same $folder; the prg, treats all files in the $folder, that has $file_motif in the filename
my $outfolder = $par{-dir}; #-dir; name of the folder to put the output files; $outfolder must exists before running the program
my $min1 = $par{-min1}; # -min1; Min number of repetitions for homopolimers; if>999 if do not want to search for homopolymers
my $min2 = $par{-min2}; #-min2; Min number of repetitions for dinucleotide motifs;
my $min3 = $par{-min3}; #-min3; Min number of repetitions for trinucleotide motifs;
my $min4 = $par{-min4}; #-min4; Min number of repetitions for tetranucleotide motifs;
my $min5 = $par{-min5}; #-min5; Min number of repetitions for pentanucleotide motifs;
my $min6 = $par{-min6}; #-min6; Min number of repetitions for hexanucleotide motifs;
my $del_prev_mask = $par{-dmask}; # if $del_prev_mask == 0 => keeps previous lower case masking in input file, else => deletes previous masking (whole sequence is put UC prior MS detection)

#Parameteres that cannot be defined from the comand line : 
my $file_motif = ''; #the prg, treats all files in the $folder, that has $file_motif in the filename (relevant only is batch submision)
my $sl_limit = 100000000; # max seq_length treated at the same time; sequences longer than $sl_limit are cut up => there is no important increase in runtime between 1 000 000 and  100 000 000 bases
my $msflank = $par{-flank_length};

my $time_beg = time;

my $log = $outfolder.'/log.txt';
unless (open(LOG, ">$log"))
{print "Cannot open $log\n";exit;}
print LOG "SCRIPT: MS_extract.pl\n\n";
print LOG "min1 = $min1\n";
print LOG "min2 = $min2\n";
print LOG "min3 = $min3\n";
print LOG "min4 = $min4\n";
print LOG "min5 = $min5\n";
print LOG "min6 = $min6\n";
print LOG "Delete previous mask = $del_prev_mask\n";
print LOG "Input file/folder = $folder\n";
print LOG "Output folder = $outfolder\n";
print LOG "Batch submission = $batch\n";
if ($batch)
{
	print LOG "File motif = $file_motif\n";
}
print LOG "Sequence length limit in one cycle = $sl_limit\n";
print LOG "\n";
print LOG "Analysed files:\n";

print "SCRIPT: MS_extract.pl\n\n";
print "min1 = $min1\n";
print "min2 = $min2\n";
print "min3 = $min3\n";
print "min4 = $min4\n";
print "min5 = $min5\n";
print "min6 = $min6\n";
print "Delete previous mask = $del_prev_mask\n";
print "Input file/folder = $folder\n";
print "Output folder = $outfolder\n";
print "Batch submission = $batch\n";
if ($batch)
{
	print "File motif = $file_motif\n";
}
print "Sequence length limit in one cycle = $sl_limit\n";
print "\n";
print "Analysed files:\n";

my @filenames = ();
if ($batch == 0) # if only one input file
{
  $folder=~ /(.*[\\\/])([^\\^\/]*)/;
  my $file = $2;
  $folder = $1;
 # print "$file, $folder\n";
  push(@filenames, $file); 
}
else # if batch submission
{
 unless ( opendir(FOLDER, $folder) )
 {
   print "Cannot access to folder $folder\n";
   exit;
 }
@filenames = grep ( !/^\.\.?$/, readdir(FOLDER) );
closedir(FOLDER);
}

foreach my $filename (sort @filenames)
{
if ($filename =~ /$file_motif/)
{
  print "$filename\n";
  print LOG "$filename\n";
 
    my $tbl_file = $par{-out_csv};
    if($tbl_file eq '' or $batch == 1)
    {
    	$tbl_file = $filename;
   	$tbl_file = $outfolder.'/'.$tbl_file; 
   	$tbl_file =~ s/\..*//;
	$tbl_file .= '.csv';
    }
    unless (open(TBL, ">$tbl_file"))
    {print "Cannot open $tbl_file\n";exit;}
    print TBL "seq_code;first_position_of_the_fragment;MS_motif;first_pos;last_pos;repeat_number\n";

    my $fas_file = $par{-out_fas};
    if($fas_file eq '' or $batch == 1)
    {
     	$fas_file = $filename;
    	$fas_file = $outfolder.'/'.$fas_file;
    	$fas_file =~ s/\..*//;
	$fas_file .= '_mask.fas';
    }

    unless (open(FAS, ">$fas_file"))
    {print "Cannot open $fas_file\n";exit;}
  
  $filename = $folder.'/'.$filename; 
  unless (open(DATA, $filename))
  {print "Cannot open $filename\n";exit;}

  my $treated_pos = 0;
  my $code = '';
  my $prev_code = '';
  my $left_over = ''; # either the def line of the next sequence, or the part of the sequences that is read but makes the fragment longer than $sl_limit
  while (my $seq = get_seq_frag($sl_limit, \$code, \$left_over)) # takes the next sequence or if the sequence is longer then $sl_limit, than the next $sl_limit bases of a sequence
  {
	unless ($code eq $prev_code)
	{
		$treated_pos = 0;
	}
	my $fr_beg = $treated_pos +1;
	my $fr_end = $treated_pos + length $seq;
    my %ms_motif = ();
    my %ms_pos = ();
  
    if ($del_prev_mask == 1)
    {
      $seq = uc $seq;
    }
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
    find_microsat_gb_motif_def_6($seq, \%ms_motif, \%ms_pos, $min1, $min2, $min3, $min4, $min5, $min6);
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
    my @ms = sort {$a <=> $b} keys %ms_pos;
    my $ms_numb = scalar @ms;

	if (scalar@ms > 0)
	{
		$seq = change_ms_into_lower_case($seq, \%ms_pos);
		# pools ms that are closer than 2*$flank and extracts the fragment with ms(s) and flankig regions on both sides.
		my %fragments = extract_ms($seq, \%ms_pos, $msflank, $par{-length_limit});
		foreach my $frag (sort {$a<=>$b} keys%fragments)
		{
			if((length $fragments{$frag})> $par{-length_limit})
			{
			my $fp_fr = $frag + $treated_pos;
			print FAS ">$code", '_', "$fp_fr\n";
			print FAS cut_up_fasta_line($fragments{$frag}, '100');
			}
		}
		my %ms_fragment = get_fragment_for_each_ms(\%ms_pos, \%fragments); # $ms_fargments{$ms_first_pos} = first position of the fragment containing the ms
		foreach my $ms (@ms)
		{
			if((length $fragments{$ms_fragment{$ms}})> $par{-length_limit}) # check if fragment containing the MS is longer than limit
			{
				my $ms_length = $ms_pos{$ms} - $ms +1;
           			my $rep_length = length $ms_motif{$ms};
           			my $rep_numb = $ms_length/$rep_length;
				my $ms_beg = $ms + $treated_pos;
		   		my $ms_end = $ms_pos{$ms} + $treated_pos;
			
				my $fragment_beg = $ms_fragment{$ms} + $treated_pos;
            			print TBL "$code;$fragment_beg;$ms_motif{$ms};$ms_beg;$ms_end;$rep_numb\n";
			}
		} 
	}

	$treated_pos += length $seq;
	$prev_code = $code;
  } # while; my $seq
  close DATA;
  close TBL;
  close FAS;
# end if file_motif
}
else
{
  print $filename.' file name does not contain "'.$file_motif.'"'.' It is not analysed'."\n";
}
#end of foreach file
}

my $run_time = time - $time_beg;
print "\nThe run took $run_time s\n";
print LOG "\nThe run took $run_time s\n";
close LOG;

exit;

##############################################################
sub get_fragment_for_each_ms
{
 my ($ms_pos, $fragments) = @_; 
# $ms_fargments{$ms_forst_pos} = first position of the fragment containing the ms,($treated positions are added)

my %ms_fragments = ();
#print Dumper(\%$fragments);
#print Dumper(\%$ms_pos);
my $fp = join(';', keys%$ms_pos).';'.join(';', keys%$fragments);
my @fp = sort {$a <=> $b} split(';', $fp);
#print "@fp\n";
for( my $i = 1; $i<scalar@fp; ++$i)
{
	if(exists $$ms_pos{$fp[$i]})
	{
		for(my $j =$i-1; $j>=0; --$j)
		{
			if (exists $$fragments{$fp[$j]})
			{
				$ms_fragments{$fp[$i]} = $fp[$j];
				$j = -1;
			}
		}
	}
} 
return %ms_fragments;
#print Dumper(\%ms_fragments);
}

############################################################

sub extract_ms
{
my ($seq, $ms_pos, $flank, $length_limit) = @_;

my %pool_ms = %$ms_pos;
# pool ms closer then 2*$flank bases
my @fp = sort {$a<=>$b} keys %pool_ms;
#print Dumper(\%pool_ms);
for(my $i = scalar@fp-1; $i>0; --$i)
{
	if (($fp[$i] - $pool_ms{$fp[$i-1]} -1)<2*$flank )
	{
		$pool_ms{$fp[$i-1]} = $pool_ms{$fp[$i]};
		delete$pool_ms{$fp[$i]};
	}
}
#print Dumper(\%pool_ms);

my %fragments = ();
foreach my $ms (sort {$a<=>$b}keys %pool_ms)
{
	my $beg = $ms -$flank;
	my $beg_missing = 0;
	if($beg <1)
	{
		$beg_missing = ($beg * -1) +1 ;
		$beg = 1;
	}
	my $end = $pool_ms{$ms} + $flank + $beg_missing;
	if($end > length $seq)
	{
		my $end_missing = $end - length $seq;
		$end = length $seq;
		$beg = $beg - $end_missing;
		if ($beg <1)
		{
			$beg = 1;
		}
	}
	my $l = $end - $beg +1;
	my $fragment = substr($seq, $beg-1, $l);
# eliminate Ns from the begigning of the fragment
	if($fragment =~ /^N+/i)
	{
		my $b = length $&;
		$beg = $beg +$b;
		$fragment =~ s/^N+//i;
	}
# eliminate Ns from the end of the fragment
	if($fragment =~ /N+$/i)
	{
		$fragment =~ s/N+$//i;
	}
# cut up fragment if there is more than 100 N in it
	while ($fragment =~ /N{100,}/gi)
	{
		my $pos = pos($fragment);
#		print "$pos\n";
		my $sub_fr = substr($fragment, 0, $pos, '');
		$sub_fr =~  s/N+$//;
		my %ms_mot = ();
		my %ms_p = ();
		find_microsat_gb_motif_def_6($sub_fr, \%ms_mot, \%ms_p, $min1, $min2, $min3, $min4, $min5, $min6);
# print out sub_fragment only if there is MS in it
		if ((scalar keys %ms_mot) > 0)
		{
			$fragments{$beg} = $sub_fr;
		}
		$beg += $pos;
	}
	$fragments{$beg} = $fragment;
}
#print Dumper(\%fragments);
return %fragments;
} 


############################################################

sub get_seq_frag
{
my ($sl_limit, $code, $left_over) = @_;
my $seq = '';
my $l = 0;

if ($$left_over =~ />/) # the last line read but not yet treated was a def line
{
	$$code = $$left_over;
	$$code =~ s/>//;
	$$code =~ s/\s.*//;
}
else
{
	$$left_over =~ s/\s//g;
	$seq = $$left_over;
	$l += length $seq;
}


while ($l < $sl_limit and my $line = <DATA>)
{
#	my $line = <DATA>;
	chomp $line;
	if ($line =~ />/)
	{
		if ($seq eq '')
		{
			$$code = $line;
			$$code =~ s/>//;
			$$code =~ s/\s.*//;
		}
		else
		{
			$$left_over = $line;
#			last;
			return $seq;
		}
	}
	else
	{
		$line =~ s/\s//g;
		$seq .= $line;
		$l += length $line;
	}

}

if ($l>= $sl_limit)
{
	$$left_over = substr($seq, $sl_limit, $l-$sl_limit, '');
}
else
{
	$$left_over = '';
}


return $seq;
}


##########################################################

sub get_lines
{
 my ($line_n, $treated_pos) = @_;
 
 my $seq = '';
 for (my $i =0; $i < $line_n; ++$i)
 {
    
    if (my $line = <DATA>)
    {
      if ($line =~ />/)
      {
        print "More than one sequences in input file\n This script is not suitable for this file\n";
        exit;
      }
      chomp $line;
      $seq .= $line;
      $$treated_pos = $$treated_pos + length $line;
    }
    else
    {
      $i = $line_n;
    }

 }
 return $seq;
}

##############################################################

sub find_microsat_gb_motif_def_6
{
    my ($seq, $ms_motif, $ms_pos, $min_single, $min2, $min3, $min4, $min5, $min6) = @_;
    my $bp = length $seq;

    for (my $i= 0; $i < $bp; ++$i)
    {
#	print "$i ";
#Check for motif of single base pair
	my $motif = substr($seq, $i, 1);
	if ( $motif =~ /N/i )
	{
	    next;
	}
	my $j = $i+1;
	if ($j > $bp-1 )
	{
	    next;
	}
	while ( substr($seq, $j, 1) eq $motif )
	{
	    ++$j;
	}
	if ($j - $i > ($min_single-1))
	{
	    $$ms_pos{$i+1}= $j;
	    $$ms_motif{$i+1}= $motif;
	    $i = $j-1;
	    next;
	}
	elsif ($j - $i > 6)
	{
		$i = $j-6;
		next;
	}

#Check for motif of two base pair
	$motif = substr($seq, $i, 2);
	if ( $motif =~ /N/i )
	{
	    next;
	}
	$j = $i+2;
	if ($j > $bp-2 )
	{
	    next;
	}
	while ( substr($seq, $j, 2) eq $motif )
	{
	    $j +=2;
	}
	if ($j - $i > (($min2*2)-1))
	{
	    $$ms_pos{$i+1}= $j;
	    $$ms_motif{$i+1}= $motif;
	    $i = $j-1;
	    next;	    
	}
	elsif ($j - $i > 6)
	{
		$i = $j-6;
		next;
	}

#Check for motif of three base pair
	$motif = substr($seq, $i, 3);
	if ( $motif =~ /N/i )
	{
	    next;
	}
	$j = $i+3;
	if ($j > $bp-3 )
	{
	    next;
	}
	while ( substr($seq, $j, 3) eq $motif )
	{
	    $j +=3;
	}
	if ($j - $i > (($min3*3)-1))
	{
	    $$ms_pos{$i+1}= $j;
	    $$ms_motif{$i+1}= $motif;
	    $i = $j-1;
	    next;   
	}
	elsif ($j - $i > 6)
	{
		$i = $j-6;
		next;
	}

#Check for motif of 4 base pair
	$motif = substr($seq, $i, 4);
	if ( $motif =~ /N/i )
	{
	    next;
	}
	$j = $i+4;
	if ($j > $bp-4 )
	{
	    next;
	}
	while ( substr($seq, $j, 4) eq $motif )
	{
	    $j +=4;
	}
	if ($j - $i > (($min4*4)-1))
	{
	    $$ms_pos{$i+1}= $j;
	    $$ms_motif{$i+1}= $motif;
	    $i = $j-1;
	    next;	    
	}

#Check for motif of 5 base pair
	$motif = substr($seq, $i, 5);
	if ( $motif =~ /N/i )
	{
	    next;
	}
	$j = $i+5;
	if ($j > $bp-5 )
	{
	    next;
	}
	while ( substr($seq, $j, 5) eq $motif )
	{
	    $j +=5;
	}
	if ($j - $i > (($min5*5)-1))
	{
	    $$ms_pos{$i+1}= $j;
	    $$ms_motif{$i+1}= $motif;
	    $i = $j-1;
	    next;	    
	}

#Check for motif of 6 base pair
	$motif = substr($seq, $i, 6);
	if ( $motif =~ /N/i )
	{
	    next;
	}
	$j = $i+6;
	if ($j > $bp-6 )
	{
	    next;
	}
	while ( substr($seq, $j, 6) eq $motif )
	{
	    $j +=6;
	}
	if ($j - $i > (($min6*6)-1))
	{
	    $$ms_pos{$i+1}= $j;
	    $$ms_motif{$i+1}= $motif;
	    $i = $j-1;	 
	    next;
	}

    }
    
#    print "\n\n";

}


###################################################################"
sub change_ms_into_lower_case
{
my ($fragment, $ms_pos) = @_;

	my @ms = sort {$a <=> $b} keys %$ms_pos;
	foreach my $ms (@ms)
	{
#	print "$ms;$$ms_pos{$ms}\n";
		$$ms_pos{$ms} =~ /[0-9]+/;
		my $l_ms = $& - $ms +1;
		my $fragment_temp = $fragment;
		my $ms_frag = substr($fragment_temp, $ms-1, $l_ms);
		$ms_frag =~ tr/ACGT/acgt/;
#	print "$ms_frag\n";
		substr($fragment, $ms-1, $l_ms, $ms_frag);
#	print	"$fragment\n";
	}
return $fragment;
}

#######################################################
sub open_file_and_hach
{
# reads each line into a hachache. First element before first ';' is the code all other in the lines are values, separates element of the line by ;
    my ($filename, $first_line) = @_;

	  unless ( open (FILE, $filename))
	  { print "Cannot open $filename\n"; exit;}
    my @data = <FILE>;
    if ($first_line == 0)
    {
      shift (@data);
    }
	  close FILE;

	  my %hach = ();
	  foreach my $ligne (@data)
    {
      chomp $ligne;
      $ligne =~ /.+?;/;
      my $code = $&;
      $ligne =~ s/$code//;
      $code =~ s/;//;
      $hach{$code} = $ligne;
    }
    return %hach;
}

############################################################
sub read_params
{
my ($arg_ref, $par) = @_;

 unless(scalar@$arg_ref%2 == 0) # if odd number of arguments
 {
	print_usage($par);
	exit;
 }
 
 for(my $i = 0; $i<@$arg_ref; $i=$i+2)
 {
	if (exists $$par{$$arg_ref[$i]})
	{
		$$par{$$arg_ref[$i]} = $$arg_ref[$i+1];
	}
	else
	{
		print_usage($par);
		exit;
	}
 }
 if ($$par{-in} eq '' or $$par{-dir} eq '')
 {
	print "\nERROR: Please, specify the input file (-in) and the output directory (-dir)\n\n";
 		print_usage($par);
		exit;
 }
 unless(-e $$par{-in})
 {
	print "The input file/directory $$par{-in} does not exist\n";
	exit;
 }
 unless(-e $$par{-dir})
 {
	print "The output directory $$par{-dir} does not exist\n";
	exit;
 }
# print Dumper(\%$par);
}

#############################################

sub	print_usage
{
my ($par) = @_;
 
print "\n\nMS_find.pl\n Search the input fasta file for pure microsatellites of 1-6 nucleotide motifs\n\nOUTPUT:\n -csv file with microsatellite positons and motifs\n -fasta file with soft masked microsatellites\n
USAGE:\n perl MS_find.pl -in filename -dir directory [OPTIONS]\n
OPTIONS:\n
	-in (compulsary): The name of the input file\n
	-dir (compulsary): Name of the directory for output files (must exist before running the script)\n  
	-min1 (optional): Minimum length of homopolymers\n\tDEFAULT: $$par{-min1}\n
	-min2 (optional): Minimum number of repeats for dinucleotide repeats\n\tDEFAULT: $$par{-min2}\n
	-min3 (optional): Minimum number of repeats for trinucleotide repeats\n\tDEFAULT: $$par{-min3}\n
	-min4 (optional): Minimum number of repeats for tertanucleotide repeats\n\tDEFAULT: $$par{-min4}\n
	-min5 (optional): Minimum number of repeats for pentanucleotide repeats\n\tDEFAULT: $$par{-min5}\n
	-min6 (optional): Minimum number of repeats for hexanucleotide repeats\n\tDEFAULT: $$par{-min6}\n
	-dmask (optional): If 1, delete already existing softmasking in the input file\n\tDEFAULT: $$par{-dmask}\n
	-batch (optional): If 1, -in is the name of folder that contains all fasta files to be analysed; if 0, -in is the name of the single input file\n\tDEFAULT: $$par{-batch}\n
	\n";

	exit;
}

#############################################

sub cut_up_fasta_line
{
 my ($seq, $l) = @_;
 my @temp = ();
 for(my $i = 0; $i<length $seq; $i = $i + $l)
 {
	my $frag = substr($seq, $i, $l);
	$frag .= "\n";
	push(@temp, $frag);
 }
 return @temp;
}
