#!/usr/bin/perl -w


  use warnings;
  use strict;
  use subprogramQDD;
  use Data::Dumper;
  
  
#INPUT: 
# fasta file with rerelatively short reads (less than ca.2000 bp)
# fasta file with adapters (if vector/adptor clipping is necessary)

# Vector/adapter clipping if adapter == 1
# selection of sequences longer than length_limit
# microsatellite search

#OUTPUT: 
#  log_file_pipe1 => input parameters + summary on output
#  file_wov => sequence after vector clipping and length selection
#  file_length => info on sequence length : 
	#[0] code; 
	#[1] orig_length (original length); 
	#[2] bases_cut_beg (number of bases cut from the beginning);  
	#[3] bases_cut_end (number of bases cut from the end); 
	#[4] length_wov (length after clipping)
#   file_ms => info on ms number and position : 
	#[0]seq_code (sequence code); 
	#[1]number_of_ms(number of microsatellite in the sequence); 
	#[2]seq_length (length of sequence); 
	#[3]ms_motif (microsatellite motif); 
	#[4]ms_first_pos (first position of microsatellite); 
	#[5]ms_last_pos (last position of microsatellite) ;
	#[6]number_of_rep (number of repetition of microsatellite)
#   file_msfas  => fasta file with sequences which contain microsatellites, and longer than length_limit after vector/adapterclipping



my $min1 = 1000000;
my $min2 = 5;
my $min3 = 5;
my $min4 = 5;
my $min5 = 5;
my $min6 = 5;
my $word_size = 7; 
my $Evector = '1E-3';
my $set_file_default = '/etc/qdd/set_qdd_default.ini';
my $blast_version = 'blast+';

$set_file_default = check_set_file($set_file_default);

print_heading_3_0('STDOUT');
print "SCRIPT pipe1.pl\n";

unless (-e $set_file_default)
{
	$set_file_default = 'set_qdd_default.ini';
	unless(-e $set_file_default)
	{
		print "set_qdd_default.ini should be either in the current working directory, or in /usr/bin/\n";
		exit;
	}
}

my %param = (
  galaxy => '', # one if run on galaxy, 0 for command line
  syst => 'linux',  
  blast_path => '',
  out_folder => '',
  adapter => '',
  length_limit => '',
  del_files => '',
  input_file => '',
  adapter_file => '',
  qdd_folder => '',
  file_log_pipe1 => '_pipe1_log.txt',
  file_wov => '_pipe1_wov.fas',
  file_length => '_pipe1_length_info.tabular',
  file_msfas => '_pipe1_for_pipe2.fas',
  file_ms => '_pipe1_ms.csv',
  outfile_string => '',
  debug => 0,
  num_threads => 1,
  qual => 0,
  qual_file => '',
  file_qual_trim => '_pipe1_trimmed.qual',
  trim_window => 10,
  mean_qs => 20,
  contig => 1,
  flank_length => 200,
  fastq => ''
);

my %param_limits = define_param_limits_pipe1(\%param); #hash of array $param_limits{syst}[type(list, real, folder, string), error message, ['list of values/lower_upper limits)']
read_set_file_galaxy(\%param, $param{qdd_folder}.$set_file_default);
modify_params_from_tags(\%param, \@ARGV);
add_slash_to_folder(\%param, \%param_limits);
my $galaxy = $param{galaxy}; # one if run on galaxy, 0 for command line

my $folder_id = time;
my $date = localtime $folder_id;
my $tempfolder = $param{out_folder}.'pipe1_'.$folder_id.'/';
make_folder($param{syst}, $tempfolder);


if($galaxy)
{
	$param{outfile_string} = 'NA';
}
elsif ($param{outfile_string} eq '')
{
	$param{outfile_string} = get_filename_root($param{input_file});
#	$param{outfile_string} =~  s/_v[0-9]+//;
}


unless($galaxy)
{
check_param_limits(\%param, \%param_limits);
	my @clefs = ('file_wov', 'file_length', 'file_msfas', 'file_ms', 'file_qual_trim');
	if($param{file_log_pipe1} =~ /^_/) # if the name of the log file has not been specified (run in pipeline)
	{
		push(@clefs, 'file_log_pipe1');
	}
	get_last_version_and_modify_1($param{out_folder}, \@clefs, $param{outfile_string}, \%param);

	if($param{file_log_pipe1} =~ /^_/) # if the name of the log file has not been specified (run in pipeline)
	{
		$param{file_log_pipe1} = $param{out_folder}.$param{outfile_string}.$param{file_log_pipe1};
	}
	$param{file_wov} = $param{out_folder}.$param{outfile_string}.$param{file_wov};
	$param{file_length} = $param{out_folder}.$param{outfile_string}.$param{file_length};
	$param{file_msfas} = $param{out_folder}.$param{outfile_string}.$param{file_msfas};
	$param{file_ms} = $param{out_folder}.$param{outfile_string}.$param{file_ms};
	$param{file_qual_trim} = $param{out_folder}.$param{outfile_string}.$param{file_qual_trim};
}

open(my $fh_log, '>>', $param{file_log_pipe1}) or die "Cannot open $param{file_log_pipe1} $!\n";
print_heading_3_0($fh_log);
print $fh_log "SCRIPT pipe1.pl\n";


unless($param{contig})
{
if($galaxy)
{
	if ($param{adapter_file} eq 'None' or $param{adapter_file} eq '')
	{
		$param{adapter} = 0;
	}
	else
	{
		$param{adapter} = 1;
	}

	if ($param{adapter} == 0)
	{
	 print "No adapter.fas file has been given or adapter search has been set to false.\nAdaptor/vector clipping step is skipped\n";
	 print $fh_log "No adapter.fas file has been given or adapter search has been set to false.\nAdaptor/vector clipping step is skipped\n\n";
	}

	if ($param{qual_file} eq 'None' or $param{qual_file} eq '')
	{
		$param{qual} = 0;
	}
	else
	{
		$param{qual} = 1;
	}
}# end if galaxy

if(0)
{
if ($param{qual} == 0)
{
 print "No quality file was provided or qual was set to false.\nThe trimming step based on quality scores is skipped.\n\n";
 print $fh_log "No quality file was provided or qual was set to false.\nThe trimming step based on quality scores is skipped.\n\n";
 if($param{adapter})
 {
	print "If you plan to use quality control in pipe3, you must include a quality file in pipe1 to clip the quality vaules of the adapters\n\n";
	print $fh_log "If you plan to use quality control in pipe3, you must include a quality file in pipe1 to clip the quality vaules of the adapters\n\n";
 }
}
}# end if(0)
}# end unless contig


if($param{debug})
{
	print $fh_log Dumper(\%param);
	print "\n";
}

print $fh_log "INPUT files:\n";
print $fh_log "\tInput file with sequences: $param{input_file}\n";
print "\tInput file with sequences: $param{input_file}\n";
unless($param{contig})
{
if($param{adapter})
{
print $fh_log "Adapter file : $param{adapter_file}\n";
print "\tAdapter file : $param{adapter_file}\n";
}
if($param{qual})
{
print $fh_log "\tQuality score file : $param{qual_file}\n\n";
print "\tQuality score file : $param{qual_file}\n\n";
}
}

print $fh_log "\nOUTPUT files:\n";
print $fh_log "String for naming output files: $param{outfile_string}\n";
print $fh_log "\t$param{file_msfas} : Fasta file with sequences ready for pipe2.pl (i) contain microsatellites (ii) longer than $param{length_limit} after vector/adapter clipping OR after being extracted from contigs\n";
print $fh_log "\t$param{file_ms} : Info on mirosatellite number and position\n";
unless($param{contig})
{
print $fh_log "\t$param{file_wov} : Sequences after vector clipping and length selection\n";
print $fh_log "\t$param{file_length} : Information on sequence length\n";
if($param{qual})
{
print $fh_log "\t$param{file_qual_trim} : Quality file after trimming\n";
}
}



print $fh_log "\nGENERAL PARAMETERS\n";
print $fh_log "System (win/linux): $param{syst}\n";
print $fh_log "Input sequences are already assembled (contigs, scaffolds or chromosomes) (YES=1/NO=0): $param{contig}\n";
print $fh_log "Input file is a fastq file (YES=1/NO=0): $param{fastq}\n";
if($param{contig})
{
print $fh_log "Length of the flanking region on both sides of the microsatellite: $param{flank_length}\n";
}
else
{
print $fh_log "Path to BLAST executables: $param{blast_path}\n";
print $fh_log "Remove adapter (YES=1/NO=0): $param{adapter}\n";
#print $fh_log "Trim sequence extremities based on quality scores (YES=1/NO=0): $param{qual}\n";
if($param{qual})
{
print $fh_log "Minimum mean quality score over 10 bases while trimming sequence extremities: $param{mean_qs}\n";
}
}
print $fh_log "Minimum sequence length: $param{length_limit}\n";
print $fh_log "Delete intermediate files (YES=1/NO=0): $param{del_files}\n";
print $fh_log "Output folder : $param{out_folder}\n";
print $fh_log "Path to QDD executables: $param{qdd_folder}\n";
print $fh_log "Minimum homopolymer length for microsatellite detection: $min1\n";
print $fh_log "Minimum repeat number for dibase motifs for microsatellite detection: $min2\n";
print $fh_log "Minimum repeat number for tribase motifs for microsatellite detection: $min3\n";
print $fh_log "Minimum repeat number for tetrabase motifs for microsatellite detection: $min4\n";
print $fh_log "Minimum repeat number for pentabase motifs for microsatellite detection: $min5\n";
print $fh_log "Minimum repeat number for hexabase motifs for microsatellite detection: $min6\n";
unless($param{contig})
{
print $fh_log "E-value for adapter detection: $Evector\n";
print $fh_log "Word size for tag and adapter detection: $word_size\n";
}

print $fh_log "************************************\n\n";



############################################## start prg

# if input file is fastq, convert it into fatsa, and put the new file into $param{input_file}
if($param{fastq})
{
	print "Converting fastq file into fasta format\n";
	my $fasta_input = $tempfolder.'input_fasta.fas';
	convert_fasq_into_fasta($param{input_file}, $fasta_input);
	$param{input_file} = $fasta_input;
}


if ($param{contig})
{
	print "\nExtracting microsatellites from input sequences with $param{flank_length} flanking region on both sides\n\n";
	my $ext_MS = 'perl "'.$param{qdd_folder}.'MS_extract.pl" -min1 '.$min1.' -min2 '.$min2.' -min3 '.$min3.' -min4 '.$min4.' -min5 '.$min5.' -min6 '.$min6;
	$ext_MS .= ' -in "'.$param{input_file}.'" -dir "'.$tempfolder.'" -dmask 1 -batch 0 -out_fas "'.$param{file_msfas}.'" -out_csv "'.$param{file_ms}.'" -length_limit "'.$param{length_limit}.'" -flank_length "'.$param{flank_length}.'" >"'.$tempfolder.'MS_extract_screen.txt"'.' 2>"'.$tempfolder.'MS_extract_error.txt"';
	system $ext_MS;
	if($param{debug})
	{
		print $fh_log "\n\nRun MS_extract:\n$ext_MS\n\n";
	}

	my $seqn = get_seq_number_from_fasta($param{file_msfas});
	print $fh_log "$seqn fragments are extracted from the input sequences\n";

}# end if contig
else
{

my $temp_qual_file = $tempfolder.'quality_file_after_vector_cut.qual';
# remove adapter
if ($param{adapter})
{
#my $syst = 'linux'; my $fas_file  my $outfolder my $blast_path my $blast_version my $Evector my $adapt_file my $outfile my $out_length 
      my $cut_vector = 'perl "'.$param{qdd_folder}.'cut_vector_galaxy.pl" "'.$param{syst}.'" "'.$param{input_file}.'" "'.$tempfolder.'" "'.$param{blast_path}.'" "'.$blast_version.'" "'.$Evector.'" "'.$param{adapter_file}.'" "'.$param{file_wov}.'" "'.$param{file_length}.'" "'.$param{length_limit}.'" "'.$param{num_threads}.'" "'.$param{qual}.'" "'.$param{qual_file}.'" "'.$temp_qual_file.'"';  
 
      if ($param{debug})
      {print $fh_log "\n$cut_vector\n\n";}
      print "Removing adapters from sequences\n";
      system $cut_vector;
}
else
{
	get_seq_length_remove_short($param{input_file}, $param{file_length}, $param{length_limit}, $param{file_wov});
}

my $seqn = get_seq_number_from_fasta($param{file_wov});
print $fh_log "$seqn sequences are longer than or equal to $param{length_limit} bp after adapter clipping\n";

# trimming low quality sequences ends
my %qual = ();
my $ms_input = $param{file_wov};
my $trimmed_qual = $tempfolder.'trimmed_qual.qual';
my $tbl_trim = $tempfolder.'trimming_info.tabular';
if ($param{qual})
{
	print "Trimming low quality sequence ends\n";
	if ($param{adapter})
	{
		%qual = read_hach_from_qual_code_wo_space($temp_qual_file);
	}
	else
	{
		%qual = read_hach_from_qual_code_wo_space($param{qual_file});
	}
	my %seq = read_hach_from_fasta_code_wo_space($param{file_wov});
	trim_extremities(\%seq, \%qual, $param{trim_window}, $param{mean_qs}, $param{length_limit}, $tbl_trim);
	$seqn = scalar(keys %seq);
	print $fh_log "$seqn sequences are longer than or equal to $param{length_limit} bp after quaility trimming\n";
	my $trimmed_fasta = $tempfolder.'trimmed_fasta.fas';
	write_hash_to_fasta(\%seq, $trimmed_fasta);
	write_hash_to_fasta(\%qual, $trimmed_qual);
	$ms_input = $trimmed_fasta;
}

my $find_ms = 'perl "'.$param{qdd_folder}.'find_and_mask_ms_galaxy.pl" '.$min1.' '.$min2.' '.$min3.' '.$min4.' '.$min5.' '.$min6.' 1 0 "'.$ms_input.'" "'.$param{file_msfas}.'" "'.$param{file_ms}.'" 1 1';
#  USAGE: $min1 $min2 $min3 $min4	$min5	$min6	$del_prev_mask $keep_non_ms $input_file  $output_fasta $output_ms  $write_seq_file $write_fas_file;
if ($param{debug})
{print $fh_log "\n$find_ms\n\n";}
print "Selecting sequences with microsatellites\n";
system $find_ms;

if($param{qual})
{
	delete_qual_if_not_in_fasta($trimmed_qual, $param{file_msfas}, $param{file_qual_trim});
}

my $seqn_ms = get_seq_number_from_fasta($param{file_msfas});
print $fh_log "$seqn_ms sequences have microsatellites\n";

}# end else contig


############################################## end prg


print $fh_log "\n";
my $time_end = time;
my $time_run = $time_end - $folder_id;
print $fh_log "\npipe1.pl started at : $date\n";
print $fh_log "The run took $time_run s.\n";


if ($param{adapter} and $param{contig} == 0)
{

  delete_tag_adapt_files($param{adapter_file}, $param{syst});
}


if ($param{del_files})
{
  delete_folder($tempfolder, $param{syst});
}

print "\nSee log_file for summary:\n$param{file_log_pipe1}\n\n";
print $fh_log "The analysis is finished\n";
close $fh_log;
exit;


#############################################################

sub delete_qual_if_not_in_fasta
{
 my ($in_qual, $fas, $out_qual) = @_;
 # input quality file, $fasta file with selected sequences, $output quality file with value onlu for selected sequnces
 my %seq = read_hach_from_fasta_code_wo_space($fas);
 my %qual = read_hach_from_qual_code_wo_space($in_qual);
 
 open(OUT, ">$out_qual") or die "Cannnot open $out_qual\n";
 foreach my $code(sort keys %seq)
 {
	print OUT ">$code\n$qual{$code}\n";
 }
 close OUT;
}

#############################################################

sub trim_extremities
{
	my  ($seq, $qual, $window, $mean_qs, $length_limit, $tbl_trim) = @_;
	
	open(TBL, ">$tbl_trim");
	foreach my $code (keys %$seq)
	{
		my $sl = length $$seq{$code};
		my ($seq_cut, $qual_cut) = trim_5_prime($$seq{$code}, $$qual{$code}, $code, $window, $mean_qs);
		my $beg_cut = $sl - length $seq_cut;
		($seq_cut, $qual_cut)= trim_3_prime($seq_cut, $qual_cut, $code, $window, $mean_qs);
		my $end_cut = $sl - $beg_cut - length $seq_cut;
		if(length $seq_cut < $length_limit)
		{
			delete$$seq{$code};
			delete$$qual{$code};
		}
		else
		{
			$$seq{$code} = $seq_cut;
			$$qual{$code} = $qual_cut;
			print TBL "$code\t$beg_cut\t$end_cut\n";
		}
	}
	close TBL;
}


##################################################################

sub define_param_limits_pipe1
{
my ($param) = @_;
#hash of array [type(list, real, folder), error message, ['list of values/lower_upper limits)']
  my %param_limits = (
      'galaxy'=> ['list', "galaxy must be 0 or 1\n", [0, 1]], 
	  'fastq'=> ['list', "fastq must be 0 or 1\n", [0, 1]], 
      'syst'=> ['list', "Operating system must be win or linux\n", ['win', 'linux']], 
      'blast_path'=> ['folder', "Cannot open $$param{blast_path} folder\n"], 
      'adapter'=> ['list', "adapter must be 0 or 1\n", [0,1]],
      'length_limit'=> ['real', "Minimum sequence length must be between 1 and 1000\n", [1,1000]],
      'del_files'=> ['list', "del_files must be 0 or 1\n", [0,1]],
      'input_file' => ['file', "Cannot open input_file ($$param{input_file})\n"],
      'qdd_folder' => ['folder', "Cannot open qdd_folder ($$param{qdd_folder})\n"],
      'out_folder' => ['folder', "Cannot open out_folder ($$param{out_folder})\n"],
      'file_wov' => ['string', "file_wov ($$param{file_wov}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_length' => ['string', "file_length ($$param{file_length}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_msfas' => ['string', "file_msfas ($$param{file_msfas}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_ms' => ['string', "file_ms ($$param{file_ms}) can contain only alphanumerical values, underscore '_' or dot\n"],
     'outfile_string' => ['string', "outfile_string ($$param{outfile_string}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'debug' => ['list', "debug must be 0 or 1\n", [0,1]],
	  'num_threads' => ['real', "The number of threads ($$param{num_threads}) must be between 1 and 100000\n", [0,100000]],
	  'qual' => ['list', "qual must be 0 or 1\n", [0,1]],
	'contig' => ['list', "contig must be 0 or 1\n", [0,1]]
      );
	  
	  if($param{qual})
	  {
		$param_limits{qual_file} = ['file', "Cannot open qual_file ($$param{qual_file})\n"];
		$param_limits{mean_qs} = ['real', "mean_qs must be between 0 and 40", [0,40]];
		$param_limits{file_qual_trim} = ['string', "file_qual_trim ($$param{file_qual_trim}) can contain only alphanumerical values, underscore '_' or dot\n"];
	  }
	  if($param{adapter})
	  {
	    $param_limits{adapter_file} = ['file', "Cannot open adapter_file ($$param{adapter_file})\n"];
	  }
	  if($param{contig})
	  {
		$param_limits{flank_length} = ['real', "The length of the flanking region ($$param{flank_length}) must be between 10 and 300\n", [0,300]];
	  }
return %param_limits;
}
##################################################################


sub copy_tag_and_adapter_file
{
my ($data1_list, $source_folder, $file, $outfolder, $bool, $syst) = @_;

my $file_new = $outfolder.'_'.$file;

if ($bool)
{
  if (exists $$data1_list{$file})
  {
    copy_one_file($source_folder, $file, $file_new, $syst);
    delete $$data1_list{$file};
  }
  else
  {
    print "The $file file is not in in the $source_folder folder\n";
    exit;
  }
}
elsif (exists $$data1_list{$file})
{
  delete $$data1_list{$file};
}

return $file_new;
}


###################################################

sub get_sorted_files
{
 my ($adapt_inp, $folder, $motif) = @_;
 my %temp = get_file_folder_list($folder);
 foreach my $file (keys %temp)
 {
  if ($file =~ /$motif/)
  {
    $$adapt_inp{$file} = $folder;
  }
 }
 %temp = ();
}

###################################################

sub delete_notag
{
my ($inp) = @_;

foreach my $file (keys %$inp)
{
  if ($file =~ /_NOTAG\./)
  {
    delete $$inp{$file};
  }
}
}

###################################################

sub get_seq_number
{
 my ($motif_list,  $seq_number, $folder) = @_;
 
unless ( opendir(FOLDER, "$folder") )
{
      print "Cannot access to folder $folder\n";
      exit;
}
my @files = grep ( !/^\.\.?$/, readdir(FOLDER) );
closedir(FOLDER);

foreach my $file (@files)
{
  foreach my $motif (keys %$motif_list)
  {
    if ($file =~ /$motif/)
    {
      open_file_in($folder.'/'.$file, 'IN');
      my $i = 0;
      while (my $line = <IN>)
      {
        if ($line =~ />/)
        {
          ++$i;
        }
      }
      $$seq_number{$motif}{$file} = $i;
      close IN;
    }
  }
}

 
}

###############################################################
sub initialise_motif_list
{
  my ($file) = @_;
 my %motif_list = ();

  $motif_list{$file} = "The number of sequences in $file";
  $motif_list{'\.woa'}= "The number of sequences without adapter at the begining (not used for microsat research)";
  $motif_list{'\.wov'} = "The number of sequences with adapter at the begining";
  $motif_list{'NOTAG\.tag'} = "The number of sequences without TAG (not used for microsat research)";
  $motif_list{'bp_mask\.fas'} = "The number of sequences with at least one microsatellite";
  $motif_list{'bp.fas'} = "The number of sequences of at least $param{length_limit} bp (without tag and adapter)";

  return %motif_list;
}


#################################################################

sub get_seq_length_remove_short
{
 my ($input, $length_file, $limit_l, $fasta_out) = @_;

open(IN, $input) or die "Cannot open $input file\n";
open(L, ">$length_file") or die "Cannot open $length_file file\n";
open(FAS, ">$fasta_out") or die "Cannot open $fasta_out file\n";

print L "code\tlength\n";

my $code = '';
my $seq = '';
my $l = 0;
while(my $line = <IN>)
{
	chomp $line;
	if ($line=~  />\s*([^\s]+)/)
	{
		unless($code eq '')
		{
			$seq =~  s/\s//g;
			my $l = length $seq;
			print L "$code\t$l\n";
			if ($l >= $limit_l)
			{
				print FAS ">$code\n";
				print FAS cut_up_fasta_line($seq, '100');
			}
			$seq = '';
		}
		$code = $1;

	}
	else
	{
		$seq .= $line;
	}
}
$seq =~  s/\s//g;
$l = length $seq;
print L "$code\t$l\n";
if ($l > $limit_l)
{
 print FAS ">$code\n";
 print FAS cut_up_fasta_line($seq, '100');
}

close IN;
close OUT;
close FAS;
}

