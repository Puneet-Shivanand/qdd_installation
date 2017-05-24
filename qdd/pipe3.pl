#!/usr/bin/perl -w


  use warnings;
  use strict;
  use subprogramQDD;
  use Data::Dumper;
  
  
#INPUT: 
#xxx_pipe2_for_pipe3.fas file with clean (no adapter, no tag) microsatellite containing non-redundant sequences (output of pipe2.pl)

# Run Priper3 in a iterative way: For each input sequences, different (more or less stringent) primer design scenarios are run, and for each of them the desired length of the PCR product is increased at each step. See documentation for more details.

# output: 
#xxx_pipe3_primers.tabular => tab separated values, info on primers, target region parameters, primer design etc.
#xxx_pipe3_targets.fas => sequences, where primer design was successful.



my $mask_mono = 5;
my $mask_di = 3;
my $mask_tri = 3;
my $mask_tetra = 3;
my $mask_penta = 3;
my $mask_hexa = 3;
my $set_file_default = '/etc/qdd/set_qdd_default.ini';
my %scoring_matrix = (
	'match' => 1,
	'mismatch' => -1,
	'N' => -0.25
);
$set_file_default = check_set_file($set_file_default);

print_heading_3_0('STDOUT');
print "SCRIPT pipe3.pl\n";

my %param = (
  galaxy => '', # one if run on galaxy, 0 for command line
  syst => '',  
  primer3_path => '',
  primer3_version => '',
  input_file => '',
  out_folder => '',
  qdd_folder => '',
  file_pipe3_inp => '_pipe2_for_pipe3.fas',
  file_log_pipe3 => '_pipe3_log.txt',  
  file_primers => '_pipe3_primers.tabular',
  file_targets => '_pipe3_targets.fas',
  outfile_string => '',
  pcr_min => '',
  pcr_max => '',
  pcr_step => '',
  PRIMER_GC_CLAMP => '',
  PRIMER_OPT_SIZE => '',
  PRIMER_MIN_SIZE => '',
  PRIMER_MAX_SIZE => '',
  PRIMER_OPT_TM => '',
  PRIMER_MIN_TM => '',
  PRIMER_MAX_TM => '',
  PRIMER_MAX_DIFF_TM => '',
  PRIMER_MIN_GC => '',
  PRIMER_OPT_GC_PERCENT => '',
  PRIMER_MAX_GC => '',
  PRIMER_SELF_ANY => '',
  PRIMER_SELF_END => '',
  PRIMER_MAX_POLY_X => '',
  PRIMER_NUM_RETURN => '',
  del_files => '',
  debug => 0,
  qual_file => '',
  min_qs => 20,
  qual => 0,
  contig => 1
);

read_set_file_galaxy(\%param, $param{qdd_folder}.$set_file_default);
modify_params_from_tags(\%param, \@ARGV);
my %param_limits = define_param_limits_pipe3(\%param); #hash of array $param_limits{syst}[type(list, real, folder), error message, ['list of values/lower_upper limits)']
add_slash_to_folder(\%param, \%param_limits);

my $galaxy = $param{galaxy}; # one if run on galaxy, 0 for command line
my $folder_id = time;
my $date = localtime $folder_id;
my $tempfolder = $param{out_folder}.'pipe3_'.$folder_id.'/';
make_folder($param{syst}, $tempfolder);

if($galaxy)
{
	
	if ($param{contig} eq 'true' or $param{contig} eq 'TRUE')
	{ 	$param{contig} = 1;	}
	else
	{	 $param{contig} = 0;	}
	
	if ($param{qual_file} eq '' or $param{qual_file} eq 'None')
	{ 	$param{qual} = 0;	}
	else
	{	 $param{qual} = 1;	}

	$param{outfile_string} = 'NA';
}
elsif ($param{outfile_string} eq '')
{
	$param{outfile_string} = get_filename_root($param{input_file});
	$param{outfile_string} =~ s/_pipe2_for_pipe3.*//;
}



unless($galaxy)
{
check_param_limits(\%param, \%param_limits);
	if ($param{outfile_string} eq '')
	{
		$param{outfile_string} = get_filename_root($param{input_file});
		$param{outfile_string} =~ s/_pipe2_for_pipe3.*//;
	}


	my @clefs = ('file_primers', 'file_targets');
	if($param{file_log_pipe3} =~ /^_/) # if the name of the log file has not been specified (run in pipeline)
	{
		push(@clefs, 'file_log_pipe3');
	}
	get_last_version_and_modify_1($param{out_folder}, \@clefs, $param{outfile_string}, \%param);

	if($param{file_log_pipe3} =~ /^_/) # if the name of the log file has not been specified (run in pipeline)
	{
		$param{file_log_pipe3} = $param{out_folder}.$param{outfile_string}.$param{file_log_pipe3};
	}

	$param{file_primers} = $param{out_folder}.$param{outfile_string}.$param{file_primers};
	$param{file_targets} = $param{out_folder}.$param{outfile_string}.$param{file_targets};
}

open(my $fh_log, '>>', $param{file_log_pipe3}) or die "Cannot open $param{file_log_pipe3} $!\n";
print_heading_3_0($fh_log);
print $fh_log "SCRIPT pipe3.pl\n";

if($param{debug})
{
	print $fh_log Dumper(\%param);
	print $fh_log "\n";
}

print $fh_log "INPUT file : $param{input_file}\n";
print "Input file : $param{input_file}\n\n";
if ($param{qual})
{
print $fh_log "Input quality file : $param{qual_file}\n\n";
}

print $fh_log "\nOUTPUT files:\n";
print $fh_log "String for naming output files: $param{outfile_string}\n";
print $fh_log "\t$param{file_primers}: Information on primers, target region parameters, primer design etc.\n";
print $fh_log "\t$param{file_targets}: Sequences, where primer design was successful\n";


print $fh_log "\nGENERAL PARAMETERS\n";
print $fh_log "System (win/linux): $param{syst}\n";
print $fh_log "Output folder: $param{out_folder}\n";
print $fh_log "Path to QDD executables: $param{qdd_folder}\n";
print $fh_log "Path to Primer3 executables: $param{primer3_path}\n";
print $fh_log "Primer3 version [1/2]: $param{primer3_version}\n";
print $fh_log "Delete intremediate files: $param{del_files}\n";
print $fh_log "Sequences are extracted from contigs: $param{contig}\n";	


#print $fh_log "Use quality file for masking low quality regions before primer design: $param{qual}\n";
if ($param{qual})
{
print $fh_log "Primers cannot include a base if its quality score is lower than $param{min_qs} (-min_qs)\n\n";
}
print $fh_log "\nPRIMER DESIGN PARAMETERS:\n";
print $fh_log "Minimum PCR product size: $param{pcr_min}\n";
print $fh_log "Maximum PCR product size: $param{pcr_max}\n";
print $fh_log "PCR size step: $param{pcr_step}\n";
print $fh_log "PRIMER_GC_CLAMP: $param{PRIMER_GC_CLAMP}\n";
print $fh_log "PRIMER_OPT_SIZE: $param{PRIMER_OPT_SIZE}\n";
print $fh_log "PRIMER_MIN_SIZE: $param{PRIMER_MIN_SIZE}\n";
print $fh_log "PRIMER_MAX_SIZE: $param{PRIMER_MAX_SIZE}\n";
print $fh_log "PRIMER_OPT_TM: $param{PRIMER_OPT_TM}\n";
print $fh_log "PRIMER_MIN_TM: $param{PRIMER_MIN_TM}\n";
print $fh_log "PRIMER_MAX_TM: $param{PRIMER_MAX_TM}\n";
print $fh_log "PRIMER_MAX_DIFF_TM: $param{PRIMER_MAX_DIFF_TM}\n";
print $fh_log "PRIMER_MIN_GC: $param{PRIMER_MIN_GC}\n";
print $fh_log "PRIMER_OPT_GC_PERCENT: $param{PRIMER_OPT_GC_PERCENT}\n";
print $fh_log "PRIMER_MAX_GC: $param{PRIMER_MAX_GC}\n";
print $fh_log "PRIMER_SELF_ANY: $param{PRIMER_SELF_ANY}\n";
print $fh_log "PRIMER_SELF_END: $param{PRIMER_SELF_END}\n";
print $fh_log "PRIMER_MAX_POLY_X: $param{PRIMER_MAX_POLY_X}\n";
print $fh_log "PRIMER_NUM_RETURN: $param{PRIMER_NUM_RETURN}\n";
print $fh_log "**********************************************************\n";


my $trans_file = $param{qdd_folder}.'motifs_transmot_hach1_6_corr.txt';
my %trans = open_file_and_hach($trans_file, '1');
  # make lc motifs in %trans
foreach my $mot (keys %trans)
{
  my $mot_lc = lc $mot;
  $trans{$mot_lc} = $trans{$mot};
}

 

my $primer3_inp_file = $tempfolder.'primer3_input.txt';
my $primer3_out_file = $tempfolder.'primer3_output.txt';

my %seq_count = ();


# %polym is a hash of tables with cons seq code as keys and a table of motif_fpos_lastpos as value
my %polym = make_polymorph_info($param{input_file});
my $seq_file = '';

# mask first microsatellites, keep mask, then mask all micro- nanosatellites and homopomlymers at the same time : TTTTTTTGTGTGTGTG: T6(TG)5

# mask only target MS in lc. In the next find MS(1) target MS are already masked and will not be taken as part of a shorter MS just before 
my $temp_ms_find0_fas = $tempfolder.'ms_find0_fas.fas';
my $temp_ms_find0_ms = $tempfolder.'ms_find0_ms.txt';
my $find_ms0 = 'perl "'.$param{qdd_folder}.'find_and_mask_ms_galaxy.pl" 1000000 5 5 5 5 5 1 1 "'.$param{input_file}.'" "'.$temp_ms_find0_fas.'" "'.$temp_ms_find0_ms.'" 1 1';
#  USAGE: $min1 $min2 $min3 $min4 $min5	$min6	$del_prev_mask $keep_non_ms $input_file  $output_fasta $output_ms  $write_seq_file $write_fas_file;
if($param{debug})
{
	print $fh_log "\n$find_ms0\n";
}

print "Soft masking sequences\n";
system $find_ms0;

# mask MS and homopolymers in the input fasta file use low limit values to avoid false BLAST hits  and eliminate CSS sequences (no autohit)
my $temp_ms_find1_fas = $tempfolder.'ms_find1_fas.fas';
my $temp_ms_find1_ms = $tempfolder.'ms_find1_ms.txt';
my $find_ms1 = 'perl "'.$param{qdd_folder}.'find_and_mask_ms_galaxy.pl" '.$mask_mono.' '.$mask_di.' '.$mask_tri.' '.$mask_tetra.' '.$mask_penta.' '.$mask_hexa.' 0 1 "'.$temp_ms_find0_fas.'" "'.$temp_ms_find1_fas.'" "'.$temp_ms_find1_ms.'" 1 1';
#  USAGE: $min1 $min2 $min3 $min4 $min5	$min6	$del_prev_mask $keep_non_ms $input_file  $output_fasta $output_ms  $write_seq_file $write_fas_file;
if($param{debug})
{
	print $fh_log "\n$find_ms1\n";
}
system $find_ms1;
 
#  $file = $outfolder_spec.'/'.$file_root.'_mask_mask.fas';
#  $seq_file = $outfolder_spec.'/'.$file_root.'_mask.seq';

  
# %ms is a hash of table of tables; hash keys are seq codes $ms{$code} that contain a table of pooled microsatellites$ms{$code}[ms]
# the first element of the ms table is the list of the first positions of target microsatellites
# Then for each microsatellite either [first_pos, last_pos, c/p, repeat_number(of the longest pure sterch), $motif(of the longest pure sterch), $repeat_length(of last strech; deleted after finishing the %ms)]
#Homopolymers are never pooled to any other MS.
# Ms are pooled if the distance between them is not greater than the longest motif length of the two neigbouring MS
my %ms =  make_ms_hash($temp_ms_find1_ms);
 
my %seq = read_hach_from_fasta_code_wo_space($temp_ms_find1_fas);
$seq_count{'Total number of input sequences:'} = scalar keys%seq;
delete_no_target_seq(\%seq, \%ms); # delete sequence from %seq that do not have taget MS (rare cases like tttttgtgtgtgtG)
$seq_count{'Total number of sequences with target MS:'} = scalar keys%seq;

my %quality_excluded = ();
if($param{qual})
{
	my %qual = read_hach_from_qual_code_wo_space($param{qual_file});
	delete_no_target_seq(\%qual, \%ms);

# for each sequence returns a list of lists [first_position, length] of regiosn with quality scores lower than $param{min_qs}
# $quality_exluded{sequence code} = [[first_pos, length], [first_pos, length]]
  %quality_excluded = select_low_quality_bases(\%qual, \%seq, $param{min_qs}); 
  # print info on low quality bases
  my $qual_excl_file = $tempfolder.'quality_excluded.txt';
  unless($param{del_files})
  {
	open(QE, ">$qual_excl_file");
	print QE Dumper(\%quality_excluded);
	close QE;
  }
}

  # primer design conditions:
  #   homo  nano_primer nano_flank  other_target_in flank target_comp
  # A no    no          no          no                    no  
  # B no    no          yes         no                    yes  
  # C no    yes         yes         no                    yes 
  # D no    no          yes         yes                   yes  
  # E no    yes         yes         yes                   yes 
  # F yes   yes         yes         yes                   yes 

  
my %primer3_base = ();
  # $$primer3_base{$code.'_'.$target_ms}{'seq'} = uc$$seq{$code};
  # $$primer3_base{$code.'_'.$target_ms}{'target'} = $target; $first_pos,length
  # $$primer3_base{$code.'_'.$target_ms}{'excluded'} = $excluded; reference to an anonymous table of $first_pos,length
my %results = (); #$results{$code_seq}{left_pos,left_l,rigth_pos,rigth_l}{left_seq/right_seq ...} => value 
my %target_ms_info = (); # $target_ms{seq_code_first_pos} = [first_pos, length_bp, length_in_rn,number_of_ms (1.5 for one compound, 1 for 1 pure),, mot_trans, ms_seq ];

my @condition_of_design = ('A','B','D','F','G');
#my @condition_of_design = ('A');
foreach my $condition (@condition_of_design)
{
    print "Identifying target and excluded regions with $condition design\n";
    select_target_excluded($condition, \%seq, \%ms, \%primer3_base, \%target_ms_info, \%trans, $param{pcr_max}); # %primer3_base is emptyed at the begining of the subroutine
 
    for(my $min= $param{pcr_min}; $min<$param{pcr_max}; $min = $min+$param{pcr_step}) # run primer3 for each PCR product length step
    {
		my $max = $min + $param{pcr_step}-1;
		if ($max > $param{pcr_max})
		{
			$max = $param{pcr_max};
		}
		print "Running Primer3 with $min - $max PCR product interval\n";
		if ($param{primer3_version} == 2)
		{
			write_primer3_input2(\%quality_excluded, $primer3_inp_file, \%primer3_base, $min, $max, $param{PRIMER_GC_CLAMP}, $param{PRIMER_OPT_SIZE}, $param{PRIMER_MIN_SIZE}, $param{PRIMER_MAX_SIZE}, $param{PRIMER_OPT_TM}, $param{PRIMER_MIN_TM}, $param{PRIMER_MAX_TM}, $param{PRIMER_MAX_DIFF_TM}, $param{PRIMER_MIN_GC}, $param{PRIMER_OPT_GC_PERCENT}, $param{PRIMER_MAX_GC}, $param{PRIMER_SELF_ANY}, $param{PRIMER_SELF_END}, $param{PRIMER_MAX_POLY_X}, $param{PRIMER_NUM_RETURN});	  
		}
		else
		{
			write_primer3_input(\%quality_excluded, $primer3_inp_file, \%primer3_base, $min, $max, $param{PRIMER_GC_CLAMP}, $param{PRIMER_OPT_SIZE}, $param{PRIMER_MIN_SIZE}, $param{PRIMER_MAX_SIZE}, $param{PRIMER_OPT_TM}, $param{PRIMER_MIN_TM}, $param{PRIMER_MAX_TM}, $param{PRIMER_MAX_DIFF_TM}, $param{PRIMER_MIN_GC}, $param{PRIMER_OPT_GC_PERCENT}, $param{PRIMER_MAX_GC}, $param{PRIMER_SELF_ANY}, $param{PRIMER_SELF_END}, $param{PRIMER_MAX_POLY_X}, $param{PRIMER_NUM_RETURN});
		}

		unless (-z $primer3_inp_file) # if file is not empty
		{
			run_primer3($param{syst}, $primer3_inp_file, $primer3_out_file, $param{primer3_path});
			if ($param{primer3_version} == 2)
			{
				read_results2($primer3_out_file, \%results, $condition, \%target_ms_info);
			}
			else
			{
				read_results($primer3_out_file, \%results, $condition, \%target_ms_info);
			}
		}
    }
  
}
  
$seq_count{'Total number of sequences with primers:'} = scalar keys%results;
make_fasta_file_with_designed_seq(\%results, \%seq, $param{file_targets});
print "\nWriting primer table\n";
write_primer_table(\%results, \%target_ms_info, $param{file_primers}, \%polym, \%seq, \%scoring_matrix);

if($param{contig})
{
# complete primer file with two columns cotaining the code of the contig and the firts position of the fragment on the contig to facilitate de selection of markers on the same contig but far away to reduce linkage
	get_contig_code_and_pos_on_contig($param{file_primers});
}


print $fh_log "\n\n";
foreach my $code (sort keys %seq_count)
{
  print $fh_log "$code $seq_count{$code}\n";
}
my $time_end = time;
my $time_run = $time_end - $folder_id;

print $fh_log "\npipe3.pl started at : $date\n";
print $fh_log "The analyses took $time_run seconds\n";
  


if ($param{del_files}==1)
{
  delete_folder($tempfolder,$param{syst});
}

print "\nSee log_file for summary:\n$param{file_log_pipe3}\n\n";
print $fh_log "The analysis is finished\n";
close $fh_log;
exit;

##########################################################


##########################################################


sub select_low_quality_bases
{
my ($qual, $seq, $min_qs) = @_;
my %quality_exluded = ();

foreach my $code (keys %$seq)
{
if(exists $$qual{$code})
{
	my @q = split(' ', $$qual{$code});
	unless(scalar@q == length $$seq{$code})
	{
		print "The number of quality scores does not equal to the number of bases in $code\n";
		exit;
	}
	my $excl = '';
	for(my $i=0; $i <scalar@q; ++$i)
	{
		if($q[$i] < $min_qs)
		{
			$excl .= 'L';
		}
		else
		{
			$excl .= 'H';
		}
	}
	while($excl =~ /L+/g)
	{
		my $end = pos$excl;
		my $l = length $&;
		my $beg = $end - $l +1;
		push(@{$quality_exluded{$code}}, [$beg, $l]);
	}
}
}
#print Dumper(\%quality_exluded);
return %quality_exluded;
}


##########################################################

sub make_ms_hash
{
my ($seq_file) = @_;
my %ms = ();
# %ms is a hash of table of tables; hash keys are seq codes $ms{$code} that contains a table of pooled microsatellites$ms{$code}[ms]
# the first element of the ms table is the list of the first positions of target microsatellites
# Then for each microsatellite either [first_pos, last_pos, c/p, repeat_number(of the longest pure sterch), $motif(of the longest pure sterch), $repeat_length(of last strech; deleted after finishing the %ms)]
#Homopolymers are never pooled to any other MS.
# Ms are pooled if the distance between them is not greater than the longest motif length of the two neigbouring MS

open_file_in($seq_file, 'SEQ');
my $title = <SEQ>;
while(my $line = <SEQ>)
{
  chomp $line;
  $line =~ s/;$//;
  my @line = split(';', $line);
  # initialize first element of the table (this will be the list of the first positions of target MS)
    $ms{$line[0]} = [[0]];
  for (my $i= 3; $i<scalar@line; $i =$i+4)
  {
    if (scalar@{${$ms{$line[0]}}[-1]} == 6 and length$line[$i]>1) # if there is a MS before and neither this nor the actual Ms is homopolymer 
    {
      my $dist_ms = $line[$i+1] - $ms{$line[0]}[-1][1] - 1;
      if ( ($dist_ms > $ms{$line[0]}[-1][5]) and ($dist_ms > length$line[$i])) # if the actual and previous MS have greater than the (longest motif length) bases between
      {
#        print "$line[0], $line[$i], $line[$i+1], $line[$i+2], $line[$i+3], $dist_ms, @{${$ms{$line[0]}}[-1]}\n";
        new_ms($ms{$line[0]}, $line[$i], $line[$i+1], $line[$i+2], $line[$i+3]);
      }
      else #distance is small between neigbouring non-homopolymer MS
      {
        pool_ms($ms{$line[0]}, $line[$i], $line[$i+1], $line[$i+2], $line[$i+3]);
      }
    }
    else # first MS in seq, or previois or actual ms is homoploymer
    {
      new_ms($ms{$line[0]}, $line[$i], $line[$i+1], $line[$i+2], $line[$i+3]);
    }
  } # end for
  # delete the motif lenght of the last strech
  for(my $j= 1; $j< scalar @{$ms{$line[0]}}; ++$j)
  {
    if (scalar @{$ms{$line[0]}[$j]} == 6)
    {
      pop@{$ms{$line[0]}[$j]};
    }
  }

  # delete seq if there is no target MS

  shift@{$ms{$line[0]}[0]};
  unless (scalar @{$ms{$line[0]}[0]})
  {
    delete $ms{$line[0]};
#    print "$line[0]\n";
  }

}
close SEQ;


return %ms; 
}
#####################################################

sub pool_ms
{
my ($code_table_ref, $motif, $beg, $end, $rep_n) = @_;

$$code_table_ref[-1][1] = $end;
$$code_table_ref[-1][2] = 'c';
if ($rep_n > $$code_table_ref[-1][3])
{
  $$code_table_ref[-1][3] = $rep_n;
  $$code_table_ref[-1][4] = $motif;
}
  $$code_table_ref[-1][5] = length$motif;

if ($rep_n > 4)
{
  unless( $$code_table_ref[0][-1] == $$code_table_ref[-1][0])
  {
    push(@{$$code_table_ref[0]}, $$code_table_ref[-1][0]);    
  }
}

}
#####################################################

sub new_ms
{
my ($code_table_ref, $motif, $beg, $end, $rep_n) = @_;

my $lm = length $motif;

if ($lm ==1)
{
  push(@{$code_table_ref}, [$beg, $end]);
} 
else
{
  push(@{$code_table_ref}, [$beg, $end, 'p', $rep_n, $motif, $lm]);
  if($rep_n > 4)
  {
    push(@{$$code_table_ref[0]}, $beg);
  }
}

}
#########################################################

sub delete_no_target_seq
{
my ($seq, $ms) = @_;
foreach my $code (keys %$seq)
{
  unless(exists $$ms{$code})
  {
    delete $$seq{$code};
  }
}
}

#########################################################

sub select_target_excluded
{
my ($condition, $seq, $ms, $primer3_base, $target_ms_info, $trans, $pcr_product_lmax) = @_;
%$primer3_base = ();
# $$primer3_base{$code.'_'.$target_ms}{'seq'} = uc$$seq{$code};
# $$primer3_base{$code.'_'.$target_ms}{'target'} = $target; $first_pos,length
# $$primer3_base{$code.'_'.$target_ms}{'excluded'} = $excluded; reference to an anonymous table of $first_pos,length
# %target_ms_info = (); # $target_ms{seqcode_first_pos} = [$first_pos, $length_bp, $length_rn, number_of_ms (1.5 for one compound, 1 for 1 pure),, $seqfr, mot_trans];


# %ms is a hash of table of tables; $ms{$code}[ms][first_pos, last_pos, c/p, repeat_number, $motif] or [first_pos, last_pos]
# hash keys are seq codes $ms{$code} that contains a table of pooled microsatellites$ms{$code}[ms]
# the first element of the ms table is the list of the first positions of target microsatellites
# Then for each microsatellite either [first_pos, last_pos, c/p, repeat_number(of the longest pure sterch), $motif(of the longest pure sterch), $repeat_length(of last strech; deleted after finishing the %ms)]
#Homopolymers are never pooled to any other MS.
# Ms are pooled if the distance between them is not greater than the longest motif length of the two neigbouring MS

  # primer design conditions:
  #   homo  nano_primer nano_flank  other_target_in flank target_comp
  # A no    no          no          no                    no  
  # B no    no          yes         no                    yes  
  # C no    yes         yes         no                    yes 
  # D no    no          yes         yes                   yes  
  # E no    yes         yes         yes                   yes 
  # F yes   yes         yes         yes                   yes 
  
foreach my $code (sort keys %$ms) # for all sequences
{
  my $target_l = '';
  my $target = '';
  my $excluded = [];
  if ($condition eq 'A' or $condition eq 'B' or $condition eq 'C' or $condition eq 'F') # only one target MS
  {
  for ( my $t = 0; $t < scalar@{$$ms{$code}[0]}; ++$t) # for all target ms
  {
    my $target_ms = $$ms{$code}[0][$t]; # first position of the target MS
    for (my $i = 1; $i< scalar@{$$ms{$code}}; ++$i) # identify target MS info in table
    {
      if ($$ms{$code}[$i][0] == $target_ms) # identify target MS info in table
      {      
        $target_l = $$ms{$code}[$i][1] - $$ms{$code}[$i][0] +1; # length of the tagter region
	if($target_l < $pcr_product_lmax) # take only target regions that are shorter than the max length of the PCR product
	{
        	$target = $$ms{$code}[$i][0].','.$target_l; # terget region ident: firtspos,length
        	$excluded = [];
        	if ($condition eq 'A')
        	{
         	 if (my $bool = design_a($code, $i, $ms, $seq, $excluded))
         	 {
           	 fill_primer3_base($primer3_base, $target_ms_info, $trans, $ms, $seq, $excluded, $target_ms, $target_l, $code, $i, $i);
         	 }
        	}
       		elsif ($condition eq 'B' or $condition eq 'C') 
        	{
          	design_bc($code, $i, $ms, $seq, $excluded, $condition);
          	fill_primer3_base($primer3_base, $target_ms_info, $trans, $ms, $seq, $excluded, $target_ms, $target_l, $code, $i, $i);
        	}
        	elsif ($condition eq 'F')
        	{
          	design_fg($code, $ms, $seq, $excluded, $i, $i);
          	fill_primer3_base($primer3_base, $target_ms_info, $trans, $ms, $seq, $excluded, $target_ms, $target_l, $code, $i, $i);
        	}
        	$i = scalar@{$$ms{$code}};
	}#end if($target_l < $pcr_product_lmax)
      } # end if  ($$ms{$code}[$i][0] == $target_ms)
    } # end for $i
  }# end for $t
  } #end if $condition ABC 
  elsif (scalar@{$$ms{$code}[0]} > 1)  # codititons DEG and more than 1 target allowed in target region
  {
    my $f_target_index = '';
    my $l_target_index = '';
    for (my $t1 = 0; $t1 <scalar@{$$ms{$code}[0]}-1; ++$t1) # determine all possible pairs of sinlge target to make a coumpound target region
    {
      for (my $t2 = $t1+1; $t2 <scalar@{$$ms{$code}[0]}; ++$t2) # determine all possible pairs of sinlge target to make a coumpound target region
      {
        for (my $i = 1; $i< scalar@{$$ms{$code}}; ++$i) # identify target MS info in table
        {
          if ($$ms{$code}[$i][0] == $$ms{$code}[0][$t1]) # identify first target MS info in table
          {
            $f_target_index = $i;
          }
          elsif ($$ms{$code}[$i][0] == $$ms{$code}[0][$t2]) # identify last target MS info in table
          {
            $l_target_index = $i;
          }
        }# end for $i
        if ($condition eq 'D' or $condition eq 'E') # next pair if there is a homopolymer between the target ms and disign is ED
        {
          my $homo_bool = check_homopolymer_between_targets($ms, $code, $f_target_index, $l_target_index);
          if ($homo_bool) # homopolymer between target mss
          {
            $t2 = scalar@{$$ms{$code}};
            ++$t1;
            next;
          }
        } # if ($condition eq 'D' or $condition eq 'E')
        
        $target_l = $$ms{$code}[$l_target_index][1] - $$ms{$code}[$f_target_index][0] +1;
        $target = $$ms{$code}[$f_target_index][0].','.$target_l;
        $excluded = [];
	if($target_l < $pcr_product_lmax) # take only target regions that are shorter than the max length of the PCR product
	{                 
       		if ($condition eq 'D' or $condition eq 'E')
        	{
          	design_de($code, $ms, $seq, $excluded, $f_target_index, $l_target_index, $condition);
          	fill_primer3_base($primer3_base, $target_ms_info, $trans, $ms, $seq, $excluded, $$ms{$code}[$f_target_index][0], $target_l, $code, $f_target_index, $l_target_index);
        	}
        	elsif ($condition eq 'G')
        	{
          	design_fg($code, $ms, $seq, $excluded, $f_target_index, $l_target_index);
          	fill_primer3_base($primer3_base, $target_ms_info, $trans, $ms, $seq, $excluded, $$ms{$code}[$f_target_index][0], $target_l, $code,  $f_target_index, $l_target_index);
       		}
	}# end if ($target_l < $pcr_product_lmax)
        
      } # end for $T2 
    }# end for $1
  }# end elsif (scalar@{$$ms{$code}[0] > 1)  # codititons DEF and more than 1 target
}# end foreach


}

##############################################

sub check_homopolymer_between_targets
{
my ($ms, $code, $f_index, $l_index) = @_;
my $bool = 0;

for (my $i = ($f_index +1); $i < $l_index; ++$i)
{
  if (scalar@{$$ms{$code}[$i]} ==2)
  {
    $bool = 1;
    last;
  }
}

return $bool;
}
################################################
sub fill_primer3_base
{
my ($primer3_base, $target_ms_info, $trans, $ms, $seq, $excluded, $target_ms, $target_l, $code, $f_target_index, $l_target_index) = @_;

#my @tagget_list = split(',', $target_ms);

 $$primer3_base{$code.'_'.$target_ms.','.$target_l}{'seq'} = uc$$seq{$code};
 $$primer3_base{$code.'_'.$target_ms.','.$target_l}{'target'} = $target_ms.','.$target_l; 
 $$primer3_base{$code.'_'.$target_ms.','.$target_l}{'excluded'} = $excluded;
 unless(exists $$target_ms_info{$code.'_'.$target_ms.','.$target_l})
 {
  my $seq_ms = substr(uc$$seq{$code}, $target_ms-1, $target_l);
  my $longest_ms = 0;
  my $longest_ms_index = 0;
  my $ms_number = 0;
  for(my $i = $f_target_index; $i<=$l_target_index; ++$i)
  {
    if (scalar@{$$ms{$code}[$i]} > 2 and $$ms{$code}[$i][3] > 4)
	{
		++$ms_number; # count the number of ms in the target region
		if($$ms{$code}[$i][3] > $longest_ms)
		{
			$longest_ms = $$ms{$code}[$i][3];
			$longest_ms_index = $i;
		}
	}
  }
  if($ms_number == 1)
  {
	if ($$ms{$code}[$longest_ms_index][2] eq 'c')
	{
		$ms_number = 1.5;
	}
  }
  my  $mot_trans = $$trans{$$ms{$code}[$longest_ms_index][4]};
  $$target_ms_info{$code.'_'.$target_ms.','.$target_l} = [$target_ms, $target_l, $$ms{$code}[$longest_ms_index][3], $ms_number, $mot_trans, $seq_ms];
 }
}
##################################################


sub design_a
{
my ($code, $i, $ms, $seq, $excluded) = @_;
my $bool = 0; 
if ($$ms{$code}[$i][2] eq 'p')
{ 
    $bool = 1;
    if ($i > 1)
    {
      push(@$excluded, '1,'.$$ms{$code}[$i-1][1])
    }
    if ($i+1 < scalar @{$$ms{$code}})
    {
      push(@$excluded, $$ms{$code}[$i+1][0].','.((length$$seq{$code}) - $$ms{$code}[$i+1][0] +1))
    }
}
return $bool;
}
###########################################

sub design_bc
{
my ($code, $i, $ms, $seq, $excluded, $condition) = @_;
# $ms{$code}[ms][first_pos, last_pos, c/p, repeat_number, $motif] or [first_pos, last_pos]

foreach (my $j = $i-1; $j > 0; $j = $j-1) 
{
  if (scalar@{$$ms{$code}[$j]} == 2 or $$ms{$code}[$j][3] > 4) # find the the closest nano of target ms before the tagret ms 
  {
      push(@$excluded, '1,'.$$ms{$code}[$j][1]);
      last;
  }
  if ($condition eq 'B')
  {
    push(@$excluded, $$ms{$code}[$j][0].','.($$ms{$code}[$j][1] - $$ms{$code}[$j][0] +1)); # exclude nanosat
  }
}

foreach (my $j = $i+1; $j < scalar @{$$ms{$code}}; ++$j)
{
  if (scalar@{$$ms{$code}[$j]} == 2 or $$ms{$code}[$j][3] > 4)  # find the the closest nano of target ms after the tagret ms 
  {
      push(@$excluded, $$ms{$code}[$j][0].','.((length$$seq{$code}) - $$ms{$code}[$j][0] +1));
      last;
  }
  if ($condition eq 'B')
  {
    push(@$excluded, $$ms{$code}[$j][0].','.($$ms{$code}[$j][1] - $$ms{$code}[$j][0] +1)); # exclude nanosat
  }
}

}

################################################

sub design_de
{
my ($code, $ms, $seq, $excluded, $f_target_index, $l_target_index, $condition) = @_;
# $ms{$code}[ms][first_pos, last_pos, c/p, repeat_number, $motif] or [first_pos, last_pos]

foreach (my $j = $f_target_index-1; $j > 0; $j = $j-1) 
{
  if (scalar@{$$ms{$code}[$j]} == 2) # find the closest homopol of target ms before the tagret ms 
  {
      push(@$excluded, '1,'.$$ms{$code}[$j][1]);
      last;
  }
  elsif ($$ms{$code}[$j][3] > 4)
  {
    push(@$excluded, $$ms{$code}[$j][0].','.($$ms{$code}[$j][1] - $$ms{$code}[$j][0] +1)); # exclude nanosat and other target ms  
  }
  elsif ($condition eq 'D')
  {
    push(@$excluded, $$ms{$code}[$j][0].','.($$ms{$code}[$j][1] - $$ms{$code}[$j][0] +1)); # exclude nanosat and other target ms
  }
}
foreach (my $j = $l_target_index+1; $j < scalar @{$$ms{$code}}; ++$j)
{
  if (scalar@{$$ms{$code}[$j]} == 2)  # find the closest homopol of target ms after the tagret ms 
  {
      push(@$excluded, $$ms{$code}[$j][0].','.((length$$seq{$code}) - $$ms{$code}[$j][0] +1));
      last;
  }
  elsif ($$ms{$code}[$j][3] > 4)
  {
    push(@$excluded, $$ms{$code}[$j][0].','.($$ms{$code}[$j][1] - $$ms{$code}[$j][0] +1)); # exclude nanosat and other target ms  
  }
  elsif ($condition eq 'D')
  {
    push(@$excluded, $$ms{$code}[$j][0].','.($$ms{$code}[$j][1] - $$ms{$code}[$j][0] +1)); # exclude nanosat and other target
  }
}
}

################################################

sub design_fg
{
my ($code, $ms, $seq, $excluded, $f_target_index, $l_target_index) = @_;
# $ms{$code}[ms][first_pos, last_pos, c/p, repeat_number, $motif] or [first_pos, last_pos]

foreach (my $j = $f_target_index-1; $j > 0; $j = $j-1) 
{
  if (scalar@{$$ms{$code}[$j]} > 2 and $$ms{$code}[$j][3] > 4) # find the closest target ms before the tagret region 
  {
      push(@$excluded, '1,'.$$ms{$code}[$j][1]);
      last;
  }
}
foreach (my $j = $l_target_index+1; $j < scalar @{$$ms{$code}}; ++$j)
{
  if (scalar@{$$ms{$code}[$j]} > 2 and $$ms{$code}[$j][3] > 4)  # find the closest target ms after the tagret region 
  {
      push(@$excluded, $$ms{$code}[$j][0].','.((length$$seq{$code}) - $$ms{$code}[$j][0] +1));
      last;
  }
}
}
################################################


sub write_primer3_input
{
my ($quality_excluded, $primer3_inp_file, $primer3_base, $min, $max, $PRIMER_GC_CLAMP, $PRIMER_OPT_SIZE, $PRIMER_MIN_SIZE, $PRIMER_MAX_SIZE, $PRIMER_OPT_TM, $PRIMER_MIN_TM, $PRIMER_MAX_TM, $PRIMER_MAX_DIFF_TM, $PRIMER_MIN_GC, $PRIMER_OPT_GC_PERCENT, $PRIMER_MAX_GC, $PRIMER_SELF_ANY, $PRIMER_SELF_END, $PRIMER_MAX_POLY_X, $PRIMER_NUM_RETURN) = @_;  

#print Dumper(\%$primer3_base);
open_file_out($primer3_inp_file, 'PR');
foreach my $code (sort keys %$primer3_base)
{
	my $seq_code = $code;
	$seq_code =~ s/(_[0-9]+,[0-9]+)+//;
  print PR "PRIMER_SEQUENCE_ID=$code\n";
  print PR "SEQUENCE=$$primer3_base{$code}{'seq'}\n";
  print PR "PRIMER_FIRST_BASE_INDEX=1\n";
  print PR "TARGET=$$primer3_base{$code}{'target'}\n";
  print PR "PRIMER_GC_CLAMP=$PRIMER_GC_CLAMP\n";
  print PR "PRIMER_OPT_SIZE=$PRIMER_OPT_SIZE\n";
  print PR "PRIMER_MIN_SIZE=$PRIMER_MIN_SIZE\n";
  print PR "PRIMER_MAX_SIZE=$PRIMER_MAX_SIZE\n";
  print PR "PRIMER_OPT_TM=$PRIMER_OPT_TM\n";
  print PR "PRIMER_MIN_TM=$PRIMER_MIN_TM\n";
  print PR "PRIMER_MAX_TM=$PRIMER_MAX_TM\n";
  print PR "PRIMER_MAX_DIFF_TM=$PRIMER_MAX_DIFF_TM\n";
  print PR "PRIMER_MIN_GC=$PRIMER_MIN_GC\n";
  print PR "PRIMER_OPT_GC_PERCENT=$PRIMER_OPT_GC_PERCENT\n";
  print PR "PRIMER_MAX_GC=$PRIMER_MAX_GC\n";
  print PR "PRIMER_SELF_ANY=$PRIMER_SELF_ANY\n";
  print PR "PRIMER_SELF_END=$PRIMER_SELF_END\n";
  print PR "PRIMER_MAX_POLY_X=$PRIMER_MAX_POLY_X\n";
  print PR "PRIMER_NUM_RETURN=$PRIMER_NUM_RETURN\n";

  if (@{$$primer3_base{$code}{'excluded'}}  or exists $$quality_excluded{$seq_code})
  {
    print PR "EXCLUDED_REGION=";
    foreach my $e(@{$$primer3_base{$code}{'excluded'}})
    {
      print PR "$e ";
    }
	if (exists$$quality_excluded{$seq_code})
	{
		foreach my $ex (@{$$quality_excluded{$seq_code}})
		{
			print PR "$$ex[0],$$ex[1] ";
		}
	}
    print PR "\n";
  }
  print PR "PRIMER_PRODUCT_SIZE_RANGE=$min-$max\n";
  print PR "=\n";

}
close PR;



}


################################################

sub write_primer3_input2
{
my ($quality_excluded, $primer3_inp_file, $primer3_base, $min, $max, $PRIMER_GC_CLAMP, $PRIMER_OPT_SIZE, $PRIMER_MIN_SIZE, $PRIMER_MAX_SIZE, $PRIMER_OPT_TM, $PRIMER_MIN_TM, $PRIMER_MAX_TM, $PRIMER_MAX_DIFF_TM, $PRIMER_MIN_GC, $PRIMER_OPT_GC_PERCENT, $PRIMER_MAX_GC, $PRIMER_SELF_ANY, $PRIMER_SELF_END, $PRIMER_MAX_POLY_X, $PRIMER_NUM_RETURN) = @_;  

open_file_out($primer3_inp_file, 'PR');
foreach my $code (sort keys %$primer3_base)
{
	my $seq_code = $code;
	$seq_code =~ s/(_[0-9]+,[0-9]+)+//;
  print PR "SEQUENCE_ID=$code\n";
  print PR "SEQUENCE_TEMPLATE=$$primer3_base{$code}{'seq'}\n";
  print PR "PRIMER_FIRST_BASE_INDEX=1\n";
  print PR "SEQUENCE_TARGET=$$primer3_base{$code}{'target'}\n";
  print PR "PRIMER_GC_CLAMP=$PRIMER_GC_CLAMP\n";
  print PR "PRIMER_OPT_SIZE=$PRIMER_OPT_SIZE\n";
  print PR "PRIMER_MIN_SIZE=$PRIMER_MIN_SIZE\n";
  print PR "PRIMER_MAX_SIZE=$PRIMER_MAX_SIZE\n";
  print PR "PRIMER_OPT_TM=$PRIMER_OPT_TM\n";
  print PR "PRIMER_MIN_TM=$PRIMER_MIN_TM\n";
  print PR "PRIMER_MAX_TM=$PRIMER_MAX_TM\n";
  print PR "PRIMER_PAIR_MAX_DIFF_TM=$PRIMER_MAX_DIFF_TM\n";
  print PR "PRIMER_MIN_GC=$PRIMER_MIN_GC\n";
  print PR "PRIMER_OPT_GC_PERCENT=$PRIMER_OPT_GC_PERCENT\n";
  print PR "PRIMER_MAX_GC=$PRIMER_MAX_GC\n";
  print PR "PRIMER_MAX_SELF_ANY=$PRIMER_SELF_ANY\n";
  print PR "PRIMER_MAX_SELF_END=$PRIMER_SELF_END\n";
  print PR "PRIMER_MAX_POLY_X=$PRIMER_MAX_POLY_X\n";
  print PR "PRIMER_NUM_RETURN=$PRIMER_NUM_RETURN\n";
  print PR 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH=',$param{primer3_path},'primer3_config/',"\n";  

  if (@{$$primer3_base{$code}{'excluded'}} or exists($$quality_excluded{$seq_code}))
  {
    print PR "SEQUENCE_EXCLUDED_REGION=";
    foreach my $e(@{$$primer3_base{$code}{'excluded'}})
    {
      print PR "$e ";
    }
	if (exists($$quality_excluded{$seq_code}))
	{
		foreach my $ex (@{$$quality_excluded{$seq_code}})
		{
			print PR "$$ex[0],$$ex[1] ";
		}
	}
    print PR "\n";
  }
  print PR "PRIMER_PRODUCT_SIZE_RANGE=$min-$max\n";
  print PR "=\n";

}
close PR;
}

#######################################################
sub run_primer3
{
my ($syst, $input, $output, $path) = @_;

my $primer3 = '';

if ($syst eq 'win')
{
	$primer3 = '"'.$path.'primer3_core.exe" <'.$input.' > '.$output;
}
elsif ($syst eq 'linux')
{
	$primer3 = $path.'primer3_core <'.$input.' > '.$output;
}
system $primer3;
}
#########################################################

sub read_results2
{
 my ($primer3_out_file, $results, $condition, $target_ms_info) = @_;
 open_file_in($primer3_out_file, 'PRO');
 
my %seq_res = (); # keeps all the results of one seq in a hach with label as a key 
while (my $l = <PRO>)
{
	chomp $l;
	my @l = split('=',$l);
	if (not($l eq '=')) # it is not the end of record (=\n) line
	{
		$seq_res{$l[0]} = $l[1];
	}
	else # end of record => extract primer_pair info for the sequence
	{
		my $primer_pair_num = 0; # the number of primer peirs that has been designed
		if(exists $seq_res{PRIMER_PAIR_NUM_RETURNED})
		{
			$primer_pair_num = $seq_res{PRIMER_PAIR_NUM_RETURNED};
		}
		if ($primer_pair_num > 0) # at least one primer pair has been designed
		{
			$seq_res{SEQUENCE_ID} =~ /(.+?)(_[0-9]+,[0-9]+$)/; # get the sequence code and the target region
			my $target_p = $2;
			my $code = $1;
			$target_p =~ s/_//;
			for(my $i = 0; $i < $primer_pair_num; ++$i)
			{
				my %pr_pair = ();
				  $pr_pair{PRIMER_PAIR_PENALTY} = $seq_res{'PRIMER_PAIR_'.$i.'_PENALTY'};
				  $pr_pair{PRIMER_LEFT_SEQUENCE} = $seq_res{'PRIMER_LEFT_'.$i.'_SEQUENCE'};
				  $pr_pair{PRIMER_RIGHT_SEQUENCE} = $seq_res{'PRIMER_RIGHT_'.$i.'_SEQUENCE'};
				  $pr_pair{PRIMER_LEFT} = $seq_res{'PRIMER_LEFT_'.$i};
				  $pr_pair{PRIMER_RIGHT} = $seq_res{'PRIMER_RIGHT_'.$i};
				  $pr_pair{PRIMER_LEFT_TM} = $seq_res{'PRIMER_LEFT_'.$i.'_TM'};
				  $pr_pair{PRIMER_RIGHT_TM} = $seq_res{'PRIMER_RIGHT_'.$i.'_TM'};
				  $pr_pair{PRIMER_LEFT_END_STABILITY} = $seq_res{'PRIMER_LEFT_'.$i.'_END_STABILITY'};
				  $pr_pair{PRIMER_RIGHT_END_STABILITY} = $seq_res{'PRIMER_RIGHT_'.$i.'_END_STABILITY'};
				  $pr_pair{PRIMER_PRODUCT_SIZE} = $seq_res{'PRIMER_PAIR_'.$i.'_PRODUCT_SIZE'};
				  $pr_pair{TARGET} = $target_p;
				add_info($results, \%pr_pair, $code, $condition, $target_ms_info);
			}
		}
		%seq_res = ();
	}
}
 close PRO;
}
#######################################################

sub read_results
{
 my ($primer3_out_file, $results, $condition, $target_ms_info) = @_;
 open_file_in($primer3_out_file, 'PRO');
 
my @seq_res = ();
while (my $l = <PRO>)
{
  push(@seq_res, $l);
 if ( $l eq "=\n")
 {
  if (scalar@seq_res >25) # at least one primer designed
  {

    my $code = shift@seq_res;
    chomp $code;
    $code =~ /(PRIMER_SEQUENCE_ID=)(.+?)(_[0-9]+,[0-9]+$)/;
    my $target_p = $3;
    $code = $2;
    $target_p =~ s/_//;
#    my @target_pos_l= split(',', $target_p);
#    print "$code $target_p\n";
    pop @seq_res;
 #   print "\n***********************\n@seq_res\n------------------------\n";

    my %pr_pair = initialise_pr_pair();
    foreach my $line (@seq_res)
    {
      chomp $line;
      my @line = split('=', $line);
      if ($line[0] =~ /PRIMER_PRODUCT_SIZE[_0-9]*$/) # end of the data on one primer pair
      {
        $line[0] =~ s/_[0-9]+$//;
        $pr_pair{$line[0]} = $line[1];
        $pr_pair{TARGET} = $target_p;
        add_info($results, \%pr_pair, $code, $condition, $target_ms_info);
        %pr_pair = initialise_pr_pair();
      }
      else
      {
        $line[0] =~ s/_[0-9]+$//;
        $line[0] =~ s/_[0-9]+_/_/;     
        if(exists $pr_pair{$line[0]})
        {
          $pr_pair{$line[0]} = $line[1];
        }
      }
    }# end foreach $line
  }# end if scalar @seq_res>
  @seq_res = ();
 }# end if $l eq "=\n"
}# end while
 
 
 close PRO;
 
}
##################################################
 
sub initialise_pr_pair
{
  my %pr_pair = ();
  $pr_pair{PRIMER_PAIR_PENALTY} = '';
  $pr_pair{PRIMER_LEFT_SEQUENCE} = '';
  $pr_pair{PRIMER_RIGHT_SEQUENCE} = '';
  $pr_pair{PRIMER_LEFT} = '';
  $pr_pair{PRIMER_RIGHT} = '';
  $pr_pair{PRIMER_LEFT_TM} = '';
  $pr_pair{PRIMER_RIGHT_TM} = '';
  $pr_pair{PRIMER_LEFT_END_STABILITY} = '';
  $pr_pair{PRIMER_RIGHT_END_STABILITY} = '';
  $pr_pair{PRIMER_PRODUCT_SIZE} = '';
  return %pr_pair;
}

##################################################

sub add_info
{
my ($results, $pr_pair, $code, $condition, $target_ms_info) = @_;
#  my %target_ms_info = (); # $target_ms_info{seq_code_first_pos} = [first_pos, length_bp, length_in_rn, number_of_ms (1.5 for one compound, 1 for 1 pure), mot_trans, ms_seq ];

my $code_pair = $$pr_pair{'PRIMER_LEFT'}.','.$$pr_pair{'PRIMER_RIGHT'};
delete $$pr_pair{'PRIMER_LEFT'};
delete $$pr_pair{'PRIMER_RIGHT'};

unless (exists $$results{$code}{$code_pair})# add info if primer pair is new to sequences
{
    $$pr_pair{'DESIGN'} = $condition;
    %{$$results{$code}{$code_pair}} = %$pr_pair;
#    print "$code $$results{$code}{$code_pair}{TARGET}\n";
}

}
##################################################################
sub write_primer_table
{
my ($results, $target_ms_info, $primer_file, $polym, $seq, $scoring_matrix) = @_;
#  my %results = (); #$results{$code_seq}{left_pos,left_l,rigth_pos,rigth_l}
#{TARGET/DESIGN/PRIMER_PAIR_PENALTY/PRIMER_PRODUCT_SIZE/PRIMER_LEFT_SEQUENCE/PRIMER_RIGHT_SEQUENCE/
#PRIMER_LEFT_TM/PRIMER_RIGHT_TM/PRIMER_LEFT_END_STABILITY/PRIMER_RIGHT_END_STABILITY}=>  

#  my %target_ms_info = (); # $target_ms{secode_first_pos} = [first_pos, length_bp, length_in_rn, number_of_ms (1.5 for one compound, 1 for 1 pure),, mot_trans, ms_seq ];


open_file_out($primer_file, 'OUT');
print OUT "SEQUENCE_CODE\tNUMBER_OF_READS\tTARGET_REGION_FIRST_POS\tTARGET_REGION_LENGTH_IN_BP\tTARGET_MS_LENGTH_IN_REPEAT_NUMBER\tNUMBER_OF_MS(1.5 for 1 compound, 1 for 1 pure)\tMOT_TRANS\tTARGET_REGION_SEQ\tPOLYMORPH\t";
print OUT "ONE_PRIMER_FOR_EACH_SEQ\tONE_PR_FOR_EACH_TARGET_REGION\tPCR_PRIMER_ALIGNSCORE\tMIN_PRIMER_TARGET_DIST\tPCR_PRODUCT_SIZE\tPCR_PRODUCT_SEQ\tPRIMER_LEFT_SEQUENCE\tPRIMER_RIGHT_SEQUENCE\tPRIMER_LEFT_DIST_FROM_MS\tPRIMER_RIGTH_DIST_FROM_MS\tPRIMER_LEFT_FIRST_POS\tPRIMER_LEFT_LENGTH\tPRIMER_RIGHT_FIRST_POS\tPRIMER_RIGHT_LENGTH\tPRIMER_LEFT_TM\tPRIMER_RIGHT_TM\tPRIMER_LEFT_END_STABILITY\tPRIMER_RIGHT_END_STABILITY\tPRIMER3_PENALTY\tDESIGN\tSEQUENE_LENGTH\tSEQUENCE\n";

foreach my $code (sort keys %$results)
{

	my $read_num = 1;
	if ($code =~ /^cons_gr[0-9]+_([0-9]+)/)
	{
		$read_num = $1;
	}
  my $b_target = '';
  my $b_design = '';
  my $b_penalty = '';
  my @target_order = (); # table that contains the order of preference of different target regions
  my %seq_pr_info = sort_primers($results, $code, $target_ms_info, $$seq{$code}, $scoring_matrix, \@target_order);
	#my %seq_pr_info= (); # $seq_pr_info{taget}{primer_pcr_alignscore_max}{min_dist}{PRIMER_PRODUCT_SIZE}{prp} = all into in tab separated form; hash of hash to determine the 'best primer pair for each target region'

  my $best_seq = 1;
  my $best_target = 1;
  
  foreach my $target (@target_order)
  {
  	$best_target = 1;
    my $polymorph_data = 'NA';
    if($code =~ /^cons/)
    {
	$polymorph_data = 'NO';
    }
    if(exists $$polym{$code})
    {
      $polymorph_data = check_overlap($$polym{$code}, $target);
    }
    foreach my $alignscore (sort {$a<=>$b} keys %{$seq_pr_info{$target}})
    {
      foreach my $min_dist (reverse sort {$a <=> $b} keys %{$seq_pr_info{$target}{$alignscore}})
      {
		foreach my $pcr_size (sort {$a <=> $b} keys %{$seq_pr_info{$target}{$alignscore}{$min_dist}})
		{
			foreach my $prp (keys %{$seq_pr_info{$target}{$alignscore}{$min_dist}{$pcr_size}})
			{
			my $line = "$code\t$read_num\t$$target_ms_info{$code.'_'.$target}[0]\t$$target_ms_info{$code.'_'.$target}[1]\t$$target_ms_info{$code.'_'.$target}[2]\t$$target_ms_info{$code.'_'.$target}[3]\t$$target_ms_info{$code.'_'.$target}[4]\t$$target_ms_info{$code.'_'.$target}[5]\t$polymorph_data\t";
			$line .= "$best_seq\t$best_target\t";
			$best_seq = 0;
			$best_target = 0;
			my $sl = length $$seq{$code};
			$line .= "$seq_pr_info{$target}{$alignscore}{$min_dist}{$pcr_size}{$prp}\t$sl\t$$seq{$code}\n";
#			$line =~  s/;/\t/g;
			print OUT $line;
			}
		}
      } # end foreach $penalty
    }  # end foreach $design
  } # end foreach $target
} # end foreach $code
close OUT;
}

##################################################################

sub sort_primers
{
my ($results, $code, $target_ms_info, $seq, $scoring_matrix, $target_order) = @_;
#  my %target_ms_info = (); # $target_ms_info{seq_code_target-first-pos} = [first_pos, length_bp, length_in_rn, number_of_ms (1.5 for one compound, 1 for 1 pure),, mot_trans, ms_seq ];


# make a hash with inro on all primers pairs of a given sequence ($code)
# get best primer pair for a target MS => 
# 1. lowest alignscore between primers and PCR product (wo primer)
# 2. highest min distance between primers and target region 
# 3. shortest PCR product

# get best primer pair for sequences by taking choosing 1 pure MS => if several=> choosing the longest (in rep number)
# if only multiply or compound MS take the target region with the longest MS (in rep number)

my %seq_pr_info= (); # $seq_pr_info{taget}{primer_pcr_alignscore_max}{min_dist}{PRIMER_PRODUCT_SIZE}{prp} = all into in tab separated form hash of hash to determine the 'best primer pair for each target region'
my %targets = (); # $targets{p/c}{length_in_rn}{length mot_trans}{target} = '';
my $revcomp_seq = reverse_complement_loc($seq);
my %primer_scores = (); # keep the aligment score between seq and primers

foreach my $prp (keys %{$$results{$code}})
{
	$targets{$$target_ms_info{$code.'_'.$$results{$code}{$prp}{TARGET}}[3]}{$$target_ms_info{$code.'_'.$$results{$code}{$prp}{TARGET}}[2]}{length($$target_ms_info{$code.'_'.$$results{$code}{$prp}{TARGET}}[4])}{$$results{$code}{$prp}{TARGET}} = '';

	my $fw = $$results{$code}{$prp}{PRIMER_LEFT_SEQUENCE};
	my $rev = $$results{$code}{$prp}{PRIMER_RIGHT_SEQUENCE};
	my @pos = split(',', $prp);   # [left_pos,left_l,rigth_pos,rigth_l] of primers		
	
	my ($pcr_prod, $pcr_prod_min)  = get_pcr_product_pr($seq, $fw, $rev); # get longest and shortest possible PCR product
	my $max_pcr_prod_length = length $pcr_prod;
	unless ($pcr_prod eq $pcr_prod_min) #if more than one PCR product possible (more than one exact match of one of the primers)
	{
		($pos[0], $pos[2]) = get_pr_pos($seq, $pcr_prod_min, $fw, $rev); # re-define primer position to correspond to the smallest PCR product
	}

    my $rigth_dist = $pos[2] - $pos[3] - $$target_ms_info{$code.'_'.$$results{$code}{$prp}{TARGET}}[0] - $$target_ms_info{$code.'_'.$$results{$code}{$prp}{TARGET}}[1] +1 ;
    my $left_dist = $$target_ms_info{$code.'_'.$$results{$code}{$prp}{TARGET}}[0] - $pos[0] - $pos[1];
	my $min_dist = min($rigth_dist, $left_dist); # min primer-MS distance in the shortest PCR product
	
# calculate the max alignment score between the sequence and the primers (taking into account the rev comp)( do not count the first perfect match, since it is a cumpulsary primer sequence match)
	unless (exists $primer_scores{$fw})
	{
	# take the second best score if the perfect match of primer-seq expected (fw-seq, rev-revcompseq) and the best score otherwise
		my ($score1_first, $score1_second) = gapless_alignment_two_best($seq, $fw, $scoring_matrix);
		my($score2_first, $score2_second) = gapless_alignment_two_best($revcomp_seq, $fw, $scoring_matrix);	
		@{$primer_scores{$fw}} = ($score1_second, $score2_first);
	}
	unless (exists $primer_scores{$rev})
	{
		my ($score1_first, $score1_second) = gapless_alignment_two_best($seq, $rev, $scoring_matrix);
		my ($score2_first, $score2_second) = gapless_alignment_two_best($revcomp_seq, $rev, $scoring_matrix);	
		@{$primer_scores{$rev}} = ($score1_first, $score2_second);
	}
	my $primer_pcr_alignscore_max = max($primer_scores{$fw}[0], $primer_scores{$fw}[1], $primer_scores{$rev}[0], $primer_scores{$rev}[1]);

	# if a primer matches the sequence mre than once, the PCR prodect and its length is given for the longest PCR product,
	# BUT the primer ositions and the primer-MS distance is given for the smallest PCR product (worst case scenarios)
	$seq_pr_info{$$results{$code}{$prp}{TARGET}}{$primer_pcr_alignscore_max}{$min_dist}{$$results{$code}{$prp}{PRIMER_PRODUCT_SIZE}}{$prp} = 
	"$primer_pcr_alignscore_max"."\t"."$min_dist"."\t".
	"$max_pcr_prod_length"."\t"."$pcr_prod"."\t".
	"$$results{$code}{$prp}{PRIMER_LEFT_SEQUENCE}"."\t"."$$results{$code}{$prp}{PRIMER_RIGHT_SEQUENCE}"."\t".
	"$left_dist"."\t"."$rigth_dist"."\t"."$pos[0]"."\t"."$pos[1]"."\t"."$pos[2]"."\t"."$pos[3]"."\t".
	"$$results{$code}{$prp}{PRIMER_LEFT_TM}"."\t"."$$results{$code}{$prp}{PRIMER_RIGHT_TM}"."\t".
	"$$results{$code}{$prp}{PRIMER_LEFT_END_STABILITY}"."\t"."$$results{$code}{$prp}{PRIMER_RIGHT_END_STABILITY}"."\t".
	"$$results{$code}{$prp}{PRIMER_PAIR_PENALTY}"."\t"."$$results{$code}{$prp}{DESIGN}";

#  $seq_pr_info{$$results{$code}{$prp}{TARGET}}{$$results{$code}{$prp}{DESIGN}}{$$results{$code}{$prp}{PRIMER_PAIR_PENALTY}} = "$primer_pcr_alignscore_max".';'."$min_dist".';'."$$results{$code}{$prp}{PRIMER_PRODUCT_SIZE}".';'."$pcr_prod".';'."$$results{$code}{$prp}{PRIMER_LEFT_SEQUENCE}".';'."$$results{$code}{$prp}{PRIMER_RIGHT_SEQUENCE}".';'."$left_dist".';'."$rigth_dist".';'."$pos[0]".';'."$pos[1]".';'."$pos[2]".';'."$pos[3]".';'."$$results{$code}{$prp}{PRIMER_LEFT_TM}".';'."$$results{$code}{$prp}{PRIMER_RIGHT_TM}".';'."$$results{$code}{$prp}{PRIMER_LEFT_END_STABILITY}".';'."$$results{$code}{$prp}{PRIMER_RIGHT_END_STABILITY}".';'."$$results{$code}{$prp}{PRIMER_PAIR_PENALTY}".';'."$$results{$code}{$prp}{DESIGN}";
}

# my %targets = (); # $targets{$number of ms}{length_in_rn}{length mot_trans}{target_first_pos,length}
# $target_order # table ref, that contains the order of preference of different target regions
foreach my $ms_number(sort {$a<=>$b} keys %targets)
{
	foreach my $rn (reverse sort {$a<=>$b} keys %{$targets{$ms_number}})
	{
		foreach my $mot_l (reverse sort {$a<=>$b} keys %{$targets{$ms_number}{$rn}})
		{
			foreach my $target (keys %{$targets{$ms_number}{$rn}{$mot_l}})
			{
				push(@$target_order, $target);
			}
		}
	}
}

return %seq_pr_info;
}
  

#####################################################################

sub get_pcr_product_pr
{
my ($seq, $fw, $rev) = @_;

$rev = reverse_complement_loc($rev);	
my $pcr_prod_wo_primer = $seq;
my $fw_mask = '';
my $rev_mask = '';
if ($pcr_prod_wo_primer =~ /($fw)(.*)/i)
{
	$fw_mask = $1;
	$pcr_prod_wo_primer = $2;
}
else
{
	print "$fw is not found in $seq\n";
}

if ($pcr_prod_wo_primer =~ /(.*)($rev)/i)
{
	$pcr_prod_wo_primer = $1;
	$rev_mask = $2;
}
else
{
	print "$rev is not found in $seq\n";
}
my $pcr_prod = $fw_mask.$pcr_prod_wo_primer.$rev_mask;


my $pcr_prod_wo_primer_min = $pcr_prod;
if ($pcr_prod_wo_primer_min =~ /.*($fw)(.*)/i)
{
	$fw_mask = $1;
	$pcr_prod_wo_primer_min = $2;
}
else
{
	print "$fw is not found in $seq\n";
}

if ($pcr_prod_wo_primer_min =~ /(.*?)($rev)/i)
{
	$pcr_prod_wo_primer_min = $1;
	$rev_mask = $2;
}
else
{
	print "$rev is not found in $seq\n";
	$pcr_prod_wo_primer_min = $pcr_prod_wo_primer
}
my $pcr_prod_min = $fw_mask.$pcr_prod_wo_primer_min.$rev_mask;

return ($pcr_prod, $pcr_prod_min);
  
}

##########################################################

sub get_pr_pos
{
my ($seq, $pcr_prod, $fw, $rev) = @_;

$seq =~ /$pcr_prod/ig;
my $right_p = (pos $seq);
my $left_p = $right_p - (length $pcr_prod) +1;
return ($left_p, $right_p);

}
###################################################################"

sub make_polymorph_info
{
my ($file) = @_;
my %polym = ();
# %polym is a hash of tables with cons seq code as keys and a table of motif_fpos_lapos as value
# reads info in the definition line of the consensus sequences

 open_file_in($file, 'IN');
 while (my $line = <IN>)
 {
    if ($line =~ />cons.[^\s]+/)
    {
      my $code = $&;
      chomp $line;
      $line =~ s/$code//;
      $line =~ s/\s//;
      $code =~ s/>//;
      unless ($line eq '')
      {
        my @line = split(';', $line);
        @{$polym{$code}} = @line;
      }
    }
 }
 close IN;
 return %polym;
}

##################################################

sub check_overlap
{
my ($polym_ref, $target) = @_;
my $poly_data = 'NO'; 
 #  $target_ms_ref if a ref to a table with [first_pos, length_bp, length_in_rn, p/c, mot_trans, ms_seq ];
  # polym_ref if a ref to a table with list of [motif_fpos_lpos]
  
  my @target_p_l = split(',', $target);
  foreach my $pol (@$polym_ref)
  {
    my @list = split('_', $pol);
    if ($target_p_l[0] >= $list[1] and $target_p_l[0] <= $list[2]) # first pos of target MS is between polymoprph pos
    {
      $poly_data = $pol;
      last;
    }
    if(($target_p_l[0] + $target_p_l[1]-1) >= $list[1] and ($target_p_l[0] + $target_p_l[1]-1)<= $list[2])# last pos of target MS is between polymoprph pos
    {
      $poly_data = $pol;
      last;    }
  }
  return $poly_data;
}

##################################################

sub make_fasta_file_with_designed_seq
{
 my ($results, $seq, $file) = @_;
 
 open_file_out($file, 'SEQ');
 my %seq_list = ();
 foreach my $code (sort keys %$results)
 {
      print SEQ ">$code\n";
	  print SEQ cut_up_fasta_line($$seq{$code}, '100');
 }
 
 close SEQ;
}

##################################################################

sub define_param_limits_pipe3
{
my ($param) = @_;
#hash of array [type(list, real, folder), error message, ['list of values/lower_upper limits)']

  my %param_limits = (
      'galaxy'=> ['list', "galaxy must be 0 or 1\n", [0,1]], 
      'syst'=> ['list', "Operating system must be win or linux\n", ['win', 'linux']], 
      'out_folder' => ['folder', "Cannot open the output_folder $$param{out_folder}\n"], 
      'primer3_path'=> ['folder', "Cannot open primer3 folder ($$param{primer3_path})\n"], 
      'primer3_version'=> ['list', "primer3_version must be 1 or 2\n", [1,2]],
      'del_files'=> ['list', "del_files must be 0 or 1\n", [0,1]],
      'input_file' => ['file', "Cannot open input_file ($$param{input_file})\n"],
      'qdd_folder' => ['folder', "Cannot open qdd_folder ($$param{qdd_folder} folder)\n"],  
	  'outfile_string' => ['string', "-outfile_string can contain only alphanumerical values, underscore '_' or dot\n"],
      'debug' => ['list', "debug must be 0 or 1\n", [0,1]],
      'file_pipe3_inp' => ['string', "file_pipe3_inp ($$param{file_pipe3_inp}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_primers' => ['string', "file_primers ($$param{file_primers}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_targets' => ['string', "file_targets ($$param{file_targets}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'pcr_min' => ['real', "Minimum PCR Product size must be between 40 and 10000\n", [40,10000]],
      'pcr_max' => ['real', "Maximum PCR Product size must be between 40 and 10000\n", [40,10000]],
      'pcr_step' => ['real', "PCR Product size interval must be between 20 and 10000\n", [20,10000]],
      'PRIMER_GC_CLAMP' => ['real', "PRIMER_GC_CLAMP must be between 0 and 20\n", [0,20]],
      'PRIMER_OPT_SIZE' => ['real', "PRIMER_OPT_SIZE must be between 1 and 50\n", [1,50]],
      'PRIMER_MIN_SIZE' => ['real', "PRIMER_MIN_SIZE must be between 1 and 50\n", [1,50]],
      'PRIMER_MAX_SIZE' => ['real', "PRIMER_MAX_SIZE must be between 1 and 50\n", [1,50]],
      'PRIMER_OPT_TM' => ['real', "PRIMER_OPT_TM must be between 1 and 100\n", [1,100]],
      'PRIMER_MIN_TM' => ['real', "PRIMER_MIN_TM must be between 1 and 100\n", [1,100]],
      'PRIMER_MAX_TM' => ['real', "PRIMER_MAX_TM must be between 1 and 100\n", [1,100]],
      'PRIMER_MAX_DIFF_TM' => ['real', "PRIMER_MAX_DIFF_TM must be between 0 and 100\n", [0,100]],
      'PRIMER_OPT_GC_PERCENT' => ['real', "PRIMER_OPT_GC_PERCENT must be between 1 and 100\n", [1,100]],
      'PRIMER_MAX_GC' => ['real', "PRIMER_MAX_GC must be between 1 and 100\n", [1,100]],
      'PRIMER_MIN_GC' => ['real', "PRIMER_MIN_GC must be between 1 and 100\n", [1,100]],
      'PRIMER_SELF_ANY' => ['real', "PRIMER_SELF_ANY must be between 1 and 100\n", [1,100]],
      'PRIMER_SELF_END' => ['real', "PRIMER_SELF_END must be between 1 and 100\n", [1,100]],
      'PRIMER_MAX_POLY_X' => ['real', "PRIMER_MAX_POLY_X must be between 1 and 10\n", [1,10]],
      'PRIMER_NUM_RETURN' => ['real', "PRIMER_NUM_RETURN must be between 1 and 10\n", [1,10]],
	  'qual' => ['list', "qual must be 0 or 1\n", [0,1]],
	'contig' => ['list', "contig must be 0 or 1\n", [0,1]]
      );
	  

	  if($param{qual})
	  {
		$param_limits{qual_file} = ['file', "Cannot open qual_file ($$param{qual_file})\n"];
		$param_limits{min_qs} = ['real', "min_qs must be between 0 and 40", [0,40]];
	  }


return %param_limits;
}
##################################################################
sub get_contig_code_and_pos_on_contig
{
 my ($pr_file) = @_;

# complete primer file with two columns cotaining the code of the contig and the firts position of the fragment on the contig to facilitate de selection of markers on the same contig but far away to reduce linkage


open(PR, $pr_file) or die "Cannot open primer table file ($pr_file)\n";
my @primers = <PR>;
close PR;

my $title = shift@primers;
chomp $title;

open(PRO, ">$pr_file") or die "Cannot open primer table output file ($pr_file)\n";
print PRO "$title\tCONTIG_CODE\tFIRST_POS_ON_CONTIG\n";
foreach my $line (@primers)
{
	chomp $line;
	my @line = split("\t", $line);
	my $scode = $line[0];
	$scode =~  s/_([0-9]+)$//;
	my $pos = $1;
	print PRO "$line\t$scode\t$pos\n";
}
close PRO;
@primers = ();
}
######################################################

sub get_distance_to_nearest_neigbour
{
 my ($pr_file) = @_;

# complete primer file with a column cotaining the distance in bases to the nearest neigpur on the contig. NA if only one marker is designed on one contig

open(PR, $pr_file) or die "Cannot open primer table file ($pr_file)\n";
my @primers = <PR>;
close PR;
my %seq_code = (); # $all_primers{seq_code}{first_pos}[list of all lignes for seq_code_fpos]
my $title = shift@primers;
chomp $title;

foreach my $line (@primers)
{
	chomp $line;
	my @line = split("\t", $line);
	my $scode = $line[0];
	$scode =~  s/_([0-9]+)$//;
	my $pos = $1;
#	print "$line[0] $scode $pos\n";
	push(@{$seq_code{$scode}{$pos}{$line[1]}}, $line);
	
}

#print Dumper(\%seq_code);
my %neighbour = ();
foreach my $scode (keys %seq_code)
{
	my @pos = sort{$a<=>$b} keys %{$seq_code{$scode}};
	if(scalar@pos ==1)
	{
		$neighbour{$scode.'_'.$pos[0]} = 'NA';
	}
	else
	{
		$neighbour{$scode.'_'.$pos[0]} = $pos[1] - $pos[0];
		for(my $i = 1; $i < (scalar@pos-1); ++$i)
		{
			my $d1 = $pos[$i] - $pos[$i-1];
			my $d2 = $pos[$i+1] - $pos[$i];
			my $dist = min($d1, $d2);
			$neighbour{$scode.'_'.$pos[$i]} = $dist;
		}
		$neighbour{$scode.'_'.$pos[-1]} = $pos[-1] - $pos[-2];
			
	}
}


open(PRO, ">$pr_file") or die "Cannot open primer table file ($pr_file)\n";
print PRO "$title\tDIST_TO_NEAREST_NEIGHBOUR\n";
foreach my $scode (sort keys%seq_code)
{
	foreach my $fpos (sort{$a<=>$b}keys %{$seq_code{$scode}})
	{
		my $seqcode_fr = $scode.'_'.$fpos;
		foreach my $taget_region_fp (sort{$a<=>$b}keys %{$seq_code{$scode}{$fpos}})
		{
			foreach my $line (@{$seq_code{$scode}{$fpos}{$taget_region_fp}})
			{
			print PRO "$line\t$neighbour{$seqcode_fr}\n";
			}
		}
	}
}
close PRO;

}

######################################################

sub gapless_alignment
{
my ($s1, $s2, $matrix) = @_;
my $score = 0;
my $l1 = length $s1;
my $l2 = length $s2;

my $l_min= min($l1, $l2); 
my $al_max_rep = abs($l1-$l2) + 1;
my @al = make_table(1, $l_min, $al_max_rep, 1); # length of the fragments in each alignement position
my @p1 = make_table_p1(($l1-1), '0', $l2); # positions of the begining of the fragment of seq 1
my @p2 = make_table_p2('0', $l1, ($l2-1)); # positions of the begining of the fragment of seq 2

#print "$s1\n$s2\n";

my $best_score = -10000000000000000;
my $best_al ='';
for(my $i=0; $i<scalar@al; ++$i)
{
	my $f1 = substr($s1, $p1[$i], $al[$i]);
	my $f2 = substr($s2, $p2[$i], $al[$i]);
	my $score = compare($f1, $f2, $matrix);
#	print "$f1\n$f2\n$score\n";
	if ($score>$best_score)
	{	
		$best_score = $score;
		$best_al = $f1."\n".$f2;
	}
}
return ($best_score, $best_al);
}

######################################################

sub gapless_alignment_two_best
{
my ($s1, $s2, $matrix) = @_;
my $score = 0;
my $l1 = length $s1;
my $l2 = length $s2;

my $l_min= min($l1, $l2); 
my $al_max_rep = abs($l1-$l2) + 1;
my @al = make_table(1, $l_min, $al_max_rep, 1); # length of the fragments in each alignement position
my @p1 = make_table_p1(($l1-1), '0', $l2); # positions of the begining of the fragment of seq 1
my @p2 = make_table_p2('0', $l1, ($l2-1)); # positions of the begining of the fragment of seq 2

#print "$s1\n$s2\n";
my @scores = ();

for(my $i=0; $i<scalar@al; ++$i)
{
	my $f1 = substr($s1, $p1[$i], $al[$i]);
	my $f2 = substr($s2, $p2[$i], $al[$i]);
	my $score = compare($f1, $f2, $matrix);
	push(@scores, $score);
}
#print "@scores\n";
@scores = sort {$a <=>$b}@scores;
my $best_score = 0;
my $second_best = 0;
if(scalar @scores >1)
{
	$best_score = $scores[-1];
	$second_best = $scores[-2];
}
#print "@scores\n";
return ($best_score, $second_best);
}

########################################################

sub make_table
{
my ($fv, $maxv, $rep, $lv) = @_;
my @al = ();

for(my $i = $fv; $i<$maxv; ++$i)
{
	push(@al, $i);
}
for(my $i = 0; $i<$rep; ++$i)
{
	push(@al, $maxv);
}
for(my $i = $maxv-1; $i>$lv-1; --$i)
{
	push(@al, $i);
}
return @al;
}
################################################

sub make_table_p1
{
my ($max, $min, $rep) = @_;
my @al = ();

for(my $i = $max; $i>$min; --$i)
{
	push(@al, $i);
}
for(my $i = 0; $i<$rep; ++$i)
{
	push(@al, $min);
}
return @al;
}

##########################################

sub make_table_p2
{
my ($min, $rep, $max) = @_;
my @al = ();

for(my $i = 0; $i<$rep; ++$i)
{
	push(@al, $min);
}
for(my $i = $min+1; $i<$max+1; ++$i)
{
	push(@al, $i);
}

return @al;
}

#######################################################
sub compare
{
my ($f1, $f2, $matrix) = @_;

my $score = 0;

$f1 = uc $f1;
$f2 = uc $f2;
my @b1 = split('', $f1);
my @b2 = split('', $f2);

for(my $i=0; $i<scalar@b1; ++$i)
{
	if($b1[$i] eq 'N' or $b2[$i] eq 'N')
	{
		$score += $$matrix{'N'};
	}
	elsif ($b1[$i] eq $b2[$i])
	{
		$score += $$matrix{'match'};	
	}
	else
	{
		$score += $$matrix{'mismatch'};			
	}
}
return $score;
}
##################################################