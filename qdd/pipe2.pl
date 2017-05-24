#!/usr/bin/perl -w


  use warnings;
  use strict;
  use subprogramQDD;
  use Data::Dumper;
  
  
#INPUT: fasta file with clean (no adapter, no tag) microsatellite containing sequences.

# output: 
# file_root_singleton.fas => the only BLASt hit is autohit; All >=4 repeat microsatellites and >=8 homoplymers are masked
# file_root_nohit_css.fas => no hit to itself; All >=4 repeat microsatellites and >=8 homoplymers are masked
# file_root_multihit_css.fas => more than one hit between a pair of sequences; All >=4 repeat microsatellites and >=8 homoplymers are masked
# file_root_gr.fas => sequences (including consensuses) that had BLAST hit to other sequences, with bellow limit identity of the overlaping region. Regions covered by BLAST hits are masked 
# file_root_consensus.fas => all unique (no hit to grouped seqs) consensus sequences, if polymorphism of ms motif_first_last pos of MS (in the aligned cons seq) is given; All >=4 repeat microsatellites are masked, but not homopolymers;
# file_root_cons_subs.fas => Consensus + aligned sequences to make a consensus ; All >=4 repeat microsatellites are masked, but not homopolymers; grouped consensus sequences are NOT removed
# file_root_pipe2_for_pipe3.fas => fasta file with all unique sequences (singletons and consensus); input for pipe3

# possibility of running RepeatMasker before de all-against-all blast, and eliminate sequences that had hits to TEs. However, this tsep #is reltively slow and one loses information on the eliminated sequences. Therefore, although this step can be called it is not #documented for users. The doc guides users to run RepeatMasker on sequences that had primers, (it is also integrated in pipe3.pl)
# and the seqence still appears in the primer table, with info on the TE hit if any.



my $set_file_default = '/etc/qdd/set_qdd_default.ini';
my $mask_mono = 8;
my $mask_di = 4;
my $mask_tri = 4;
my $mask_tetra = 4;
my $mask_penta = 4;
my $mask_hexa = 4;
my $e = 1e-40;
my $blast_version = 'blast+';

$set_file_default = check_set_file($set_file_default);

print_heading_3_0('STDOUT');
print "SCRIPT pipe2.pl\n";

my %param = (
  galaxy => '', # one if run on galaxy, 0 for command line
  syst => '',  
  blast_path => '',
  make_cons => '', # if 0 make only a file with singleton sequences
  ident_limit => '',
  clust_path => '',
  prop_maj => '',
  input_file => '', # input file with full path
  out_folder => '',
  del_files => '',
#  file_msfas => '_pipe1_for_pipe2.fas', # motif of the output pipe1 that will be the input pipe2
  qdd_folder => '',
  file_log_pipe2 => '_pipe2_log.txt',
  file_singleton => '_pipe2_singleton.fas',
  file_nohit_css => '_pipe2_nohit_css.fas',
  file_multihit_css => '_pipe2_multihit_css.fas',
  file_gr => '_pipe2_groupped.fas',
  file_consensus => '_pipe2_consensus.fas',
  file_cons_subs => '_pipe2_cons_subs.fas',
  file_pipe3_inp => '_pipe2_for_pipe3.fas',
  outfile_string => '',
  debug => 0,
  num_threads => 1,
  file_rm_hit => '_pipe2_rm_hit.fas',
  rm_pipe2 => 0,
  rm_path => '/usr/local/RepeatMasker/',
  rm_lib => 'metazoa'
);

my $clustal_ex = '';
my %param_limits = define_param_limits_pipe2(\%param); #hash of array $param_limits{syst}[type(list, real, folder), error message, ['list of values/lower_upper limits)']
read_set_file_galaxy(\%param, $param{qdd_folder}.$set_file_default);
modify_params_from_tags(\%param, \@ARGV);
add_slash_to_folder(\%param, \%param_limits);
my $galaxy = $param{galaxy}; # one if run on galaxy, 0 for command line

if ($param{make_cons} eq 'false' or $param{make_cons} eq 'FALSE' or $param{make_cons} eq '0')
{
 $param{make_cons} = 0;
}
else
{
 $param{make_cons} = 1;
}

my $folder_id = time;
my $date = localtime $folder_id;
my $tempfolder = $param{out_folder}.'pipe2_'.$folder_id.'/';
make_folder($param{syst}, $tempfolder);

if($param{rm_pipe2})
{
print "Checking the validity of $param{rm_lib} for selecting elements from the RepeatMasker database\n\n";
check_repeatmasker_db($param{rm_lib}, $param{rm_path}, $tempfolder);
}

if($galaxy)
{
	$param{outfile_string} = 'NA';
}
elsif ($param{outfile_string} eq '')
{
	$param{outfile_string} = get_filename_root($param{input_file});
	$param{outfile_string} =~ s/_pipe1_for_pipe2.*//;
}


unless($galaxy)
{
check_param_limits(\%param, \%param_limits);
#print "$clustal_ex\n";

	if ($param{outfile_string} eq '')
	{
		$param{outfile_string} = get_filename_root($param{input_file});
		$param{outfile_string} =~ s/_pipe1_for_pipe2.*//;
	}
	my @clefs = ('file_singleton', 'file_nohit_css', 'file_multihit_css', 'file_gr', 'file_consensus', 'file_cons_subs', 'file_pipe3_inp', 'file_rm_hit');
	if($param{file_log_pipe2} =~ /^_/) # if the name of the log file has not been specified (run in pipeline)
	{
		push(@clefs, 'file_log_pipe2');
	}

	get_last_version_and_modify_1($param{out_folder}, \@clefs, $param{outfile_string}, \%param);

	if($param{file_log_pipe2} =~ /^_/) # if the name of the log file has not been specified (run in pipeline)
	{
		$param{file_log_pipe2} = $param{out_folder}.$param{outfile_string}.$param{file_log_pipe2};
	}


	$param{file_singleton} = $param{out_folder}.$param{outfile_string}.$param{file_singleton};
	$param{file_nohit_css} = $param{out_folder}.$param{outfile_string}.$param{file_nohit_css};
	$param{file_multihit_css} = $param{out_folder}.$param{outfile_string}.$param{file_multihit_css};
	$param{file_gr} = $param{out_folder}.$param{outfile_string}.$param{file_gr};
	$param{file_consensus} = $param{out_folder}.$param{outfile_string}.$param{file_consensus};
	$param{file_cons_subs} = $param{out_folder}.$param{outfile_string}.$param{file_cons_subs};
	$param{file_pipe3_inp} = $param{out_folder}.$param{outfile_string}.$param{file_pipe3_inp};
	$param{file_rm_hit} = $param{out_folder}.$param{outfile_string}.$param{file_rm_hit};
}
$clustal_ex = get_clustal_executable_name($param{clust_path}, $param{syst});

open(my $fh_log, '>>', $param{file_log_pipe2}) or die "Cannot open $param{file_log_pipe2} $!\n";
print_heading_3_0($fh_log);
print $fh_log "SCRIPT pipe2.pl\n";

if($param{debug})
{
	print $fh_log Dumper(\%param);
	print "\n";
}


print $fh_log "INPUT file : $param{input_file}\n";
print "Input file : $param{input_file}\n\n";

print $fh_log "\nOUTPUT files:\n";
print $fh_log "String for naming output files: $param{outfile_string}\n";
if($param{rm_pipe2})
{
print $fh_log "\t$param{file_rm_hit}: Sequences that had a hit in RepeatMasker\n"
}

print $fh_log "\t$param{file_nohit_css}: Low complexity sequences\n";
print $fh_log "\t$param{file_multihit_css}: Putatif minisatellites\n";
if($param{make_cons})
{
print $fh_log "\t$param{file_singleton}: Singleton sequences (the only hit is an autohit)\n";
print $fh_log "\t$param{file_gr}: Sequences similar to others but bellow $param{ident_limit} % identity\n";
print $fh_log "\t$param{file_consensus}: All unique consensus sequences\n";
print $fh_log "\t$param{file_cons_subs}: Consensus + aligned reads\n";
}
print $fh_log "\t$param{file_pipe3_inp}: All unique sequences (including consensus if make_cons =1)\n";



print $fh_log "\nGENERAL PARAMETERS\n";
print $fh_log "System (win/linux): $param{syst}\n";
print $fh_log "Output folder: $param{out_folder}\n";
print $fh_log "Path to QDD executables: $param{qdd_folder}\n";
print $fh_log "Path to BLAST executables: $param{blast_path}\n";
print $fh_log "Number of threads for BLAST: $param{num_threads}\n";
print $fh_log "Path to Clustal2 executables: $param{clust_path}\n";
print $fh_log "Make consensus sequences(0/1): $param{make_cons}\n";
if($param{make_cons})
{
print $fh_log "Minimum identity between two sequences to make a consensus: $param{ident_limit}\n";
print $fh_log "Proportion of sequences that must have the same base on the aligned site to accept it as a consensus: $param{prop_maj}\n";
}
print $fh_log "Delete intremediate files: $param{del_files}\n";
print $fh_log "E-value for all-against-all BLAST: $e\n";

print $fh_log "************************************\n";

my %seq_count = ();

# run RepatMasker
my $input_for_ms_find1 = $param{input_file};
if($param{rm_pipe2})
{
	print "Running RepeatMasker\n";
	print $fh_log "Run RepeatMasker: $param{rm_pipe2}\n";	
	print $fh_log "RepeatMasker library: $param{rm_lib}\n";
	print $fh_log "Path to RepeatMasker: $param{rm_path}\n";
	my $rm_tbl = run_RepeatMasker($param{input_file}, $param{rm_path}, $param{rm_lib}, $param{num_threads}, $tempfolder, $param{debug}, $fh_log);
	$input_for_ms_find1 = $tempfolder.'input_for_ms_find1.fas';
# select sequences that did not have hit to TE elements
	select_seq_wo_rm_hit($param{input_file}, $rm_tbl, $input_for_ms_find1, \%seq_count, $param{file_rm_hit});
}

# mask MS and homopolymers in the input fasta file use low limit values to avoid false BLAST hits  and eliminate CSS sequences (no autohit)
my $temp_ms_find1_fas = $tempfolder.'ms_find1_fas.fas';
my $temp_ms_find1_ms = $tempfolder.'ms_find1_ms.txt';
my $find_ms1 = 'perl "'.$param{qdd_folder}.'find_and_mask_ms_galaxy.pl" '.$mask_mono.' '.$mask_di.' '.$mask_tri.' '.$mask_tetra.' '.$mask_penta.' '.$mask_hexa.' 1 0 "'.$input_for_ms_find1.'" "'.$temp_ms_find1_fas.'" "'.$temp_ms_find1_ms.'" 1 1';
#  USAGE: $min1 $min2 $min3 $min4	$min5	$min6	$del_prev_mask $keep_non_ms $input_file  $output_fasta $output_ms  $write_seq_file $write_fas_file;

if ($param{debug})
{print $fh_log "\n$find_ms1\n\n";}

print "Soft masking sequences\n";
system $find_ms1;

my $blast_out = $tempfolder.'ms_find1_fas_blastres.txt';
print "All-against-all BLAST\n";
call_blast_query_to_db_soft_mask_1($temp_ms_find1_fas, $temp_ms_find1_fas, $blast_out, $e, $param{blast_path}, $blast_version, $param{syst}, $param{num_threads});

# delete_blast_temp_files_win_linux($outfolder_spec, $param{syst});
delete_tag_adapt_files($temp_ms_find1_fas, $param{syst});

if ($param{make_cons} == 0)
{
   $param{file_singleton} = $param{file_pipe3_inp};
}

print "Select singleton sequences\n";
my %seq = read_hach_from_fasta_code_wo_space($temp_ms_find1_fas);

unless($param{rm_pipe2})
{
$seq_count{'Total number of input sequences:'} = scalar keys%seq;
}
my %css =(); # seq code as key, '' as value
make_singleton_only(\%seq, $param{file_singleton}, $param{file_nohit_css}, $blast_out, \%seq_count, \%css);


  my %blast = (); # hash of anonymous hashes $blast{code1}{code2}=>[blast ident, qfirts, qlast, sfirst, slast]
      # does not contain autohits; only one hit is recorded in case of reciproc hits; 
      # several hits between the two same sequences is marked as css

  print "Sorting BLAST results\n";
  sort_blast_results(\%seq, \%blast, \%css, $blast_out);

  # delete_multihit_css from the %blast hash and from %seq
  # print multihit_css_file;

  print "Remove multihit sequences\n";
  treat_multihit_css(\%blast, \%css, \%seq, $param{file_multihit_css}, \%seq_count);
 
if ($param{make_cons} == 1)
{ 

  # for each remaining blast line align the two sequence in the rigth orientation, and determine %ident for the flanking region
  # if ident is bigger than limit, pool sequences into, one group

  print "Identifying sequence groups for consensus sequences\n";
  my %gr_seq = (); # group_id as key, [seq_code;orientation (1/+1)] as value
  my %seq_gr = (); # seq_id as key, group_id as value

  calc_perc_ident_flanking_make_big_group(\%blast, \%seq, \%gr_seq, \%seq_gr, $param{ident_limit}, $tempfolder, $clustal_ex);

  print "Making consensus sequences\n";
  # conensus file sequences and subcons sequences are with gaps (aligned) 
  # subcons sequences are deleted from %seq
  make_cons_file_1(\%gr_seq, \%seq_gr, \%seq, $param{prop_maj}, $param{file_cons_subs}, $clustal_ex, $tempfolder, \%seq_count);

  my $temp_ms_find2_fas = $tempfolder.'ms_find2_fas.fas';
  my $temp_ms_find2_ms = $tempfolder.'ms_find2_ms.txt';
  my $find_ms2 = 'perl "'.$param{qdd_folder}.'find_and_mask_ms_galaxy.pl" 1000000 5 5 5 5 5 1 1 "'.$param{file_cons_subs}.'" "'.$temp_ms_find2_fas.'" "'.$temp_ms_find2_ms.'" 1 1';
 #  USAGE: $min1 $min2 $min3 $min4	$min5	$min6	$del_prev_mask $keep_non_ms $input_file  $output_fasta $output_ms  $write_seq_file $write_fas_file;

 if($param{debug})
  {  
	print $fh_log "\n$find_ms2\n";
  }

  print "Check polymorphism in consensus sequences\n";
  system $find_ms2;

  my %cons_poly = ();
  check_polymorphism(\%cons_poly, $temp_ms_find2_ms, $temp_ms_find2_fas);

  # replace non_masked sequences of cons subs file by MS masked (4 repeat, but not homopolymer)
  # write consensus file with ms masked and polymorphism info (motif_first_last pas in the aligned cons seq)
  write_consensus_file_with_polymorphism_info(\%cons_poly, $param{file_consensus}, $param{file_cons_subs}, $temp_ms_find2_fas);


  # identify regions of %gr_sequences that have hit to other sequences from %blast; regions of hit printed in lc;
  # make gr_list.fas
  my $gr_file_temp = $tempfolder.'gr_file_temp.fas'; # MS and homopolymer masked; $gr_file hit to other sequeces are masked
  treat_group_seq(\%seq, \%blast, $param{file_gr}, $gr_file_temp, \%seq_count);

  # delete %blast OR make consensus more precise... TO DO
  
 
# BLast consensus.fas againts gr_list.fas
unless(-z $param{file_consensus} or -z $param{file_gr}) #execute blast only if there are groupped or concensus files
{
 	 print "Blasting consensus sequences against grouped sequences\n";
 	my $blast_out_cons_gr = $tempfolder.'cons_gr_blast.txt';


  	call_blast_query_to_db_soft_mask_1($param{file_consensus}, $param{file_gr}, $blast_out_cons_gr, $e, $param{blast_path}, $blast_version, $param{syst}, $param{num_threads});
  	delete_tag_adapt_files($param{file_gr}, $param{syst});


  # if hit of cons, print it as grouped in class file, else print it as unique cons

  # remove cons sequences from consensus file if they had BLAST hit to grouped sequences; and add them to gr_seq

 	 print "Moving consensus sequences with partial hits to grouped sequences into grouped file\n";
 	identify_grouped_cons_seq($blast_out_cons_gr, $param{file_cons_subs}, $param{file_consensus}, $param{file_gr}, \%seq_count);
}

  pool_files($param{file_singleton}, $param{file_consensus}, $param{file_pipe3_inp});

#  print_hachage(\%seq);
} # fin if $make_cons

  print $fh_log "\n\n";
  foreach my $code (sort keys %seq_count)
  {
    print $fh_log "$code $seq_count{$code}\n";
  }
  print $fh_log "\n\n";

  my $time_end = time;
  my $time_run = $time_end - $folder_id;
  print $fh_log "pipe2.pl started at : $date\n";
  print $fh_log "The run took $time_run s.\n";
  
#  move_files_local($outfolder_spec, $pipe2_folder, $param{syst}); # move files from working folder to pipe2 folder
 

delete_blast_temp(\%param);
if ($param{del_files} ==1)
{
  delete_folder($tempfolder, $param{syst});
}


print "\nSee log_file for summary:\n$param{file_log_pipe2}\n\n";
print $fh_log "The analysis is finished\n";

close $fh_log;
exit;

#################################################################

sub select_seq_wo_rm_hit
{
my ($input_fas, $rm_tbl, $output_fas, $seq_count, $rm_outfas) = @_;

my %rm = ();
open(TBL, $rm_tbl) or die "Cannot open rm tbl file ($rm_tbl)\n";
my $title = <TBL>;
$title = <TBL>;
$title = <TBL>;
while(my $line = <TBL>)
{
	chomp $line;
	my @line = split(' ', $line);
	$rm{$line[4]} = '';
}
close TBL;

my %seq = read_hach_from_fasta_code_wo_space($input_fas);
$$seq_count{'Total number of input sequences:'} = scalar keys %seq;
$$seq_count{'Number of sequences with NO hit to trasposable elements:'} = (scalar keys %seq) - (scalar keys %rm);
$$seq_count{'Number of sequences with hit to trasposable elements:'} = scalar keys %rm;

open(OUT, ">$output_fas") or die "Cannot open output fasta file ($output_fas)\n";
open(RM, ">$rm_outfas") or die "Cannot open output fasta file ($rm_outfas)\n";
foreach my $code (sort keys %seq)
{
	if(exists $rm{$code})
	{
		print RM ">$code\n";
		print RM cut_up_fasta_line($seq{$code}, '100');
	}
	else
	{
		print OUT ">$code\n";
		print OUT cut_up_fasta_line($seq{$code}, '100');
	}
}
close OUT;
close RM;
%seq = ();
%rm = ();
}


############################################################

sub delete_blast_temp
{
 my ($param) = @_;
foreach my $code (keys %$param)
{
	my $temp = $$param{$code}
}
}

############################################################


sub sort_blast_results
{
my ($seq, $blast, $css, $blast_out) =@_;

# %seq = (); # sequences that had autohit and to other sequences
# %blast = (); # hash of anonumous hash $blast{code1}{code2}=>[blast line]
# does not contain autohits; only one hit is recorded in case of reciproc hits; 
#several hits between the two same sequences is marked as css
# %css =(); # seq code as key, '' as value

  open_file_in($blast_out, 'BL');
  while (my $line = <BL>)
  {
    my @line = split(' ', $line);
    
    unless ($line[0] eq $line[1])
    {
    if (exists $$seq{$line[0]} and exists $$seq{$line[1]}) 
    {
      if (exists $$blast{$line[0]} and exists $$blast{$line[0]}{$line[1]}) # hit between the two same sequences already tested => css_multihit
      {
        $$blast{$line[0]}{$line[1]} = 'css';
        $$css{$line[0]} = '';
        $$css{$line[1]} = '';     
      }
      elsif (exists $$blast{$line[1]} and exists $$blast{$line[1]}{$line[0]}) # reciproc hit already tested
      {
        next;
      }
      else 
      {
        $$blast{$line[0]}{$line[1]} = [$line[2], $line[6], $line[7], $line[8], $line[9]];
      }
    }
    }

    
  } # end while
  close BL;
} # end sub

########################################################

sub delete_multihit_css_make_gr_seqlist
{
# delete_multihit_css from the %blast hash
# make a hash with the list of seq code that have genuine hit to other sequences and are not css.
# the case of very few sequences if their only valid hit is towards a css_multihit sequence, they are classed as unique, since their css pair is deleted
my ($blast, $css, $autohit) = @_;
my %gr_list = ();

foreach my $code1 (keys %$blast)
{
  if (exists $$css{$code1}) # if code1 css delete all blast hits with code1
  {
    delete$$blast{$code1};
  }
  elsif(exists $$autohit{$code1})
  {
    foreach my $code2 (keys %{$$blast{$code1}})
    {
      if (exists $$css{$code2}) # if code2 css but not code1 delete specific code2 blast hits
      {
        delete$$blast{$code1}{$code2};   
      }
      elsif (exists $$autohit{$code2})
      {
        $gr_list{$code2} = '';
      }
      else # code2 is css no_hit
      {
        delete$$blast{$code1}{$code2};         
      }
    }
    
    if (scalar(keys %{$$blast{$code1}}) == 0) # if all hits to code1 are deleted
    {
      delete$$blast{$code1};
    }
    else
    {
        $gr_list{$code1} = '';
    }
  }
  else # code1 is css no_hit
  {
    delete$$blast{$code1};    
  }
}

return %gr_list;  
} # end sub


##############################################################

sub calc_perc_ident_flanking_make_big_group
{
my ($blast, $seq, $gr_seq, $seq_gr, $ident_limit, $temp_folder, $clustal_ex) =@_;

 #  %gr_seq = (); # group_id (first seq in group) as key, {seq_code =>orientation (1/+1, relative to the group name seq)} as value
 #  %seq_gr = (); # seq_id as key, group_id as value

foreach my $code1 (keys %$blast)
{
  foreach my $code2 (keys %{$$blast{$code1}})
  {
 #   print "$code1; $code2;\n";

   my $orient = determine_orient($$blast{$code1}{$code2});
    # if identity is calculated, BLASt resukts are replaced by flanking ident, oriens seq1 align, seq2 align
  
  unless (exists $$seq{$code2})
  {
    print STDERR "$code2\n";
    print STDERR "$code1 $code2\n";
#    print STDERR Dumper "$$blast{$code1}\n";
  }
   unless (exists $$seq{$code1})
  {
    print STERR "$code1\n";

    print STERR "$code1 $code2\n";
#    print STERR Dumper "$$blast{$code1}\n";
  }
   my $ident = determine_ident($$blast{$code1}{$code2}, $orient, $ident_limit, $$seq{$code1}, $$seq{$code2}, $temp_folder, $clustal_ex);
   
 #   print "$$blast{$code1}{$code2}[0], $$blast{$code1}{$code2}[1]\n";
 #   print "$$blast{$code1}{$code2}[2]\n$$blast{$code1}{$code2}[3]\n\n";

    if ($ident >= $ident_limit)
    {
      if (exists $$seq_gr{$code1} and $$seq_gr{$code2} )
      {
        if (not ($$seq_gr{$code1} eq $$seq_gr{$code2}))
        {
          merge_groups1($gr_seq, $seq_gr, $code1, $code2, $orient);
        }
      }
      elsif (exists $$seq_gr{$code1})
      {
        add_to_group($gr_seq, $seq_gr, $code2, $code1, $orient);
      }
      elsif (exists $$seq_gr{$code2})
      {
        add_to_group($gr_seq, $seq_gr, $code1, $code2, $orient);
      }
      else
      {
        make_new_group($gr_seq, $seq_gr, $code1, $code2, $orient);
      }
    } # end if ident > limit
  } # end foreach $code2
} # end foreach $code1

#print Dumper $gr_seq;
} # end sub

#############################################

sub determine_ident
{
my ($blast, $orient, $ident_limit, $seq1, $seq2, $temp_folder, $clustal_ex) = @_;
# $blast is a ref to an anonymous table with [blast_ident, qfirst, qlast, sfirst, slast]
my %seq_cut = (); # keys: s1before, s1overlap, s1after, s2before, s2overlap, s2after

my $ident = 0;
my $overhang = 50;

if ($$blast[0] >= ($ident_limit - 5)) # examine hits only if the balst identity is at most 5% smaller than limit
{
    # length of the overlaping sequence that is not alingned by BLAST
    # determine %$seq_cut
  my $non_align = determine_non_aligned_length($blast, $orient, $seq1, $seq2, \%seq_cut); # length of the overlaping sequence that is not alingned by BLAST
  if ($non_align < $overhang) # examine hits only if the non_blast_aligned is smaller than 50 bases 
  {
    my %seq1_mspos = determine_ms_pos($seq_cut{s1overlap});
    my %seq2_mspos = determine_ms_pos($seq_cut{s2overlap});
    my %seq_align = align_2seq($seq_cut{s1overlap}, $seq_cut{s2overlap}, $clustal_ex, $temp_folder);
#    print_hachage(\%seq_align);
    $ident= calc_ident(\%seq_align, \%seq1_mspos, \%seq2_mspos);
#    print_hachage(\%seq_align);
#    print TMP "$overhang;$$blast[0];$ident\n";
      # replace the BLAST results by the %identity on the flanking region, $orientation, seq1 align; seq2 align
    if ($ident >= $ident_limit)
    {
      add_seq_cut_to_blast($blast, \%seq_align, \%seq_cut, $ident, $orient); 
    }
  }
}

return $ident;   
}
#############################################

sub determine_orient
{
my ($blast) = @_;
# $blast is a ref to an anonymous table with [blast_ident, qfirst, qlast, sfirst, slast]
my $orient = 1;

if ($$blast[3] > $$blast[4])
{
  $orient = -1;
}

return $orient;
}
#############################################
sub determine_non_aligned_length
{
my ($blast, $orient, $seq1, $seq2, $seq_cut) = @_;
# $blast is a ref to an anonymous table with [blast_ident, qfirst, qlast, sfirst, slast]


my $l1 = length$seq1;
my $l2 = length$seq2;
my $beg = $$blast[1] -1;
my $end = ($l1) - $$blast[2];
if ($orient == 1)
{
  if ($$blast[3] <= $beg)
  {
     $beg = $$blast[3] -1;
  }
  my $temp = ($l2) - $$blast[4];
  if ( $temp < $end)
  {
     $end = $temp;
  }

}
else
{
  my $temp1 = ($l2) - $$blast[3];

  if ($temp1 < $beg)
  {
     $beg = $temp1;
  }
  if ($$blast[4] <= $end)
  {
     $end = $$blast[4] -1;
  }
#   print "@$blast; $l1; $l2; $beg; $end\n";

}
my $non = $beg + $end;

my $s1before = $$blast[1] - $beg;
my $s1overlap = $$blast[2] - $$blast[1] +1 + $non;
my $s1after = $l1 - $$blast[2] - $end;
my $s2before = '';
my $s2overlap = '';
my $s2after = '';

if ($orient == 1)
{
  $s2before = $$blast[3] - $beg;
  $s2overlap = $$blast[4] - $$blast[3] +1 + $non;
  $s2after = $l1 - $$blast[4] - $end;
}
else
{
  $s2before = $l2 - $$blast[3] - $beg;
  $s2overlap = $$blast[3] - $$blast[4] +1 + $non;
  $s2after = $$blast[4] - $end;
  $seq2 = reverse_complement_loc($seq2);
}

#print "$seq1\n";
$$seq_cut{s1before} = substr($seq1, 0, $s1before, '');
$$seq_cut{s1overlap} = substr($seq1, 0, $s1overlap, '');
$$seq_cut{s1after} = $seq1;
$$seq_cut{s2before} = substr($seq2, 0, $s2before, '');
$$seq_cut{s2overlap} = substr($seq2, 0, $s2overlap, '');
$$seq_cut{s2after} = $seq2;

#unless ($orient ==1)
#{
#print "$beg;$end;\n1:$$seq_cut{s1overlap}\n2:$$seq_cut{s2overlap}\n\n";
#}

return $non;
}

######################################################

sub determine_ms_pos
{
my ($seq) = @_;
my %seq_mspos = ();
my @seq = split('', $seq);
for(my $i=0; $i< scalar@seq; ++$i)
{
  if ($seq[$i] =~ /[atgc]/)
  {
    $seq_mspos{$i} = '';
  }
}

return %seq_mspos;
}
######################################################

sub align_2seq
{
my ($seq1, $seq2, $clustal_ex, $temp_folder) = @_;
my $file = $temp_folder.'clustal.fas';
my $out_cl = $temp_folder.'clustal.gde';
my $log = $temp_folder.'logfile';
open_file_out($file, 'IN');
print IN ">seq1\n$seq1\n>seq2\n$seq2\n";
close 'IN';

system '"'.$clustal_ex.'" -INFILE='.$file.' -QUIET -OUTPUT=GDE >'.$log;

open_file_in($out_cl, 'AL');
my $code = '';
my %seqa = ();
while (my $line = <AL>)
{
  chomp $line;
  if ($line =~ /\#/)
  {
     $code = $line;
     $code =~ s/\#//;
  }
  else
  {
    $seqa{$code} .= $line;
  }
}
close AL;
return %seqa;

}

###################################################################
sub calc_ident
{
my ($seq_align, $seq1_ms, $seq2_ms) = @_;
my @seq1 = split('', $$seq_align{seq1});
my @seq2 = split('', $$seq_align{seq2});
my $s1 = -1; 
my $s2 = -1; 
my $flank_l = 0;
my $mism = 0;
for (my $i = 0; $i<scalar@seq1; ++$i)
{
  unless($seq1[$i] eq '-')
  {++$s1;}
  unless($seq2[$i] eq '-')
  {++$s2;}
  if (exists $$seq1_ms{$s1} or exists $$seq2_ms{$s2})
  {next;}
  
  ++$flank_l;
  unless($seq1[$i] eq $seq2[$i])
  { ++$mism;}
}
my $ident = 0;
if ($flank_l >20)
{
$ident = int((($flank_l-$mism)/$flank_l)*100);
}
return $ident;
}
###################################################################
sub add_seq_cut_to_blast
{
my ($blast, $seq_align, $seq_cut, $ident, $orient) = @_;
# $blast is a ref to an anonymous table with [blast_ident, qfirst, qlast, sfirst, slast]
# This is changed to [flank_ident, orient, qalig, salig]
my $seq1 = '';
my $seq2 = '';

# write sequences before the overlap
if ($$seq_cut{s1before} eq '')
{
  unless($$seq_cut{s2before} eq '')
  {
    $seq2 = $$seq_cut{s2before};
    $seq1 = '-' x (length $$seq_cut{s2before});
  }
}
else
{
    $seq1 = $$seq_cut{s1before};
    $seq2 = '-' x (length $$seq_cut{s1before});
}

# write aligned sequences in the overlap
    $seq2 .= $$seq_align{seq2};
    $seq1 .= $$seq_align{seq1};

# write sequences after the overlap
if ($$seq_cut{s1after} eq '')
{
  unless($$seq_cut{s2after} eq '')
  {
    $seq2 .= $$seq_cut{s2after};
    $seq1 .= '-' x (length $$seq_cut{s2after});
  }
}
else
{
    $seq1 .= $$seq_cut{s1after};
    $seq2 .= '-' x (length $$seq_cut{s1after});
}

@$blast = ($ident, $orient, $seq1, $seq2);
}
#################################################################
sub merge_groups1
{
 my ($gr_seq, $seq_gr, $code1, $code2, $orient) =@_;
 
 my $gr = $$seq_gr{$code1};
 my $gr_to_del = $$seq_gr{$code2};
 my $code1_orient = $$gr_seq{$gr}{$code1};

# print Dumper $gr_seq;
 foreach my $code (keys %{$$gr_seq{$gr_to_del}})
 {
  my $o = $$gr_seq{$gr_to_del}{$code};
  $$gr_seq{$gr}{$code} = $o * $code1_orient;
  $$seq_gr{$code} = $gr;
 }
  delete $$gr_seq{$gr_to_del};
 #  print Dumper $gr_seq;

}
#############################################
sub add_to_group
{
 my ($gr_seq, $seq_gr, $new_code, $old_code, $orient) =@_;
# old_code is already in a group
 
 $$seq_gr{$new_code} = $$seq_gr{$old_code};

# print "$new_code, $old_code, $orient\n";
# print Dumper $$gr_seq{$$seq_gr{$old_code}};
 
  my $old_code_orient = $$gr_seq{$$seq_gr{$old_code}}{$old_code};
 $orient = $orient * $old_code_orient;
 $$gr_seq{$$seq_gr{$old_code}}{$new_code} = $orient; 
# print Dumper $$gr_seq{$$seq_gr{$old_code}};
# print "*********************\n";
}

#############################################

sub make_new_group
{
 my ($gr_seq, $seq_gr, $code1, $code2, $orient) =@_;
 #  %gr_seq = (); # group_id (first seq in group) as key, {seq_code =>orientation (1/+1, relative to the group name seq)} as value
 #  %seq_gr = (); # seq_id as key, group_id as value

#print "$code1\n$code2\n";
#print Dumper $gr_seq;
 $$seq_gr{$code1} = $code1;
 $$seq_gr{$code2} = $code1;
 $$gr_seq{$code1}{$code2} = $orient;
 $$gr_seq{$code1}{$code1} = 1;
# print Dumper $gr_seq;

}

#################################################
sub make_cons_file_1
{
my ($gr_seq, $seq_gr, $seq, $prop_maj, $cons_file_subseq, $clustal_ex, $temp_folder, $seq_count_hash) = @_;

my $file = $temp_folder.'clustal.fas';
my $out_cl = $temp_folder.'clustal.gde';
my $log = $temp_folder.'logfile';

open_file_out($cons_file_subseq, 'CONSUB');

#system '"'.$clustal_ex.'" -INFILE='.$file.' -QUIET -OUTPUT=GDE >'.$log;

$$seq_count_hash{'Number of sequences pooled into consensus sequences:'} = 0;

my $gr_count = 0;
foreach my $gr (sort keys %$gr_seq)
{
  ++$gr_count;
  open_file_out($file, 'IN');
  my $seq_count = 1; # code 'seqx' as a key; real code as a value; for clustal short coe is necessary (seqx), but keep track of real codes
  my %code_corr = ();
  foreach my $code (sort keys %{$$gr_seq{$gr}})
  {
    if ($$gr_seq{$gr}{$code} == -1)
    {
      $$seq{$code} = reverse_complement_loc($$seq{$code});
    }
    print IN ">seq_$seq_count\n";
	print IN cut_up_fasta_line($$seq{$code}, '100');
    $code_corr{'seq_'.$seq_count} = $code;
    ++$seq_count;
  }
  close 'IN';

  system '"'.$clustal_ex.'" -INFILE='.$file.' -QUIET -OUTPUT=GDE >'.$log;
  open_file_in($out_cl, 'AL');
  my $code = '';
  my %seqa = ();
  while (my $line = <AL>)
  {
    chomp $line;
    if ($line =~ /\#/)
    {
      $code = $line;
      $code =~ s/\#//;
      $code = $code_corr{$code};
    }
    else
    {
      $seqa{$code} .= $line;
    }
  }
  close AL;
  my $cons = make_cons(\%seqa, $prop_maj);
  
  --$seq_count;
  print CONSUB ">cons_gr$gr_count";
  print CONSUB "_$seq_count\n$cons\n";
  foreach my $code (sort keys %seqa)
  {
    print CONSUB ">$code orientation:$$gr_seq{$gr}{$code}\n$seqa{$code}\n";
    delete $$seq{$code};
    ++$$seq_count_hash{'Number of sequences pooled into consensus sequences:'};
  }

}
close CONSUB;

}

#######################################
sub make_cons
{
my ($seq_align, $prop_maj) = @_;
my @codes = sort keys %$seq_align;
my $cons = '';
my %seq_temp = %$seq_align;

# replace end gaps by '*';
foreach my $code (@codes)
{
  my $l = '';
  my $temp = '';
  if ($seq_temp{$code} =~  s/^-+//)
  {
    $l = length $&;
    $temp = '*' x $l;
    $seq_temp{$code} = $temp.$seq_temp{$code};
  }
  if ($seq_temp{$code} =~  s/-+$//)
  {
    $l = length $&;
    $temp = '*' x $l;
    $seq_temp{$code} .= $temp;
  }
}


for(my $i= 0; $i < length($seq_temp{$codes[0]}); ++$i)
{
  my %bases = ('A' => 0,
               'T' => 0,
                'C' => 0,
                'G' => 0,
                '-' => 0,
                '*' => 0);
                
  foreach my $code (@codes)
  {
    my $base = substr($seq_temp{$code}, $i, 1);
    ++$bases{uc$base};
  }

  my $n = scalar@codes - $bases{'*'};
  delete $bases{'*'};
  my @count = sort values %bases;
  if ($count[-1]/$n >= $prop_maj)
  {
      foreach my $b (keys %bases)
      {
        if ($bases{$b} == $count[-1])
        {
          $cons .= $b;
          last;
        }
      }
  }
  else 
  {$cons .= 'N';}
}
return $cons;
}

#################################################

sub treat_group_seq
{
my ($seq, $blast, $gr_file, $gr_file_temp, $seq_count) = @_;
  # identify regions of %gr_sequences that have hit to other sequences from %blast
  # print %gr_list into class file with regions of hit prined in lc; delete gr_list from %seq
  # make gr_list.fas
  
  # mask hit regions in all grouped sequences; delete %blast lines with grouped sequences
  my %seq_mask = ();
  foreach my $code1 (keys %$blast)
  {
    foreach my $code2 (keys %{$$blast{$code1}})
    {
      my $bool = 0;
      my $seq_temp = '';
      if (exists $$seq{$code1})
      {
        if (exists $seq_mask{$code1})
        {$seq_temp = $seq_mask{$code1};}
        else
        {$seq_temp = $$seq{$code1};}
        
        $seq_mask{$code1} = mask_interval($seq_temp, $$blast{$code1}{$code2}[1], $$blast{$code1}{$code2}[2]);
        $bool = 1;
      }
      if (exists $$seq{$code2})
      {
        if (exists $seq_mask{$code2})
        {$seq_temp = $seq_mask{$code2};}
        else
        {$seq_temp = $$seq{$code2};}
 
         $seq_mask{$code2} = mask_interval($seq_temp, $$blast{$code1}{$code2}[3], $$blast{$code1}{$code2}[4]);   
         $bool = 1;    
      }
      if ($bool)
      {
        delete $$blast{$code1}{$code2};
      }
    }
  }
  
    # print %gr_list into class file with regions of hit prined in lc; delete gr_list from %seq
     # make gr_list.fas

  open_file_out($gr_file, 'GR');
  open_file_out($gr_file_temp, 'GRT');
  $$seq_count{'Number of grouped sequenes (without consensus):'} = 0;
  foreach my $code (sort keys %$seq)
  {
    print GRT ">$code\n";
	print GRT cut_up_fasta_line($$seq{$code}, '100');
    ++$$seq_count{'Number of grouped sequenes (without consensus):'};
    if (exists $seq_mask{$code})
    {
      print GR ">$code\n";
	  print GR cut_up_fasta_line($seq_mask{$code}, '100');
    }
    else # the only blast hit of the sequence was to nohit_css sequence and thus elimnated earlier 
    {
      print GR ">$code\n";
	  print GR cut_up_fasta_line($$seq{$code}, '100');    
    }

    delete $$seq{$code};
  }
  close GR;
  close GRT;
  
}
 
######################################################

sub mask_interval
{
my ($seq, $beg, $end) = @_;
if ($beg > $end)
{
  my $temp = $beg;
  $beg = $end;
  $end = $temp;
}
#print "$beg, $end\n $seq\n";
my $frag = lc substr($seq, $beg-1, $end-$beg +1);
substr($seq, $beg-1, $end-$beg +1, $frag);
#print "$seq\n";

return $seq;
}

######################################################

sub identify_grouped_cons_seq
{
my ($blast_out_cons, $cons_subs_file, $cons_file, $gr_file, $seq_count) = @_;

open_file_in($blast_out_cons, 'BL');
my %blast = ();

# identify hit region for each grouped consensus sequences
while (my $line = <BL>)
{
  my @line = split(' ', $line);

  unless(exists $blast{$line[0]})
  {
    $blast{$line[0]} = [$line[6], $line[7]];
    next;      
  }
  
    if ($blast{$line[0]}[0] > $line[6])
    {
      $blast{$line[0]}[0] = $line[6];
    }
    if ($blast{$line[0]}[1] < $line[7])
    {
      $blast{$line[0]}[1] = $line[7]; 
    }   
}
close BL;


# delete grouped consensus sequences from cons_file
# add grouped cons seq to gr_file;
$$seq_count{'Number of unique consensus sequences:'} = 0;
$$seq_count{'Number of grouped consensus sequences:'} = 0;
open_file_in($cons_file, 'CONS');
my %cons_seq = (); # consensus sequences
my %cons_line = (); # keep the definition line with polymorphism info
my $code = '';
while (my $line = <CONS>)
{
  chomp $line;
  if ($line =~ />/)
  {
	$code = $line;
	$code =~ s/>//;
	$code =~ s/\s.*//;
	$cons_line{$code} = $line;
  }
  else
  {
	$cons_seq{$code} .= $line;
  }
}
close CONS;
open_file_out($cons_file, 'COOUT');
unless(open('GR', ">>$gr_file"))
{print "Cannot open $gr_file";}

foreach my $code (sort keys %cons_seq)
{
    if (exists $blast{$code})
    {
      my $fr = substr($cons_seq{$code}, $blast{$code}[0]-1, $blast{$code}[1]-$blast{$code}[0] +1);
      substr($cons_seq{$code}, $blast{$code}[0]-1, $blast{$code}[1]-$blast{$code}[0] +1, lc$fr);
      print GR "$cons_line{$code}\n";
	  print GR cut_up_fasta_line($cons_seq{$code}, '100');
      ++$$seq_count{'Number of grouped consensus sequences:'};
    }
    else
    {
      print COOUT "$cons_line{$code}\n";
	  print COOUT cut_up_fasta_line($cons_seq{$code}, '100');
      ++$$seq_count{'Number of unique consensus sequences:'};
    }
  
}
close COOUT;
close GR;

if (0)
{
open_file_in($cons_subs_file, 'COS');
my $cons_code= '';
my $seq_code = '';

while (my $line = <COS>)
{
  chomp $line;
  if($line =~ />cons/)
  {
    $cons_code = $line;
    $cons_code =~ s/>//;
    $seq_code = '';
  }
  elsif ($line =~ />/)
  {
    $seq_code = $line;
    $seq_code =~ s/>//; 
  }
  else
  {
    if (exists $blast{$cons_code})
    {
      if ($seq_code eq '')
      {
        $line =~ s/-//g;
        my $fr = substr($line, $blast{$cons_code}[0]-1, $blast{$cons_code}[1]-$blast{$cons_code}[0] +1);
        substr($line, $blast{$cons_code}[0]-1, $blast{$cons_code}[1]-$blast{$cons_code}[0] +1, lc$fr);
        print GR ">$cons_code\n";
		print GR cut_up_fasta_line($line, '100');
      }
      else
      {
        $line = uc $line;
        print GR ">$cons_code";
        print GR "_$seq_code\n";
		print GR cut_up_fasta_line($line, '100');
      }
    }
  }
}
close COS;
close GR;
}

}
######################################################

sub check_polymorphism
{
my ($cons_poly, $ms_file, $masked_file) = @_;

my %seq = read_hach_from_fasta_code_wo_space($masked_file);
open_file_in($ms_file, 'IN');
my $title = <IN>;
my %cons = (); # hach of hach of tables $cons{$cons_code}{$ms_code}[ms_mot, beg, end, rep_number_seq1...]
my $cons_code = '';
while (my $line = <IN>)
{
  chomp $line;
  my @line = split(';', $line);
  if($line[0] =~ /cons/) # collect info on Ms position in the aligned cons seq
  {
    for (my $i = 3; $i< scalar@line; $i = $i + 4)
    {
      $cons_code = $line[0];
      if (($line[$i+1] < length$line[$i]) or (($line[2]- $line[$i+2]) <length$line[$i])) # if ms starts at the begining or end at the end of the cons seq.
      { next; }
      $cons{$line[0]}{($i+1)/4} = [$line[$i], $line[$i+1], $line[$i+2]];
    }
  }
  else
  {
     for (my $i = 3; $i< scalar@line; $i = $i + 4)
    {
      my $gap_beg = 0;
      my $gap_end = 0;
      if ($seq{$line[0]} =~ /^-+/)
	{$gap_beg = length $&}
      if ($seq{$line[0]} =~ /-+$/)
       {$gap_end = length $&;}
      if ((($line[$i+1]-$gap_beg) < length$line[$i]) or $line[2]-$line[$i+2]-$gap_end < length$line[$i]) # if ms is within one motif length of the begining or the end of the seqeunce (after eliminating end gaps)
      {next;}
      check_overlap($line[$i+1], $line[$i+2], $line[$i+3], $cons{$cons_code}); # if ms overlaps with one of the MS in cons seq, add repeat number to $cons{$cons_code}{$ms_code}[ms_mot, beg, end, rep_number_seq1...]
    }   
  }
}
close IN;

foreach my $code (keys %cons)
{
  foreach my $ms (keys %{$cons{$code}})
  {
    my $polymorph = 0;
    for (my $j = 4; $j < scalar@{$cons{$code}{$ms}}; ++$j)
    {
      if ($cons{$code}{$ms}[3] != $cons{$code}{$ms}[$j])
      {
        $polymorph = 1;
        last;
      }
    }
    if ($polymorph)
    {
      # adjust MS beg and Ms end positions in the consensus sequence WO gaps (tha actual values are with gaps)
      my $adjusted_pos = ''; # msmotif_beg_end
      adjust_ms_pos(\$adjusted_pos, $cons{$code}{$ms}[0], $cons{$code}{$ms}[1], $cons{$code}{$ms}[2], $seq{$code});
#      $$cons_poly{$code} .= $cons{$code}{$ms}[0].'_'.$cons{$code}{$ms}[1].'_'.$cons{$code}{$ms}[2].';';
      $$cons_poly{$code} .= $adjusted_pos.';';
    }
  }
}

%cons = ();
%seq = ();
#print Dumper %cons;

}
######################################################

sub adjust_ms_pos
{
 my ($adjusted_pos_ref, $mot, $ms_beg, $ms_end, $seq_align) =@_;
 
 my $beg = substr($seq_align, 0, $ms_beg-1);
 $beg =~ s/-//g;
 $ms_beg = (length $beg) +1;
 
 my $end = substr($seq_align, 0, $ms_end-1);
 $end =~ s/-//g;
 $ms_end = (length $end) +1;
 $$adjusted_pos_ref = $mot.'_'.$ms_beg.'_'.$ms_end;
 
}
#########################################################
sub check_overlap
{
my ($beg, $end, $rep_n, $cons_ref) =@_;
foreach my $ms (keys %$cons_ref)
{
  if (($$cons_ref{$ms}[1] <= $beg and $beg <= $$cons_ref{$ms}[2]) or ($$cons_ref{$ms}[1] <= $end and $end <= $$cons_ref{$ms}[2]) or ($$cons_ref{$ms}[1] >= $beg and $end >= $$cons_ref{$ms}[2]))
  {
    push(@{$$cons_ref{$ms}}, $rep_n);
  }
}
}

######################################################

sub write_consensus_file_with_polymorphism_info
{
my ($cons_poly, $consensus_file, $cons_subs_file, $masked_file) = @_;
open_file_out($consensus_file, 'CON');
open_file_out($cons_subs_file, 'SUB');
open_file_in($masked_file, 'MASK');
my $code = '';
my %cons_local = ();
while(my $line = <MASK>)
{
  print SUB "$line";
  if ($line =~ />cons/)
  {
    $code = $line;
    chomp $code;
    $code =~ s/>//;
  }
  elsif($line =~ />/)
  {$code = '';}
  elsif($code)
  {
	chomp $line;
	$cons_local{$code} .= $line;
  }
}

foreach my $code (sort keys %cons_local)
{
	print CON ">$code";
	if (exists $$cons_poly{$code})
	{
		print CON " $$cons_poly{$code}";
	}
	print CON "\n";
    
	$cons_local{$code} =~ s/-//g;
	print CON cut_up_fasta_line($cons_local{$code}, '100');		
}
close MASK;
close CON;
close SUB;
}

#####################################################

sub make_singleton_only
{
my ($seq, $singleton_file, $nohit_css_file, $blast_out, $seq_count, $multihit_css_hash) =@_;
#make file with sequenes that had only an autohit, but nothing else
#make file with sequenes that had did not had autohit
# delete all these sequences from %$seq

my %hit = (); #codes of all sequences that had non-autohit
my %autohit = (); # codes of all sequences that had autohit
open_file_in($blast_out, 'BL');
while(my $line = <BL>)
{
  my @line = split(' ', $line);
  if($line[0] eq $line[1])
  {
	if(exists $autohit{$line[0]})
	{
		++$autohit{$line[0]};
	}
	else
	{
    		$autohit{$line[0]} = 1;
	}
  }
  else
  {
    $hit{$line[0]} = '';
    $hit{$line[1]} = '';
  }
}
close BL;

open_file_out($singleton_file, 'OUT');
open_file_out($nohit_css_file, 'NOHIT');

$$seq_count{'Number of singleton sequences:'} = 0;
$$seq_count{'Number of sequences without autohit (nohit_css):'} = 0;
foreach my $code (sort keys %$seq)
{
  if (exists $autohit{$code})
  {
	if($autohit{$code} > 1)
	{
		$$multihit_css_hash{$code} = '';
	}
	else
	{
    		unless(exists $hit{$code})
    		{
     		 print OUT ">$code\n";
			 print OUT cut_up_fasta_line($$seq{$code}, '100');
     		 delete $$seq{$code};
     		 ++$$seq_count{'Number of singleton sequences:'};
    		}
  	}
   }
  else
  {
   	 print NOHIT ">$code\n";
	 print NOHIT cut_up_fasta_line($$seq{$code}, '100');  
   	 delete $$seq{$code};
 	 ++$$seq_count{'Number of sequences without autohit (nohit_css):'};
  }
}
close OUT;
close NOHIT;

%hit = ();
%autohit = ();

}

########################################################

sub treat_multihit_css
{
my ($blast, $css, $seq, $multihit_css_file, $seq_count) = @_;
# delete_multihit_css from the %blast hash and from %seq

foreach my $code1 (keys %$blast)
{
  if (exists $$css{$code1}) # if code1 css delete all blast hits with code1
  {
    delete$$blast{$code1};
  }
  else
  {
    foreach my $code2 (keys %{$$blast{$code1}})
    {
      if (exists $$css{$code2}) # if code2 css but not code1 delete specific code2 blast hits
      {
        delete$$blast{$code1}{$code2};   
      }
    }
    if (scalar(keys %{$$blast{$code1}}) == 0) # if all hits to code1 are deleted
    {
      delete$$blast{$code1};
    }
  }
}


# delete_multihit_css from %seq
# print multihit_css_file;
open_file_out($multihit_css_file, 'MULT');
$$seq_count{'Number of multihit css sequences:'} = 0;
foreach my $code (keys %$seq)
{
  if (exists $$css{$code})
  {
    print MULT ">$code\n";
	print MULT cut_up_fasta_line($$seq{$code}, '100');
    delete $$seq{$code};
    ++$$seq_count{'Number of multihit css sequences:'};
  } 
}

close MULT;

}

########################################################

sub define_param_number_pipe2
{
  my %param_list = (
      1=> ['syst', 'Operating system (win/linux):'], 
      2=> ['project_folder', 'Project folder:'], 
      3=> ['make_cons', 'Make consensus sequences (YES=1/NO=0):'],
      4=> ['ident_limit', 'Minimum % of identity between sequences of a contig (80-100):'], 
      5=> ['prop_maj', 'Proportion of sequences that must have the same base at a site to accept it as a consensus (0.5-1):' ], 
      6=> ['clust_path', 'Pathway to CLUSTALW2 executables:'], 
      7=> ['blast_path', 'Pathway to BLAST executables:'], 
#      8=> ['blast_version', 'BLAST version (blast/blast+):'],
      8=> ['del_files', 'Delete intermediate files (YES=1/NO=0):']
      );
  return %param_list;
}

##################################################################

sub define_param_limits_pipe2
{
my ($param) = @_;
#hash of array [type(list, real, folder), error message, ['list of values/lower_upper limits)']
  my %param_limits = (
      'galaxy'=> ['list', "galaxy must be 0 or 1\n", [0, 1]], 
      'syst'=> ['list', "Operating system must be win or linux\n", ['win', 'linux']], 
      'blast_path'=> ['folder', "Cannot open blast folder ($$param{blast_path} folder)\n"], 
	  'num_threads' => ['real', "The number of threads ($$param{num_threads}) must be between 1 and 100000\n", [0,100000]],
      'make_cons'=> ['list', "make_cons must be 1 or 0\n", [0,1]],
      'ident_limit'=> ['real', "ident_limit must be between 80-100\n", [80,100]], 
      'prop_maj'=> ['real', "prop_maj must be between 0.5-1\n", [0.5,1]], 
      'input_file' => ['file', "Cannot open input_file ($$param{input_file}\n"],
      'del_files'=> ['list', "del_files must be 0 or 1\n", [0,1]],
      'out_folder' => ['folder', "Cannot open out_folder ($$param{out_folder})\n"], 
      'qdd_folder' => ['folder', "Cannot open qdd_folder ($$param{qdd_folder})\n"],
      'clust_path'=> ['folder', "Cannot open clustal folder ($$param{clust_path})\n"], 
      'outfile_string' => ['string', "-outfile_string can contain only alphanumerical values, underscore '_' or dot\n"],
      'debug' => ['list', "debug must be 0 or 1\n", [0,1]],
#      'file_msfas' => ['string', "file_msfas ($$param{file_msfas}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_singleton' => ['string', "file_singleton ($$param{file_singleton}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_nohit_css' => ['string', "file_nohit_css ($$param{file_nohit_css}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_multihit_css' => ['string', "file_multihit_css ($$param{file_multihit_css}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_gr' => ['string', "file_gr ($$param{file_gr}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_consensus' => ['string', "file_consensus ($$param{file_consensus}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_cons_subs' => ['string', "file_cons_subs ($$param{file_cons_subs}) can contain only alphanumerical values, underscore '_' or dot\n"],
	'file_pipe3_inp' => ['string', "file_pipe3_inp ($$param{file_pipe3_inp}) can contain only alphanumerical values, underscore '_' or dot\n"]
      );
return %param_limits;
}

