#!/usr/bin/perl -w


  use warnings;
  use strict;
  use subprogramQDD;
  use Data::Dumper;


# Run all three parts of QDD + option of sorting sequences by tag before the run.
# Runs QDD on all sequences in the input folder  
# Pipe1=> sequence checks and microsatellite detection
# Pipe2=> redundancy elimination
# Pipe3=> primer design
# Pipe4=> Contamination check and RepeatMasker run
  
#INPUT: 
# -Input folder name with full path. It must contain all the inut fasta files (with the sequences of interest) and nothing else. Inut fasta file must have relatively short reads (less than ca.2000 bp)
# -Fasta file with adapters (if vector/adptor clipping is necessary); should not in the input folder
# -Fasta file with tags (if sorting by tags is necessary); should not in the input folder


#OUTPUT: 
# PIPE1
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

#PIPE2
# file_root_singleton.fas => the only BLASt hit is autohit; All >=4 repeat microsatellites and >=8 homoplymers are masked
# file_root_nohit_css.fas => no hit to itself; All >=4 repeat microsatellites and >=8 homoplymers are masked
# file_root_multihit_css.fas => more than one hit between a pair of sequences; All >=4 repeat microsatellites and >=8 homoplymers are masked
# file_root_gr.fas => sequences (including consensuses) that had BLAST hit to other sequences, with bellow limit identity of the overlaping region. Regions covered by BLAST hits are masked 
# file_root_consensus.fas => all unique (no hit to grouped seqs) consensus sequences, if polymorphism of ms motif_first_last pos of MS (in the aligned cons seq) is given; All >=4 repeat microsatellites are masked, but not homopolymers;
# file_root_cons_subs.fas => Consensus + aligned sequences to make a consensus ; All >=4 repeat microsatellites are masked, but not homopolymers; grouped consensus sequences are NOT removed
# file_root_pipe2_for_pipe3.fas => fatasta file with all unique sequences (singletons and consensus); input for pipe3

#PIPE3
#xxx_pipe3_primers.csv => semicolon separated values, info on primers, target region parameters, primer design etc.
#xxx_pipe3_targets.fas => sequences, where primer design was successful.

#PIPE4
#xxx_pipe4_primers.csv => xxx_pipe3_primers.csv completed with information on BLAST hits againts NCBI nt and RepeatMasker hits.

my $set_file_default = '/etc/qdd/set_qdd_default.ini';
my $blast_version = 'blast+';
my $e_tag = '1E-1';
my $word_size = 7; 

$set_file_default = check_set_file($set_file_default);

print_heading_3_0('STDOUT');

my %param = (
  galaxy => '',
  tag => '',
  tag_file => '',
  run_all => '', #[0/1] 1 for running pipe1 2 3 at once
  input_folder => '',
  out_folder => '',
  syst => '',  
  debug => 0,
  del_files => '',
  qdd_folder => '',  
  blast_path => '',
  clust_path => '',
  primer3_path => '',
  primer3_version => '', 
  num_threads => '',
  qual => 0,


# pipe1
  fastq => '',
  adapter => '',
  adapter_file => '',
  length_limit => '',
  mean_qs => '',

  file_log_pipe1 => '_pipe1_log.txt',
  file_wov => '_pipe1_wov.fas',
  file_length => '_pipe1_length_info.tabular',
  file_msfas => '_pipe1_for_pipe2.fas',
  file_ms => '_pipe1_ms.csv',
  file_qual_trim => '_pipe1_trimmed.qual',

  mean_qs => 20,
  contig => 1,
  flank_length => 200,


# pipe2
  make_cons => '', # if 0 make only a file with unique autohit sequences
  ident_limit => '',
  prop_maj => '',

  file_log_pipe2 => '_pipe2_log.txt',
  file_singleton => '_pipe2_singleton.fas',
  file_nohit_css => '_pipe2_nohit_css.fas',
  file_multihit_css => '_pipe2_multihit_css.fas',
  file_gr => '_pipe2_groupped.fas',
  file_consensus => '_pipe2_consensus.fas',
  file_cons_subs => '_pipe2_cons_subs.fas',
  file_pipe3_inp => '_pipe2_for_pipe3.fas',

# pipe3
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
  min_qs => '',

  file_log_pipe3 => '_pipe3_log_file.txt',  
  file_primers => '_pipe3_primers.tabular',
  file_targets => '_pipe3_targets.fas',

# pipe4
  rm => 1,
  rm_path => '',
  rm_lib => '',
  blastdb => '',
  db_evalue => '1e-20',
  local_blast => '',
  check_contamination => '',
  file_primers_completed => '_pipe4_primers.tabular'
);

read_set_file_galaxy(\%param, $param{qdd_folder}.$set_file_default);
modify_params_from_tags(\%param, \@ARGV);
#print Dumper(\%param);
unless($param{galaxy} eq '0')
{
	print "QDD.pl cannot be run from Galaxy.\n-galaxy must be set to 0, while using QDD.pl\n";
	exit;
}
my %param_limits_pipe123 = define_param_limits_pipe123(\%param);
my %param_limits_sort = define_param_limits_sort(\%param);
add_slash_to_folder(\%param, \%param_limits_pipe123);
add_slash_to_folder(\%param, \%param_limits_sort);

my $folder_id = time;
my $date = localtime $folder_id;
my $tempfolder = $param{out_folder}.$folder_id.'/';

make_folder($param{syst}, $tempfolder);

my $log_file = $param{out_folder}.'QDD_log_'.$folder_id.'.txt';
open(my $fh_log, '>>', $log_file) or die "Cannot open $log_file $!\n";
print_heading_3_0($fh_log);
print $fh_log "Running QDD.pl\n\n";

my %file_list_pipe = ();

# print out an error message if some of the fasta files do not have a corresponding qual file (xxx.fas=> xxx.qual) and stops the program
my %file_list_tag = get_file_list($param{input_folder}); # make hash with files (including path); does not include directories
if($param{qual})
{
	check_fasta_qual_pairs(\%file_list_tag);
}

# sort sequences according to tags
if($param{tag}) 
{
	check_param_limits(\%param, \%param_limits_sort);


	my $sorted_files_folder = $param{out_folder}.'sorted_files_'.$folder_id.'/';
	make_folder($param{syst}, $sorted_files_folder);
		
	foreach my $file (sort keys %file_list_tag)
	{
	unless($file =~ /\.qual$/)
	{
		my $qual_file_sort = $file;
		$qual_file_sort =~ s/\.[^\.]*$/\.qual/;
		
		if($param{fastq})
		{
			print "Converting fastq file ($file) into fasta format\n";
			my $fasta = $file;
			$fasta =~ s/.*[\\\/]//;
			$fasta =~ s/\..*/.fas/;
			$fasta = $tempfolder.$fasta;
			convert_fasq_into_fasta($file, $fasta);
			$file = $fasta;	
		}
	
		
		my $sort_tag = 'perl "'.$param{qdd_folder}.'sort_tag_galaxy.pl" "'.$param{syst}.'" "'.$file.'" "'.$sorted_files_folder.'" "'.$param{blast_path}.'" "'.$blast_version.'" "'.$e_tag.'" "'.$word_size.'" "'.$param{tag_file}.'" "'.$tempfolder.'" "'.$param{num_threads}.'" "'.$qual_file_sort.'" "'.$param{qual}.'"';
	#$syst,	$fas_file (with path),	$outfolder, $blast_path , $blast_version $Evector, $word_size$tag_file, $tempfolder, $num_threads;

		if($param{debug})
		{
		print $fh_log "$sort_tag\n";
		}
		print "Sorting sequences in $file by tags\n";
		system $sort_tag;
	}
	}# end foreach
	$param{fastq} = 0;
	%file_list_pipe = get_file_list($sorted_files_folder);
	print_seq_numbers_delete_notag(\%file_list_pipe, $fh_log, '_NOTAG.fas');
}# end if tag
else # get file list 
{
	%file_list_pipe = get_file_list($param{input_folder});
}

# run pipe1,2,3
if ($param{run_all})
{
	check_param_limits(\%param, \%param_limits_pipe123);
	if($param{rm})
	{
		print "Checking the validity of $param{rm_lib} for selecting elements from the RepeatMasker database\n\n";
		check_repeatmasker_db($param{rm_lib}, $param{rm_path}, $tempfolder, $param{syst});
	}
	foreach my $file (sort keys %file_list_pipe)
	{
	unless($file=~ /\.qual/)
	{
		my $qual_file = $file;
		$qual_file =~ s/\.[^\.]*$/\.qual/;
		my $pipe1 = 'perl "'.$param{qdd_folder}.'pipe1.pl" -input_file "'.$file.
				'" -galaxy "'.$param{galaxy}.
				'" -fastq "'.$param{fastq}.
				'" -syst "'.$param{syst}.
				'" -blast_path "'.$param{blast_path}.
				'" -out_folder "'.$param{out_folder}.
				'" -adapter "'.$param{adapter}.
				'" -length_limit "'.$param{length_limit}.
				'" -del_files "'.$param{del_files}.
				'" -adapter_file "'.$param{adapter_file}.
				'" -qdd_folder "'.$param{qdd_folder}.
				'" -file_wov "'.$param{file_wov}.
				'" -file_length "'.$param{file_length}.
				'" -file_msfas "'.$param{file_msfas}.
				'" -file_ms "'.$param{file_ms}.
				'" -debug "'.$param{debug}.
				'" -file_log_pipe1 "'.$log_file.
				'" -num_threads "'.$param{num_threads}.
				'" -qual "'.$param{qual}.
				'" -qual_file "'.$qual_file.
				'" -file_qual_trim "'.$param{file_qual_trim}.
				'" -mean_qs "'.$param{mean_qs}.
				'" -contig "'.$param{contig}.
				'" -flank_length "'.$param{flank_length}.'"';

		if($param{debug})
		{print $fh_log "PIPE1\n$pipe1\n";}
		system $pipe1;

		$param{file_msfas} = get_last_version($param{out_folder}, $param{file_msfas});
		my $filename_root = get_filename_root($file);
		my $pipe2_inp = $param{out_folder}.$filename_root.$param{file_msfas};
		my $pipe2 = 'perl "'.$param{qdd_folder}.'pipe2.pl" -input_file "'.$pipe2_inp.
			'" -galaxy "'.$param{galaxy}.
			'" -syst "'.$param{syst}.
			'" -blast_path "'.$param{blast_path}.
			'" -out_folder "'.$param{out_folder}.
			'" -make_cons "'.$param{make_cons}.
			'" -ident_limit "'.$param{ident_limit}.
			'" -clust_path "'.$param{clust_path}.
			'" -prop_maj "'.$param{prop_maj}.
			'" -del_files "'.$param{del_files}.
			'" -qdd_folder "'.$param{qdd_folder}.
#			'" -file_msfas "'.$param{file_msfas}.
			'" -file_singleton "'.$param{file_singleton}.
			'" -file_nohit_css "'.$param{file_nohit_css}.
			'" -file_multihit_css "'.$param{file_multihit_css}.
			'" -file_gr "'.$param{file_gr}.
			'" -file_consensus "'.$param{file_consensus}.
			'" -file_cons_subs "'.$param{file_cons_subs}.
			'" -file_pipe3_inp "'.$param{file_pipe3_inp}. 
			'" -debug "'.$param{debug}.'" -file_log_pipe2 "'.$log_file.'" -num_threads "'.$param{num_threads}.'"';

		if($param{debug})
		{print $fh_log "PIPE2 \n$pipe2\n";}
		system $pipe2;

		$param{file_pipe3_inp} = get_last_version($param{out_folder}, $param{file_pipe3_inp});
		$param{file_qual_trim} = get_last_version($param{out_folder}, $param{file_qual_trim});
		my $pipe3_inp = $param{out_folder}.$filename_root.$param{file_pipe3_inp};
		my $pipe3_inp_qual = $param{out_folder}.$filename_root.$param{file_qual_trim};
		my $pipe3 ='perl "'.$param{qdd_folder}.'pipe3.pl" -input_file "'.$pipe3_inp.
			'" -galaxy "'.$param{galaxy}.
			'" -syst "'.$param{syst}.
			'" -out_folder "'.$param{out_folder}.
			'" -del_files "'.$param{del_files}.
			'" -qdd_folder "'.$param{qdd_folder}.
			'" -primer3_path "'.$param{primer3_path}.
			'" -primer3_version "'.$param{primer3_version}.
			'" -file_pipe3_inp "'.$param{file_pipe3_inp}.
			'" -file_primers "'.$param{file_primers}.
			'" -file_targets "'.$param{file_targets}.
			'" -pcr_min "'.$param{pcr_min}.
			'" -pcr_max "'.$param{pcr_max}.
			'" -pcr_step "'.$param{pcr_step}.
			'" -PRIMER_GC_CLAMP "'.$param{PRIMER_GC_CLAMP}.
			'" -PRIMER_OPT_SIZE "'.$param{PRIMER_OPT_SIZE}.
			'" -PRIMER_MIN_SIZE "'.$param{PRIMER_MIN_SIZE}.
			'" -PRIMER_MAX_SIZE "'.$param{PRIMER_MAX_SIZE}.
			'" -PRIMER_OPT_TM "'.$param{PRIMER_OPT_TM}.
			'" -PRIMER_MIN_TM "'.$param{PRIMER_MIN_TM}.
			'" -PRIMER_MAX_TM "'.$param{PRIMER_MAX_TM}.
			'" -PRIMER_MAX_DIFF_TM "'.$param{PRIMER_MAX_DIFF_TM}.
			'" -PRIMER_MIN_GC "'.$param{PRIMER_MIN_GC}.
			'" -PRIMER_OPT_GC_PERCENT "'.$param{PRIMER_OPT_GC_PERCENT}.
			'" -PRIMER_MAX_GC "'.$param{PRIMER_MAX_GC}.
			'" -PRIMER_SELF_ANY "'.$param{PRIMER_SELF_ANY}.
			'" -PRIMER_SELF_END "'.$param{PRIMER_SELF_END}.
			'" -PRIMER_MAX_POLY_X "'.$param{PRIMER_MAX_POLY_X}.
			'" -PRIMER_NUM_RETURN "'.$param{PRIMER_NUM_RETURN}.
			'" -debug "'.$param{debug}.
			'" -qual_file "'.$pipe3_inp_qual.
			'" -min_qs "'.$param{min_qs}.
			'" -qual "'.$param{qual}.
			'" -file_log_pipe3 "'.$log_file.
			'" -contig "'.$param{contig}.'"';

	
  
		if($param{debug})
		{print $fh_log "PIPE3\n$pipe3\n";}
		system $pipe3;

		$param{file_primers} = get_last_version($param{out_folder}, $param{file_primers});
		my $pipe4_inp = $param{out_folder}.$filename_root.$param{file_primers};
		my $pipe4 ='perl "'.$param{qdd_folder}.'pipe4.pl" -input_file "'.$pipe4_inp.
			'" -galaxy "'.$param{galaxy}.
			'" -syst "'.$param{syst}.
			'" -out_folder "'.$param{out_folder}.
			'" -del_files "'.$param{del_files}.
			'" -qdd_folder "'.$param{qdd_folder}.
			'" -blast_path "'.$param{blast_path}.
			'" -blastdb "'.$param{blastdb}.
			'" -db_evalue "'.$param{db_evalue}.
			'" -num_threads "'.$param{num_threads}.
			'" -local_blast "'.$param{local_blast}.
			'" -check_contamination "'.$param{check_contamination}.
			'" -debug "'.$param{debug}.
			'" -rm "'.$param{rm}.
			'" -rm_path "'.$param{rm_path}.
			'" -rm_lib "'.$param{rm_lib}.
			'" -file_log_pipe4 "'.$log_file.
			'" -file_primers_completed "'.$param{file_primers_completed}.'"';

 
		if($param{debug})
		{print $fh_log "PIPE4\n$pipe4\n";}
		system $pipe4;

	}# end unless qual
	}
}



if ($param{del_files} ==1)
{
  delete_folder($tempfolder, $param{syst});
}

my $time_end = time;
my $time_run = $time_end - $folder_id;
print $fh_log "\nQDD.pl started at : $date\n";
print $fh_log "The run took $time_run s.\n";
print "\nSee log_file for summary:\n$log_file\n\n";
print $fh_log "QDD.pl is finished\n";
close $fh_log;
exit;

#########################################################

sub print_seq_numbers_delete_notag
{
my ($file_list_hr, $fh_log, $motif) = @_;

print $fh_log "\nThe number of sequences in each of the sorted files:\n";
foreach my $file (sort keys %$file_list_hr)
{
	my $seqn = get_seq_number_from_fasta($file);
	print $fh_log "$file:$seqn\n";
	if($file =~  /$motif/)
	{
		delete($$file_list_hr{$file});
	}

}
print $fh_log "\n";

}


###################################################
sub define_param_limits_sort
{
my ($param) = @_;
#hash of array [type(list, real, folder), error message, ['list of values/lower_upper limits)']

  my %param_limits = (
      'syst'=> ['list', "Operating system must be win or linux\n", ['win', 'linux']], 
      'blast_path'=> ['folder', "Cannot open $$param{blast_path} folder\n"], 
      'tag'=> ['list', "tag must be 0 or 1\n", [0,1]],
      'del_files'=> ['list', "del_files must be 0 or 1\n", [0,1]],
      'input_folder' => ['folder', "Cannot open input_folder ($$param{input_folder})\n"],
      'tag_file' => ['file', "Cannot open tag_file ($$param{tag_file})\n"],
      'qdd_folder' => ['folder', "Cannot open qdd_folder ($$param{qdd_folder})\n"],
      'out_folder' => ['folder', "Cannot open out_folder ($$param{out_folder})\n"],
      'debug' => ['list', "debug must be 0 or 1\n", [0,1]],
	  'num_threads' => ['real', "The number of threads ($$param{num_threads}) must be between 1 and 100000\n", [0,100000]],
	  'qual' => ['list', "qual must be 0 or 1\n", [0,1]],
      );
	  if($param{qual})
	  {
		$param_limits{file_qual_trim} = ['string', "file_qual_trim ($$param{file_qual_trim}) can contain only alphanumerical values, underscore '_' or dot\n"];
	  }


return %param_limits;
}

###################################################

sub define_param_limits_pipe123
{
my ($param) = @_;
#hash of array [type(list, real, folder), error message, ['list of values/lower_upper limits)']

  my %param_limits = (
      'run_all'=> ['list', "run_all must be 0 or 1\n", [0,1]],
      'input_folder' => ['folder', "Cannot open input_folder ($$param{input_folder})\n"],
      'out_folder' => ['folder', "Cannot open out_folder ($$param{out_folder})\n"],
      'syst'=> ['list', "Operating system must be win or linux\n", ['win', 'linux']],
	  'fastq'=> ['list', "fastq must be 0 or 1\n", [0, 1]], 
      'debug' => ['list', "debug must be 0 or 1\n", [0,1]],
      'del_files'=> ['list', "del_files must be 0 or 1\n", [0,1]],
      'qdd_folder' => ['folder', "Cannot open qdd_folder ($$param{qdd_folder})\n"],
      'blast_path'=> ['folder', "Cannot open $$param{blast_path} folder\n"], 
      'clust_path'=> ['folder', "Cannot open clustal folder ($$param{clust_path})\n"], 
      'primer3_path'=> ['folder', "Cannot open primer3 folder ($$param{primer3_path})\n"], 
      'primer3_version'=> ['list', "primer3_version must be 1 or 2\n", [1,2]],
	  'num_threads' => ['real', "The number of threads ($$param{num_threads}) must be between 1 and 100000\n", [0,100000]],
	  'local_blast' => ['list', "local_blast must be 1 (run local blast for contamination check) or 0 (run remote blast for contamination check)\n", [0,1]],
	  'check_contamination' => ['list', "check_contamination must be 1 (check contamination by blast against nt of NCBI) or 0 (skip contamination check)\n", [0,1]],
	  'qual' => ['list', "qual must be 0 or 1\n", [0,1]],
	  
      'adapter'=> ['list', "adapter must be 0 or 1\n", [0,1]],
      'length_limit'=> ['real', "Minimum sequence length must be between 1 and 1000\n", [1,1000]],
      'file_wov' => ['string', "file_wov ($$param{file_wov}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_length' => ['string', "file_length ($$param{file_length}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_msfas' => ['string', "file_msfas ($$param{file_msfas}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_ms' => ['string', "file_ms ($$param{file_ms}) can contain only alphanumerical values, underscore '_' or dot\n"],
	'contig' => ['list', "contig must be 0 or 1\n", [0,1]],

      'make_cons'=> ['list', "make_cons must be 1 or 0\n", [0,1]],
      'ident_limit'=> ['real', "ident_limit must be between 80-100\n", [80,100]], 
      'prop_maj'=> ['real', "prop_maj must be between 0.5-1\n", [0.5,1]], 
      'file_msfas' => ['string', "file_msfas ($$param{file_msfas}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_singleton' => ['string', "file_singleton ($$param{file_singleton}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_nohit_css' => ['string', "file_nohit_css ($$param{file_nohit_css}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_multihit_css' => ['string', "file_multihit_css ($$param{file_multihit_css}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_gr' => ['string', "file_gr ($$param{file_gr}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_consensus' => ['string', "file_consensus ($$param{file_consensus}) can contain only alphanumerical values, underscore '_' or dot\n"],      
	'file_cons_subs' => ['string', "file_cons_subs ($$param{file_cons_subs}) can contain only alphanumerical values, underscore '_' or dot\n"],
	'file_pipe3_inp' => ['string', "file_pipe3_inp ($$param{file_pipe3_inp}) can contain only alphanumerical values, underscore '_' or dot\n"],

      'rm' => ['list', "rm must be 0 or 1\n", [0,1]],
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
      'file_primers' => ['string', "file_primers ($$param{file_primers}) can contain only alphanumerical values, underscore '_' or dot\n"],
      'file_targets' => ['string', "file_targets ($$param{file_targets}) can contain only alphanumerical values, underscore '_' or dot\n"]

      );
	  
	  if ($$param{check_contamination} and $$param{local_blast})
	  {
		$param_limits{blast_path} = ['folder', "Cannot open blast folder ($$param{blast_path} folder)\n"];
#		$param_limits{blastdb} = ['file', "BLAST db file ($$param{blastdb}) does not exists\n"];
		$param_limits{db_evalue} = ['real', "db_evalue ($$param{db_evalue}) must be between 0 and 1e-5\n", [0,1e-5]];
	  }
	  
	  if($param{qual})
	  {
		$param_limits{mean_qs} = ['real', "mean_qs must be between 0 and 40", [0,40]];
		$param_limits{file_qual_trim} = ['string', "file_qual_trim ($$param{file_qual_trim}) can contain only alphanumerical values, underscore '_' or dot\n"];
	  	$param_limits{min_qs} = ['real', "min_qs must be between 0 and 40", [0,40]];
	  }
	  if($param{adapter})
	  {
	    $param_limits{adapter_file} = ['file', "Cannot open adapter_file ($$param{adapter_file})\n"];
	  }
	  if($param{contig})
	  {
		$param_limits{flank_length} = ['real', "The length of the flanking region ($$param{flank_length}) must be between 10 and 300\n", [0,300]];
	  }

	  if($param{rm})
	  {
		$param_limits{rm_path} = ['folder', "Cannot open RepeatMasker folder ($$param{rm_path})\n"],;
	  }



return %param_limits;
}

######################################################
	# print out an error message if some of the fasta files do not have a corresponding qual file (xxx.fas=> xxx.qual) and stops the program
sub check_fasta_qual_pairs
{
 my ($file_hr) = @_;
		
	foreach my $file (sort keys %$file_hr)
	{
	unless($file =~ /\.qual$/)
	{
		my $qual_file = $file;
		$qual_file =~ s/\.[^\.]*$/\.qual/;
		unless(-e $qual_file)
		{
			print "The quality file '$qual_file' is missing from the input folder\n
Since you have selected the check quality option (-qual 1),
all input fasta file must have a corresponding quality file. (xxx.fas => xxx.qual)\n";
			exit;
		}
	}
	}
}

######################################################

sub get_last_version
{
 my ($folder, $motif) = @_;
my %files = get_file_folder_list($folder);
my $version = 0;

my $motif_temp = $motif;
my $end = '';
if ($motif_temp =~ /(\..*)/)
{
	$end = $1;
	$motif_temp =~  s/$end//;
}
foreach my $file (keys %files)
{
	if ($file =~  /$motif_temp[_]v([0-9]*)/)
	{
		my $v = $1;
		if ($v > $version)
		{
			$version = $v;
		}
	}
}

if($version)
{
	$motif = $motif_temp.'_v'.$version.$end;
}
return $motif;
}
