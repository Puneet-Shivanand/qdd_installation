#!/usr/bin/perl -w


  use warnings;
  use strict;
  use subprogramQDD;
  use Data::Dumper;
  
  
#INPUT: 
#xxx_pipe3_primers.tabular => tab separated values, info on primers, target region parameters, primer design etc.

# BLAST sequences againts nt database of ncbi (remote or local blast) and completes primer table with info on the best hit (lowest e-value)
# runs RepeatMasker on sequences and completes primer table on the best hit (highest score)

# output: 
#xxx_pipe3_primers_rm_nt.tabular => semicolon separated values, info on primers, target region parameters, primer design etc. with rm and/or nt hit info



my $set_file_default = '/etc/qdd/set_qdd_default.ini';
$set_file_default = check_set_file($set_file_default);

print_heading_3_0('STDOUT');
print "SCRIPT pipe4.pl\n";

my %param = (
  galaxy => '', # one if run on galaxy, 0 for command line
  syst => '',  
  input_file => '',
  out_folder => '',
  qdd_folder => '',
  blast_path => '',
  blastdb => '',
  db_evalue => '1e-20',
  num_threads => '',
  local_blast => '',
  check_contamination => '',
  outfile_string => '',
  rm => '',
  rm_path => '',
  rm_lib => '',
  del_files => '',
  debug => 0,
  file_log_pipe4 => '_pipe4_log.txt',
  file_primers_completed => '_pipe4_primers.tabular'
);

read_set_file_galaxy(\%param, $param{qdd_folder}.$set_file_default);
modify_params_from_tags(\%param, \@ARGV);
my %param_limits = define_param_limits_pipe4(\%param); #hash of array $param_limits{syst}[type(list, real, folder), error message, ['list of values/lower_upper limits)']
add_slash_to_folder(\%param, \%param_limits);

my $galaxy = $param{galaxy}; # one if run on galaxy, 0 for command line
my $folder_id = time;
my $date = localtime $folder_id;
my $tempfolder = $param{out_folder}.'pipe4_'.$folder_id.'/';
make_folder($param{syst}, $tempfolder);

if($param{rm} and $param{file_log_pipe4} =~ /pipe4_log.txt/)
{
print "Checking the validity of $param{rm_lib} for selecting elements from the RepeatMasker database\n\n";
check_repeatmasker_db($param{rm_lib}, $param{rm_path}, $tempfolder, $param{syst});
}

if($galaxy)
{
	if ($param{check_contamination} eq 'true' or $param{check_contamination} eq 'TRUE')
	{ 	
		$param{check_contamination} = 1;
		if ($param{local_blast} eq 'true' or $param{check_contamination} eq 'TRUE')
		{$param{local_blast} = 1}
		else
		{$param{local_blast} = 0}	
	}
	else
	{	 $param{check_contamination} = 0;	}
	
	$param{outfile_string} = 'NA';
}
elsif ($param{outfile_string} eq '')
{
	$param{outfile_string} = get_filename_root($param{input_file});
	$param{outfile_string} =~ s/_pipe3_primers.*//;
}



unless($galaxy)
{
check_param_limits(\%param, \%param_limits);
	if ($param{outfile_string} eq '')
	{
		$param{outfile_string} = get_filename_root($param{input_file});
		$param{outfile_string} =~ s/_pipe3_primers.*//;
	}


	my @clefs = ('file_primers_completed');
	if($param{file_log_pipe4} =~ /^_/) # if the name of the log file has not been specified (run in pipeline)
	{
		push(@clefs, 'file_log_pipe4');
	}
	get_last_version_and_modify_1($param{out_folder}, \@clefs, $param{outfile_string}, \%param);

	if($param{file_log_pipe4} =~ /^_/) # if the name of the log file has not been specified (run in pipeline)
	{
		$param{file_log_pipe4} = $param{out_folder}.$param{outfile_string}.$param{file_log_pipe4};
	}

	$param{file_primers_completed} = $param{out_folder}.$param{outfile_string}.$param{file_primers_completed};
}

open(my $fh_log, '>>', $param{file_log_pipe4}) or die "Cannot open $param{file_log_pipe4} $!\n";
print_heading_3_0($fh_log);
print $fh_log "SCRIPT pipe4.pl\n";

if($param{debug})
{
	print $fh_log Dumper(\%param);
	print $fh_log "\n";
}

print $fh_log "INPUT file : $param{input_file}\n";
print "Input file : $param{input_file}\n\n";
print $fh_log "\nOUTPUT files:\n";
print $fh_log "String for naming output files: $param{outfile_string}\n";
print $fh_log "\t$param{file_primers_completed}: Information on primers, target region parameters, primer design etc., completed with information of the best ncbi and RepeatMasker hits\n";

print $fh_log "\nGENERAL PARAMETERS\n";
print $fh_log "System (win/linux): $param{syst}\n";
print $fh_log "Output folder: $param{out_folder}\n";
print $fh_log "Path to QDD executables: $param{qdd_folder}\n";
print $fh_log "Delete intremediate files: $param{del_files}\n";
print $fh_log "Check contamination by BLASTing sequences against nt db of NCBI: $param{check_contamination}\n";
print $fh_log "Run RepeatMasker: $param{rm}\n";	
if($param{rm})
{
	print $fh_log "RepeatMasker library: $param{rm_lib}\n";
	print $fh_log "Path to RepeatMasker: $param{rm_path}\n";
}


if ($param{check_contamination})
{
 print $fh_log "Check contamination by local BLAST : $param{local_blast}\n";
 if ($param{local_blast})
 {
	print $fh_log "Path to BLAST executables: $param{blast_path}\n";
	print $fh_log "Name of database for BLAST: $param{blastdb}\n";
	print $fh_log "Evalue for the BLAST agaisnt nt: $param{db_evalue}\n";
	print $fh_log "The number of threads for BLAST: $param{num_threads}\n";
 }
}

print $fh_log "**********************************************************\n";

my $fasta = $tempfolder.'seq_to_blast.fas';


print "Preparing sequence file for BLAST/RepeatMasker\n";
write_fasta_from_table($param{input_file}, $fasta);

if($param{rm})
{
	print "Running RepeatMasker\nThis step can be long\n";

	my $rm_tbl = run_RepeatMasker($fasta, $param{rm_path}, $param{rm_lib}, $param{num_threads}, $tempfolder, $param{debug}, $fh_log);
	add_rm_info($rm_tbl, $param{input_file}, $param{file_primers_completed});
}
	

if ($param{check_contamination})
{
	my $input_primer_table = $param{input_file};
	if($param{rm})
	{
		$input_primer_table = $param{file_primers_completed};
	}

	if($param{local_blast})
	{
		my $nt_blast_out = $tempfolder.'nt_blast.tabular';
		blast_against_nt_local($fasta, $param{blast_path}, $param{blastdb}, $param{db_evalue}, $param{num_threads}, $param{syst}, $nt_blast_out, $param{debug});
		complete_primer_tbl($nt_blast_out, $input_primer_table, $param{file_primers_completed}, $param{blast_path}, $param{blastdb}, $param{syst}, $param{debug});
	}
	else
	{
		my $tax_file = $tempfolder.'tax.txt';
		my $blast_summary_file = $tempfolder.'BLAST_summary.csv';
		my $remote_blast = 'perl "'.$param{qdd_folder}.'remote_blast.pl" "'.$fasta.'" "'.$blast_summary_file.'" "'.$tax_file.'" "'.$param{db_evalue}.'" 2>'.$tempfolder.'remote_blast_messages.txt';
		if($param{debug})
		{	print $fh_log "$remote_blast\n";}
		print "\nBLASTing sequences against GenBank (remote BLAST)\nThis step can be very long\n";
		system $remote_blast;

		my %taxonomy = ();
		read_taxonomy_data($tax_file, \%taxonomy); 
		#  print Dumper(\%taxonomy);
		pool_tax_blast_primer_data(\%taxonomy, $blast_summary_file, $input_primer_table, $param{file_primers_completed});
	}
}
print $fh_log "\n\n";
my $time_end = time;
my $time_run = $time_end - $folder_id;

print $fh_log "\npipe4.pl started at : $date\n";
print $fh_log "The analyses took $time_run seconds\n";
  


if ($param{del_files}==1)
{
  delete_folder($tempfolder,$param{syst});
}

print "\nSee log_file for summary:\n$param{file_log_pipe4}\n\n";
print $fh_log "The analysis is finished\n";
close $fh_log;
exit;

##########################################################

sub add_rm_info
{
 my ($rm_file, $primer_file_in, $primer_file_out) = @_;

# make hash with the best rm hit
my %rm = ();
if (open(IN, $rm_file))
{

for(my $i = 0; $i< 3; ++$i)
{
	my $title = <IN>;
}
while(my $line = <IN>)
{
	chomp $line;
	my @line = split(' ', $line);
	if(exists $rm{$line[4]})
	{
		if($rm{$line[4]}[0]<$line[0])
		{
			$rm{$line[4]} = [$line[0], $line[1],$line[2],$line[3],$line[5],$line[6],$line[9],$line[10]];
		}
	}
	else
	{
		$rm{$line[4]} = [$line[0], $line[1],$line[2],$line[3],$line[5],$line[6],$line[9],$line[10]];
	}
}
close IN;
#print Dumper(\%rm);
# complete primer table
open(PR, $primer_file_in) or die "Cannot open primer file ($primer_file_in)\n";
my @prim = <PR>;
close PR;
open(PRO, ">$primer_file_out") or die "Cannot open primer file ($primer_file_out)\n";
my $title = shift @prim;
chomp $title;
print PRO "$title\tRM_score\tRM_div\tRM_del\tRM_ins\tRM_hit_beg\tRM_hit_end\tRM_hit_repeat\tRM_hit_class\n";
foreach my $line (@prim)
{
	chomp $line;
	my @line = split('\t', $line);
	print PRO "$line";
	if(exists $rm{$line[0]})
	{
		print PRO "\t$rm{$line[0]}[0]\t$rm{$line[0]}[1]\t$rm{$line[0]}[2]\t$rm{$line[0]}[3]\t$rm{$line[0]}[4]\t$rm{$line[0]}[5]\t$rm{$line[0]}[6]\t$rm{$line[0]}[7]\n";
	}
	else
	{
		print PRO "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
	}
}
}# end if open rm file
else
{
	print "rm file ($rm_file) does not exists\n";
	open(PR, $primer_file_in) or die "Cannot open primer file ($primer_file_in)\n";
	my @prim = <PR>;
	close PR;
	open(PRO, ">$primer_file_out") or die "Cannot open primer file ($primer_file_out)\n";
	my $title = shift @prim;
	chomp $title;
	print PRO "$title\tRM_run\n";
	foreach my $line (@prim)
	{
		chomp $line;
		print PRO "$line\tError in RM run\n";
	}
	close PRO;
}

}


##########################################################


sub blast_against_nt_local
{
 my ($query, $blast_path, $nt, $e, $num_threads, $syst, $outfile, $debug) = @_;
 
 my $blast = '';
 
 if($syst eq 'win')
 {
	$blast = '"'.$blast_path.'blastn.exe" -task blastn -db "'.$nt.'" -query "'.$query.'" -evalue '.$e.' -out "'.$outfile.'" -outfmt 6 -max_target_seqs 1 -num_threads '.$num_threads;
 }
 else
 {
	$blast = '"'.$blast_path.'blastn" -task blastn -db "'.$nt.'" -query "'.$query.'" -evalue '.$e.' -out "'.$outfile.'" -outfmt 6 -max_target_seqs 1 -num_threads '.$num_threads; 
 }
if($debug)
{
 print $fh_log "\n$blast\n\n";
} 
print "\nBLASTing sequences against $nt\nThis step can be long\n";
system $blast;
}

############################################################

sub complete_primer_tbl
{
my ($nt_blast_out, $primer_file_in, $primer_file_out, $blast_path, $blastdb, $syst, $debug) = @_;

my %code = (); # $code{seqid} = gi;
my %gi = (); # $gi{gi} = tax_id;
my %taxid = (); # $taxid{taxid} = lineage;

# read blast results file and make %code and start %gi

if (-e $nt_blast_out)
{
	if ( -z $nt_blast_out)
	{
		print $fh_log "No hit to nt database\n";
		print "No hit to nt database\n";
  		open(IN, $primer_file_in) or print "Cannot open input primer file ($primer_file_in)\n";
		my @data = <IN>;
		close IN;	
		
		open(OUT, ">$primer_file_out") or print "Cannot open output primer file ($primer_file_out)\n";
		my $title = shift @data;
		chomp $title;
		print OUT "$title\tBLAST_nt\n";  
		foreach my $line (@data)
		{
			chomp $line;
		        print OUT "$line\tNO_HIT\n";
		}
		close OUT;
	}
	else
	{
		open(IN, $nt_blast_out);
		while ( my $line = <IN>)
		{
			chomp $line;
			my @line = split("\t", $line);
			$line[1] =~ /gi\|([0-9]+)\|/;
			my $gi = $1;
			my $seq_code = $line[0];
			shift(@line);
			$line = join("\t", @line);
			$code{$seq_code} = [$gi, $line];
			$gi{$gi} = '';
		}
		close IN;

# write unique gis to file
		my $gi_file = $nt_blast_out.'.gi';
		open(GI, ">$gi_file") or die "Cannot open tempfile ($gi_file)\n";
		foreach my $gi (sort keys %gi)
		{
			print GI "$gi\n";
		}
		close GI;

# run blastdbcmd to get taxid from gis
		my $gitax = $nt_blast_out.'.gitax';
		my $blastdbcmd = '';
		if ($syst eq 'win')
		{
			$blastdbcmd = '"'.$blast_path.'blastdbcmd.exe" -db "'.$blastdb.'" -entry_batch "'.$gi_file.'" -outfmt "%g %T" -out "'.$gitax.'"';
		}
		else
		{
			$blastdbcmd = $blast_path.'blastdbcmd -db '.$blastdb.' -entry_batch '.$gi_file.' -outfmt "%g %T" -out '.$gitax;
		}
		if($debug)
		{
			print $fh_log "\n$blastdbcmd\n\n";
		}
		system $blastdbcmd;

# read output to finish %gi and start %taxid
		open(IN, $gitax) or die "Cannot open gitax file ($gitax)\n";
		while(my $line = <IN>)
		{
			my @line = split(' ', $line);
			$gi{$line[0]} = $line[1];
			$taxid{$line[1]} = '';
		}
		close IN;

# get lineage for all taxid
		print "Extract lineage for the best hits\n";
		my $taxlin = $nt_blast_out.'.taxlin';
		foreach my $tax(sort keys %taxid)
		{
			my $taxcmd = 'perl "'.$param{qdd_folder}.'ncbi_tax.pl" '.$tax.' '.$param{qdd_folder}.' >>"'.$taxlin.'"';
			system $taxcmd;
		}
# re_define %taxid $taxid{taxid}{rank} = taxon_name
		%taxid = get_lineages($taxlin);
#print Dumper(\%taxid);
		add_lineage_to_primer(\%taxid, \%code, \%gi, $primer_file_in, $primer_file_out);
	}# fin unless file empty
}# fin if exitst file
else
{
  print  "Cannot open blast output file ($nt_blast_out)\n";
}

}


###################################################

sub get_lineages
{
my ($file) = @_;
my  %tax_test= ();
open(IN, $file) or die "Cannot open $file\n";

my @data = <IN>;
close IN;
for(my $i = 0; $i<scalar@data; $i=$i+15)
{
	my $taxid = $data[$i+2];
	chomp $taxid; 
	my $lineage = $data[$i+10];
	chomp $lineage;
	$lineage =~ s/;\s*/;/g;
	my @lineage = split(';', $lineage);

	foreach my $tax (@lineage)
	{
		$tax =~ /\[(.+)\](.+)/;
		my $rank = $1;
		my $taxon = $2;

		unless($rank eq 'no rank')
		{
			$tax_test{$taxid}{$rank} = $taxon;
		}
	}
	my $spec = $data[$i+14];
	my @spec = split(';', $spec);
	$spec = $spec[-2];
	$spec =~ s/^\s*//;
	$tax_test{$taxid}{species} = $spec;
}
@data = ();
#print Dumper(\%tax_test);
return %tax_test;
}
############################################################

sub add_lineage_to_primer
{
 my ($taxid, $code, $gihash, $primer_file_in, $primer_file_out) = @_;
 # $$code{seqid} = [gi, blast_line], $$gihash{gi} = taxid, $$taxid{taxid}{rank} = taxon_name
 open(IN, $primer_file_in) or die "Cannot open the primer file ($primer_file_in)\n";
 my @data = <IN>;
 close IN;
 my @taxlevel = ('superkingdom', 'kingdom', 'phylum', 'subphylum', 'class', 'order', 'superfamily', 'family', 'genus', 'species');

open(OUT, ">$primer_file_out") or die "Cannot open the primer file for output ($primer_file_out)\n";

my $title = shift@data;
chomp $title;
print OUT "$title";
print OUT "\t", 'subject_id_best_hit', "\t", '%_identity', "\t", 'alignment_length', "\t", 'mismatches', "\t", 'gap_opens', "\t", 'q_start', "\t", 'q_end', "\t", 's_start', "\t", 's_end', "\t", 'evalue', "\t", 'bit_score'; 

foreach my $tl (@taxlevel)
{
	print OUT "\t$tl";
}
print OUT "\n";

foreach my $line (@data)
{
	chomp $line;
	my @line = split("\t", $line);
	print OUT $line;
	my $acc = $line[0];
	if(exists $$code{$acc})
	{
		print OUT "\t$$code{$acc}[1]";
		my $gi = $$code{$acc}[0];
		my $ti = ''; # taxid 
		if (exists $$gihash{$gi})
		{$ti = $$gihash{$gi}};

		if (exists $$taxid{$ti})
		{
			foreach my $tl (@taxlevel)
			{
				if (exists $$taxid{$ti}{$tl})
				{print OUT "\t$$taxid{$ti}{$tl}";}
				else
				{print OUT "\tNA";}
			}
		}
	}
	else
	{
		my $temp = "\tNO_NT_HIT" x 21;
		print OUT $temp;
	}
	print OUT "\n";

}

 close OUT;
 
}

##################################################################

sub define_param_limits_pipe4
{
my ($param) = @_;
#hash of array [type(list, real, folder), error message, ['list of values/lower_upper limits)']

  my %param_limits = (
      'galaxy'=> ['list', "galaxy must be 0 or 1\n", [0,1]], 
      'syst'=> ['list', "Operating system must be win or linux\n", ['win', 'linux']], 
      'out_folder' => ['folder', "Cannot open the output_folder $$param{out_folder}\n"], 
      'del_files'=> ['list', "del_files must be 0 or 1\n", [0,1]],
      'input_file' => ['file', "Cannot open input_file ($$param{input_file})\n"],
      'qdd_folder' => ['folder', "Cannot open qdd_folder ($$param{qdd_folder} folder)\n"],  
	  'local_blast' => ['list', "local_blast must be 1 (run local blast for contamination check) or 0 (run remote blast for contamination check)\n", [0,1]],
	  'check_contamination' => ['list', "check_contamination must be 1 (check contamination by blast against nt of NCBI) or 0 (skip contamination check)\n", [0,1]],
	  'outfile_string' => ['string', "-outfile_string can contain only alphanumerical values, underscore '_' or dot\n"],
      'debug' => ['list', "debug must be 0 or 1\n", [0,1]],
      'rm' => ['list', "rm must be 0 or 1\n", [0,1]]
      );
	  
	  if ($$param{check_contamination} and $$param{local_blast})
	  {
		$param_limits{blast_path} = ['folder', "Cannot open blast folder ($$param{blast_path} folder)\n"];
#		$param_limits{blastdb} = ['file', "BLAST db file ($$param{blastdb}) does not exists\n"];
		$param_limits{db_evalue} = ['real', "db_evalue ($$param{db_evalue}) must be between 0 and 1e-5\n", [0,1e-5]];
		$param_limits{num_threads} = ['real', "The number of threads ($$param{num_threads}) must be between 1 and 100000\n", [0,100000]];
	  }

	  if($param{rm})
	  {
		$param_limits{rm_path} = ['folder', "Cannot open RepeatMasker folder ($$param{rm_path})\n"],;
	  }

return %param_limits;
}
##################################################################

sub read_taxonomy_data
{
my ($tax_file, $taxonomy) = @_;
#my %taxonomy = (); # hash of hash $taxonomy{$TaxId}{species=> , superkingdom=> , kingdom=>, phylum=>, class=>, family=>, genus=>}

unless(open(TAX, $tax_file))
{
  print $fh_log "Cannot open $tax_file\n";
}

my $id_data = '';
while(my $line = <TAX>)
{
  if ($line =~ /^<\/Taxon>$/)
  {
    $id_data =~ /(<TaxId>)(.+?)(<\/TaxId>)/m;
    my $taxId = $2;
    $id_data =~ /(<ScientificName>)(.+?)(<\/ScientificName>)/m;
    my $spec = $2;
    initialize_hash($spec, $taxId, $taxonomy);
    while ($id_data =~ /(<Taxon>)(.+?)(<\/Taxon>)/sg)
    {
      my $taxon = $2;
      $taxon =~ /(<Rank>)(.+?)(<\/Rank>)/m;
      {
        my $rank = $2;
        if(exists $$taxonomy{$taxId}{$rank})
        {
          $taxon =~ /(<ScientificName>)(.+?)(<\/ScientificName>)/m;
          $$taxonomy{$taxId}{$rank} = $2;
        }
      }
    }
    $id_data = '';
  }
  else
  {
    $id_data .= $line;
  }
}
close TAX;
}

#################################################

sub initialize_hash
{
my ($spec, $TaxId, $taxonomy) = @_;

$$taxonomy{$TaxId}{species} = $spec;
$$taxonomy{$TaxId}{superkingdom} ='';
$$taxonomy{$TaxId}{kingdom} ='';
$$taxonomy{$TaxId}{phylum} ='';
$$taxonomy{$TaxId}{class} ='';
$$taxonomy{$TaxId}{family} ='';
$$taxonomy{$TaxId}{genus} ='';
}

#################################################

sub pool_tax_blast_primer_data
{
my ($taxonomy, $blast_summary_file, $primer_file3, $primer_file4) = @_;

if (open(BL, $blast_summary_file))
{
my %bl = ();
my $title_bl = <BL>;
chomp $title_bl;
$title_bl =~ s/;/\t/g;
$title_bl =~ s/[^\t]+\t//;
while (my $line = <BL>)
{
  chomp $line;
  my @line = split(';', $line);
  my $code = shift@line;
  $line = join("\t", @line);
  $bl{$code}[0] = $line[5];
  $bl{$code}[1] = $line;

}
close BL;

open_file_in($primer_file3, 'PRIN');
my @prin = <PRIN>;
close PRIN;
open_file_out($primer_file4, 'OUT');
my $title = shift @prin;
chomp $title;
print OUT "$title\t$title_bl\tsuperkingdom\tkingdom\tphylum\tclass\tfamily\tgenus\tspecies\n";
foreach my $line (@prin)
{
  chomp $line;
  my @line = split('\t', $line);
  if (exists $bl{$line[0]})
  {
    print OUT "$line\t$bl{$line[0]}[1]\t";
    print OUT "$$taxonomy{$bl{$line[0]}[0]}{superkingdom}\t$$taxonomy{$bl{$line[0]}[0]}{kingdom}\t$$taxonomy{$bl{$line[0]}[0]}{phylum}\t$$taxonomy{$bl{$line[0]}[0]}{class}\t$$taxonomy{$bl{$line[0]}[0]}{family}\t$$taxonomy{$bl{$line[0]}[0]}{genus}\t$$taxonomy{$bl{$line[0]}[0]}{species}\n";
  }
  else
  {
	my $temp = "\tNO_NT_HIT" x 13;

    print OUT "$line";
    print OUT "$temp\n";
  }
}
close PRIN;
close OUT;
}# end if open blast summary file
else
{print $fh_log "Cannot open blast summary file ($blast_summary_file)\n"}
}


#################################################################

sub write_fasta_from_table
{
 my ($primer_file, $fasta) = @_;

open(IN, $primer_file) or die "Cannot open primer_file ($primer_file)\n";
my %seq = ();
my $title = <IN>;
while ( my $line = <IN>)
{
	chomp $line;
	my @line = split("\t", $line);
	unless (exists $seq{$line[0]})
	{
		$seq{$line[0]} = uc $line[30];
	}
}
close IN;

write_hash_to_fasta(\%seq, $fasta);

}
###################################################################
