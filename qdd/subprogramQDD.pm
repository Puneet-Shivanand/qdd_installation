

sub add_slash_to_folder
{
 my ($param, $param_limits) = @_;

foreach my $par (keys %$param_limits)
{
	if ($$param_limits{$par}[0] eq 'folder' and not($$param{$par} eq ''))
	{
		$$param{$par} =~  s/[\/\\]+$//;
		$$param{$par} .= '/';
	}
}
}

################################################
sub call_blast_adapter_no_filter
{
    my ($filename, $filename_out, $adaptor_file, $Evector, $path_blast, $blast_version, $syst) = @_;
    my $formatdb = '';
    my $blast = '';

    if ($blast_version eq 'blast')
    {
      $formatdb = '"'.$path_blast.'formatdb" -i '.$adaptor_file.' -p F -o T';
      $blast = '"'.$path_blast.'blastall" -p blastn -d "'.$adaptor_file.'" -i "'.$filename.'" -o "'.$filename_out.'" -T F -e '.$Evector.' -m 8 -F F' ;
    }
    elsif ($syst eq 'win') 
    {
      $formatdb = '"'.$path_blast.'makeblastdb.exe" -dbtype nucl -in "'.$adaptor_file.'" -parse_seqids >blast_message.txt';
      $blast = '"'.$path_blast.'blastn.exe" -task blastn -db "'.$adaptor_file.'" -query "'.$filename.'" -evalue '.$Evector.' -out "'.$filename_out.'" -outfmt 6 -dust no';      
    }
    else
    {
      $formatdb = '"'.$path_blast.'makeblastdb" -dbtype nucl -in "'.$adaptor_file.'" -parse_seqids >blast_message.txt';
      $blast = '"'.$path_blast.'blastn" -task blastn -db "'.$adaptor_file.'" -query "'.$filename.'" -evalue '.$Evector.' -out "'.$filename_out.'" -outfmt 6 -dust no';      
    }
    
    system $formatdb;
#    print "$formatdb\n";

    system $blast;
#    print "$blast\n"; 

}
##########################################################

sub call_blast_adapter_no_filter_1
{
    my ($filename, $filename_out, $adaptor_file, $Evector, $path_blast, $blast_version, $syst, $num_threads) = @_;
    my $formatdb = '';
    my $blast = '';

    if ($blast_version eq 'blast')
    {
      $formatdb = '"'.$path_blast.'formatdb" -i '.$adaptor_file.' -p F -o T';
      $blast = '"'.$path_blast.'blastall" -p blastn -d "'.$adaptor_file.'" -i "'.$filename.'" -o "'.$filename_out.'" -T F -e '.$Evector.' -m 8 -F F' ;
    system $formatdb;
    system $blast;
    }
    elsif ($syst eq 'win') 
    {
      $formatdb = '"'.$path_blast.'makeblastdb.exe" -dbtype nucl -in "'.$adaptor_file.'" -parse_seqids >blast_message.txt';
      $blast = '"'.$path_blast.'blastn.exe" -task blastn -db "'.$adaptor_file.'" -query "'.$filename.'" -evalue '.$Evector.' -out "'.$filename_out.'" -outfmt 6 -dust no -num_threads '.$num_threads;      
    system $formatdb;
    system $blast;
	system 'del blast_message.txt';
    }
    else
    {
      $formatdb = '"'.$path_blast.'makeblastdb" -dbtype nucl -in "'.$adaptor_file.'" -parse_seqids >/dev/null';
      $blast = '"'.$path_blast.'blastn" -task blastn -db "'.$adaptor_file.'" -query "'.$filename.'" -evalue '.$Evector.' -out "'.$filename_out.'" -outfmt 6 -dust no -num_threads '.$num_threads;      
    system $formatdb;
    system $blast;
    }

}

################################################
sub call_blast_adapter_no_filter_word_size
{
    my ($filename, $filename_out, $adaptor_file, $Evector, $path_blast, $blast_version, $word_size, $syst) = @_;

    my $formatdb = '';
    my $blast = '';
    
    if ($blast_version eq 'blast')
    {
      $formatdb = '"'.$path_blast.'formatdb" -i "'.$adaptor_file.'" -p F -o T';
      $blast = '"'.$path_blast.'blastall" -p blastn -d "'.$adaptor_file.'" -i "'.$filename.'" -o "'.$filename_out.'" -T F -e '.$Evector.' -m 8 -F F -W '.$word_size;
    }
    elsif ($syst eq 'win') 
    {
      $formatdb = '"'.$path_blast.'makeblastdb.exe" -dbtype nucl -in "'.$adaptor_file.'" -parse_seqids >blast_message.txt';
      $blast = '"'.$path_blast.'blastn.exe" -task blastn -db "'.$adaptor_file.'" -query "'.$filename.'" -evalue '.$Evector.' -word_size '.$word_size.' -out "'.$filename_out.'" -outfmt 6 -dust no';
    }
    else
    {
      $formatdb = '"'.$path_blast.'makeblastdb" -dbtype nucl -in "'.$adaptor_file.'" -parse_seqids >blast_message.txt';
      $blast = '"'.$path_blast.'blastn" -task blastn -db "'.$adaptor_file.'" -query "'.$filename.'" -evalue '.$Evector.' -word_size '.$word_size.' -out "'.$filename_out.'" -outfmt 6 -dust no';
    }
    
#    print "$formatdb\n";
    system $formatdb;

#    print "$blast\n"; 
    system $blast;

}
################################################
sub call_blast_adapter_no_filter_word_size_1
{
    my ($filename, $filename_out, $adaptor_file, $Evector, $path_blast, $blast_version, $word_size, $syst, $num_threads) = @_;

    my $formatdb = '';
    my $blast = '';
    
    if ($blast_version eq 'blast')
    {
      $formatdb = '"'.$path_blast.'formatdb" -i "'.$adaptor_file.'" -p F -o T';
      $blast = '"'.$path_blast.'blastall" -p blastn -d "'.$adaptor_file.'" -i "'.$filename.'" -o "'.$filename_out.'" -T F -e '.$Evector.' -m 8 -F F -W '.$word_size;
    system $formatdb;
    system $blast;
    }
    elsif ($syst eq 'win') 
    {
      $formatdb = '"'.$path_blast.'makeblastdb.exe" -dbtype nucl -in "'.$adaptor_file.'" -parse_seqids >blast_message.txt';
      $blast = '"'.$path_blast.'blastn.exe" -task blastn -db "'.$adaptor_file.'" -query "'.$filename.'" -evalue '.$Evector.' -word_size '.$word_size.' -out "'.$filename_out.'" -outfmt 6 -dust no -num_threads '.$num_threads;
    system $formatdb;
    system $blast;
	system 'del blast_message.txt';
    }
    else
    {
      $formatdb = '"'.$path_blast.'makeblastdb" -dbtype nucl -in "'.$adaptor_file.'" -parse_seqids >/dev/null';
      $blast = '"'.$path_blast.'blastn" -task blastn -db "'.$adaptor_file.'" -query "'.$filename.'" -evalue '.$Evector.' -word_size '.$word_size.' -out "'.$filename_out.'" -outfmt 6 -dust no -num_threads '.$num_threads;
    system $formatdb;
    system $blast;
    }
    


}

#################################################################
sub call_blast_query_to_db_soft_mask
{
my ($query_file, $db_file, $filename_out, $Eblast, $blast_path, $blast_version, $syst) = @_;

my $formatdb = '';
my $blast = '';
if ($blast_version eq 'blast')
{
    $formatdb = '"'.$blast_path.'formatdb" -i "'.$db_file.'" -p F -o T';
    $blast = '"'.$blast_path.'blastall" -p blastn -d "'.$db_file.'" -i "'.$query_file.'" -o "'.$filename_out.'" -T F -e '.$Eblast.' -m 8 -F m -U T';
}
elsif ($syst eq 'win') 
{
    $formatdb = '"'.$blast_path.'makeblastdb.exe" -dbtype nucl -in "'.$db_file.'" -parse_seqids >blast_message.txt';
    $blast = '"'.$blast_path.'blastn.exe" -task blastn -db "'.$db_file.'" -query "'.$query_file.'" -out "'.$filename_out.'" -evalue '.$Eblast.' -outfmt 6 -lcase_masking -soft_masking true';
}
else
{
    $formatdb = '"'.$blast_path.'makeblastdb" -dbtype nucl -in "'.$db_file.'" -parse_seqids >blast_message.txt';
    $blast = '"'.$blast_path.'blastn" -task blastn -db "'.$db_file.'" -query "'.$query_file.'" -out "'.$filename_out.'" -evalue '.$Eblast.' -outfmt 6 -lcase_masking -soft_masking true';

}
    system $formatdb;
 #   print "$formatdb\n";

    system $blast;
 #   print "$blast\n"; 
}

############################################################"
sub call_blast_query_to_db_soft_mask_1
{
my ($query_file, $db_file, $filename_out, $Eblast, $blast_path, $blast_version, $syst, $num_threads) = @_;

my $formatdb = '';
my $blast = '';
if ($blast_version eq 'blast')
{
    $formatdb = '"'.$blast_path.'formatdb" -i "'.$db_file.'" -p F -o T';
    $blast = '"'.$blast_path.'blastall" -p blastn -d "'.$db_file.'" -i "'.$query_file.'" -o "'.$filename_out.'" -T F -e '.$Eblast.' -m 8 -F m -U T';
   system $formatdb;
    system $blast;
}
elsif ($syst eq 'win') 
{
    $formatdb = '"'.$blast_path.'makeblastdb.exe" -dbtype nucl -in "'.$db_file.'" -parse_seqids >blast_message.txt';
    $blast = '"'.$blast_path.'blastn.exe" -task blastn -db "'.$db_file.'" -query "'.$query_file.'" -out "'.$filename_out.'" -evalue '.$Eblast.' -outfmt 6 -lcase_masking -soft_masking true -num_threads '.$num_threads;
   system $formatdb;
   system $blast;
	system 'del blast_message.txt';
}
else
{
    $formatdb = '"'.$blast_path.'makeblastdb" -dbtype nucl -in "'.$db_file.'" -parse_seqids >/dev/null';
    $blast = '"'.$blast_path.'blastn" -task blastn -db "'.$db_file.'" -query "'.$query_file.'" -out "'.$filename_out.'" -evalue '.$Eblast.' -outfmt 6 -lcase_masking -soft_masking true -num_threads '.$num_threads;
   system $formatdb;
    system $blast;

}
 
}

############################################################"
sub call_blast_tbl_soft_filter
{
    my ($filename, $filename_out, $Eblast, $blast_path, $blast_version, $syst) = @_;

my $formatdb = '';
my $blast = '';
if ($blast_version eq 'blast')
{
    $formatdb = '"'.$blast_path.'formatdb" -i "'.$filename.'" -p F -o T';
    $blast = '"'.$blast_path.'blastall" -p blastn -d "'.$filename.'" -i "'.$filename.'" -o "'.$filename_out.'" -T F -e '.$Eblast.' -m 8 -F m -U T';
}
elsif ($syst eq 'win') 
{
    $formatdb = '"'.$blast_path.'makeblastdb.exe" -dbtype nucl -in "'.$filename.'" -parse_seqids';
    $blast = '"'.$blast_path.'blastn.exe" -task blastn -db "'.$filename.'" -query "'.$filename.'" -out "'.$filename_out.'" -evalue '.$Eblast.' -outfmt 6 -lcase_masking -soft_masking true';
}
else
{
    $formatdb = '"'.$blast_path.'makeblastdb" -dbtype nucl -in "'.$filename.'" -parse_seqids';
    $blast = '"'.$blast_path.'blastn" -task blastn -db "'.$filename.'" -query "'.$filename.'" -out "'.$filename_out.'" -evalue '.$Eblast.' -outfmt 6 -lcase_masking -soft_masking true';
}
    system $formatdb;
 #   print "$formatdb\n";

    system $blast;
 #   print "$blast\n"; 
}

###################################################################

sub check_repeatmasker_db
{
	my ($spec, $path, $tempfolder, $syst) = @_;
	
my $res_file = $tempfolder.'/queryRepeatDatabase_stat.txt';
my $rm = $path.'/util/queryRepeatDatabase.pl -species '.$spec.' -stat >'.$res_file.' 2>&1';
system $rm;

open(IN, $res_file) or die "Cannot open queryRepeatDatabase_stat.txt\n";
my @res = <IN>;
close IN;
my $res = join('', @res);
if($res =~ /Species $spec is not in the database/)
{
	print $res;
	delete_folder($tempfolder, $syst);
	exit;
}
else
{
	my $sum = $res[-3].$res[-2].$res[-1];
	print "In the RepeaMasker Database for $spec there are \n$sum\n";
}

}

######################################################################
sub check_set_file
{
my ($set_file) = @_;

unless (-e $set_file)
{
	$set_file = 'set_qdd_default.ini';
	unless(-e $set_file)
	{
		print "set_qdd_default.ini should be either in the current working directory, or in /etc/qdd/\n";
		exit;
	}
}
return $set_file;
}

########################################################
sub copy_files
{
my ($infolder, $outfolder, $syst, $motif) = @_;


unless ( opendir(FOLDER, "$infolder") )
{
      print "Cannot access to folder $infolder\n";
      exit;
}
my @files = grep ( !/^\.\.?$/, readdir(FOLDER) );
closedir(FOLDER);

foreach my $file (@files)
{
  if ($file =~ /$motif/)
  {
    my $copy = '';
    if ($syst eq 'win')
    {
  #    $copy = 'copy "'.$infolder.'\\'.$file.'" "'.$outfolder.'\\'.$file.'"';
     $copy = 'copy "'.$infolder.'\\'.$file.'" "'.$outfolder.'"';

      system $copy; 
    }
    if ($syst eq 'linux')
    {
      $copy = 'cp '.$infolder.'/'.$file.' '.$outfolder;
      system $copy; 
    }
  }
}

}

####################################################

sub copy_one_file
{
my ($source_folder, $file, $dest_folder, $syst) = @_;

    my $copy = '';
    if ($syst eq 'win')
    {
      $file = $source_folder.'\\'.$file;
      $copy = 'copy "'.$file.'" "'.$dest_folder.'"';
      system $copy; 
    }
    if ($syst eq 'linux')
    {
      $file = $source_folder.'/'.$file;
      $copy = 'cp '.$file.' '.$dest_folder;
      system $copy; 
    }

}

####################################################
sub copy_rename_file
{
($input, $copied, $syst) = @_;


    my $copy = '';
    if ($syst eq 'win')
    {
      $copy = 'copy "'.$input.'" "'.$copied.'"';
      system $copy; 
    }
    if ($syst eq 'linux')
    {
      $copy = 'cp '.$input.' '.$copied;
      system $copy; 
    }


}
#####################################################

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

##################################################################

sub check_param_limits
{
my ($param, $param_limits) = @_;

my $bool = 1;
foreach my $p (keys %$param_limits)
{
  if ($$param_limits{$p}[0] eq 'list')
  {
    my $bool_loc = 0;
    foreach my $v (@{$$param_limits{$p}[2]})
    {
      if ($$param{$p} eq $v)
      {
        $bool_loc = 1;
      }
    }
    if ($bool_loc == 0)
    {
      print "$$param_limits{$p}[1]\n";
      $bool = 0;
    }
  }
  elsif ($$param_limits{$p}[0] eq 'real')
  {
    unless ($$param{$p}>= $$param_limits{$p}[2][0] and $$param{$p}<= $$param_limits{$p}[2][1])
    {
      print "$$param_limits{$p}[1]\n";
      $bool = 0;
    }
  }
  elsif ($$param_limits{$p}[0] eq 'folder')
  {
  unless($$param{$p} eq '')
  {
    if (opendir(F, $$param{$p}))
    {
      closedir (F);    
    }
    else
    {
      print "$$param_limits{$p}[1]\n";
      $bool = 0;
    }
  }
  }
  elsif ($$param_limits{$p}[0] eq 'file')
  {
    unless (-e $$param{$p})
    {
      print "$$param_limits{$p}[1]\n";
      $bool = 0;
    }
  }
  elsif ($$param_limits{$p}[0] eq 'string')
  {
    if ($$param{$p} =~  /\w/)
    {
	if ($$param{$p} =~  /[^a-z^A-Z^0-9^_^\.^\s]/)
    	{

     	 print "$$param_limits{$p}[1]\n";
      	$bool = 0;
    	}
    }
    else
    {
      print "$$param_limits{$p}[1]\n";
      $bool = 0;
    }
  }
}

if ($bool == 0)
{
	print Dumper(\%$param);
  	exit;
}
}

####################################################

sub convert_fasq_into_fasta
{
 my ($fastq, $fas) = @_;
 open(IN, $fastq) or die "Cannnot open fastq file ($fastq)\n";
 open(OUT, ">$fas") or die "Cannnot open output fasta file ($fas)\n";
 
 my $count = 0;
 while(my $line = <IN>)
 {
	++$count;
	if($count == 1)
	{
		if ($line =~ s/^@/>/)
		{
			print OUT $line;
		}
		else
		{
			print "Problem of input fastq file format. Each sequence should be printed on one line.\n";
		}
	}
	elsif ($count == 2)
	{
		print OUT $line;
	}
	elsif($count == 4)
	{
		$count = 0;
	}
 }
 close IN;
 close OUT;
}

###################################################

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


####################################################
sub delete_folder
{
my ($folder, $syst) =@_;

if ( opendir(FOLDER, $folder) )
{
	closedir(FOLDER);
	my $rmdir = '';
	if ($syst eq 'win')
	{
		$folder =~ s/\//\\/g;
		$rmdir = 'rmdir /s /q "'.$folder.'"';
	}
	else
	{
		$rmdir = 'rm -r '.$folder;
	}
	system $rmdir;

}

}

######################################################


sub delete_tag_adapt_files
{
my ($file, $syst) = @_;

if ($syst eq 'win')
{
	$file =~ s/\//\\/g;
	my $del = 'del '.$file.'.n*';
  system $del;
}
else
{
  system 'rm '.$file.'.*';
}
}
#################################################################


sub determine_last_folder
{
my ($folder, $motif) = @_;
my $sub_folder = '';
# supposes that folder names are in the format of [a-z0-9]+_date
 
my %list = get_file_folder_list($folder);

foreach my $folder (sort keys %list)
{
  if ($folder =~ /$motif/)
  {
    unless($folder =~ /\./)
    {
      $sub_folder = $folder;
    }
  }
}

$sub_folder = $folder.'/'.$sub_folder;
return $sub_folder;
}


#############################################################

sub file_versioning_output_old
{
my ($param) = @_;

my $c_max = 0;
foreach my $par (keys %$param)
{
    if ($par =~ /^file_/)
    {
	if( -e $$param{$par})
	{
		print "$$param{$par}\n";
		my $c = 1;
		if($$param{$par} =~  /\.([^\.]*)$/)
		{
			my $ext = $1;
			$$param{$par} =~ s/\.$ext/\(1\)\.$ext/;
		}
		print "$$param{$par}\n";
		while (-e $$param{$par} and $$param{$par} =~ /\(([0-9]+)\)\./)
		{
			$c=$1;
			my $new_c = $c+1;
			$$param{$par} =~  s/\($c\)\./\($new_c\)\./;
		}
		print "$$param{$par}\n";
		if($c_max < $c)
		{$c_max=$c;}
	}
    }
}
print "c_max = $c_max\n";
if ($c_max > 0)
{
	++$c_max;
	foreach my $par (keys %$param)
	{
		if ($par =~ /^file_/)
		{
			if($$param{$par} =~ /\([0-9]+\)\./)
			{
				$$param{$par} =~ s/\([0-9]+\)\./\($c_max)\./;
			}
			else
			{
				$$param{$par} =~ s/\./\($c_max)\./;
			}
		}

	}
}
}

#############################################################

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

######################################################

sub get_clustal_executable_name
{# determines the exact executable name of clustal path/
	my ($path, $syst) = @_;
	my $ex_name = '';


	if($syst eq 'linux')
	{
		$ex_name = `which clustalw`;
		chomp $ex_name;
		if ($ex_name eq '')
		{
			$ex_name = `which clustalw2`;
			chomp $ex_name;
			if ($ex_name eq '')
			{
				print "Clustal excutable is not found in $path\n";
				exit;
			}
		}
	}
	else
	{
		if(-e $path.'clustalw.exe')
		{ $ex_name = $path.'clustalw.exe'}
		elsif (-e $path.'clustalw2.exe')
		{ $ex_name = $path.'clustalw2.exe'}
		else
		{print "Clustal excutable is not found in $path\n";
		exit;
		}
	}
return $ex_name;
}

#########################################################
sub get_file_folder_list
{
my ($folder) = @_;
my %list = ();

unless ( opendir(FOLDER, "$folder") )
{
      print "Cannot access to folder $folder\n";
      exit;
}
my @files = grep ( !/^\.\.?$/, readdir(FOLDER) );
closedir(FOLDER);
foreach my $file (@files)
{
  $list{$file} = '';
}
@files = ();
return %list;
}

##############################################################

sub get_file_list
{
my ($folder) = @_;
my %list = ();

unless ( opendir(FOLDER, "$folder") )
{
      print "Cannot access to folder $folder\n";
      exit;
}
my @files = grep ( !/^\.\.?$/, readdir(FOLDER) );
closedir(FOLDER);
foreach my $file (@files) # adds path and keeps only files not directories
{
	my $file = $folder.$file;
	if (-f $file)
	{
		$list{$file} = '';
	}

}
@files = ();
return %list;
}
##############################################################
sub get_filename_root
{
my ($file) = @_;
my $root = $file;

$root =~ s/.*[\\\/]//;
$root =~ s/\..*//;
return $root;
}

#################################################################

sub get_last_version_and_modify_1
{
my ($folder, $motifs, $froot, $param) = @_;

# returns the highest version to be added to new files e.g._v2
#if no file with any of the motifs => 0

my %files = get_file_folder_list($folder);
my $c_max = 0;
foreach my $file (keys %files) # get the higest version among all files in the folder
{
	if ($file =~ /_v([0-9]+)(\.[^\.]*)*$/)
	{
		if($1> $c_max)
		{
			$c_max = $1;
		}
	}
}


my $f_max = -1;
foreach my $motif (@$motifs)
{
   my $f = -1; 
   my $file = $folder.$froot.$$param{$motif};
   if (-e $file)
   {
	$f = 0;
   }

   for(my $i=1; $i<=$c_max; ++$i)
   {
	$motif_temp = $$param{$motif};
	$motif_temp =~ s/\./_v$i./; 
	$file = $folder.$froot.$motif_temp;
	if ( -e $file)
	{$f = $i}
   }

   if($f>$f_max)
   {
	$f_max = $f;
   }

}

++$f_max;

if ($f_max > 0)
{
	modify_motif($motifs, $f_max, $param);
}
}

#################################################################

sub get_last_version_and_modify
{
my ($folder, $motifs, $froot, $param, $inout) = @_;

# returns the highest version to be aded to new files e.g.(1)
#if no file with any of the motifs => 0

my %files = get_file_folder_list($folder);
my $c_max = 0;
foreach my $file (keys %files) # get the higest version among all files in the folder
{
	if ($file =~ /\(([0-9]+)\)/)
	{
		if($1> $c_max)
		{
			$c_max = $1;
		}
	}
}


my $f_max = -1;
foreach my $motif (@$motifs)
{
   my $f = -1; 
   my $file = $folder.$froot.$$param{$motif};
   if (-e $file)
   {
	$f = 0;
   }

   for(my $i=1; $i<=$c_max; ++$i)
   {
	$motif_temp = $$param{$motif};
	$motif_temp =~ s/\./($i)./; 
	$file = $folder.$froot.$motif_temp;
	if ( -e $file)
	{$f = $i}
   }

   if($f>$f_max)
   {
	$f_max = $f;
   }

}

if($inout eq 'output')
{
	++$f_max;
}

if ($f_max > 0)
{
	modify_motif($motifs, $f_max, $param);
}
}
##################################################################
sub get_seq_number_from_fasta
{
 my ($file) = @_;
my $seqn = 0;
open(IN, $file) or die "Cannot open $file $!\n";
while(my $line = <IN>)
{
	if ($line =~ />/)
	{
		++$seqn;
	}
}
close IN;

return $seqn;
}
##################################################################
sub length_sequences
{
    my ($hach) = @_;
    my @code = keys %$hach;
    my %orig_length  = ();

    foreach my $code (@code)
    {
	$orig_length{$code} = length($$hach{$code});
    }
    return %orig_length;
}

##################################################

sub make_folder
{
my ($syst, $folder_name) = @_;

my $make_folder = ''; 
if ($syst eq 'win')
{
	$make_folder = 'md "'.$folder_name.'"';
#	print "$make_folder\n";
	system $make_folder;
}
elsif ($syst eq 'linux')
{
	$make_folder = 'mkdir '.$folder_name;
	system $make_folder;
}
}

#####################################################

sub min
{
	my @sorted = sort{$a<=>$b}@_;
	my $min = shift@sorted;
	return $min;  
}

#####################################################

sub max
{
	my @sorted = sort{$a<=>$b}@_;
	my $max = pop@sorted;
	return $max;  
}

######################################################
sub modify_params_from_tags
{
  my ($param, $inp) = @_;
  
my @bad_tags = ();
  for(my $i = 0; $i<scalar@$inp; $i=$i+2)
  {
	$$inp[$i] =~ s/^-*//;
	if(exists $$param{$$inp[$i]})
	{
		$$param{$$inp[$i]} = $$inp[$i+1];
	}
	else
	{
		push(@bad_tags, $$inp[$i]);

	}

  } 
if(scalar @bad_tags > 0)
{
	print "The following tags are not accepted: @bad_tags\n";
	print_usage();
	exit;
}
}
####################################################

sub modify_motif
{
 my ($motifs, $version, $param) = @_;
foreach my $mot (@$motifs)
{
	my $temp = $$param{$mot};
	$temp =~  s/\./_v$version\./;
	$$param{$mot} = $temp;
}
}

#######################################################"
sub move_files
{
my ($infolder, $outfolder, $syst, $motif) = @_;


unless ( opendir(FOLDER, "$infolder") )
{
      print "Cannot access to folder $infolder\n";
      exit;
}
my @files = grep ( !/^\.\.?$/, readdir(FOLDER) );
closedir(FOLDER);

foreach my $file (@files)
{
  if ($file =~ /$motif/)
  {
    my $move = '';
    if ($syst eq 'win')
    {
     $move = 'move "'.$infolder.'\\'.$file.'" "'.$outfolder.'"';

      system $move; 
    }
    if ($syst eq 'linux')
    {
      $move = 'mv '.$infolder.'/'.$file.' '.$outfolder;
      system $move; 
    }
  }
}

}


###################################################################"
sub move_files_local
{
my ($source_folder, $dest_folder, $syst) =@_;

my $copy = '';
if ($syst eq 'win')
{
   $copy = 'MOVE "'.$source_folder.'\*.*" "'.$dest_folder.'" >nul';
}
else
{
   $copy = 'mv '.$source_folder.'/*.* '.$dest_folder;
#print "$copy\n";
}

system $copy;
}

#####################################################################
sub open_file_in
{
 my ($filename, $handle) = @_;
 
 unless (open($handle, "$filename"))
 {
     print "Cannot open $filename\n";
     exit;
 }

}
######################################################

sub open_file_and_hach
{
    my ($filename, $first_line) = @_;
#   print "$first_line"; 
	  unless ( open (FILE, "$filename"))
	  { print "Cannot open $filename\n"; exit;}
	  my @data = <FILE>;
    if ($first_line == 0)
    {
	shift (@data);
    }
#	  print "@data";
	  my %hach = txt_into_hachage_1(@data);
	  close FILE;

    return %hach;
}

#############################################################

sub open_file_out
{
 my ($filename, $handle) = @_;
 
 unless (open($handle, ">$filename"))
 {
     print "Cannot open $filename\n";
     exit;
 }

}

###########################################################
sub pool_files
{
my ($file1, $file2, $pool_file) = @_;

open_file_out($pool_file, 'OUT');
open_file_in($file1, 'F1');
while (my $line = <F1>)
{
print OUT "$line";
}
close F1;

open_file_in($file2, 'F2');
while (my $line = <F2>)
{
print OUT "$line";
}
close OUT;
close F2;
}

##############################################################

sub print_hachage
 {
 # print out a hachage in a human readable format
  use warnings;
  use strict;   
   my ($hach) = @_;
   my @cles = keys %$hach;
   my @sorted_keys = sort @cles;

   foreach my $sorted_key (@sorted_keys)
     {
       print "$sorted_key     $$hach{$sorted_key}\n";
     }
  print "\n\n";
  }
  
##################################################################

sub print_heading
{
  print "**********************************************************\n";
  print "QDD version2.1(beta) 28 June 2011\n";
  print "Emese Meglecz, Aix-Marseille University, Marseille, France\n";
  print 'emese.meglecz@univ-provence.fr';
  print "\n";
  print "Plesae, read the Documentation_QDD2.pdf\n";
  print "**********************************************************\n";
}

##################################################################

sub print_heading_3_0
{
my ($fh) = @_;
  print $fh "\n**********************************************************\n";
  print $fh "QDD version3.1.2 July 2014\n";
  print $fh "Emese Meglecz, Aix-Marseille University, Marseille, France\n";
  print $fh 'emese.meglecz@imbe.fr';
  print $fh "\n";
  print $fh "http://net.imbe.fr/~emeglecz/qdd.html\n";
  print $fh "**********************************************************\n\n";
}


######################################################

sub read_hach_from_fasta_code_wo_space
{
my ($filename) = @_;
my %seq = ();
	my $i = 0;
	open_file_in($filename, 'IN');
	$/ = ">";
	while (my $seq = <IN>)
	{
		$seq =~ s/>//;
		unless ($seq eq '')
		{
			$seq =~ /.*\n/;
			my $code = $&;
			$seq =~ s/$code//;
			$seq =~ s/\s//g;
			$code =~ s/\s.*//;
			$code =~ s/\s//;
			if (exists $seq{$code})
			{
				++$i;
			}
			$seq{$code} = $seq;
		}
	}
close IN;
	if ($i>0)
	{
		print "number of sequences with code already used by other sequence: $i\n";
	}
	$/ = "\n";

return %seq;
}

######################################################

sub read_hach_from_qual_code_wo_space
{
my ($filename) = @_;

open_file_in($filename, 'IN');
my $code = '';
my $qual = '';
my %qual = ();
while (my $line = <IN>)
{
	chomp $line;
	if ($line =~ />/)
	{
		$code = $line;
		$code =~ s/\s.*//;
		$code =~ s/>//;
	}
	else
	{
		$line =~ s/\s+$//;
		$qual{$code} .= $line.' ';
	}
}
close IN;

return %qual;
}

####################################################

sub read_hash_from_tbl
{
	my ($file, $sep, $key_col, $title_line_numb) = @_;
	my %hash = ();
	open(IN, $file) or die "Cannor open file to be hashed ($file)\n";
	for(my $i = 0; $i< $title_line_numb; ++$i)
	{
		my $title = <IN>;
	}
	while(my $line = <IN>)
	{
		chomp $line;
		my @line = split($sep, $line);
		if(exists $hash{$line[$key_col]})
		{
			print "$line[$key_col] is present at more than once in $file\n";
		}
		else
		{
			$hash{$line[$key_col]} = $line;
		}
	}
return %hash;
}
######################################################################

sub rename_files
{
my ($infolder, $outfolder, $syst, $motif_source, $motif_dest) = @_;


unless ( opendir(FOLDER, "$infolder") )
{
      print "Cannot access to folder $infolder\n";
      exit;
}
my @files = grep ( !/^\.\.?$/, readdir(FOLDER) );
closedir(FOLDER);

foreach my $file (@files)
{
  if ($file =~ /$motif_source/)
  {
    my $dest_file = $file;
    $dest_file =~ s/$motif_source/$motif_dest/;
    my $move = '';
    if ($syst eq 'win')
    {
     $move = 'move "'.$infolder.'\\'.$file.'" "'.$outfolder.'\\'.$dest_file.'"';

      system $move; 
    }
    if ($syst eq 'linux')
    {
      $move = 'mv '.$infolder.'/'.$file.' '.$outfolder.'/'.$dest_file;
      system $move; 
    }
  }
}

}

##################################################################

sub print_new_set_file
{
my ($param, $set_file) = @_;
open_file_out($set_file, 'SETOUT');

foreach my $p (keys %$param)
{
  print SETOUT "$p=$$param{$p}\n";
}

close SETOUT;
}

##############################################################

sub print_param_menu
{
  my ($param_list, $param) = @_;
  print "\n";
  foreach my $numb (sort {$a<=>$b} keys %$param_list)
  {
    print "$numb : $$param_list{$numb}[1] $$param{$$param_list{$numb}[0]}\n";
  }
  print "\n";
  my $reply = '';
  print "Press enter if all of the settings are correct, or the number of the parameter if you whish to change the settings!\n";
  $reply = <STDIN>;
  return $reply;
}

#######################################################

sub print_usage
{

print "\nThe following tags are accepted:
The parameters in the set_qdd_default.txt are used if they are not scecified in the command line

GENERAL PARAMETERS:

\t-input_file [fasta file with input sequences to be analysed]
\t-galaxy [0/1] (run scrip from galaxy) 
\t-syst [linux/win] (operating system; default= linux)
\t-blast_path [Full path to BLAST+ executables]
\t-clust_path [Full path to custalw executables]
\t-primer3_path [Full path to Primer3 executables]
\t-primer3_version [1/2] (Primer3 version 1.xx => 1, 2.xx =>2)
\t-qdd_folder [Full path to QDD scripts]
\t-out_folder[folder for output files]
\t-del_files [0/1] (1 for deleting temporary files after the run)
\t-outfile_string [string] string to specify the begining of the names of the output files
\t-debug [0/1] (1 for printing out more details in log file for debugging)
\t-blastdb [name, inluding full path of a local ncbi database] (Only needed if local BLAST is used for contamination check)
\t-num_threads [integer] (number of threads for BLAST)
\t-local_blast [0/1] (1=>run local blast for contamination check; 0=>run remote blast for contamination check)

 

PIPE1 SPECIFIC PARAMETERS
\t-fastq [0/1] (1 if infor file is fastq format, 0 if fasta)
\t-contig [0/1] (1 if sequences has been assembled (contigs, scaffolds, chromosomes), 0 if they are short sequencing reads)
\t-flank_length [integer] (length of the flanking region, if extracting microsatellites from contigs)
\t-adapter [0/1] (1 for running adaptor/vector clipping step)
\t-length_limit [integer] (sequences shorter then length_limit are eliminated) 
\t-adapter_file [fasta file with adaters] (optional)


PIPE2 SPECIFIC PARAMETERS
\t-make_cons [0/1] (Make consensus sequences (YES=1/NO=0; default 1))
\t-ident_limit [integer] Minimum % of pirwise identity between sequences of a contig (80-100; default 95)
\t-prop_maj [floating] Proportion of sequences that must have the same base at a site to accept it as a consensus (0.5-1, default 0.66)


PIPE3 SPECIFIC PARAMETERS
\t-rm [0/1] (1: Run RepeatMasker on the sequences with primers (not available for windows))
\t-rm_path = [Full path to RepeatMasker executables]
\t-rm_lib [RepeatMasker library]
\t-check_contamination [0/1] (1: check for contamination by blasting sequences against the nt database)
\t-db_evalue [0, 1e-5] (e-value for testing contamination; local blast only)
\t-pcr_min [40,10000] Minimum PCR Product size
\t-pcr_max [40,10000] Maximum PCR Product size
\t-pcr_step [20,10000] PCR Product size interval
\t-PRIMER_GC_CLAMP [0,20]
\t-PRIMER_OPT_SIZE [1,50]
\t-PRIMER_MIN_SIZE [1,50]
\t-PRIMER_MAX_SIZE [1,50]
\t-PRIMER_OPT_TM [1,100]
\t-PRIMER_MIN_TM [1,100]
\t-PRIMER_MAX_TM [1,100]
\t-PRIMER_MAX_DIFF_TM [1,100]
\t-PRIMER_MIN_GC [1,100]
\t-PRIMER_OPT_GC_PERCENT [1,100]
\t-PRIMER_MAX_GC [1,100]
\t-PRIMER_SELF_ANY [1,100]
\t-PRIMER_SELF_END [1,100]
\t-PRIMER_MAX_POLY_X [1,10]
\t-PRIMER_NUM_RETURN [1,10]\n\n";


#\t-qual [0/1] 1 if sequence quality check should be done 
#\t-qual_file [name of the quality file including full path]
#\t-mean_qs  [integer betwee 0 and 40] (trimming sequence extremities if the mean quality score of a 10 bp window is lower than mean_qs)
#\t-min_qs [0,40] (Primers cannot include a base if its quality score is lower than min_qs (only necessary if qual is set to 1 ))

}

#########################################################
sub read_set_file_for_command_line
{
my ($param, $set_file, $param_list, $set_file_default) = @_;

read_set_file_for_GUI($param, $set_file, $set_file_default);

my $reply = print_param_menu($param_list, $param);

while ( not ($reply eq "\n"))
{
    chomp $reply;
    {
      if (exists $$param_list{$reply})
      {
        print "$$param_list{$reply}[1]\n";
        $$param{$$param_list{$reply}[0]} = <STDIN>;
        chomp $$param{$$param_list{$reply}[0]};
	    }
	    else
	    {
        print "$reply is not a menu option\n";
	    }
    }
    $reply = print_param_menu($param_list, $param);
}
print_new_set_file($param, $set_file);

}

#######################################################
sub read_set_file_galaxy
{
  my ($param, $set_file) = @_;

  unless(open('SET', $set_file))
  {
    print "Cannot open $set_file\n";
  }
 
my @bad_tag = ();
  while(my $line = <SET>)
  {
    $line =~  s/#.*//;
    if ($line =~  /\w/) # other than comment line
    {
	$line =~ s/\s*$//; # eliminate any spaces from the end of the line
	$line =~ s/^\s*//;
	$line =~ s/\s*=\s*/=/;
	my @line = split('=', $line);
	if(exists $$param{$line[0]})
	{
	    if (exists $line[1])
	    {
		$$param{$line[0]} = $line[1];
	    }
	    else
	    {
		$$param{$line[0]} = '';
	    }
	}
    }
  }
  close SET;

}
#######################################################################

sub run_RepeatMasker
{
 my ($input_file, $rm_path, $rm_lib, $num_threads, $tempfolder, $debug, $fh) = @_;
 my $rm_tbl = $input_file;
 $rm_tbl =~ s/.*[\\\/]//;
 $rm_tbl = $tempfolder.$rm_tbl.'.out';
my $screen = $tempfolder.'RM_message_on_screen.out';
 my $rmc = $rm_path.'RepeatMasker -nolow -q -pa '.$num_threads.' -spec "'.$rm_lib.'" -dir "'.$tempfolder.'" "'.$input_file.'" >'.$screen;
 if ($debug)
 {
	print $fh "$rmc\n";
 }
system $rmc;
return $rm_tbl;
}


########################################################################
sub read_set_file_for_GUI
{
  my ($param, $set_file, $set_file_default) = @_;

  unless(open('SET', $set_file))
  {
    print "Cannot open $set_file\n";
  }
  while(my $line = <SET>)
  {
    chomp $line;
    $line =~ s/\s*=\s*/=/;
    $line =~ s/\s*$//;
    $line =~ s/^\s*//;

    my @line = split('=', $line);
    if(exists $$param{$line[0]})
    {
      $$param{$line[0]} = $line[1];
    }
  }
  close SET;

  my $bool = 1;
  foreach my $p (keys %$param)
  {
    if ($$param{$p} eq '')
    {
      print "$p is not defined in $set_file. $set_file_default is used for reading parameters\n";
      $bool = 0;
    }
  }
  unless($bool)
  {
    read_set_file_for_GUI($param, $set_file_default);
  }
}
#######################################################################

sub reverse_complement_loc
{
    my ($temp) = @_;
    my $revcomp = reverse $temp;
    $revcomp =~ tr/ATCGatcg/TAGCtagc/;
    return $revcomp;
}


#########################################################

sub trim_5_prime
{
my ($seq, $qual, $code, $window, $mean_qs) = @_;
my @q = split(' ', $qual);
unless(scalar @q == length $seq)
{
	print "The number of quality score and the number of bases does not match in $code\n";
exit;
}

for(my $i = 0; $i<= (scalar @q - $window); ++$i)
{
	my $sum = 0;
	for(my $j = $i; $j<$i+$window; ++$j)
	{
		$sum += $q[$j];
	}
	my $mean = $sum/$window;
	if($mean_qs <= $mean)
	{
#		print "$code $i\n$seq\n$qual\n";
		$qual =~ s/^([0-9]+\s+){$i}//;
		substr($seq, 0, $i, '');
#		print "$code $i\n$seq\n$qual\n";
		$i = scalar @q;
	}
}
return ($seq, $qual);
}

############################################
sub trim_3_prime
{
my ($seq, $qual, $code, $window, $mean_qs) = @_;
my @q = split(' ', $qual);
unless(scalar @q == length $seq)
{
	print "The number of quality score and the number of bases does not match in $code\n";
	exit;
}

for(my $i = scalar@q-1; $i>= $window-1; --$i)
{
	my $sum = 0;
	for(my $j = $i; $j>$i-$window; --$j)
	{
		$sum += $q[$j];
	}
	my $mean = $sum/$window;
	if($mean_qs <= $mean)
	{
		my $cut = scalar@q - $i -1;
#		print "$code $cut\n$$seqref\n$$qualref\n";
		$qual =~ s/([0-9]+\s+){$cut}$//;
		substr($seq, $i+1, $cut, '');
#		print "$code $cut\n$$seqref\n$$qualref\n";
		$i = -1;
	}
}
return ($seq, $qual);
}

###########################################################

sub txt_into_hachage_1
{
# reads each line into a hachache. First element before first ';' is the code all other in the lines are values, separates element of the line by ;
    my (@data) = @_;
    my %hach = ();

foreach my $ligne (@data)
{
    
#    $ligne =~ s/\s+/\s/g;
    chomp $ligne;
	$ligne =~ s/\s*$//;
    $ligne =~ /.+?;/;
    my $code = $&;
    $ligne =~ s/$code//;
    $code =~ s/;//;
    $hach{$code} = $ligne;
}

return %hach;

}


###############################################

sub write_hash_to_fasta
{
	my ($hashref, $outfile) = @_;
	
	open(OUT, ">$outfile") or die "Cannot open $outfile file\n";
	foreach my $code (sort keys %$hashref)
	{
		print OUT ">$code\n";
		print OUT cut_up_fasta_line($$hashref{$code}, '100');
	}
	close OUT;
}


1;
