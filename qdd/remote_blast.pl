  use warnings;
  use strict;
#  use subprogramQDD;
  use Bio::Tools::Run::RemoteBlast;
  use Bio::DB::EUtilities;
  use Bio::Perl;


my $input_fas = "d:/test with space/pipe3_1308054947/nuc_pipe3_targets.fas";
#my $input_fas = "d:/test with space/pipe3_1308054947/sample_pipe3_target.fas";
my $summary_file = 'd:/pipe_clean/QDD2/summ.csv';
my $tax_file = 'd:/pipe_clean/QDD2/tax_lin.txt';
my $prog = 'blastn';
my $db = 'nr';
my $e_val = '1e-20';
my @taxids = ();

if ( scalar@ARGV == 4 )
{
  $input_fas = $ARGV[0];
  $summary_file = $ARGV[1];
  $tax_file = $ARGV[2];
  $e_val = $ARGV[3];
}
else
{
  print "Direct parameters are used in remote_blast.pl\n";
}

remote_blast($input_fas, $summary_file, $prog, $db, $e_val, \@taxids);
get_tax_lineage($tax_file, \@taxids);


exit;
###############################################################
sub remote_blast
{
my ($input_fas, $summary_file, $prog, $db, $e_val, $taxids) = @_;



 unless (open(SUM, ">$summary_file"))
  {
    print "Cannot open $summary_file\n";
  }
print SUM "Sequence_code;Accession(best_hit);Name(best_hit);Description(best_hit);E_value(best_hit);Score(best_hit);TaxId(best_hit)\n";

  my @params = ( '-prog' => $prog,
         '-data' => $db,
         '-expect' => $e_val,
         '-readmethod' => 'SearchIO');

  my $factory = Bio::Tools::Run::RemoteBlast->new(@params);
  my $v = 1; #$v is just to turn on and off the messages
  
#    print STDERR "Waiting for NCBI BLAST..." if( $v > 0 );
  print "waiting for NCBI BLAST..." if( $v > 0 );
  my $r = $factory->submit_blast($input_fas);

    my %taxid_list = ();
#    my %query_spec = (); # Query code as a key and taxId of the best hit as a value
    while ( my @rids = $factory->each_rid ) 
    {
      foreach my $rid ( @rids ) 
      {
        my $rc = $factory->retrieve_blast($rid);
        if( !ref($rc) ) 
        {
          if( $rc < 0 ) 
          {
            $factory->remove_rid($rid);
          }
#          print STDERR "." if ( $v > 0 );
          print "." if ( $v > 0 );
          sleep 5;
        } 
        else 
        {
          my $result = $rc->next_result();
          #save the output
          if (0)
          {
            my $filename = $result->query_name()."\.out";
            $factory->save_output($filename);
          }
          $factory->remove_rid($rid);
          print "\nQuery Name: ", $result->query_name();
   # take only the best hit and only one HSP
          my $hit = $result->next_hit;
          if ($hit)
          {
       # print info into file
            print SUM $result->query_name(), ";";
            print SUM $hit->accession, ";";
            print SUM $hit->name, ";";
            my $description = $hit->description;
            $description =~ s/;//g;
            print SUM $description, ";";
            print SUM $hit->significance, ";";
            print SUM $hit->raw_score, ";";
            
        # get gb file for sbj and extract tax id 
            my $seq_object_sbj = get_sequence('genbank', $hit->accession);
            my $taxId = $seq_object_sbj->species->ncbi_taxid();
            $taxid_list{$taxId} = '';
#            $query_spec{$result->query_name()} = $taxId;
            print SUM $taxId, "\n";
           }
        } # end else
      } # end foreach my $rid ( @rids ) 
    } # end while  (my @rids = $factory->each_rid)

#print Dumper(\%taxid_list);
close SUM;

@$taxids = keys %taxid_list;

}

######################################################

sub get_tax_lineage
{
my ($tax_file, $taxids) = @_;

my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                       -db      => 'taxonomy',
                             #          -rettype => 'genbank',
                             #          -email   => 'your@mail',
                                       -id      => $taxids);
 
$factory->get_Response(-file => $tax_file);
}
