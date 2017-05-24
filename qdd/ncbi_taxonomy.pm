 
  package ncbi_taxonomy;
 
  use strict;
  use warnings;
  use DBI;

  our $VERSION = '1.00';
  
  use base 'Exporter';
  
  our @EXPORT = qw(init_taxo approx_taxname taxname_from_taxid taxid_from_taxname tax_info tax_lineage tax_lineage_ranked tax_lineage_full tax_division tax_rank);

my $dbd;
my $database;
our $debug = 0;
our $ROOT = '';
our $approx = 0;

sub init_taxo{
	$database = shift;
}
sub approx_taxname{
	return $approx;
}
sub connect_to_database{
	$database = 'ncbi_tax.db' unless defined ($database);
	eval {
		$dbd=DBI->connect("DBI:SQLite:dbname=$database",'','',{RaiseError=>1}) or die "Connection impossible : $@";
	};
	if ($@) {
		die "Connection à ncbi_tax.db impossible : $@";
	}
}

sub dbQuery {
	my $sql = shift;
	if (!defined $dbd) {connect_to_database()}
	my $sth = $dbd->prepare($sql);
  	$sth->execute();
  	return $sth;
	#while (my $ref = $sth->fetchrow_hashref()) {
}

sub taxname_from_taxid{
	my $id = shift;
	if (!$id>0) {return '';}
	$debug && print "taxname_from_taxid Taxid=$id\n"; 
	my $result = dbQuery( "select * from tax_names where tax_id=$id and class='scientific name'");
	my $t = $result->fetchrow_hashref();
	$debug && print "taxname_from_taxid t->{'name'}=$t->{'name'}\n"; 
	my $res = ($t->{'name'}?$t->{'name'}:'no tax name');
	$debug && !defined($res) && print "NOT DEFINED RES\n";
	return $res;
}

sub taxid_from_taxname{
	my $name = shift;
	$debug && print "...searching for '$name'\n";
	if ($name eq '') {
		return 'NULL';}
	$name =~ s/\'/\\'/g;
	my $result = dbQuery( "select * from tax_names where name='$name'");
	my $t = $result->fetchrow_hashref();
	if (defined $t->{'tax_id'} && $t->{'tax_id'} > 1) {
		return $t->{'tax_id'};}
	$result = dbQuery( "select * from tax_names where name like '$name %'");
	$t = $result->fetchrow_hashref();
	if (defined $t->{'tax_id'} && $t->{'tax_id'} > 1) {
		$approx = "$name *";
		return $t->{'tax_id'};}
	$result = dbQuery( "select * from tax_names where name like '% $name'");
	$t = $result->fetchrow_hashref();
	if (defined $t->{'tax_id'} && $t->{'tax_id'} > 1) {
		$approx = "* $name";
		return $t->{'tax_id'};}
	$result = dbQuery( "select * from tax_names where name like '% $name %'");
	$t = $result->fetchrow_hashref();
	if (defined $t->{'tax_id'} && $t->{'tax_id'} > 1) {
		$approx = "* $name *";
		return $t->{'tax_id'};}
	$result = dbQuery( "select * from tax_names where name like '$name%'");
	$t = $result->fetchrow_hashref();
	if (defined $t->{'tax_id'} && $t->{'tax_id'} > 1) {
		$approx = "$name*";
		return $t->{'tax_id'};}
	$result = dbQuery( "select * from tax_names where name like '%".$name."'");
	$t = $result->fetchrow_hashref();
	if (defined $t->{'tax_id'} && $t->{'tax_id'} > 1) {
		$approx = "*$name";
		return $t->{'tax_id'};}
	$result = dbQuery( "select * from tax_names where name like '%".$name."%'");
	$t = $result->fetchrow_hashref();
	if (defined $t->{'tax_id'} && $t->{'tax_id'} > 1) {
		$approx = "*$name*";
		return $t->{'tax_id'};}
	$name =~ s/\s+\S+$//;
	if ($name) {
		$approx = "$name";
		return taxid_from_taxname($name);
	}
	return 0;
}
sub tax_info{
	my $id = shift;
	if (!$id>0) {return '';}
	my $result = dbQuery( "select * from tax_nodes join tax_genetic_codes using(gen_code) where tax_id=$id");
	my $t = $result->fetchrow_hashref();
	return "Rank: ".$t->{'rank'}." - Genetic Code: ".$t->{'name'}." - taxid: $id\n".
	'kingdom: '.tax_rank('superkingdom',$id)." - phylum: ".tax_rank('phylum',$id).
	" - class: ".tax_rank('class',$id)." - order: ".tax_rank('order',$id)."\n".
	tax_lineage($id)."\n";
}

sub tax_lineage{
	my $id = shift;
	if ($id <= 1) {return $ROOT;}
	$debug && print "tax_lineage Taxid=$id\n"; 
	my $result = dbQuery( "select * from tax_nodes where tax_id=$id");
	my $t = $result->fetchrow_hashref();
	if ($id == $t->{'parent_id'}) {return;}
	my $print = ($t->{'hidden'}?'':taxname_from_taxid($id).'; ');
	$debug && print "tax_lineage print=$print\n"; 
	$debug && !defined($print) && print "tax_lineage NOT DEFINED print!!!!!!\n";
	my $res = (defined($t->{'parent_id'}) && $t->{'parent_id'} > 0?tax_lineage($t->{'parent_id'}):'');
	$debug && !defined($res) && print "tax_lineage NOT DEFINED res t->{'parent_id'}=$t->{'parent_id'}!!!!!!\n";
	return $res.$print;
}
sub tax_lineage_ranked{
	my $id = shift;
	if ($id <= 1) {return $ROOT;}
	$debug && print "tax_lineage Taxid=$id\n"; 
	my $result = dbQuery( "select * from tax_nodes where tax_id=$id");
	my $t = $result->fetchrow_hashref();
	if ($id == $t->{'parent_id'}) {return;}
	my $print = ($t->{'hidden'}?'':'['.$t->{'rank'}.']'.taxname_from_taxid($id).'; ');
	$debug && print "tax_lineage print=$print\n"; 
	$debug && !defined($print) && print "tax_lineage NOT DEFINED print!!!!!!\n";
	my $res = (defined($t->{'parent_id'}) && $t->{'parent_id'} > 0?tax_lineage_ranked($t->{'parent_id'}):'');
	$debug && !defined($res) && print "tax_lineage NOT DEFINED res t->{'parent_id'}=$t->{'parent_id'}!!!!!!\n";
	return $res.$print;
}

sub tax_lineage_full{
	my $id = shift;
	if ($id <= 1) {return $ROOT;}
	$debug && print "tax_lineage Taxid=$id\n"; 
	my $result = dbQuery( "select * from tax_nodes where tax_id=$id");
	my $t = $result->fetchrow_hashref();
	if ($id == $t->{'parent_id'}) {return;}
	my $print = taxname_from_taxid($id).'; ';
	$debug && print "tax_lineage print=$print\n"; 
	$debug && !defined($print) && print "tax_lineage NOT DEFINED print!!!!!!\n";
	my $res = (defined($t->{'parent_id'}) && $t->{'parent_id'} > 0?tax_lineage_full($t->{'parent_id'}):'');
	$debug && !defined($res) && print "tax_lineage NOT DEFINED res t->{'parent_id'}=$t->{'parent_id'}!!!!!!\n";
	return $res.$print;
}

sub tax_division{
	my $id = shift;
	if ($id <= 1) {return;}
	my $result = dbQuery( "select * from tax_nodes join tax_division using(div_id) where tax_id=$id");
	my $t = $result->fetchrow_hashref();
	if ($id == $t->{'parent_id'}) {return;}
	my $print = ($t->{'inh_div_id'}?'':$t->{'name'});
	if ($print eq '') {
		return tax_division($t->{'parent_id'});}
	else {
		return $print;
	}
}
sub tax_rank{
	my $rank = shift;
	my $id = shift;
	if ($id eq 'NULL' || $id <= 1) {return;}
	my $result = dbQuery( "select * from tax_nodes where tax_id=$id");
	my $t = $result->fetchrow_hashref();
	if ($id == $t->{'parent_id'}) {return;}
	my $print = ($t->{'rank'} eq $rank?taxname_from_taxid($id):'');
	if ($print eq '') {
		return tax_rank($rank,$t->{'parent_id'});}
	else {
		return $print;
	}
}


  1;

