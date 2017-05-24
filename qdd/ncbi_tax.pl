#!/usr/bin/perl -w
# Pascal Hingamp 2013 pascal.hingamp@univ-amu.fr

use lib './';
use strict;
use warnings;
use ncbi_taxonomy;


my $taxid = $ARGV[0] if defined $ARGV[0];
my $db_path = $ARGV[1] if defined $ARGV[1];

my $USAGE = "USAGE:\n./ncbi_tax 9606\n./ncbi_tax \"homo sapiens\"\n";

unless (defined $taxid && $taxid =~ /\S/){
	die $USAGE;
}

#if ncbi_tax.db not in usual place/different name, use following optional line:
my $db_name = $db_path.'ncbi_tax.db';
init_taxo($db_name);

unless ($taxid =~ /^\d+$/){
	$taxid = taxid_from_taxname($taxid);
}

unless ($taxid =~ /^\d+$/){
	die "Did not understand taxonomy term to search: '$taxid'\n$USAGE";
}
print "NCBI taxid:\n===========\n$taxid\n\n".(approx_taxname()?"WARNING: exact term not found, found approximate match against '".approx_taxname()."'!!!\n\n":'')."NCBI abbreviated:\n=================\n".tax_lineage($taxid)."\n\nNCBI abbreviated with NCBI ranks:\n=================================\n".tax_lineage_ranked($taxid)."\n\nNCBI full lineage (with hidden nodes):\n======================================\n".tax_lineage_full($taxid)."\n";

