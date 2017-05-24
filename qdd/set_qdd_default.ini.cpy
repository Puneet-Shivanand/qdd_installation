# These parameters are read as a default, and can be overwriten by command line tag-value pairs (e.g. perl pipe1.pl -adapter 0)

#GENEAL PARAMETERS:
#[0/1] 1 for funning QDD from Galaxy, 0 for running it from terminal
galaxy =0
#operating system [linux/win]
syst = linux
# Full path to blast executables (including the bin folder). If the folder is in your path it can be left empty. (e.g. C:\Program Files\NCBI\blast-2.2.25+\bin or /home/pri/galaxy-dist/tools/qdd/ncbi-blast-2.2.27+/bin/) If the folder is in your path it can be left empty.
blast_path=/home/puneet/qdd/tools/ncbi-blast-2.6.0+/bin/
# Full path to clustalw executables. If the folder is in your path it can be left empty. (e.g. C:\CLUSTALW2 or /home/pri/galaxy-dist/tools/qdd/clustalw-2.0.10-linux-i386-libcppstatic/)
clust_path =/home/puneet/qdd/tools/clustalw-2.1/src/ 
#  Full path to Primer3_core executable. (e.g. C:\Primer3-release-2.3.6\ or C:\primer3-1.1.4-WINXP\bin or /usr/bin/ /usr/local/primer3-2.3.6/src/)
primer3_path = /home/puneet/qdd/tools/primer3-2.3.7/src/
# Primer3 version [1/2] 1 for primer3-1.xxx, 2 for primer3-2.xxx
primer3_version = 2
# Full path to qdd scripts. (e.g. /home/qdd/galaxy-dist/tools/qdd/ OR D:\QDD2.3\)
qdd_folder = /home/puneet/qdd/qdd/
# Output folder name with full path. Must be created before running qdd. If not specified output files are written in the current working directory
out_folder = /home/puneet/qdd/qdd/qdd_output/
#[0/1] (1 for deleting temporary files after the run)
del_files =1
# string for naming outfiles. All outfile name starts with this string. If empty, the input filename is used for naming outfiles 
outfile_string =
#[0/1] (1 for printing out supplementary information, only needed for debugging)
debug = 0
# name of the local database including full path  (e.g. /usr/local/nt/nt or D:\blastdb\nt) Only needed if local BLAST is used for contamination check 
blastdb = /home/puneet/qdd/tools/NCBI_nt/nt
#number of threads for BLAST (if unsure, use 1)
num_threads = 10
# [0/1] 1:run local blast for contamination check, 0:run remote blast for contamination check
local_blast = 1
# Input sequences type for pipe1. 1 if sequences has been assembled (contigs, scaffolds, chromosomes), 0 if they are short sequencing reads. Same parameter is used for pipe3 to check the distance between nearest neigbours
contig = 0


#PIPE1 SPECIFIC PARAMETERS
#input file is in fastq format; [0/1] (1 for fatsq, 0 for fasta)
fastq = 0 
# if extracting microsatellites from contigs, get flank_length bp of flanking region on both sides of the microsatellite
flank_length = 200
#[0/1] (1 for running adapter/vector clipping step); Not applicable if contig =1
adapter= 0
#[integer] (sequences shorter then length_limit are eliminated)
length_limit= 80
#[fasta file with adaters] (can be empty if adapter=0) ; Not applicable if contig =1
adapter_file =


#PIPE2 SPECIFIC PARAMETERS
# [0/1] Make consensus sequences (YES=1/NO=0)
make_cons=1
# [integer] Minimum % of pirwise identity between sequences of a contig (80-100)
ident_limit =95
# [floating] Proportion of sequences that must have the same base at a site to accept it as a consensus (0.5-1)
prop_maj =0.66


#PIPE3 SPECIFIC PARAMETERS
#Minimum size of PCR product
pcr_min = 90
#Maximum size of PCR product
pcr_max = 300
#PCR Product size interval
pcr_step = 50
PRIMER_GC_CLAMP = 0
PRIMER_MIN_SIZE = 18
PRIMER_MAX_SIZE = 27
PRIMER_OPT_SIZE = 20
PRIMER_OPT_TM = 60.0
PRIMER_MIN_TM = 57.0
PRIMER_MAX_TM = 63.0
PRIMER_MAX_DIFF_TM = 10.0
PRIMER_MIN_GC = 20.0
PRIMER_OPT_GC_PERCENT = 50.0
PRIMER_MAX_GC = 80.0
PRIMER_SELF_ANY = 8.0
PRIMER_SELF_END = 3.0
PRIMER_MAX_POLY_X = 3
PRIMER_NUM_RETURN = 3

#PIPE4 SPECIFIC PARAMETERS
# [0/1] 1: Run RepeatMasker on the sequences with primers (not available for windows). Primer table is completed with info on RM hits
rm = 0
# Full path to RepeatMasker executables. (e.g. /usr/local/RepeatMasker/) If the folder is in your path it can be left empty. 
rm_path = /usr/local/RepeatMasker/
# RepeatMasker library
rm_lib = eukaryota
# [0/1] 1: check for contamination by blasting sequences against the nt database of NCBI (option remote or local)
check_contamination = 0
# e-value for local blast db for testing contamination
db_evalue = 1e-20

#QDD.pl SPECIFIC PARAMETERS
#[0/1] (1 for sorting sequences according to tags)
tag = 0
# fasta file with tag sequences 
tag_file =
#[0/1] (1 for running pipe1, pipe2, pipe3 in one go)
run_all = 1
# folder that contains all input files for batch submission; must not contain other files
input_folder =
