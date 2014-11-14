#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Copy;
use File::Basename;

# perl uneak_analysis_pipeline.pl -i ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/GQ03122013_5_fastq.txt.gz -n MPB_MALE_UNEAK_GBS -k ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/mpb_barcodes.txt -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_GBS_UNEAK_ANALYSIS

my ($fastq_infile, $project_name, $reference_genome_infile, $fastq_barcodes_infile, $restriction_enzymes, $max_num_barcode_reads, 
$min_num_present_ufastq_to_tag, $min_num_present_merge_multiple_tag, $max_num_tags_tbt, $min_minor_allele_freq_tag_to_snp, $min_minor_allele_freq_gbs_hap_map, 
$min_minor_allele_count, $max_num_sites, $start_chromosome, $end_chromosome, $genotypic_mismat_thres, $min_site_coverage, $tassel_num_ram, $bwa_num_cpu, $uneak_project_dir);

GetOptions(
	'i=s'    => \$fastq_infile,
	'n=s'    => \$project_name,
	'g=s'    => \$reference_genome_infile,
	'k=s'    => \$fastq_barcodes_infile,
	'e=s'    => \$restriction_enzymes,
	's=s'    => \$max_num_barcode_reads,
	'c=s'    => \$min_num_present_ufastq_to_tag,
	'u=s'    => \$min_num_present_merge_multiple_tag,
	't=s'    => \$max_num_tags_tbt,
	'f=s'    => \$min_minor_allele_freq_tag_to_snp,
	'a=s'    => \$min_minor_allele_freq_gbs_hap_map,
	'r=s'    => \$min_minor_allele_count,
	'x=s'    => \$max_num_sites,
	'h=s'    => \$start_chromosome,
	'v=s'    => \$end_chromosome,
	'y=s'    => \$genotypic_mismat_thres,
	'z=s'    => \$min_site_coverage,
	'm=s'    => \$tassel_num_ram,
	'p=s'    => \$bwa_num_cpu,
	'o=s'    => \$uneak_project_dir,
);

usage() unless (
	defined $fastq_infile
# 	and defined $project_name
# 	and defined $reference_genome_infile
#  	and defined $fastq_barcodes_infile
	and defined $uneak_project_dir
);

$restriction_enzymes = 'PstI-MspI' unless defined $restriction_enzymes;
$max_num_barcode_reads = 250000000 unless defined $max_num_barcode_reads;
# $max_num_tags_tbt = 200000000 unless defined $max_num_tags_tbt;
$min_num_present_ufastq_to_tag = 1 unless defined $min_num_present_ufastq_to_tag;
# $min_num_present_merge_multiple_tag = 5 unless defined $min_num_present_merge_multiple_tag;
$tassel_num_ram = 12 unless defined $tassel_num_ram;
# $bwa_num_cpu = 7 unless defined $bwa_num_cpu;
# $max_num_sites = 500000 unless defined $max_num_sites;
# $min_minor_allele_freq_tag_to_snp = 0.02 unless defined $min_minor_allele_freq_tag_to_snp;
# $min_minor_allele_freq_gbs_hap_map = 0.01 unless defined $min_minor_allele_freq_gbs_hap_map;
# $min_minor_allele_count = 100000 unless defined $min_minor_allele_count;
# $start_chromosome = 1 unless defined $start_chromosome;
# $end_chromosome = 8188 unless defined $end_chromosome;
# $genotypic_mismat_thres = 0.1 unless defined $genotypic_mismat_thres;
# $min_site_coverage = 0.2 unless defined $min_site_coverage;

my ($run_pipeline, $reformat_refgen_to_bwatassel, $bwa);
$run_pipeline		= '/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl';
$bwa 			= '/usr/local/bin/bwa';

sub usage {

die <<"USAGE";


Usage: $0 -i fastq_infile -k fastq_barcodes_infile -e restriction_enzymes -s max_num_barcode_reads -c min_num_present_tag -m tassel_num_ram -p bwa_num_cpu -o output_dir

Description - 

OPTIONS:

	-i fastq_infile - 
	
	-k fastq_barcodes_infile -
	
	-e restriction_enzymes -

	-s max_num_barcode_reads - 

	-c min_num_present_tag -

	-m tassel_num_ram - 
	
	-p bwa_num_cpu -

	-o output_dir -

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $uneak_project_dir){
      mkdir($uneak_project_dir, 0777) or die "Can't make directory: $!";
}


# ./Illumina/ (original raw data files, one file per flowcell lane)
# Create the Illumina output directory if it doesn't already exist.
my $illumina_output_dir = join('/', $uneak_project_dir, "Illumina");
unless(-d $illumina_output_dir){
      mkdir($illumina_output_dir, 0777) or die "Can't make directory: $!";
}

# Copying fastq_filename to the Illumina directory 
my $fastq_filename = basename($fastq_infile);
warn "Copying $fastq_filename to the Illumina directory.....\n";
my $uneak_fastq_infile = join('/',$illumina_output_dir, $fastq_filename);
unless(-s $uneak_fastq_infile){
	copy($fastq_infile, $uneak_fastq_infile) or die "Copying $fastq_infile to $uneak_fastq_infile failed: $!";
}
# ./key/ (Barcode key file of original raw data files)
# Create the fastq barcodes key output directory if it doesn't already exist.
my $barcodes_key_output_dir = join('/', $uneak_project_dir, "key");
unless(-d $barcodes_key_output_dir){
      mkdir($barcodes_key_output_dir, 0777) or die "Can't make directory: $!";
}

# Copying fastq_barcodes_infile to the barcodes key directory
my $fastq_barcodes_filename = basename($fastq_barcodes_infile);
warn "Copying $fastq_barcodes_infile to the barcodes key directory.....\n";
my $uneak_fastq_barcodes_infile = join('/', $barcodes_key_output_dir, $fastq_barcodes_filename);
unless(-s $uneak_fastq_barcodes_infile){
	copy($fastq_barcodes_infile, $uneak_fastq_barcodes_infile) or die "Copying $fastq_barcodes_infile to $uneak_fastq_barcodes_infile failed: $!";
}

# ./tagCounts/ (for output from UQseqToTagCountPlugin OR UFastqToTagCountPlugin OR UMergeTaxaTagCountPlugin)
# Create the tag counts output directory if it doesn't already exist.
my $tag_counts_output_dir = join('/', $uneak_project_dir, "tagCounts");
unless(-d $tag_counts_output_dir){
      mkdir($tag_counts_output_dir, 0777) or die "Can't make directory: $!";
}

ufastq_to_tag_counts($uneak_project_dir, $project_name, $fastq_barcodes_infile, $restriction_enzymes, $max_num_barcode_reads, $min_num_present_ufastq_to_tag, 
	$tassel_num_ram, $tag_counts_output_dir);
	
# ./mergedTagCounts/ (for output fromUMergeTaxaTagCountPlugin)
# Create the merged tag counts output directory if it doesn't already exist.
my $merged_tag_counts_output_dir = join('/', $uneak_project_dir, "mergedTagCounts");
unless(-d $merged_tag_counts_output_dir){
      mkdir($merged_tag_counts_output_dir, 0777) or die "Can't make directory: $!";
}

# ./tagPair/ (for output fromUTagCountToTagPairPlugin)
# Create the tag pair output directory if it doesn't already exist.
my $tag_pair_output_dir = join('/', $uneak_project_dir, "tagPair");
unless(-d $tag_pair_output_dir){
      mkdir($tag_pair_output_dir, 0777) or die "Can't make directory: $!";
}

# ./tagsByTaxa/ (for output fromUTBTToMapInfoPlugin)
# Create the tags by taxa output directory if it doesn't already exist.
my $tags_by_taxa_output_dir = join('/', $uneak_project_dir, "tagsByTaxa");
unless(-d $tags_by_taxa_output_dir){
      mkdir($tags_by_taxa_output_dir, 0777) or die "Can't make directory: $!";
}

# ./mapInfo/ (for output from UTBTToMapInfoPlugin)
# Create the map info output directory if it doesn't already exist.
my $map_info_output_dir = join('/', $uneak_project_dir, "mapInfo");
unless(-d $map_info_output_dir){
      mkdir($map_info_output_dir, 0777) or die "Can't make directory: $!";
} 

# ./hapMap/ (for output from TagsToSNPByAlignmentPlugin)
# Create the hap map output directory if it doesn't already exist.
my $hap_map_output_dir = join('/', $uneak_project_dir, "hapMap");
unless(-d $hap_map_output_dir){
      mkdir($hap_map_output_dir, 0777) or die "Can't make directory: $!";
}


# uneak_create_working_dir($uneak_project_dir);
# /home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -fork1 -UCreatWorkingDirPlugin -w M:/UNEAK/ -endPlugin -runfork1
# sub uneak_create_working_dir{
# 
#         my $uneak_working_dir = shift;
#         die "Error lost the UNEAK working directory" unless defined $uneak_working_dir;
#         
# 	warn "Generating UNEAK sub-directories using the UCreatWorkingDirPlugin.....\n";
# 	my $uneakCreateWorkingDirCmd  = "$run_pipeline -fork1 -UCreatWorkingDirPlugin -w $uneak_working_dir -endPlugin -runfork1";
# 	warn $uneakCreateWorkingDirCmd . "\n\n";
# 	system($uneakCreateWorkingDirCmd) == 0 or die "Error calling $uneakCreateWorkingDirCmd: $?";
# }

#run ufastqtotagcount plugin from the current directory (-w) using PstI-MspI combo (-e), retaining all reads with atleast 1 copy (-c). Expected number of good barcoded reads to be less than 250 million (-s). 
#don't change the -c option setting here since you'll be able to filter by read depth at the umergetaxatagcount step later on.
# run_pipeline.pl -Xmx4g -fork1 -UFastqToTagCountPlugin -w ./ -e PstI-MspI -s 250000000 -c 1 -endPlugin -runfork1 > ./UtagCount.log
sub ufastq_to_tag_counts{

        my $uneak_working_dir = shift;
        die "Error lost the UNEAK working directory" unless defined $uneak_working_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $fastq_barcodes_infile = shift;
        die "Error lost the key file listing barcodes for each sample and plate layout" unless defined $fastq_barcodes_infile;

        my $restriction_enzymes = shift;
        die "Error lost the enzyme(s) used to create the GBS library" unless defined $restriction_enzymes;
        
        my $max_num_barcode_reads = shift;
        die "Error lost the maximum number of good barcoded reads per lane" unless defined $max_num_barcode_reads;
        
        my $min_num_present_tag = shift;
        die "Error lost the minimum number of times a tag must be present to be output" unless defined $min_num_present_tag;
        
        my $tassel_num_ram = shift;
        die "Error lost the amount of tassel memory" unless defined $tassel_num_ram;
        
        my $tag_count_output_dir = shift;
        die "Error lost the output directory to contain output .cnt (tag count) files, one per input FASTQ file" unless defined $tag_count_output_dir;
        
        my $uneak_fastq_to_tag_counts_log_outfile = join('/', $tag_count_output_dir, join("_", $project_name, "UFastqToTagCounts.log"));
	warn "Generating the UNEAK tag counts list file using the UFastqToTagCountPlugin.....\n";
	my $fastqToTagCountsCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -UFastqToTagCountPlugin -w $uneak_working_dir -e $restriction_enzymes -s $max_num_barcode_reads -c $min_num_present_tag -endPlugin -runfork1 > $uneak_fastq_to_tag_counts_log_outfile");
	warn $fastqToTagCountsCmd . "\n\n";
	system($fastqToTagCountsCmd) == 0 or die "Error calling $fastqToTagCountsCmd: $?";
}

#run binarytotext plugin to make the ufastqtotagcount result from above human readable
#change GQ25032013_5 in the following line to the actual name of your files (should be GQ25032013_yourlanenumber 
# run_pipeline.pl -Xmx4g -fork1 -BinaryToTextPlugin -i ./tagCounts/GQ25032013_7.cnt -o ./tagCounts/GQ25032013_7_cnt.txt -t TagCounts -endPlugin -runfork1 > ./bintotex.log

#run umergetaxatagcounts plugin from the current directory (-w), using a read depth cutoff of 5 (-c). You can change the cutoff to whatever you like.
# run_pipeline.pl -Xmx4g -fork1 -UMergeTaxaTagCountPlugin -w ./ -c 5 -endPlugin -runfork1 > UmergeTagCount.log
sub umerge_taxa_tag_counts{

        my $uneak_working_dir = shift;
        die "Error lost the UNEAK working directory" unless defined $uneak_working_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $fastq_barcodes_infile = shift;
        die "Error lost the key file listing barcodes for each sample and plate layout" unless defined $fastq_barcodes_infile;

        my $restriction_enzymes = shift;
        die "Error lost the enzyme(s) used to create the GBS library" unless defined $restriction_enzymes;
        
        my $max_num_barcode_reads = shift;
        die "Error lost the maximum number of good barcoded reads per lane" unless defined $max_num_barcode_reads;
        
        my $min_num_present_tag = shift;
        die "Error lost the minimum number of times a tag must be present to be output" unless defined $min_num_present_tag;
        
        my $tassel_num_ram = shift;
        die "Error lost the amount of tassel memory" unless defined $tassel_num_ram;
        
        my $tag_count_output_dir = shift;
        die "Error lost the output directory to contain output .cnt (tag count) files, one per input FASTQ file" unless defined $tag_count_output_dir;
        
        my $uneak_fastq_to_tag_counts_log_outfile = join('/', $tag_count_output_dir, join("_", $project_name, "UFastqToTagCounts.log"));
	warn "Generating the UNEAK tag counts list file using the UFastqToTagCountPlugin.....\n";
	my $fastqToTagCountsCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -UFastqToTagCountPlugin -w $uneak_working_dir -e $restriction_enzymes -s $max_num_barcode_reads -c $min_num_present_tag -endPlugin -runfork1 > $uneak_fastq_to_tag_counts_log_outfile");
	warn $fastqToTagCountsCmd . "\n\n";
	system($fastqToTagCountsCmd) == 0 or die "Error calling $fastqToTagCountsCmd: $?";
}

#run binarytotext plugin to make the umergetaxatagcount result from above human readable
# run_pipeline.pl -Xmx4g -fork1 -BinaryToTextPlugin -i ./mergedTagCounts/mergeAll.cnt -o ./mergedTagCounts/mergeAll_cnt.txt -t TagCounts -endPlugin -runfork1 >> ./bintotex.log

#run utagcounttotagpair plugin from current directory (-w), with a sequencing error threshold of 3% (this should not be higher than 5% based on known Illumina 
#sequencing error rates.
# run_pipeline.pl -Xmx4g -fork1 -UTagCountToTagPairPlugin -w ./ -e 0.03 -endPlugin -runfork1 > utagtotp.log

#run utagpairtoTBTplugin from current directory (-w).
# run_pipeline.pl -Xmx4g -fork1 -UTagPairToTBTPlugin -w ./ -endPlugin -runfork1 > utagtotbt.log

#run uTBTtomapinfoplugin from current directory (-w).
# run_pipeline.pl -Xmx4g -fork1 -UTBTToMapInfoPlugin -w ./ -endPlugin -runfork1 > utbttomapinfo.log

#run umapinfotohapmapplugin from current directory (-w), with a minor allele frequency of 1%, a major allele frequency of 99%, minimum call rate of 0, maximum 
#callrate of 1.
# run_pipeline.pl -Xmx4g -fork1 -UMapInfoToHapMapPlugin -w ./ -mnMAF 0.01 -mxMAF 0.99 -mnC 0 -mxC 1 -endPlugin -runfork1 > umapinfotohapmap.log

