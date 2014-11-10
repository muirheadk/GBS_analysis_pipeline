#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

# perl tassel_analysis_pipeline.pl -i ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_GBS_ANALYSIS/GBS_FASTQ_DIR -n MPB_MALE_GBS -k ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_GBS_ANALYSIS/mpb_barcodes.txt -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_GBS_ANALYSIS

my ($fastq_infile_dir, $project_name, $reference_genome_infile, $fastq_barcodes_infile, $restriction_enzymes, $max_num_barcode_reads, 
$min_num_present_fastq_to_tag, $min_num_present_merge_multiple_tag, $max_num_tags_tbt, $min_minor_allele_freq_tag_to_snp, $min_minor_allele_freq_gbs_hap_map, 
$min_minor_allele_count, $max_num_sites, $start_chromosome, $end_chromosome, $genotypic_mismat_thres, $min_site_coverage, $tassel_num_ram, $bwa_num_cpu, $output_dir);

GetOptions(
	'i=s'    => \$fastq_infile_dir,
	'n=s'    => \$project_name,
	'g=s'    => \$reference_genome_infile,
	'k=s'    => \$fastq_barcodes_infile,
	'e=s'    => \$restriction_enzymes,
	's=s'    => \$max_num_barcode_reads,
	'c=s'    => \$min_num_present_fastq_to_tag,
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
	'o=s'    => \$output_dir,
);

usage() unless (
	defined $fastq_infile_dir
	and defined $project_name
	and defined $reference_genome_infile
 	and defined $fastq_barcodes_infile
	and defined $output_dir
);

$restriction_enzymes = 'PstI-MspI' unless defined $restriction_enzymes;
$max_num_barcode_reads = 300000000 unless defined $max_num_barcode_reads;
$max_num_tags_tbt = 200000000 unless defined $max_num_tags_tbt;
$min_num_present_fastq_to_tag = 1 unless defined $min_num_present_fastq_to_tag;
$min_num_present_merge_multiple_tag = 5 unless defined $min_num_present_merge_multiple_tag;
$tassel_num_ram = 12 unless defined $tassel_num_ram;
$bwa_num_cpu = 7 unless defined $bwa_num_cpu;
$max_num_sites = 500000 unless defined $max_num_sites;
$min_minor_allele_freq_tag_to_snp = 0.02 unless defined $min_minor_allele_freq_tag_to_snp;
$min_minor_allele_freq_gbs_hap_map = 0.01 unless defined $min_minor_allele_freq_gbs_hap_map;
$min_minor_allele_count = 100000 unless defined $min_minor_allele_count;
$start_chromosome = 1 unless defined $start_chromosome;
$end_chromosome = 8188 unless defined $end_chromosome;
$genotypic_mismat_thres = 0.1 unless defined $genotypic_mismat_thres;
$min_site_coverage = 0.2 unless defined $min_site_coverage;

my ($run_pipeline, $reformat_refgen_to_bwatassel, $bwa);
$run_pipeline		= '/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl';
$bwa 			= '/usr/local/bin/bwa';

sub usage {

die <<"USAGE";


Usage: $0 -i fastq_file_dir -k fastq_barcodes_infile -e restriction_enzymes -s max_num_barcode_reads -c min_num_present_tag -m tassel_num_ram -p bwa_num_cpu -o output_dir

Description - 

OPTIONS:

	-i fastq_file_dir - 
	
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
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create the tag counts output directory if it doesn't already exist.
# my $tag_count_output_dir = join('/', $output_dir, "tagCounts");
my $tag_count_output_dir = join('/', $output_dir, "TAG_COUNTS_DIR");
unless(-d $tag_count_output_dir){
      mkdir($tag_count_output_dir, 0777) or die "Can't make directory: $!";
}

fastq_to_tag_counts($fastq_infile_dir, $project_name, $fastq_barcodes_infile, $restriction_enzymes, $max_num_barcode_reads, $min_num_present_fastq_to_tag, 
	$tassel_num_ram, $tag_count_output_dir);

# Create the merged tag counts output directory if it doesn't already exist.
# my $merged_tag_count_output_dir = join('/', $output_dir, "mergedTagCounts");
my $merged_multiple_tag_count_output_dir = join('/', $output_dir, "MERGED_TAG_COUNTS_DIR");
unless(-d $merged_multiple_tag_count_output_dir){
      mkdir($merged_multiple_tag_count_output_dir, 0777) or die "Can't make directory: $!";
}


my ($merged_multiple_tag_counts_outfile, $merged_multiple_tag_counts_fasta_outfile) = merge_multiple_tag_counts($tag_count_output_dir, $project_name, 
	$min_num_present_merge_multiple_tag, $tassel_num_ram, $merged_multiple_tag_count_output_dir);

my $merged_multiple_tag_counts_to_fastq_outfile = merge_multiple_tag_counts_to_fastq($merged_multiple_tag_counts_outfile, $project_name, $min_num_present_merge_multiple_tag, 
	$tassel_num_ram, $merged_multiple_tag_count_output_dir);

# Create the reference genome output directory if it doesn't already exist.
my $reference_genome_output_dir = join('/', $output_dir, "REFERENCE_GENOME");
unless(-d $reference_genome_output_dir){
      mkdir($reference_genome_output_dir, 0777) or die "Can't make directory: $!";
}

my $reference_genome_fasta_outfile = convert_reference_genome_to_bwa_tassel($reference_genome_infile, $project_name, $reference_genome_output_dir);

bwa_index($reference_genome_fasta_outfile, $project_name, $reference_genome_output_dir);

my $bwa_alignment_outfile = bwa_aln($reference_genome_fasta_outfile, $project_name, $bwa_num_cpu, $merged_multiple_tag_counts_to_fastq_outfile, 
	$merged_multiple_tag_count_output_dir);

my $bwa_aligned_master_outfile = bwa_samse($reference_genome_fasta_outfile, $project_name, $bwa_alignment_outfile, $merged_multiple_tag_counts_to_fastq_outfile, 
	$merged_multiple_tag_count_output_dir);

# Create the topm output directory if it doesn't already exist.
my $topm_output_dir = join('/', $output_dir, "TAGS_ON_PHYSICAL_MAP_DIR");
unless(-d $topm_output_dir){
      mkdir($topm_output_dir, 0777) or die "Can't make directory: $!";
}

my $topm_outfile = sam_converter($bwa_aligned_master_outfile, $project_name, $tassel_num_ram, $topm_output_dir);

# Create the tbt output directory if it doesn't already exist.
my $tbt_output_dir = join('/', $output_dir, "TAGS_BY_TAXA_DIR");
unless(-d $tbt_output_dir){
      mkdir($tbt_output_dir, 0777) or die "Can't make directory: $!";
}

fastq_to_tags_by_taxa($fastq_infile_dir, $project_name, $fastq_barcodes_infile, $restriction_enzymes, $topm_outfile, $tassel_num_ram, $tbt_output_dir);

# Create the mergedtbt output directory if it doesn't already exist.
my $merged_tbt_output_dir = join('/', $output_dir, "MERGED_TAGS_BY_TAXA_DIR");
unless(-d $merged_tbt_output_dir){
      mkdir($merged_tbt_output_dir, 0777) or die "Can't make directory: $!";
}

my $merged_tags_by_taxa_outfile = merge_tags_by_taxa_files($tbt_output_dir, $project_name, $max_num_tags_tbt, $tassel_num_ram, $merged_tbt_output_dir);

# Create the hap map output directory if it doesn't already exist.
my $hap_map_output_dir = join('/', $output_dir, "HAP_MAP_DIR");
unless(-d $hap_map_output_dir){
      mkdir($hap_map_output_dir, 0777) or die "Can't make directory: $!";
}

# Create the hap map raw output directory if it doesn't already exist.
my $hap_map_raw_output_dir = join('/', $hap_map_output_dir, "HAP_MAP_RAW_DIR");
unless(-d $hap_map_raw_output_dir){
      mkdir($hap_map_raw_output_dir, 0777) or die "Can't make directory: $!";
}


my $hap_map_genotype_outfile = tags_to_snp_by_alignment($merged_tags_by_taxa_outfile, $topm_outfile, $project_name, $max_num_sites, 
	$min_minor_allele_freq_tag_to_snp, $min_minor_allele_count, $reference_genome_fasta_outfile, $start_chromosome, $end_chromosome, 
	$tassel_num_ram, $topm_output_dir, $hap_map_raw_output_dir);
        
# Create the hap map merged SNP output directory if it doesn't already exist.
my $hap_map_merged_snps_output_dir = join('/', $hap_map_output_dir, "HAP_MAP_MERGED_SNPS_DIR");
unless(-d $hap_map_merged_snps_output_dir){
      mkdir($hap_map_merged_snps_output_dir, 0777) or die "Can't make directory: $!";
}

my $merged_hap_map_genotype_outfile = merge_duplicate_snps($hap_map_genotype_outfile, $project_name, $genotypic_mismat_thres, $start_chromosome, $end_chromosome, 
	$tassel_num_ram, $hap_map_merged_snps_output_dir);

# Create the hap map filtered SNP output directory if it doesn't already exist.
#my $hap_map_filtered_snps_output_dir = join('/', $hap_map_output_dir, "HAP_MAP_FILTERED_SNPS_DIR");
#unless(-d $hap_map_filtered_snps_output_dir){
#      mkdir($hap_map_filtered_snps_output_dir, 0777) or die "Can't make directory: $!";
#}

#my $filtered_hap_map_snp_outfile = gbs_hap_map_filters($merged_hap_map_genotype_outfile, $project_name, $min_minor_allele_freq_gbs_hap_map, $min_site_coverage, 
#	$start_chromosome, $end_chromosome, $tassel_num_ram, $hap_map_filtered_snps_output_dir);

#/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx11g -fork1 -FastqToTagCountPlugin -i ./fastq -k mpb_barcodes.txt -e PstI-MspI -s 300000000 -c 1 -o ./tagCounts -endPlugin -runfork1 > ./tagCounts.log
sub fastq_to_tag_counts{

        my $fastq_infile_dir = shift;
        die "Error lost the input directory containing FASTQ text (fastq.txt) or gzipped FASTQ (fastq.gz) text files" unless defined $fastq_infile_dir;
        
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
        
        my $fastq_to_tag_counts_log_outfile = join('/', $tag_count_output_dir, join("_", $project_name, "FastqToTagCounts.log"));
	warn "Generating tag counts list file using the FastqToTagCountPlugin.....\n";
	my $fastqToTagCountsCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -FastqToTagCountPlugin -i $fastq_infile_dir -k $fastq_barcodes_infile -e $restriction_enzymes -s $max_num_barcode_reads -c $min_num_present_tag -o $tag_count_output_dir -endPlugin -runfork1 > $fastq_to_tag_counts_log_outfile");
	warn $fastqToTagCountsCmd . "\n\n";
	system($fastqToTagCountsCmd) == 0 or die "Error calling $fastqToTagCountsCmd: $?";
}

#/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx11g -fork1 -MergeMultipleTagCountPlugin -i ./tagCounts -o ./mergedTagCounts/mpbGBSTags.cnt -c 5 -t -endPlugin -runfork1 > ./MTC.log
sub merge_multiple_tag_counts{

        my $tag_count_input_dir = shift;
        die "Error lost the input directory containing tag count (.cnt) files." unless defined $tag_count_input_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $min_num_present_merge_multiple_tag = shift;
        die "minimum number of times a tag must be present to be output" unless defined $min_num_present_merge_multiple_tag;
        
        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $merged_multiple_tag_count_output_dir = shift;
        die "Error lost the output directory to contain output .cnt (merged tag count) files" unless defined $merged_multiple_tag_count_output_dir;
        
        # The merged tag counts output file name.
        my $merged_multiple_tag_counts_outfile = join('/', $merged_multiple_tag_count_output_dir, join("_", $project_name, "Tags.cnt"));
        
        my $merge_multiple_tag_counts_log_outfile = join('/', $merged_multiple_tag_count_output_dir, join("_", $project_name, "MergeMultipleTagCounts.log"));
	warn "Generating merged multiple tag counts list file using the MergeMultipleTagCountPlugin.....\n";
	my $mergeMultipleTagCountCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -MergeMultipleTagCountPlugin -i $tag_count_input_dir -o $merged_multiple_tag_counts_outfile -c $min_num_present_merge_multiple_tag -endPlugin -runfork1 > $merge_multiple_tag_counts_log_outfile");
	warn $mergeMultipleTagCountCmd . "\n\n";
	system($mergeMultipleTagCountCmd) == 0 or die "Error calling $mergeMultipleTagCountCmd: $?";
	
	my $merge_multiple_tag_fastq_counts_log_outfile = join('/', $merged_multiple_tag_count_output_dir, join("_", $project_name, "MergeMultipleTagFastqCounts.log"));
	warn "Generating merged multiple tag counts fasta list file using the MergeMultipleTagCountPlugin.....\n";
	my $mergeMultipleTagCountFastaCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -MergeMultipleTagCountPlugin -i $tag_count_input_dir -o $merged_multiple_tag_counts_outfile -c $min_num_present_merge_multiple_tag -t -endPlugin -runfork1 > $merge_multiple_tag_fastq_counts_log_outfile");
	warn $mergeMultipleTagCountFastaCmd . "\n\n";
	system($mergeMultipleTagCountFastaCmd) == 0 or die "Error calling $mergeMultipleTagCountFastaCmd: $?";
	
	my $merged_multiple_tag_counts_fasta_outfile = $merged_multiple_tag_counts_outfile . ".fq";
	return ($merged_multiple_tag_counts_outfile, $merged_multiple_tag_counts_fasta_outfile);
}

#/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx11g -fork1 -TagCountToFastqPlugin -i ./mergedTagCounts/mpbGBSTags.cnt -o ./mergedTagCounts/mpbGBSTags.fq -c 5 -endPlugin -runfork1 > ./MTCtoFQ.log
sub merge_multiple_tag_counts_to_fastq{

        my $merged_multiple_tag_counts_infile = shift;
        die "Error lost the merged tag counts input file" unless defined $merged_multiple_tag_counts_infile;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $min_num_present_merge_multiple_tag = shift;
        die "minimum number of times a tag must be present to be output" unless defined $min_num_present_merge_multiple_tag;
        
        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $merged_multiple_tag_count_output_dir = shift;
        die "Error lost the output directory to contain output .cnt (merged tag count) files" unless defined $merged_multiple_tag_count_output_dir;
        
        # The merged tag counts output file name.
        my $merged_multiple_tag_counts_to_fastq_outfile = join('/', $merged_multiple_tag_count_output_dir, join("_", $project_name, "Tags.fq"));
        
        my $merge_multiple_tag_counts_to_fastq_log_outfile = join('/', $merged_multiple_tag_count_output_dir, join("_", $project_name, "MergeMultipleTagCountsToFastq.log"));
	warn "Generating tag counts to fastq list file using the TagCountToFastqPlugin.....\n";
	my $mergeMultipleTagCountToFastqCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -TagCountToFastqPlugin -i $merged_multiple_tag_counts_infile -o $merged_multiple_tag_counts_to_fastq_outfile -c $min_num_present_merge_multiple_tag -endPlugin -runfork1 > $merge_multiple_tag_counts_to_fastq_log_outfile");
	warn $mergeMultipleTagCountToFastqCmd . "\n\n";
	system($mergeMultipleTagCountToFastqCmd) == 0 or die "Error calling $mergeMultipleTagCountToFastqCmd: $?";
	
	return $merged_multiple_tag_counts_to_fastq_outfile;
}

sub convert_reference_genome_to_bwa_tassel{

	my $reference_genome_fasta_infile = shift;
	die "Error lost the tab-delimited adaptor blastn file" unless defined $reference_genome_fasta_infile;
	
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
	my $reference_genome_output_dir = shift;
	die "Error lost reference genome to bwa-tassel output directory" unless defined $reference_genome_output_dir;

	warn "Calling reformat_refgen_to_bwatassel for $reference_genome_fasta_infile....\n";
	my $reference_genome_fasta_filename = fileparse($reference_genome_fasta_infile);

	my $reference_genome_fasta_outfile = join('/', $reference_genome_output_dir, join("_", $project_name, $reference_genome_fasta_filename . ".fasta"));

	my $reference_genome_toc_outfile = join('/', $reference_genome_output_dir, join("_", $project_name, $reference_genome_fasta_filename, "toc.txt"));
	
	unless(-s $reference_genome_fasta_outfile and -s $reference_genome_toc_outfile){
		open(OUTFILE1, ">$reference_genome_fasta_outfile") or die "Couldn't open file $reference_genome_fasta_outfile for writting, $!";
		open(OUTFILE2, ">$reference_genome_toc_outfile") or die "Couldn't open file $reference_genome_toc_outfile for writting, $!"; 
		print OUTFILE2 join("\t", "bwatassel_fasta_header", "original_fasta_header") . "\n";
		warn "Converting $reference_genome_fasta_filename to bwatassel format....\n\n";
		my $seqio = Bio::SeqIO->new(-file => $reference_genome_fasta_infile, '-format' => 'Fasta');
		my $i = 1;
		while(my $seq_entry = $seqio->next_seq) {

			my $seq_id = $seq_entry->id;
			my $sequence = $seq_entry->seq;
			my $seq_desc = $seq_entry->desc;
		
			my $fasta_header = join(" ", $seq_id, $seq_desc);
			warn join("\t", $i, $fasta_header) . "\n";
			print OUTFILE1 join("\n", ">$i", $sequence) . "\n";
			print OUTFILE2 join("\t", $i, $fasta_header) . "\n";
			
			$i++;
		
		}
		close(OUTFILE1) or die "Couldn't close file $reference_genome_fasta_outfile"; 
		close(OUTFILE2) or die "Couldn't close file $reference_genome_toc_outfile";

		# Clean out the sequence I/O object.
		$seqio = ();
	}
	return $reference_genome_fasta_outfile;
}

#bwa index -a bwtsw ./refgen/renumbered_Msequence.fasta
sub bwa_index{

	my $renumbered_reference_genome_fasta_infile = shift;
	die "Error lost the renumbered reference genome input file" unless defined $renumbered_reference_genome_fasta_infile;
	
	my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $reference_genome_output_dir = shift;
	die "Error lost reference genome to bwa-tassel output directory" unless defined $reference_genome_output_dir;
	
	# format the index file into .amb .ann .bwt .pac .sa files.
	my ($renum_ref_gen_fastadbAMB, $renum_ref_gen_fastadbANN, $renum_ref_gen_fastadbBWT, $renum_ref_gen_fastadbPAC, $renum_ref_gen_fastadbSA);
	$renum_ref_gen_fastadbAMB = $renumbered_reference_genome_fasta_infile . '.amb';
	$renum_ref_gen_fastadbANN = $renumbered_reference_genome_fasta_infile . '.ann';
	$renum_ref_gen_fastadbBWT = $renumbered_reference_genome_fasta_infile . '.bwt';
	$renum_ref_gen_fastadbPAC = $renumbered_reference_genome_fasta_infile . '.pac';
	$renum_ref_gen_fastadbSA = $renumbered_reference_genome_fasta_infile . '.sa';
	unless(-s $renum_ref_gen_fastadbAMB and -s $renum_ref_gen_fastadbANN and -s $renum_ref_gen_fastadbBWT 
		and -s $renum_ref_gen_fastadbPAC and -s $renum_ref_gen_fastadbSA){
		warn "Generating the bwa index file.....\n";
		my $bwaIndexCmd  = "$bwa index -a bwtsw $renumbered_reference_genome_fasta_infile";
		warn $bwaIndexCmd . "\n\n"; 
		system($bwaIndexCmd) == 0 or die "Error calling $bwaIndexCmd: $?";
	}
}

#bwa aln -t 4 ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/mpbGBSTags.cnt.fq > ./mergedTagCounts/AlignedGBSTags1.sai
sub bwa_aln{

	my $renumbered_reference_genome_fasta_infile = shift;
	die "Error lost the renumbered reference genome input file" unless defined $renumbered_reference_genome_fasta_infile;
	
	my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $bwa_num_cpu = shift;
        die "Error lost the number of cores for BWA" unless defined $bwa_num_cpu;
        
        my $merged_multiple_tag_counts_fasta_infile = shift;
        die "Error lost the tag count file in FASTQ format (.cnt.fq) produced by the MergeMultipleTagCountPlugin" unless defined $merged_multiple_tag_counts_fasta_infile;
        
	my $merged_multiple_tag_count_output_dir = shift;
        die "Error lost the output directory to contain output .cnt (merged tag count) files" unless defined $merged_multiple_tag_count_output_dir;
	
        my $bwa_alignment_outfile = join('/', $merged_multiple_tag_count_output_dir, join("_", "Aligned", $project_name, "Tags.sai"));
        
        unless(-s $bwa_alignment_outfile){
		warn "Generating the bwa alignment file.....\n";
		my $bwaAlnCmd  = "$bwa aln -t $bwa_num_cpu $renumbered_reference_genome_fasta_infile $merged_multiple_tag_counts_fasta_infile > $bwa_alignment_outfile";
		warn $bwaAlnCmd . "\n\n";
		system($bwaAlnCmd) == 0 or die "Error calling $bwaAlnCmd: $?";
	}
	return $bwa_alignment_outfile;
}

#bwa samse ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/AlignedGBSTags1.sai ./mergedTagCounts/mpbGBSTags.cnt.fq > mergedTagCounts/AlignedMasterTagsMPB.sam
sub bwa_samse{

	my $renumbered_reference_genome_fasta_infile = shift;
	die "Error lost the renumbered reference genome input file" unless defined $renumbered_reference_genome_fasta_infile;
	
	my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $bwa_alignment_infile = shift;
        die "Error lost the BWA SAI formatted alignment (.sai) file" unless defined $bwa_alignment_infile;
        
        my $merged_multiple_tag_counts_fasta_infile = shift;
        die "Error lost the tag count file in FASTQ format (.fq) produced by the TagCountToFastqPlugin" unless defined $merged_multiple_tag_counts_fasta_infile;
        
	my $merged_multiple_tag_count_output_dir = shift;
        die "Error lost the output directory to contain output .cnt (merged tag count) files" unless defined $merged_multiple_tag_count_output_dir;
	
        my $bwa_aligned_master_outfile = join('/', $merged_multiple_tag_count_output_dir, join("_", "AlignedMasterTags", $project_name . ".sam"));
	unless(-s $bwa_aligned_master_outfile){
		warn "Generating the bwa sam file.....\n";
		my $bwaSamseCmd  = "$bwa samse $renumbered_reference_genome_fasta_infile $bwa_alignment_infile $merged_multiple_tag_counts_fasta_infile > $bwa_aligned_master_outfile";
		warn $bwaSamseCmd . "\n\n";
		system($bwaSamseCmd) == 0 or die "Error calling $bwaSamseCmd: $?";
	}
	return $bwa_aligned_master_outfile;
}

#/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx12g -fork1 -SAMConverterPlugin -i ./mergedTagCounts/AlignedMasterTagsMPB.sam -o ./topm/MasterTagsMPB.topm -endPlugin -runfork1 > ./SAMconvert.log
sub sam_converter{

        my $bwa_aligned_master_infile = shift;
        die "Error lost the bwa alignment input file in SAM format (.sam)" unless defined $bwa_aligned_master_infile;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $topm_output_dir = shift;
        die "Error lost the output directory to contain output TagsOnPhysicalMap (TOPM) file" unless defined $topm_output_dir;
        
        # The master tags topm output file name.
        my $topm_outfile = join('/', $topm_output_dir, join("_", "MasterTags", $project_name . ".topm"));
        
        my $sam_convert_log_outfile = join('/', $topm_output_dir, join("_", $project_name, "SAMConvert.log"));
	warn "Generating tags on physical map file using the SAMConverterPlugin.....\n";
	my $SAMConverterCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -SAMConverterPlugin -i $bwa_aligned_master_infile -o $topm_outfile -endPlugin -runfork1 > $sam_convert_log_outfile");
	warn $SAMConverterCmd . "\n\n";
	system($SAMConverterCmd) == 0 or die "Error calling $SAMConverterCmd: $?";
	
	return $topm_outfile;
}

#/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx12g -fork1 -FastqToTBTPlugin -i ./fastq -k mpb_barcodes.txt -e PstI-MspI -o ./tbt -y -m ./topm/MasterTagsMPB.topm -endPlugin -runfork1 > ./TBT.log
sub fastq_to_tags_by_taxa{

        my $fastq_infile_dir = shift;
        die "Error lost the input directory containing FASTQ text (fastq.txt) or gzipped FASTQ (fastq.gz) text files" unless defined $fastq_infile_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $fastq_barcodes_infile = shift;
        die "Error lost the key file listing barcodes for each sample and plate layout" unless defined $fastq_barcodes_infile;

        my $restriction_enzymes = shift;
        die "Error lost the enzyme(s) used to create the GBS library" unless defined $restriction_enzymes;
        
        my $topm_infile = shift;
        die "Error lost the TagsOnPhysicalMap (TOPM) input file" unless defined $topm_infile;
        
        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $tbt_output_dir = shift;
        die "Error lost the output directory to contain output TagsByTaxa (TBT) file" unless defined $tbt_output_dir;
       
        my $tags_by_taxa_log_outfile = join('/', $tbt_output_dir, join("_", $project_name, "TagsByTaxa.log"));
	warn "Generating tags by taxa file using the FastqToTBTPlugin.....\n";
	my $FastqToTBTCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -FastqToTBTPlugin -i $fastq_infile_dir -k $fastq_barcodes_infile -e $restriction_enzymes -o $tbt_output_dir -y -m $topm_infile -endPlugin -runfork1 > $tags_by_taxa_log_outfile");
	warn $FastqToTBTCmd . "\n\n";
	system($FastqToTBTCmd) == 0 or die "Error calling $FastqToTBTCmd: $?";
	
}

#/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx12g -fork1 -MergeTagsByTaxaFilesPlugin -i ./tbt -o ./mergedTBT/mpbGBSstudy.tbt.byte -s 200000000 -endPlugin -runfork1 > ./mergeTBT.log
sub merge_tags_by_taxa_files{
        my $tbt_input_dir = shift;
        die "Error lost the input directory containing the TagsByTaxa (TBT) file" unless defined $tbt_input_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $max_num_tags_tbt = shift;
        die "Error lost maximum number of tags the TBT can hold while merging" unless defined $max_num_tags_tbt;
        
        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $merged_tbt_output_dir = shift;
        die "Error lost the output directory to contain output MergeTagsByTaxa (mergedTBT) file" unless defined $merged_tbt_output_dir;
        
        # The merged tags by taxa output file name.
        my $merged_tags_by_taxa_outfile = join('/', $merged_tbt_output_dir, join("_", $project_name, "Study.tbt.byte"));
        
        my $tags_by_taxa_log_outfile = join('/', $merged_tbt_output_dir, join("_", $project_name, "MergeTagsByTaxa.log"));
	warn "Generating merged tags by taxa file using the MergeTagsByTaxaFilesPlugin.....\n";
	my $mergeTagsByTaxaFilesCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -MergeTagsByTaxaFilesPlugin -i $tbt_input_dir -o $merged_tags_by_taxa_outfile -s $max_num_tags_tbt -endPlugin -runfork1 > $tags_by_taxa_log_outfile");
	warn $mergeTagsByTaxaFilesCmd . "\n\n";
	system($mergeTagsByTaxaFilesCmd) == 0 or die "Error calling $mergeTagsByTaxaFilesCmd: $?";
	
	return $merged_tags_by_taxa_outfile;
	
}

#### START CHECKING SCRIPT HERE BECAUSE THIS IS WHERE YOU STARTED ON NOVEMBER 06, 2014 ####
#/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx12g -fork1 -TagsToSNPByAlignmentPlugin -i ./mergedTBT/mpbGBSstudy.tbt.byte -y -m ./topm/MasterTagsMPB.topm -mUpd ./topm/MasterTagsMPBwVariants.topm -o ./hapmap/raw/mpbGBSGenos_chr+.hmp.txt -mxSites 500000 -mnMAF 0.02 -mnMAC 100000 -ref ./refgen/renumbered_Msequence.fasta -sC 1 -eC 8188 -endPlugin -runfork1 > ./TagstoSNPAlign.log
sub tags_to_snp_by_alignment{

        my $merged_tags_by_taxa_infile = shift;
        die "Error lost the merged tags by taxa input file" unless defined $merged_tags_by_taxa_infile;
        
        my $topm_infile = shift;
        die "Error lost the TagsOnPhysicalMap (TOPM) input file" unless defined $topm_infile;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $max_num_sites = shift;
	die "Error lost maximum number of sites (SNPs) output per chromosome" unless defined $max_num_sites;
	
	my $min_minor_allele_freq_tag_to_snp = shift;
	die "Error lost the minimum minor allele frequency" unless defined $min_minor_allele_freq_tag_to_snp;
	
	my $min_minor_allele_count = shift;
	die "Error lost the minimum minor allele count " unless defined $min_minor_allele_count;
	
        my $renumbered_reference_genome_fasta_infile = shift;
	die "Error lost the renumbered reference genome input file" unless defined $renumbered_reference_genome_fasta_infile;
	
	my $start_chromosome = shift;
	die "Error lost the start position of the chromosome" unless defined $start_chromosome;
	
	my $end_chromosome = shift;
	die "Error lost the end position of the chromosome" unless defined $end_chromosome;
	
        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $topm_output_dir = shift;
        die "Error lost the output directory to contain output TagsOnPhysicalMap (TOPM) file" unless defined $topm_output_dir;

        my $hap_map_raw_dir = shift;
        die "Error lost the output directory to contain output HapMap genotype file." unless defined $hap_map_raw_dir;
        
        # The updated TOPM file with variants called during SNP calling
        my $updated_topm_with_variants_outfile = join('/', $topm_output_dir, join("_", "MasterTags", $project_name, "WithVariants.topm"));
        
        # The HapMap genotype file
        my $hap_map_genotype_outfile = join('/', $hap_map_raw_dir, join("_", $project_name, "Genos_Chr+.hmp.txt"));
        
        my $tags_to_snp_align_log_outfile = join('/', $hap_map_raw_dir, join("_", $project_name, "TagsToSNPByAlign.log"));
	warn "Generating tags to SNP by alignment file using the TagsToSNPByAlignmentPlugin.....\n";
	my $tagsToSNPByAlignmentCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -TagsToSNPByAlignmentPlugin -i $merged_tags_by_taxa_infile -y -m $topm_infile -mUpd $updated_topm_with_variants_outfile -o $hap_map_genotype_outfile -mxSites $max_num_sites -mnMAF $min_minor_allele_freq_tag_to_snp -mnMAC $min_minor_allele_count -ref $renumbered_reference_genome_fasta_infile -sC $start_chromosome -eC $end_chromosome -endPlugin -runfork1 > $tags_to_snp_align_log_outfile");
	warn $tagsToSNPByAlignmentCmd . "\n\n";
	system($tagsToSNPByAlignmentCmd) == 0 or die "Error calling $tagsToSNPByAlignmentCmd: $?";
	
	return $hap_map_genotype_outfile;
}

#/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx12g -fork1 -MergeDuplicateSNPsPlugin -hmp ./hapmap/raw/mpbGBSGenos_chr+.hmp.txt -o ./hapmap/mergedSNPs/mpbGBSGenos_mergedSNPs_chr+.hmp.txt -misMat 0.1 -callHets -sC 1 -eC 8188 -endPlugin -runfork1 > MergeDupSNP.log
sub merge_duplicate_snps{

        my $hap_map_genotype_infile = shift;
        die "Error lost the hap map genotype input file" unless defined $hap_map_genotype_infile;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
	
	my $genotypic_mismat_thres = shift;
	die "Error lost the threshold genotypic mismatch rate" unless defined $genotypic_mismat_thres;
	
	my $start_chromosome = shift;
	die "Error lost the start position of the chromosome" unless defined $start_chromosome;
	
	my $end_chromosome = shift;
	die "Error lost the end position of the chromosome" unless defined $end_chromosome;
	
        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $hap_map_merged_snps_output_dir = shift;
        die "Error lost the output directory to contain output HapMap merged SNPs file." unless defined $hap_map_merged_snps_output_dir;
        
        # The merged HapMap genotype file
        my $merged_hap_map_genotype_outfile = join('/', $hap_map_merged_snps_output_dir, join("_", $project_name, "GenosMergedSNPs_Chr+.hmp.txt"));
        
        my $merge_duplicate_snps_log_outfile = join('/', $hap_map_merged_snps_output_dir, join("_", $project_name, "MergeDuplicateSNPs.log"));
	warn "Generating merged duplicate SNPs file using the MergeDuplicateSNPsPlugin.....\n";
	my $mergeDuplicateSNPsCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -MergeDuplicateSNPsPlugin -hmp $hap_map_genotype_infile -o $merged_hap_map_genotype_outfile -misMat $genotypic_mismat_thres -callHets -sC $start_chromosome -eC $end_chromosome -endPlugin -runfork1 > $merge_duplicate_snps_log_outfile");
	warn $mergeDuplicateSNPsCmd . "\n\n";
	system($mergeDuplicateSNPsCmd) == 0 or die "Error calling $mergeDuplicateSNPsCmd: $?";
	
	return $merged_hap_map_genotype_outfile;
}

# /home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl -Xmx12g -fork1 -GBSHapMapFiltersPlugin -hmp ./hapmap/mergedSNPs/mpbGBSGenos_mergedSNPs_chr+.hmp.txt -o ./hapmap/filt/mpbGBSGenos_mergedSNPsFilt_chr+.hmp.txt -mnSCov 0.2 -MnMAF 0.01 -sC 1 -eC 8188 -endPlugin -runfork1 > ./hmpFilt.log
sub gbs_hap_map_filters{

        my $merged_hap_map_genotype_infile = shift;
        die "Error lost the merged hap map snps input file" unless $merged_hap_map_genotype_infile;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
	
	my $min_minor_allele_freq_gbs_hap_map = shift;
	die "Error lost the minimum minor allele frequency" unless defined $min_minor_allele_freq_gbs_hap_map;
	
	my $min_site_coverage = shift;
	die "Error lost the minimum site coverage" unless defined $min_site_coverage;
	
	my $start_chromosome = shift;
	die "Error lost the start position of the chromosome" unless defined $start_chromosome;
	
	my $end_chromosome = shift;
	die "Error lost the end position of the chromosome" unless defined $end_chromosome;
	
        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $hap_map_filtered_snps_output_dir = shift;
        die "Error lost the output directory to contain output HapMap filtered SNPs file." unless defined $hap_map_filtered_snps_output_dir;
        
        # The filtered HapMap snp file
        my $filtered_hap_map_snp_outfile = join('/', $hap_map_filtered_snps_output_dir, join("_", $project_name, "GenosMergedSNPsFilters_Chr+.hmp.txt"));
        
        my $gbs_hap_map_filters_log_outfile = join('/', $hap_map_filtered_snps_output_dir, join("_", $project_name, "GBSHapMapFilters.log"));
        
	warn "Generating filtered SNPs file using the GBSHapMapFiltersPlugin.....\n";
	my $GBSHapMapFiltersCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -GBSHapMapFiltersPlugin -hmp $merged_hap_map_genotype_infile -o $filtered_hap_map_snp_outfile -mnSCov $min_site_coverage -MnMAF $min_minor_allele_freq_gbs_hap_map -sC $start_chromosome -eC $end_chromosome -endPlugin -runfork1 > $gbs_hap_map_filters_log_outfile");
	warn $GBSHapMapFiltersCmd . "\n\n";
	system($GBSHapMapFiltersCmd) == 0 or die "Error calling $GBSHapMapFiltersCmd: $?";
	
	return $filtered_hap_map_snp_outfile;
}
