#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use IPC::Open2;
use File::Basename;
use File::Copy;

# VERSION 1.21

#### PROGRAM NAME ####
# refgen_stacks_analysis_pipeline.pl - Program that aligns quality filtered, demultiplexed, and adapter trimmed GBS data sequences (unpadded) against a reference genome using the BWA alignment program. It then runs the Stacks reference genome pipeline using the pstacks, cstacks, and sstacks programs of the Stacks Software Suite.

#### DESCRIPTION ####
# This program takes the quality filtered, demultiplexed, and adapter trimmed *.fastq input files (unpadded) and a reference genome fasta input file as input. Converts a reference genome fasta file to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers. It then performs a BWA alignment to align the GBS fastq sequences to a reference genome. It then executes the pstacks program, which extracts sequence stacks that were aligned to a reference genome using the BWA alignment program and identifies SNPs. These sequence stacks are then processed using cstacks and sstacks to obtain the filtered SNP stacks output files.

#### SAMPLE COMMANDS ####

## Standard example
# perl refgen_stacks_analysis_pipeline.pl -i ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR_UNPADDED/STEPHEN_TREVOY/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_sequence_data/DendPond_male_1.0/Primary_Assembly/unplaced_scaffolds/FASTA/DendPond_male_1.0_unplaced.scaf.fa -c 7 -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3

## Using the cstacks catalog prefix path.
# perl ~/refgen_stacks_analysis_pipeline.pl -i ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_TRIMMED_GBS/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_MALE_GBS_STACKS/REFERENCE_GENOME/DendPond_male_1.0_unplaced.scaf.fa.fasta -p ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_MALE_GBS_STACKS/STACKS_OUTFILES/batch_1 -t fastq -l 92 -b 1 -d 5 -a 0.05 -s 1 -m 20 -c 24 -o ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/test_output_dir/
my ($gbs_fastq_dir, $gbs_fastq_file_type, $refgen_infile, $cstacks_catalog_prefix, $use_existing_catalog, $gbs_sequence_length, $stacks_sql_id, $min_depth_coverage_pstacks, $alpha_value_pstacks, $num_mismatches_tag, $sam_mapq_threshold, $num_threads, $output_dir);
GetOptions(
    'i=s'    => \$gbs_fastq_dir, # The absolute path to the quality filtered, demultiplexed, and adapter trimmed *.fastq input file (unpadded) directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
    't=s'    => \$gbs_fastq_file_type, # The fastq input file type. Can be either fastq or gzfastq. Default: gzfastq
    'g=s'    => \$refgen_infile, # The absolute path to the reference genome input fasta file to align GBS fastq sequences.
    'p=s'    => \$cstacks_catalog_prefix, # The cstacks catalog prefix file path consisting of a set of consensus loci built from a set of samples processed by the pstacks program. i.e. /path/to/catalog/batch_1
    'e=s'    => \$use_existing_catalog, # The use existing catalog bitwise flag. If true add data to the existing catalog, Otherwise run the default cstacks program. Default: false
    'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
    'b=s'    => \$stacks_sql_id, # The SQL ID to insert into the output to identify this sample. Default: 1
    'd=s'    => \$min_depth_coverage_pstacks, # The minimum depth of coverage to report a stack. Default: 5
    'a=s'    => \$alpha_value_pstacks, # The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
    's=s'    => \$num_mismatches_tag, # The number of mismatches allowed between sample tags when generating the catalog. Default: 1
    'm=s'    => \$sam_mapq_threshold, # The mapping quality score threshold filter for sam alignments. Any mapping quality score below this value is filtered out of the sam alignment files. Default: 20
    'c=s'    => \$num_threads, # The number of cpu threads to use for the stacks programs. Default: 2
    'o=s'    => \$output_dir, # The absolute path to the output directory to contain the renumerated reference genome, BWA sam, padded sam, and Stacks output files and directories.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
    defined $gbs_fastq_dir
    and defined $refgen_infile
    and defined $output_dir
);

# The fastq input file type. Can be either fastq or gzfastq. Default: gzfastq
$gbs_fastq_file_type = 'gzfastq' unless defined $gbs_fastq_file_type;

# The cstacks catalog prefix file path consisting of a set of consensus loci built from a set of samples processed by the pstacks program. i.e. /path/to/catalog/batch_1
$cstacks_catalog_prefix = "" unless defined $cstacks_catalog_prefix;

# The use existing catalog bitwise flag. If true add data to the existing catalog, Otherwise run the default cstacks program. Default: false
$use_existing_catalog = "false" unless defined $use_existing_catalog;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

# The SQL ID to insert into the output to identify this sample. Default: 1
$stacks_sql_id = 1 unless defined $stacks_sql_id;

# The minimum depth of coverage to report a stack. Default: 5
$min_depth_coverage_pstacks = 5 unless defined $min_depth_coverage_pstacks;

# The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
$alpha_value_pstacks = 0.05 unless defined $alpha_value_pstacks;

# The number of mismatches allowed between sample tags when generating the catalog. Default: 1
$num_mismatches_tag = 1 unless defined $num_mismatches_tag;

# The mapping quality score threshold filter for sam alignments. Any mapq score below this value is filtered out of the sam alignment files. Default: 20
$sam_mapq_threshold = 20 unless defined $sam_mapq_threshold;

# The number of cpu threads to use for the stacks programs. Default: 2
$num_threads = 2 unless defined $num_threads;

# Program dependencies - The absolute paths to gunzip to uncompress fastq.gz input files, BWA alignment program, and Stacks related programs.
my ($gunzip, $bwa, $pstacks, $cstacks, $sstacks);
$gunzip				= '/bin/gunzip';
$bwa				= '/global/software/bwa/bwa0712/bin/bwa';
$pstacks			= '/global/software/stacks/stacks131/bin/pstacks';
$cstacks			= '/global/software/stacks/stacks131/bin/cstacks';
$sstacks			= '/global/software/stacks/stacks131/bin/sstacks';

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i gbs_fastq_dir -t gbs_fastq_file_type -g refgen_infile -p cstacks_catalog_prefix -e use_existing_catalog -l gbs_sequence_length -b stacks_sql_id -d min_depth_coverage_pstacks -a alpha_value_pstacks -s num_mismatches_tag -m sam_mapq_threshold -c num_threads -o output_dir
   
    
VERSION 1.21
    
DESCRIPTION - This program takes the quality filtered, demultiplexed, and adapter trimmed *.fastq input files (unpadded) and reference genome fasta input files as input. Converts the reference genome fasta file to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers. It then performs a BWA alignment to align the GBS fastq sequences to the reference genome. It then executes the pstacks program, which extracts sequence stacks that were aligned to the reference genome using the BWA alignment program and identifies SNPs. These sequence stacks are then processed using cstacks and sstacks to obtain the filtered SNP stacks output files.
    
OPTIONS:
    
-i gbs_fastq_dir - The absolute path to the quality filtered, demultiplexed, and adapter trimmed *.fastq input file (unpadded) directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.

-t gbs_fastq_file_type - The fastq input file type. Can be either fastq or gzfastq. Default: gzfastq

-g refgen_infile - The absolute path to the reference genome input fasta file to align GBS fastq sequences using the BWA alignment program.

-p cstacks_catalog_prefix - The cstacks catalog prefix file path consisting of a set of consensus loci built from a set of samples processed by the pstacks program. i.e. /path/to/catalog/batch_1

-e use_existing_catalog - The use existing catalog bitwise flag. If true add data to the existing catalog, Otherwise run the default cstacks program. Default: false

-l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92

-b stacks_sql_id - The SQL ID to insert into the output to identify this sample. Default: 1

-d min_depth_coverage_pstacks - The minimum depth of coverage to report a stack. Default: 5

-a alpha_value_pstacks - The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05

-s num_mismatches_tag - The number of mismatches allowed between sample tags when generating the catalog. Default: 1

-m sam_mapq_threshold - The mapping quality score threshold filter for sam alignments. Any mapq score below this value is filtered out of the sam alignment files. Default: 20

-c num_threads - The number of cpu cores to use for the stacks programs. Default: 2

-o output_dir - The absolute path to the output directory to contain the renumerated reference genome, BWA, and Stacks output files and directories.

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create the reference genome output directory if it doesn't already exist.
my $refgen_output_dir = join('/', $output_dir, "REFERENCE_GENOME");
unless(-d $refgen_output_dir){
      mkdir($refgen_output_dir, 0777) or die "Can't make directory: $!";
}

# Converts the reference genome to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers.
my ($refgen_fasta_outfile, $refgen_toc_outfile, $num_chromosomes) = convert_refgen_bwa_input_format($refgen_infile, $refgen_output_dir);

# Creates the BWA reference genome index file
bwa_index($refgen_fasta_outfile, $refgen_output_dir);

# Create the BWA alignment file output directory if it doesn't already exist.
my $bwa_output_dir = join('/', $output_dir, "BWA_ALIGNMENT_FILES");
unless(-d $bwa_output_dir){
      mkdir($bwa_output_dir, 0777) or die "Can't make directory: $!";
}

# Find all files in the processed GBS fastq file directory with the extension *.fastq or *.fastq.qz.
my ($gbs_fastq_files, $gbs_fastq_file_count);
if($gbs_fastq_file_type eq "gzfastq"){
	($gbs_fastq_files, $gbs_fastq_file_count) = find_files($gbs_fastq_dir, "fastq.gz");
}elsif($gbs_fastq_file_type eq "fastq"){
	($gbs_fastq_files, $gbs_fastq_file_count) = find_files($gbs_fastq_dir, "fastq");
}else{
	die "Error: $gbs_fastq_file_type is not a recognized input file type. Use the -t option and specify gzfastq for *.fastq.gz files or fastq for *.fastq files.";
}

die "Error: The GBS fastq file count was $gbs_fastq_file_count. The input file type was specified as $gbs_fastq_file_type. Either you specified the wrong directory or the files in $gbs_fastq_dir are a different format than $gbs_fastq_file_type. Use the -t option and specify gzfastq for *.fastq.gz files or fastq for *fastq files." unless($gbs_fastq_file_count > 0);

# Iterate through each GBS fastq file and perform an alignment to the reference genome using the BWA alignment program.
foreach my $file_name (sort keys %{$gbs_fastq_files}){
	warn "Processing " . $file_name . ".....\n";
	my $gbs_fastq_infile = $gbs_fastq_files->{$file_name};
    
	# If the fastq file is compressed, uncompress the file and set the resulting fastq filename to be the fastq infile.
	if($gbs_fastq_file_type eq "gzfastq"){
		my $uncompressed_fastq_file = gunzip_fastq_file($gbs_fastq_infile);
		$gbs_fastq_infile = $uncompressed_fastq_file;
	}elsif(($gbs_fastq_file_type eq "fastq") and ($file_name =~ m/trimmed_offset_\d+/)){ # If the fastq file is not compressed set the resulting fastq filename to be the fastq infile.
	        my ($fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/_trimmed_offset_\d+\.fastq/);
        	my $fastq_infile = join('/', $fastq_dir, $fastq_filename . ".fastq");
        	warn "Copying $gbs_fastq_infile to $fastq_infile.....";
        	copy($gbs_fastq_infile, $fastq_infile) or die "Copy failed: $!";
        	$gbs_fastq_infile = $fastq_infile;
   	}

	# Creates the BWA alignment file using the GBS fastq input file to align the fastq sequence reads to the reference genome.
	my $bwa_alignment_outfile = bwa_aln($refgen_fasta_outfile, $gbs_fastq_infile, $num_threads, $bwa_output_dir);
	
	# Creates the BWA single-ended alignment file in sam format using the GBS fastq input file to align the fastq sequence reads to the reference genome and the BWA aln format file.
	my $bwa_aligned_master_outfile = bwa_samse($refgen_fasta_outfile, $gbs_fastq_infile, $bwa_alignment_outfile, $sam_mapq_threshold, $bwa_output_dir);
	
	# Remove the uncompressed fastq file to save space as we do not need file after this point.
	unlink($gbs_fastq_infile) or die "Could not unlink $gbs_fastq_infile: $!" if(($gbs_fastq_file_type eq "gzfastq") and ($gbs_fastq_infile =~ m/\.fastq$/));
	unlink($gbs_fastq_infile) or die "Could not unlink $gbs_fastq_infile: $!" if(($gbs_fastq_file_type eq "fastq") and ($file_name =~ m/trimmed_offset_\d+/));
}

# Create the padded sam alignment file output directory if it doesn't already exist.
my $padded_sam_output_dir = join('/', $output_dir, "PADDED_SAM_OUTFILES");
unless(-d $padded_sam_output_dir){
	mkdir($padded_sam_output_dir, 0777) or die "Can't make directory: $!";
}

# Find all the BWA sam alignment output files from the BWA alignment directory with the extension *.sam.
my ($bwa_align_files, $bwa_align_file_count) = find_files($bwa_output_dir, "sam");

# Iterate through each BWA sam alignment output file with extension *.sam and padd sequences with poly-Ns. If the sequence is aligned to the reference genome in the sense strand orientation padd the sequence on the 3' end. Otherwise the sequence is aligned to the reference genome in the antisense strand orientation padd the sequence on the 5' end.
foreach my $file_name (sort keys %{$bwa_align_files}){
	warn "Processing " . $file_name . ".....\n";
	my $bwa_align_infile = $bwa_align_files->{$file_name};
	bwa_pad_sam_files($bwa_align_infile, $gbs_sequence_length, $sam_mapq_threshold, $padded_sam_output_dir);
	
	# Remove the BWA sam file to save space as we do not need the file after this point.
	#unlink($bwa_align_infile) or die "Could not unlink $bwa_align_infile: $!";
}

# Create the stacks output directory if it doesn't already exist.
my $stacks_output_dir = join('/', $output_dir, "STACKS_OUTFILES");
unless(-d $stacks_output_dir){
	mkdir($stacks_output_dir, 0777) or die "Can't make directory: $!";
}

# Find all padded sam alignment output files from the padded sam output output directory with the extension *.sam.
my ($padded_sam_files, $padded_sam_file_count) = find_files($padded_sam_output_dir, "sam");

# Iterate through each padded sam alignment output file with extension *.sam and execute the pstacks program.
my $pstacks_sql_id = 1;
foreach my $file_name (sort keys %{$padded_sam_files}){
	warn "Processing " . $file_name . ".....\n";
	my $padded_sam_infile = $padded_sam_files->{$file_name};
    
	# Execute the pstacks program, which extract stacks that have been aligned to a reference genome and identify SNPs. These stacks can then be processed with cstacks and sstacks.
	pstacks($padded_sam_infile, $pstacks_sql_id, $min_depth_coverage_pstacks, $num_threads, $alpha_value_pstacks, $stacks_output_dir);
	
	$pstacks_sql_id++;
}

# Execute the cstacks program to build a catalog from a set of samples processed by the pstacks program if a catalog tags file was not specified. The cstacks program creates a set of consensus loci, merging alleles together. If a catalog tags file was specified and the use_existing_catalog bitwise flag is set to true, then execute the cstacks program to add data to an existing catalog. If a catalog tags file was specified and the use_existing_catalog bitwise flag is set to false, then skip the execution of the cstacks program and use the catalog when executing sstacks.
my $cstacks_file = "";
if($cstacks_catalog_prefix eq ""){
    $cstacks_file = cstacks($stacks_output_dir, $cstacks_catalog_prefix, $stacks_sql_id, $num_mismatches_tag, $num_threads);
    
}elsif($cstacks_catalog_prefix =~ m/batch_\d+$/){
    
    if($use_existing_catalog eq "true"){
        $cstacks_file = cstacks($stacks_output_dir, $cstacks_catalog_prefix, $stacks_sql_id, $num_mismatches_tag, $num_threads);
    }elsif($use_existing_catalog eq "false"){
    	
	# The existing cstacks catalog alleles, snps, and tags file paths.
    	my ($existing_cstacks_alleles_file, $existing_cstacks_snps_file, $existing_cstacks_tags_file);
    	$existing_cstacks_alleles_file = $cstacks_catalog_prefix . '.catalog.alleles.tsv';
    	$existing_cstacks_snps_file = $cstacks_catalog_prefix . '.catalog.snps.tsv';
    	$existing_cstacks_tags_file = $cstacks_catalog_prefix . '.catalog.tags.tsv';
    	
    	# The cstacks batch file prefix file path.
    	$cstacks_file = join('/', $stacks_output_dir, join("_", "batch", $stacks_sql_id));
    	
    	# The cstacks catalog alleles, snps, and tags file paths.
    	my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file);
    	$cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
    	$cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
    	$cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';    
	
	warn "Copying $existing_cstacks_alleles_file to $cstacks_alleles_file.....";
        copy($existing_cstacks_alleles_file, $cstacks_alleles_file) or die "Copy failed: $!";

        warn "Copying $existing_cstacks_snps_file to $cstacks_snps_file.....";
        copy($existing_cstacks_snps_file, $cstacks_snps_file) or die "Copy failed: $!";

        warn "Copying $existing_cstacks_tags_file to $cstacks_tags_file.....";
        copy($existing_cstacks_tags_file, $cstacks_tags_file) or die "Copy failed: $!";
    }
}

# Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.
my ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");

# Iterate through each tags output file with extension *.tags.tsv and execute the sstacks program.
foreach my $file_name (sort keys %{$pstacks_tags_files}){
	if($file_name !~ m/batch_\d+\.catalog\.tags\.tsv/){ # If *.tags.tsv file does not match batch_*.catalog.tags.tsv.
		my $pstacks_tags_infile = $pstacks_tags_files->{$file_name};
		
		# Get the basename of the tags filename without the .tags.tsv extension.
		my $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv/);
		warn "Processing " . $pstacks_filename . ".....\n";
		my $sstacks_infile = join('/', $stacks_output_dir, $pstacks_filename);
		
		# Execute the sstacks program. Sets of stacks constructed by the pstacks program is searched against the catalog produced by cstacks.
		sstacks($cstacks_file, $sstacks_infile, $stacks_sql_id, $num_threads, $stacks_output_dir);
	}
}

# convert_refgen_bwa_input_format($renum_refgen_fasta_infile, $refgen_output_dir) - Converts the reference genome to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers.
#
# Input paramater(s):
#
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
#
# $refgen_output_dir - The reference genome output directory.
sub convert_refgen_bwa_input_format{
    
	# The reference genome input fasta file.
	my $refgen_fasta_infile = shift;
	die "Error lost the reference genome input fasta file." unless defined $refgen_fasta_infile;
	
	# The reference genome output directory.
	my $refgen_output_dir = shift;
	die "Error lost reference genome output directory" unless defined $refgen_output_dir;
    
	# Format the renumerated fasta and table of contents .toc files.
	my $refgen_fasta_filename = fileparse($refgen_fasta_infile);
	my $refgen_fasta_outfile = join('/', $refgen_output_dir, $refgen_fasta_filename . ".fasta");
	my $refgen_toc_outfile = join('/', $refgen_output_dir, $refgen_fasta_filename . ".toc.txt");
	
	# Generate the renumerated fasta and table of contents .toc files if the files are not already generated.
	unless(-s $refgen_fasta_outfile and -s $refgen_toc_outfile){
		warn "Converting $refgen_fasta_filename to BWA input format....\n\n";
		my $seq_counter = 0;
        my %fasta_header = ();
		my %fasta_seqs = ();
        
        # Parse the contents of the reference genome fasta file to filter.
        open(INFILE, "<$refgen_fasta_infile") or die "Couldn't open file $refgen_fasta_infile for reading, $!";
        while(<INFILE>){
            chomp $_;
            
            if($_ !~ /^$/){
                if($_ =~ /^>/){
                    # Parse stacks reference genome fasta input file for the reference genome.
                    $seq_counter++;
                    my $refgen_seq_id = substr($_, 1, length $_);
                    $fasta_header{$seq_counter}{"HEADER"} = $refgen_seq_id;
                    print $fasta_header{$seq_counter}{"HEADER"} . "\n";
                }elsif($_ =~ m/^[ACGTNRYSWKM]+/i){
                    my $sequence = $_;
                    push(@{$fasta_seqs{$seq_counter}{"SEQUENCE"}}, $sequence);
                }else{
                    die "Error: $refgen_fasta_infile is not in the correct format!\n$_";
                }
            }
        }
        close(INFILE) or die "Couldn't close file $refgen_fasta_infile";
        
 		my $num_digits = length($seq_counter - 1);
		open(REFGEN_FASTA_OUTFILE, ">$refgen_fasta_outfile") or die "Couldn't open file $refgen_fasta_outfile for writting, $!";
		open(REFGEN_TOC_OUTFILE, ">$refgen_toc_outfile") or die "Couldn't open file $refgen_toc_outfile for writting, $!";
		print REFGEN_TOC_OUTFILE join("\t", "bwa_input_fasta_header", "original_fasta_header") . "\n";
		foreach my $seq_num (sort {$a <=> $b} keys %fasta_seqs){
            #warn "Printing $seq_num";
 			my $fasta_header = $fasta_header{$seq_num}{"HEADER"};
            
 			my $sequence = join("", @{$fasta_seqs{$seq_num}{"SEQUENCE"}});
 			my $padded_zeros_length = ($num_digits - length($seq_num));
 			my $padded_seq_num = '0' x $padded_zeros_length . $seq_num;
			
            # 			warn join("\t", $padded_seq_num, $fasta_header) . "\n";
			print REFGEN_FASTA_OUTFILE join("\n", ">$padded_seq_num", $sequence) . "\n";
			print REFGEN_TOC_OUTFILE join("\t", $padded_seq_num, $fasta_header) . "\n";
            
		}
		close(REFGEN_FASTA_OUTFILE) or die "Couldn't close file $refgen_fasta_outfile";
		close(REFGEN_TOC_OUTFILE) or die "Couldn't close file $refgen_toc_outfile";
        
	}
	
	# Need to get the actual number of chromosomes or scaffolds for the -eC ($num_chromosomes option).
	open(REFGEN_TOC_INFILE, "<$refgen_toc_outfile") or die "Couldn't open file $refgen_toc_outfile for reading, $!";
	my $i = 0;
	my $num_chromosomes = 0;
	while(<REFGEN_TOC_INFILE>){
		chomp $_;
        # 		warn $_ . "\n";
		if($i ne 0){
			$num_chromosomes++;
		}
		$i++;
	}
	close(REFGEN_TOC_INFILE) or die "Couldn't close file $refgen_toc_outfile";
	
	return ($refgen_fasta_outfile, $refgen_toc_outfile, $num_chromosomes);
}

# bwa_index($renum_refgen_fasta_infile, $refgen_output_dir) - Executes the BWA alignment program using the index option to generate BWA index .amb .ann .bwt .pac .sa files from a renumerated reference genome file.
#
# Input paramater(s):
#
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
#
# $refgen_output_dir - The reference genome output directory that contains the .amb .ann .bwt .pac .sa index files.

#bwa index -a bwtsw ./refgen/renumbered_Msequence.fasta
sub bwa_index{

	# The renumerated reference genome input fasta file.
	my $renum_refgen_fasta_infile = shift;
	die "Error lost the renumbered reference genome input fasta file" unless defined $renum_refgen_fasta_infile;
    
	# The reference genome output directory that contains the .amb .ann .bwt .pac .sa index files.
	my $refgen_output_dir = shift;
	die "Error lost reference genome output directory" unless defined $refgen_output_dir;

	# Format the index file into .amb .ann .bwt .pac .sa files.
	my ($renum_refgen_fastadbAMB, $renum_refgen_fastadbANN, $renum_refgen_fastadbBWT, $renum_refgen_fastadbPAC, $renum_refgen_fastadbSA);
	$renum_refgen_fastadbAMB = $renum_refgen_fasta_infile . '.amb';
	$renum_refgen_fastadbANN = $renum_refgen_fasta_infile . '.ann';
	$renum_refgen_fastadbBWT = $renum_refgen_fasta_infile . '.bwt';
	$renum_refgen_fastadbPAC = $renum_refgen_fasta_infile . '.pac';
	$renum_refgen_fastadbSA = $renum_refgen_fasta_infile . '.sa';
	unless(-s $renum_refgen_fastadbAMB and -s $renum_refgen_fastadbANN and -s $renum_refgen_fastadbBWT 
		and -s $renum_refgen_fastadbPAC and -s $renum_refgen_fastadbSA){
		warn "Generating the bwa index file.....\n\n";
		my $bwaIndexCmd  = "$bwa index -a bwtsw $renum_refgen_fasta_infile";
		warn $bwaIndexCmd . "\n\n"; 
		system($bwaIndexCmd) == 0 or die "Error calling $bwaIndexCmd: $?";
	}
}

# $bwa_alignment_outfile = bwa_samse($renum_refgen_fasta_infile, $gbs_fastq_infile, $num_threads, $bwa_output_dir) - Executes the BWA alignment program using the aln option to generate BWA sai files from a renumerated reference genome and quality filtered and trimmed adapter GBS fastq input files.
#
# Input paramater(s):
#
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
#
# $gbs_fastq_infile - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
#
# $num_threads - The number of threads to use for BWA.
#
# $bwa_output_dir - The BWA alignment output directory that contains the sai alignment files.
#
# Output paramater(s):
#
# $bwa_alignment_outfile - The BWA sai alignment output file.

#bwa aln -t 4 ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/mpbGBSTags.cnt.fq > ./mergedTagCounts/AlignedGBSTags1.sai
sub bwa_aln{

    # The renumerated reference genome input fasta file.
    my $renum_refgen_fasta_infile = shift;
    die "Error lost the renumbered reference genome input fasta file" unless defined $renum_refgen_fasta_infile;

    # The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
    my $gbs_fastq_infile = shift;
    die "Error lost the GBS fastq input file" unless defined $gbs_fastq_infile;

    # The number of threads to use for BWA.
    my $num_threads = shift;
    die "Error lost the number of cores for BWA" unless defined $num_threads;

    # The BWA alignment output directory that contains the sai alignment files.
    my $bwa_output_dir = shift;
    die "Error lost the bwa output directory" unless defined $bwa_output_dir;

    my $individual_id = "";
    if($gbs_fastq_infile =~ m/trimmed_offset_\d+\.fastq/){
        # Get the basename of the fastq filename without the _trimmed_offset_\d+\.fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/_trimmed_offset_\d+\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }else{
        # Get the basename of the fastq filename without the .fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }

    # Format the sai file.
    my $bwa_alignment_outfile = join('/', $bwa_output_dir, $individual_id . ".sai");

    # Execute the BWA alignment program if the SAI file is not already generated.
    unless(-s $bwa_alignment_outfile){
        warn "Generating the bwa alignment file.....\n\n";
        my $bwaAlnCmd  = "$bwa aln -t $num_threads $renum_refgen_fasta_infile $gbs_fastq_infile > $bwa_alignment_outfile";
        warn $bwaAlnCmd . "\n\n";
        system($bwaAlnCmd) == 0 or die "Error calling $bwaAlnCmd: $?";
    }

    # Return the BWA sai alignment output file.
    return $bwa_alignment_outfile;
}

# $bwa_aligned_master_outfile = bwa_samse($renum_refgen_fasta_infile, $gbs_fastq_infile, $bwa_alignment_infile, $sam_mapq_threshold, $bwa_output_dir) - Executes the BWA alignment program using the samse option to generate sam alignment files from a renumerated reference genome, BWA sai, and quality filtered and trimmed adapter GBS fastq input files.
# 
# Input paramater(s):
# 
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
# 
# $gbs_fastq_infile - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
# 
# $bwa_alignment_infile - The BWA alignment index input file.
#
# $sam_mapq_threshold - The sam alignment file mapping quality score threshold.
#
# $bwa_output_dir - The BWA alignment output directory that contains the sam alignment files.
# 
# Output paramater(s):
# 
# $bwa_aligned_master_outfile - The BWA sam alignment output file.

#bwa samse ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/AlignedGBSTags1.sai ./mergedTagCounts/mpbGBSTags.cnt.fq > mergedTagCounts/AlignedMasterTagsMPB.sam
sub bwa_samse{

    # The renumerated reference genome input fasta file.
    my $renum_refgen_fasta_infile = shift;
    die "Error lost the renumbered reference genome input file" unless defined $renum_refgen_fasta_infile;

    # The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
    my $gbs_fastq_infile = shift;
    die "Error lost the GBS fastq input file" unless defined $gbs_fastq_infile;

    # The BWA alignment index input file.
    my $bwa_alignment_infile = shift;
    die "Error lost the BWA SAI formatted alignment (.sai) file" unless defined $bwa_alignment_infile;

    # The sam alignment file mapping quality score threshold.
    my $sam_mapq_threshold = shift;
    die "Error lost the sam alignment file mapping quality score threshold" unless defined $sam_mapq_threshold;

    # The BWA alignment output directory that contains the sam alignment files.
    my $bwa_output_dir = shift;
    die "Error lost the bwa output directory" unless defined $bwa_output_dir;

    my $individual_id = "";
    if($gbs_fastq_infile =~ m/trimmed_offset_\d+\.fastq/){
        # Get the basename of the fastq filename without the _trimmed_offset_\d+\.fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/_trimmed_offset_\d+\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }else{
        # Get the basename of the fastq filename without the .fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }

	# Execute the BWA alignment program if the sam alignment file is not already generated.
	my $bwa_aligned_master_outfile = join('/', $bwa_output_dir, $individual_id . ".sam");
	unless(-s $bwa_aligned_master_outfile){
		warn "Generating the bwa sam file.....\n\n";
        my $bwaSamseCmd  = "$bwa samse $renum_refgen_fasta_infile $bwa_alignment_infile $gbs_fastq_infile";
		warn $bwaSamseCmd . "\n\n";
        
        open(SAM_OUTFILE, ">$bwa_aligned_master_outfile") or die "Couldn't open file $bwa_aligned_master_outfile for writting, $!";
        local (*BWA_OUT, *BWA_IN);
        my $pid = open2(\*BWA_OUT,\*BWA_IN, $bwaSamseCmd) or die "Error calling open2: $!";
        close BWA_IN or die "Error closing STDIN to bwa process: $!";
        while(<BWA_OUT>){
            chomp $_;
            if($_ =~ m/^\@HD|^\@SQ|^\@RG|^\@PG|^\@CO/){ # Parse header lines starting with @HD, @SQ, @RG, @PG, or @CO.
				print SAM_OUTFILE $_ . "\n";
			}else{ #8_1308_8038_19954_1LL-06_TTCTG_5_A3length=66    16      7376    784845  37      66M     *       0       0       AATAATGCGCCAGCCAAGAGCTGTTTGGTAGCATTCTGCTGACGCTCGCACTCCCGTACACCTGCA      DDCB@>>>AFEFFFHC@GHGGGEGEGGIGCDEHBGGCIGIGGIIHGHEFCE:F<GIGIHFHHDFDD       XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:66
                
				my @split_sam_entry = split(/\t/, $_);
				my ($fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext,
                $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, @optional_fields) = @split_sam_entry;
				
				my $optional_fields = join("\t", @optional_fields);
				
                # die join("\t", $fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext, $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, $optional_fields);
                
				# Filter alignment entry if the read is unmapped in the sense and antisense orientation.
				next if(($bam_bitwise_flag eq 4) or ($bam_bitwise_flag eq 20));
				
				# Filter alignment entry if the mapping quality of the read is below this threshold.
				#next if($mapq < $sam_mapq_threshold);
				
				# Filter alignment entry if the cigar string contains an insertion/deletion in the alignment.
				# next if($cigar !~ m/^\d+M$/);
				
				# Filter alignment entry if the read is not uniquely mapped to the reference.
				next if($optional_fields !~ m/XT:A:U/);
				
				print SAM_OUTFILE $_ . "\n";
            }

        }
        close BWA_OUT or die "Error closing STDOUT from bwa process: $!";
        wait;
        close(SAM_OUTFILE) or die "Couldn't close file $bwa_aligned_master_outfile";
	}
    
	# The BWA sam alignment output file.
	return $bwa_aligned_master_outfile;
}


# bwa_pad_sam_files($sam_infile, $gbs_sequence_length, $sam_mapq_threshold, $padded_sam_output_dir) - Pads the BWA sam alignment files with poly-Ns so that all sequences are of uniform length.
#
# Input paramater(s):
#
# $sam_infile - The unpadded BWA sam alignment input file.
#
# $gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences.
#
# $sam_mapq_threshold - The sam alignment file mapping quality score threshold.
#
# $padded_sam_output_dir - The output directory that contains all the BWA sam alignment files.
sub bwa_pad_sam_files{

	# The unpadded BWA sam alignment input file.
	my $sam_infile = shift;
	die "Error lost the sam input file" unless defined $sam_infile;
	
	# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences.
	my $gbs_sequence_length = shift;
	die "Error lost the GBS fastq sequence length in base pairs (bps)" unless defined $gbs_sequence_length;

    # The sam alignment file mapping quality score threshold.
	my $sam_mapq_threshold = shift;
	die "Error lost the sam alignment file mapping quality score threshold" unless defined $sam_mapq_threshold;
    
	# The output directory that contains all the BWA sam alignment files.
	my $padded_sam_output_dir = shift;
	die "Error lost the padded sam output file directory" unless defined $padded_sam_output_dir;
	
	# Get the basename of the sam filename without the .sam extension.
	my $sam_filename = fileparse($sam_infile, qr/\.sam/);
	
	# Generate the padded sam output file if the file is not already generated.
	my $padded_sam_outfile = join('/', $padded_sam_output_dir, $sam_filename . ".sam");
	unless(-s $padded_sam_outfile){
		open(PADDED_SAM_OUTFILE, ">$padded_sam_outfile") or die "Couldn't open file $padded_sam_outfile for writting, $!";
		open(SAM_INFILE, "<$sam_infile") or die "Couldn't open file $sam_infile for reading, $!";
		while(<SAM_INFILE>){
			chomp $_;
	#  		warn $_ . "\n";
			if($_ =~ m/^\@HD|^\@SQ|^\@RG|^\@PG|^\@CO/){ # Parse header lines starting with @HD, @SQ, @RG, @PG, or @CO.
				print PADDED_SAM_OUTFILE $_ . "\n";
			}else{ #8_1308_8038_19954_1LL-06_TTCTG_5_A3length=66    16      7376    784845  37      66M     *       0       0       AATAATGCGCCAGCCAAGAGCTGTTTGGTAGCATTCTGCTGACGCTCGCACTCCCGTACACCTGCA      DDCB@>>>AFEFFFHC@GHGGGEGEGGIGCDEHBGGCIGIGGIIHGHEFCE:F<GIGIHFHHDFDD       XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:66

				my @split_sam_entry = split(/\t/, $_);
				my ($fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext, 
					$pnext, $tlen, $fastq_sequence, $fastq_quality_scores, @optional_fields) = @split_sam_entry;
				
				my $optional_fields = join("\t", @optional_fields);
				
				my $fastq_sequence_length = length($fastq_sequence);
				my $fastq_quality_scores_length = length($fastq_quality_scores);
# 				die join("\t", $fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext, $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, $optional_fields);

				# Filter alignment entry if the read is unmapped in the sense and antisense orientation.
				next if(($bam_bitwise_flag eq 4) or ($bam_bitwise_flag eq 20));
				
				# Filter alignment entry if the mapping quality of the read is below this threshold.
				next if($mapq < $sam_mapq_threshold);
				
				# Filter alignment entry if the cigar string contains an insertion/deletion in the alignment.
				next if($cigar !~ m/^\d+M$/);
				
				# Filter alignment entry if the read is not uniquely mapped to the reference.
				next if($optional_fields !~ m/XT:A:U/);
				
				if(($fastq_sequence_length ne $gbs_sequence_length) and ($fastq_quality_scores_length ne $gbs_sequence_length)){
                    
					# Get the number of poly-Ns to pad.
					my $padded_nucleotide_length = ($gbs_sequence_length - $fastq_sequence_length);
					my $padded_nucleotide_seq = 'N' x $padded_nucleotide_length;
					
					# Pad sequence based on alignment type.
					my $padded_fastq_sequence = "";
					
					# Pad sequence if aligned to sense strand.
					$padded_fastq_sequence = join("", $fastq_sequence, $padded_nucleotide_seq) if($bam_bitwise_flag eq 0);
					
					# Pad sequence if aligned to antisense strand.
					$padded_fastq_sequence = join("", $padded_nucleotide_seq, $fastq_sequence) if($bam_bitwise_flag eq 16);
					
					# Get the number of quality scores to pad.
					my $padded_score_length =  ($gbs_sequence_length - $fastq_quality_scores_length);
					my $padded_score_seq = '#' x $padded_score_length;
					
					# Pad quality scores based on alignment type.
					my $padded_fastq_quality_scores = "";
					
					# Pad quality scores if aligned to sense strand sequence in sense orientation.
					$padded_fastq_quality_scores = join("", $fastq_quality_scores, $padded_score_seq) if($bam_bitwise_flag eq 0);
					
					# Pad quality scores if aligned to antisense strand in antisense orientation.
					$padded_fastq_quality_scores = join("", $padded_score_seq, $fastq_quality_scores) if($bam_bitwise_flag eq 16);
					
					my $padded_fastq_sequence_length = length($padded_fastq_sequence);
					my $padded_fastq_quality_scores_length = length($padded_fastq_quality_scores);
					
					die "Error: $fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($padded_fastq_sequence_length ne $gbs_sequence_length);
					die "Error: $fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length ne padded_fastq_quality_scores_length=$padded_fastq_quality_scores_length" if($padded_fastq_sequence_length ne $padded_fastq_quality_scores_length);

					$split_sam_entry[9] = $padded_fastq_sequence;
					$split_sam_entry[10] = $padded_fastq_quality_scores;
					
					print PADDED_SAM_OUTFILE join("\t", @split_sam_entry) . "\n";
				}elsif(($fastq_sequence_length eq $gbs_sequence_length) and ($fastq_quality_scores_length eq $gbs_sequence_length)){
				
					print PADDED_SAM_OUTFILE $_ . "\n";
				}
	# 			die $_;
			}
		}
		close(SAM_INFILE) or die "Couldn't close file $sam_infile";
		close(PADDED_SAM_OUTFILE) or die "Couldn't close file $padded_sam_outfile";
	}

}

# pstacks($sam_infile, $sql_id, $min_depth_coverage, $num_threads, $alpha_value) - Executes the pstacks program in the Stacks Software Suite. Extracts stacks that have been aligned to a reference genome and identify SNPs. These stacks can then be processed with cstacks and sstacks.
# 
# Input paramater(s):
# 
# $sam_infile - The padded sam alignment input file.
#
# $sql_id - The SQL ID to insert into the output to identify this sample.
#
# $min_depth_coverage - The minimum depth of coverage to report a stack.
#
# $num_threads - The number of threads to use for pstacks.
#
# $alpha_value - The chi square significance level required to call a heterozygote or homozygote.

# pstacks -t sam -f ../PADDED_SAM_OUTFILES/LL-06_MPB-MALE-GBS.sam -o .. -i 1 -m 1 -p 7 --model_type snp --alpha 0.05
sub pstacks{

	# The padded sam alignment input file.
	my $sam_infile = shift;
	die "Error lost the sam alignment input file" unless defined $sam_infile;
	
	# The SQL ID to insert into the output to identify this sample.
	my $sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $sql_id;
	
	# The minimum depth of coverage to report a stack.
	my $min_depth_coverage = shift;
	die "Error lost the minimum depth of coverage to report a stack" unless defined $min_depth_coverage;
	
	# The number of threads to use for pstacks.
	my $num_threads = shift;
	die "Error lost the number of cores for stacks" unless defined $num_threads;
	
	# The chi square significance level required to call a heterozygote or homozygote.
	my $alpha_value = shift;
	die "Error lost the chi square significance level required to call a heterozygote or homozygote" unless defined $alpha_value; 
	
	# The stacks output directory that contains the results from the pstacks program.
	my $stacks_output_dir = shift;
	die "Error lost the stacks output file directory" unless defined $stacks_output_dir;
	
	# Get the basename of the sam filename without the .sam extension.
	my $sam_filename = fileparse($sam_infile, qr/\.sam/);
	
	# Format the pstacks individual alleles, snps, and tags output files.
	my ($pstacks_alleles_file, $pstacks_snps_file, $pstacks_tags_file);
	$pstacks_alleles_file = join('/', $stacks_output_dir, $sam_filename . '.alleles.tsv');
	$pstacks_snps_file = join('/', $stacks_output_dir, $sam_filename . '.snps.tsv');
	$pstacks_tags_file = join('/', $stacks_output_dir, $sam_filename . '.tags.tsv');
	
	#### USED PARAMETERS ####
	# t — input file Type. Supported types: bowtie, sam, or bam.
	# f — input file path.
	# o — output path to write results.
	# i — SQL ID to insert into the output to identify this sample.
	# m — minimum depth of coverage to report a stack (default 1).
	# p — enable parallel execution with num_threads threads.

	# Model options:
	# --model_type [type] — either 'snp' (default), 'bounded', or 'fixed'

	# For the SNP or Bounded SNP model:
	# --alpha [num] — chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.

	#### UNUSED PARAMETERS ####
	# h — display this help messsage.

	# For the Bounded SNP model:
	# --bound_low [num] — lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).
	# --bound_high [num] — upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).

	# For the Fixed model:
	# --bc_err_freq [num] — specify the barcode error frequency, between 0 and 1.0.

	# Create the stacks log output directory if it doesn't already exist.
	my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
	unless(-d $stacks_log_output_dir){
		mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
	}

	# The standard output file for the pstacks program.
	my $pstacks_log_outfile = join('/', $stacks_log_output_dir, join("_", $sam_filename, "pstacks.log"));
        
	# Execute the pstacks program if the pstacks alleles, snps, and tags output files are not already generated.
	unless(-s $pstacks_alleles_file and -s $pstacks_snps_file and -s $pstacks_tags_file){
		warn "Executing pstacks.....\n\n";
		my $pstacksCmd  = "$pstacks -t sam -f $sam_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -p $num_threads --model_type snp --alpha $alpha_value 2> $pstacks_log_outfile";
		warn $pstacksCmd . "\n\n";
		system($pstacksCmd) == 0 or die "Error calling $pstacksCmd: $?";
	}
}

# $cstacks_file = cstacks($stacks_output_dir, $cstacks_catalog_prefix, $stacks_sql_id, $num_mismatches_tag, $num_threads) - Executes the cstacks program in the Stacks Software Suite. Build a catalog from a set of samples processed by the ustacks program. Creates a set of consensus loci, merging alleles together.
# 
# Input paramater(s):
# 
# $stacks_output_dir - The stacks directory that contains the results from the pstacks program.
#
# $cstacks_catalog_prefix - The absolute path to an existing catalog. If path is not defined execute the default cstacks program, otherwise add data to an existing catalog.
#
# $stacks_sql_id - The SQL ID to insert into the output to identify this sample.
#
# $num_mismatches_tag - The number of mismatches allowed between sample tags when generating the catalog.
#
# $num_threads - The number of threads to use for sstacks.
# 
# Output paramater(s):
# 
# $cstacks_file - The cstacks output file prefix for the tab-delmited catalog alleles, snps, and tags files.

# cstacks -b 1 -o /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -g -n 1 -p 7 \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/LL-06_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M004-13-01-1D-G14-DA01_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M004-13-01-2D-G05-DA01_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M004-13-01-2F-G24-DA01_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M004-13-01-5ST-G36-DA02_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M024-07-01-02-DA07_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M033-06-01-05-DA126_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M033-07-01-09-DA36_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M033-11-02-05-DA41_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/NG-38_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/QU-48_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/RR-45_MPB-MALE-GBS
sub cstacks{

	# The stacks directory that contains the results from the pstacks program.
	my $stacks_output_dir = shift;
	die "Error lost the stacks output directory" unless defined $stacks_output_dir;
	
    # The absolute path to an existing catalog. If path is not defined execute the default cstacks program, otherwise add data to an existing catalog.
    my $cstacks_catalog_prefix = shift;
    die "Error lost he absolute path to an existing catalog. If path is not defined execute the default cstacks program, otherwise add data to an existing catalog" unless defined $cstacks_catalog_prefix;
    
	# The SQL ID to insert into the output to identify this sample.
	my $stacks_sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $stacks_sql_id;
	
	# The number of mismatches allowed between sample tags when generating the catalog.
	my $num_mismatches_tag = shift;
	die "Error lost the number of mismatches allowed between sample tags when generating the catalog" unless defined $num_mismatches_tag;
	
	# The number of threads to use for cstacks.
	my $num_threads = shift;
	die "Error lost the number of threads to use for cstacks" unless defined $num_threads;
	
    # Create the stacks log output directory if it doesn't already exist.
    my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
    unless(-d $stacks_log_output_dir){
        mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
    }
    
    # The standard out log file for the cstacks program.
    my $cstacks_log_outfile = join('/', $stacks_log_output_dir, "cstacks.log");
    
    if($cstacks_catalog_prefix eq ""){
        # The cstacks batch file prefix file path.
        my $cstacks_file = join('/', $stacks_output_dir, join("_", "batch", $stacks_sql_id));
        
        # The cstacks catalog alleles, snps, and tags file paths.
        my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file);
        $cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
        $cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
        $cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
        
        # Execute the cstacks program if the following files are not already generated.
        unless(-s $cstacks_alleles_file and -s $cstacks_snps_file and -s $cstacks_tags_file){
        
            # Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.
            my ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");
            my @cstacks_soptions = ();
            # Iterate through each pstacks tags output file with extension *.tags.tsv and execute the cstacks program.
            foreach my $file_name (sort keys %{$pstacks_tags_files}){
                my $pstacks_tags_infile = $pstacks_tags_files->{$file_name};
                
                # Get the basename of the tags filename without the .tags.tsv extension.
                my $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv/);
                
                # Obtain a list of all the file prefix paths specifed by the -s option for cstacks.
                warn "Processing " . $pstacks_filename . ".....\n";
                my $cstacks_infile = join('/', $stacks_output_dir, $pstacks_filename);
                push(@cstacks_soptions, "-s $cstacks_infile ");
            }
            
            #### USED PARAMETERS ####
            # b — MySQL ID of this batch.
            # o — output path to write results.
            # s — TSV file from which to load radtags.
            # g — base catalog matching on genomic location, not sequence identity.
            # n — number of mismatches allowed between sample tags when generating the catalog.
            # p — enable parallel execution with num_threads threads.

            #### UNUSED PARAMETERS ####
            # m — include tags in the catalog that match to more than one entry.
            # h — display this help messsage.

            # Catalog editing:
            # --catalog [path] — provide the path to an existing catalog. cstacks will add data to this existing catalog.

            # Advanced options:
            # --report_mmatches — report query loci that match more than one catalog locus.

            warn "Executing cstacks.....\n\n";
            my $cstacks_joined_soptions = join("\\\n", @cstacks_soptions);
            my $cstacksCmd  = "$cstacks -b $stacks_sql_id -o $stacks_output_dir -g -n $num_mismatches_tag -p $num_threads \\\n $cstacks_joined_soptions 2> $cstacks_log_outfile";
            warn $cstacksCmd . "\n\n";
            system($cstacksCmd) == 0 or die "Error calling $cstacksCmd: $?";
        }
    }elsif($cstacks_catalog_prefix =~ m/batch_\d+$/){
        
        
        # The existing cstacks catalog alleles, snps, and tags file paths.
        my ($existing_cstacks_alleles_file, $existing_cstacks_snps_file, $existing_cstacks_tags_file);
        $existing_cstacks_alleles_file = $cstacks_catalog_prefix . '.catalog.alleles.tsv';
        $existing_cstacks_snps_file = $cstacks_catalog_prefix . '.catalog.snps.tsv';
        $existing_cstacks_tags_file = $cstacks_catalog_prefix . '.catalog.tags.tsv';
        
        # The cstacks batch file prefix file path.
        my $cstacks_file = join('/', $stacks_output_dir, join("_", "batch", $stacks_sql_id));
        
        # The cstacks catalog alleles, snps, and tags file paths.
        my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file);
        $cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
        $cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
        $cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
        
        warn "Copying $existing_cstacks_alleles_file to $cstacks_alleles_file.....";
        copy($existing_cstacks_alleles_file, $cstacks_alleles_file) or die "Copy failed: $!";
        
        warn "Copying $existing_cstacks_snps_file to $cstacks_snps_file.....";
        copy($existing_cstacks_snps_file, $cstacks_snps_file) or die "Copy failed: $!";
        
        warn "Copying $existing_cstacks_tags_file to $cstacks_tags_file.....";
        copy($existing_cstacks_tags_file, $cstacks_tags_file) or die "Copy failed: $!";
        
        # Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.
        my ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");
        my @cstacks_soptions = ();
        # Iterate through each pstacks tags output file with extension *.tags.tsv and execute the cstacks program.
        foreach my $file_name (sort keys %{$pstacks_tags_files}){
            my $pstacks_tags_infile = $pstacks_tags_files->{$file_name};
            
            # Get the basename of the tags filename without the .tags.tsv extension.
            my $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv/);
            
            # Obtain a list of all the file prefix paths specifed by the -s option for cstacks.
            warn "Processing " . $pstacks_filename . ".....\n";
            my $cstacks_infile = join('/', $stacks_output_dir, $pstacks_filename);
            push(@cstacks_soptions, "-s $cstacks_infile ");
        }
        
        #### USED PARAMETERS ####
        # b — MySQL ID of this batch.
        # o — output path to write results.
        # s — TSV file from which to load radtags.
        # g — base catalog matching on genomic location, not sequence identity.
        # n — number of mismatches allowed between sample tags when generating the catalog.
        # p — enable parallel execution with num_threads threads.
        
        #### UNUSED PARAMETERS ####
        # m — include tags in the catalog that match to more than one entry.
        # h — display this help messsage.
        
        # Catalog editing:
        # --catalog [path] — provide the path to an existing catalog. cstacks will add data to this existing catalog.
        
        # Advanced options:
        # --report_mmatches — report query loci that match more than one catalog locus.
        
        warn "Executing cstacks.....\n\n";
        my $cstacks_joined_soptions = join("\\\n", @cstacks_soptions);
        my $cstacksCmd  = "$cstacks --catalog $cstacks_file -b $stacks_sql_id -o $stacks_output_dir -g -n $num_mismatches_tag -p $num_threads \\\n $cstacks_joined_soptions 2> $cstacks_log_outfile";
        warn $cstacksCmd . "\n\n";
        system($cstacksCmd) == 0 or die "Error calling $cstacksCmd: $?";
    }
	# Returns the cstacks batch file prefix file path.
	return $cstacks_file;
}

# sstacks($stacks_catalog_infile, $stacks_sample_infile, $stacks_sql_id, $num_threads, $stacks_output_dir) - Executes the sstacks program in the Stacks Software Suite. Sets of stacks constructed by the pstacks program is searched against the catalog produced by cstacks.
# 
# Input paramater(s):
# 
# $stacks_catalog_infile - The stacks catalog input batch_*.catalog.*.tsv file from which to load the catalog GBS-Tags.
# 
# $stacks_sample_infile - The stacks sample input *.tags.tsv file from which to load sample GBS-Tags.
#
# $stacks_sql_id - The SQL ID to insert into the output to identify this sample.
# 
# $num_threads - The number of threads to use for sstacks.
# 
# $stacks_output_dir - The stacks output file directory that contains all the results files generated by sstacks.
sub sstacks{
    
	# The stacks catalog input batch_*.catalog.*.tsv file from which to load the catalog GBS-Tags.
	my $stacks_catalog_infile = shift;
	die "Error lost the stacks catalog input file" unless defined $stacks_catalog_infile;
    
	# The stacks sample input *.tags.tsv file from which to load sample GBS-Tags.
	my $stacks_sample_infile = shift;
	die "Error lost the stacks sample input file" unless defined $stacks_sample_infile;
	
	# The SQL ID to insert into the output to identify this sample.
	my $stacks_sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $stacks_sql_id;
	
	# The number of threads to use for sstacks.
	my $num_threads = shift;
	die "Error lost the number of threads to use for sstacks" unless defined $num_threads;
	
	# The stacks output file directory that contains all the results files generated by sstacks.
	my $stacks_output_dir = shift;
	die "Error lost the stacks output file directory that contains all files generated by sstacks" unless defined $stacks_output_dir;
	
	# Get the basename of the sam filename without the .sam extension.
	my $stacks_sample_filename = fileparse($stacks_sample_infile, qr//);

	# The sstacks output matches file.
	my $sstacks_matches_file = join('/', $stacks_output_dir, $stacks_sample_filename . '.matches.tsv');
	
	# Create the stacks log output directory if it doesn't already exist.
	my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
	unless(-d $stacks_log_output_dir){
		mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
	}

	# The standard out log file for the sstacks program.
	my $sstacks_log_outfile = join('/', $stacks_log_output_dir, join("_", $stacks_sample_filename, "sstacks.log"));
        
	#### USED PARAMETERS ####
	# b — MySQL ID of this batch.
	# c — TSV file from which to load the catalog RAD-Tags.
	# s — TSV file from which to load sample RAD-Tags.
	# o — output path to write results.
	# g — base matching on genomic location, not sequence identity.
	# p — enable parallel execution with num_threads threads.

	#### UNUSED PARAMETERS ####
	# r — Load the TSV file of a single sample instead of a catalog.
	# x — don’t verify haplotype of matching locus.
	# v — print program version.
	# h — display this help messsage.

	# Execute the sstacks program if the matches sstacks results files are not already generated.
	unless(-s $sstacks_matches_file){
		warn "Executing sstacks.....\n\n";
		my $sstacksCmd  = "$sstacks -b $stacks_sql_id -c $stacks_catalog_infile -s $stacks_sample_infile -o $stacks_output_dir -g -p $num_threads 2> $sstacks_log_outfile";
		warn $sstacksCmd . "\n\n";
		system($sstacksCmd) == 0 or die "Error calling $sstacksCmd: $?";
	}
}

# (\%files, $file_counter) = find_files($infile_dir) - Find all files in the specified input file directory with the file extension *.suffix.
# 
# Input paramater(s):
# 
# $infile_dir - The input file directory.
# 
# $suffix - The file extension suffix.
# 
# Output paramater(s):
# 
# \%files - A hash reference containing all the files with file extension *.suffix in key/value pairs.
# 
# key => filename ( e.g. filename.suffix )
# value => absolue filepath ( e.g. /path/to/filename.suffix )
# 
# $file_count - The number of files stored with file extension *.suffix.
sub find_files{
    
	# The input file directory.
	my $infile_dir = shift;
	die "Error lost input file directory" unless defined $infile_dir;
	
	# The file extension suffix.
	my $suffix = shift;
	die "Error lost file extension suffix directory" unless defined $suffix;
	
	if(-d $infile_dir){ # Check if $infile_dir is a directory.
		my %files = ();
		my $file_counter = 0;
		opendir(DIR, $infile_dir) || die "Error in opening dir $infile_dir\n";
		while( my $file_name = readdir(DIR)){
			my $infile_name = join('/', $infile_dir, $file_name) if ($file_name =~ m/\.$suffix$/);
			warn "$infile_name\n" if ($file_name =~ m/\.$suffix$/);
			$files{$file_name} = $infile_name if ($file_name =~ m/\.$suffix$/);
			$file_counter++ if ($file_name =~ m/\.$suffix$/);
		}
		closedir(DIR);
		
		return (\%files, $file_counter);
	}else{
		die "Error $infile_dir does not exist!\n";
	}
}

# $output_dir = gunzip_fastq_file($fastq_file) - Execute the gunzip program to uncompress the compressed fastq file.
#
# Input paramater(s):
#
# $fastq_file - The fastq file to uncompress using gunzip.
#
# Output paramater(s):
#
# $uncompressed_fastq_file - The uncompressed fastq file path.
sub gunzip_fastq_file{
	
	# The fastq file to uncompress using gunzip.
	my $fastq_file = shift;
	die "Error lost the fastq file to uncompress using gunzip" unless defined $fastq_file;

	my ($fastq_filename, $fastq_dir, $individual_id) = "";
	if($fastq_file =~ m/trimmed_offset_\d+\.fastq\.gz/){
        	# get the basename of the fastq filename without the _trimmed_offset_\d+\.fastq extension.
        	($fastq_filename, $fastq_dir) = fileparse($fastq_file, qr/_trimmed_offset_\d+\.fastq\.gz/);
        	$individual_id = $fastq_filename;
        }else{
		# get the basename of the fastq filename without the \.fastq extension.
		($fastq_filename, $fastq_dir) = fileparse($fastq_file, qr/\.fastq/);
		$individual_id = $fastq_filename;
	}         	
	my $uncompressed_fastq_file = join('/', $fastq_dir, $individual_id . ".fastq");
	
	unless(-s $uncompressed_fastq_file){
		warn "calling gunzip for $fastq_file....\n";
		my $gunzipCmd  = "$gunzip -c $fastq_file > $uncompressed_fastq_file";
		warn $gunzipCmd . "\n\n";
		system($gunzipCmd) == 0 or die "Error calling $gunzipCmd: $?";
	}
	return $uncompressed_fastq_file;
}
