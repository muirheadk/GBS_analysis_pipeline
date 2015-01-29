#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

#### PROGRAM NAME ####
# refgen_stacks_analysis_pipeline.pl - Program that aligns quality filtered, demultiplexed, and adapter trimmed GBS data sequences (unpadded) against a reference genome using the BWA alignment program. It then runs the Stacks reference genome pipeline using the pstacks, cstacks, and sstacks programs of the Stacks Software Suite.

#### DESCRIPTION ####
# This program takes the quality filtered, demultiplexed, and adapter trimmed *.fastq input files (unpadded) and a reference genome fasta input file as input. Converts a reference genome fasta file to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers. It then performs a BWA alignment to align the GBS fastq sequences to a reference genome. It then executes the pstacks program, which extracts sequence stacks that were aligned to a reference genome using the BWA alignment program and identifies SNPs. These sequence stacks are then processed using cstacks and sstacks to obtain the filtered SNP stacks output files.

#### SAMPLE COMMAND ####
# perl refgen_stacks_analysis_pipeline.pl -i ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR_UNPADDED/STEPHEN_TREVOY/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_sequence_data/DendPond_male_1.0/Primary_Assembly/unplaced_scaffolds/FASTA/DendPond_male_1.0_unplaced.scaf.fa -c 7 -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3
my ($gbs_fastq_dir, $gbs_fastq_file_type, $refgen_infile, $gbs_sequence_length, $stacks_sql_id, $min_depth_coverage_pstacks, $alpha_value_pstacks, $num_threads, $output_dir);
GetOptions(
	'i=s'    => \$gbs_fastq_dir, # The absolute path to the quality filtered, demultiplexed, and adapter trimmed *.fastq input file (unpadded) directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	't=s'    => \$gbs_fastq_file_type, # The fastq input file type. Default: gzfastq
	'g=s'    => \$refgen_infile, # The absolute path to the reference genome input fasta file to align GBS fastq sequences.
	'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
	'b=s'    => \$stacks_sql_id, # The SQL ID to insert into the output to identify this sample. Default: 1
	'd=s'    => \$min_depth_coverage_pstacks, # The minimum depth of coverage to report a stack. Default: 1
	'a=s'    => \$alpha_value_pstacks, # The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
	'c=s'    => \$num_threads, # The number of cpu threads to use for the stacks programs. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2
	'o=s'    => \$output_dir, # The absolute path to the output directory to contain the renumerated reference genome, BWA sam, padded sam, and Stacks output files and directories.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
	defined $gbs_fastq_dir
	and defined $refgen_infile
	and defined $output_dir
);

# The fastq input file type. Default: gzfastq
$gbs_fastq_file_type = 'gzfastq' unless defined $gbs_fastq_file_type;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

# The SQL ID to insert into the output to identify this sample.
$stacks_sql_id = 1 unless defined $stacks_sql_id;

# The minimum depth of coverage to report a stack. Default: 1
$min_depth_coverage_pstacks = 1 unless defined $min_depth_coverage_pstacks;

# The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
$alpha_value_pstacks = 0.05 unless defined $alpha_value_pstacks;

# The number of cpu threads to use for the stacks programs. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2
$num_threads = 2 unless defined $num_threads;

# Program dependencies - The absolute paths to gunzip to uncompress fastq.gz input files, BWA alignment program, and Stacks related programs.
my ($gunzip, $bwa, $pstacks, $cstacks, $sstacks);
$gunzip				= '/bin/gunzip';
$bwa				= '/usr/local/bin/bwa';
$pstacks			= '/usr/local/bin/pstacks';
$cstacks			= '/usr/local/bin/cstacks';
$sstacks			= '/usr/local/bin/sstacks';

sub usage {

die <<"USAGE";

Usage: $0 -i gbs_fastq_dir -t gbs_fastq_file_type -g refgen_infile -l gbs_sequence_length -b stacks_sql_id -d min_depth_coverage_pstacks -a alpha_value_pstacks -c num_cpu_cores -o output_dir

DESCRIPTION - This program takes the quality filtered, demultiplexed, and adapter trimmed *.fastq input files (unpadded) and reference genome fasta input files as input. Converts the reference genome fasta file to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers. It then performs a BWA alignment to align the GBS fastq sequences to the reference genome. It then executes the pstacks program, which extracts sequence stacks that were aligned to the reference genome using the BWA alignment program and identifies SNPs. These sequence stacks are then processed using cstacks and sstacks to obtain the filtered SNP stacks output files.

OPTIONS:

-i gbs_fastq_dir - The absolute path to the quality filtered, demultiplexed, and adapter trimmed *.fastq input file (unpadded) directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.

-t gbs_fastq_file_type - The fastq input file type. Default: gzfastq

-g refgen_infile - The absolute path to the reference genome input fasta file to align GBS fastq sequences using the BWA alignment program.

-l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92

-b stacks_sql_id - The SQL ID to insert into the output to identify this sample.

-d min_depth_coverage_pstacks - The minimum depth of coverage to report a stack. Default: 1

-a alpha_value_pstacks - The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05

-c num_cpu_cores - The number of cpu cores to use for the stacks programs. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2

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
    
	# If the bulk fastq file is compressed, uncompress the file and set the resulting fastq filename to be the fastq infile.
	if($gbs_fastq_file_type eq "gzfastq"){
		my $uncompressed_fastq_file = gunzip_fastq_file($gbs_fastq_infile);
		$gbs_fastq_infile = $uncompressed_fastq_file;
	}

	# Creates the BWA alignment file using the GBS fastq input file to align the fastq sequence reads to the reference genome.
	my $bwa_alignment_outfile = bwa_aln($refgen_fasta_outfile, $gbs_fastq_infile, $num_threads, $bwa_output_dir);

	# Creates the BWA single-ended alignment file in sam format using the GBS fastq input file to align the fastq sequence reads to the reference genome and the BWA aln format file.
	my $bwa_aligned_master_outfile = bwa_samse($refgen_fasta_outfile, $gbs_fastq_infile, $bwa_alignment_outfile, $bwa_output_dir);
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
	bwa_pad_sam_files($bwa_align_infile, $gbs_sequence_length, $padded_sam_output_dir);
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

# Execute the cstacks program to build a catalog from a set of samples processed by the pstacks program. The cstacks program creates a set of consensus loci, merging alleles together.
my $cstacks_file = cstacks($stacks_output_dir, $stacks_sql_id, $num_threads);

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
		my $seqio = Bio::SeqIO->new(-file => $refgen_fasta_infile, '-format' => 'Fasta');
		my $seq_counter = 1;
		my %fasta_seqs = ();
		while(my $seq_entry = $seqio->next_seq) {

			my $seq_id = $seq_entry->id;
			my $sequence = $seq_entry->seq;
			my $seq_desc = $seq_entry->desc;
		
			my $fasta_header = join(" ", $seq_id, $seq_desc);
			$fasta_seqs{$seq_counter} = join("\t", $fasta_header, $sequence);
			$seq_counter++;
		
		}
		
 		my $num_digits = length($seq_counter - 1);
		open(REFGEN_FASTA_OUTFILE, ">$refgen_fasta_outfile") or die "Couldn't open file $refgen_fasta_outfile for writting, $!";
		open(REFGEN_TOC_OUTFILE, ">$refgen_toc_outfile") or die "Couldn't open file $refgen_toc_outfile for writting, $!"; 
		print REFGEN_TOC_OUTFILE join("\t", "bwa_input_fasta_header", "original_fasta_header") . "\n";
		foreach my $seq_num (sort {$a <=> $b} keys %fasta_seqs){
 			my ($fasta_header, $sequence) = split(/\t/, $fasta_seqs{$seq_num});
 			my $padded_zeros_length = ($num_digits - length($seq_num));
 			my $padded_seq_num = '0' x $padded_zeros_length . $seq_num;
			
# 			warn join("\t", $padded_seq_num, $fasta_header) . "\n";
			print REFGEN_FASTA_OUTFILE join("\n", ">$padded_seq_num", $sequence) . "\n";
			print REFGEN_TOC_OUTFILE join("\t", $padded_seq_num, $fasta_header) . "\n";

		}
		close(REFGEN_FASTA_OUTFILE) or die "Couldn't close file $refgen_fasta_outfile";
		close(REFGEN_TOC_OUTFILE) or die "Couldn't close file $refgen_toc_outfile";

		# Clean out the sequence I/O object.
		$seqio = ();
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

	# Get the basename of the fastq filename without the .fastq extension.
	my $gbs_fastq_filename = fileparse($gbs_fastq_infile, qr/\.fastq/);
	
	# Split the GBS fastq filename so that we can get the individual id.
	my @split_gbs_fastq_filename = split(/_/, $gbs_fastq_filename);
	my $individual_id = $split_gbs_fastq_filename[0];
    
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

# $bwa_aligned_master_outfile = bwa_samse($renum_refgen_fasta_infile, $gbs_fastq_infile, $bwa_alignment_infile, $bwa_output_dir) - Executes the BWA alignment program using the samse option to generate sam alignment files from a renumerated reference genome, BWA sai, and quality filtered and trimmed adapter GBS fastq input files.
# 
# Input paramater(s):
# 
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
# 
# $gbs_fastq_infile - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
# 
# $bwa_alignment_infile - The BWA alignment index input file.
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

	# The BWA alignment output directory that contains the sam alignment files.
	my $bwa_output_dir = shift;
	die "Error lost the bwa output directory" unless defined $bwa_output_dir;

	# Get the basename of the fastq filename without the .fastq extension.
	my $gbs_fastq_filename = fileparse($gbs_fastq_infile, qr/\.fastq/);
	
	# Split the GBS fastq filename so that we can get the individual id.
	my @split_gbs_fastq_filename = split(/_/, $gbs_fastq_filename);
	my $individual_id = $split_gbs_fastq_filename[0];

	# Execute the BWA alignment program if the sam alignment file is not already generated.
	my $bwa_aligned_master_outfile = join('/', $bwa_output_dir, $individual_id . ".sam");
	unless(-s $bwa_aligned_master_outfile){
		warn "Generating the bwa sam file.....\n\n";
		my $bwaSamseCmd  = "$bwa samse $renum_refgen_fasta_infile $bwa_alignment_infile $gbs_fastq_infile > $bwa_aligned_master_outfile";
		warn $bwaSamseCmd . "\n\n";
		system($bwaSamseCmd) == 0 or die "Error calling $bwaSamseCmd: $?";
	}
    
	# The BWA sam alignment output file.
	return $bwa_aligned_master_outfile;
}


# bwa_pad_sam_files($sam_infile, $gbs_sequence_length, $padded_sam_output_dir) - Pads the BWA sam alignment files with poly-Ns so that all sequences are of uniform length.
# 
# Input paramater(s):
# 
# $sam_infile - The unpadded BWA sam alignment input file.
#
# $gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences.
#
# $padded_sam_output_dir - The output directory that contains all the BWA sam alignment files.
sub bwa_pad_sam_files{

	# The unpadded BWA sam alignment input file.
	my $sam_infile = shift;
	die "Error lost the sam input file" unless defined $sam_infile;
	
	# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences.
	my $gbs_sequence_length = shift;
	die "Error lost the GBS fastq sequence length in base pairs (bps)" unless defined $gbs_sequence_length;

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
				my ($fastq_header, $bam_bitwise_flag, $fastq_sequence, $fastq_quality_scores) = ($split_sam_entry[0], $split_sam_entry[1], $split_sam_entry[9], $split_sam_entry[10]);
				
				my $fastq_sequence_length = length($fastq_sequence);
				my $fastq_quality_scores_length = length($fastq_quality_scores);
# 				die join("\t", $fastq_sequence_length, $fastq_header, $bam_bitwise_flag, $fastq_sequence, $fastq_quality_scores);
				if(($fastq_sequence_length ne $gbs_sequence_length) and ($fastq_quality_scores_length ne $gbs_sequence_length)){
                    
					# Get the number of poly-Ns to pad.
					my $padded_nucleotide_length =  ($gbs_sequence_length - $fastq_sequence_length);
					my $padded_nucleotide_seq = 'N' x $padded_nucleotide_length;
					
					# Pad sequence based on alignment type.
					my $padded_fastq_sequence = "";
					
					# Pad sequence if aligned to sense strand or unmatched sequence.
					$padded_fastq_sequence = join("", $fastq_sequence, $padded_nucleotide_seq) if(($bam_bitwise_flag eq 0) or ($bam_bitwise_flag eq 4));
					
					# Pad sequence if aligned to antisense strand.
					$padded_fastq_sequence = join("", $padded_nucleotide_seq, $fastq_sequence) if(($bam_bitwise_flag eq 16) or ($bam_bitwise_flag eq 20));
					
					# Get the number of quality scores to pad.
					my $padded_score_length =  ($gbs_sequence_length - $fastq_quality_scores_length);
					my $padded_score_seq = '#' x $padded_score_length;
					
					# Pad quality scores based on alignment type.
					my $padded_fastq_quality_scores = "";
					
					# Pad quality scores if aligned to sense strand or unmatched sequence in sense orientation.
					$padded_fastq_quality_scores = join("", $fastq_quality_scores, $padded_score_seq) if(($bam_bitwise_flag eq 0) or ($bam_bitwise_flag eq 4));
					
					# Pad quality scores if aligned to antisense strand or unmatched sequence in antisense orientation.
					$padded_fastq_quality_scores = join("", $padded_score_seq, $fastq_quality_scores) if(($bam_bitwise_flag eq 16) or ($bam_bitwise_flag eq 20));
					
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

	# The standard out file for the pstacks program.
	my $pstacks_log_outfile = join('/', $stacks_log_output_dir, "pstacks.log");
        
	# Execute the pstacks program if the pstacks alleles, snps, and tags output files are not already generated.
	unless(-s $pstacks_alleles_file and -s $pstacks_snps_file and -s $pstacks_tags_file){
		warn "Executing pstacks.....\n\n";
		my $pstacksCmd  = "$pstacks -t sam -f $sam_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -p $num_threads --model_type snp --alpha $alpha_value 2> $pstacks_log_outfile";
		warn $pstacksCmd . "\n\n";
		system($pstacksCmd) == 0 or die "Error calling $pstacksCmd: $?";
	}
}

# $cstacks_file = cstacks($stacks_output_dir, $stacks_sql_id, $num_threads) - Executes the cstacks program in the Stacks Software Suite. Build a catalog from a set of samples processed by the ustacks program. Creates a set of consensus loci, merging alleles together.
# 
# Input paramater(s):
# 
# $stacks_output_dir - The stacks directory that contains the results from the pstacks program.
#
# $stacks_sql_id - The SQL ID to insert into the output to identify this sample.
#
# $num_threads - The number of threads to use for sstacks.
# 
# Output paramater(s):
# 
# $cstacks_file - The cstacks output file prefix for the tab-delmited catalog alleles, snps, and tags files.

# cstacks -b 1 -o /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -g -p 7 \
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
	
	# The SQL ID to insert into the output to identify this sample.
	my $stacks_sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $stacks_sql_id;
	
	# The number of threads to use for cstacks.
	my $num_threads = shift;
	die "Error lost the number of threads to use for cstacks" unless defined $num_threads;
	
	# The cstacks batch file prefix file path.
	my $cstacks_file = join('/', $stacks_output_dir, join("_", "batch", $stacks_sql_id));
	
	# The cstacks catalog alleles, snps, and tags file paths.
	my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file);
	$cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
	$cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
	$cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
	
	# Create the stacks log output directory if it doesn't already exist.
	my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
	unless(-d $stacks_log_output_dir){
		mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
	}

	# The standard out log file for the cstacks program.
	my $cstacks_log_outfile = join('/', $stacks_log_output_dir, "cstacks.log");
	
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
		# p — enable parallel execution with num_threads threads.

		#### UNUSED PARAMETERS ####
		# m — include tags in the catalog that match to more than one entry.
		# n — number of mismatches allowed between sample tags when generating the catalog.
		# h — display this help messsage.

		# Catalog editing:
		# --catalog [path] — provide the path to an existing catalog. cstacks will add data to this existing catalog.

		# Advanced options:
		# --report_mmatches — report query loci that match more than one catalog locus.

		warn "Executing cstacks.....\n\n";
		my $cstacks_joined_soptions = join("\\\n", @cstacks_soptions);
		my $cstacksCmd  = "$cstacks -b $stacks_sql_id -o $stacks_output_dir -g -p $num_threads \\\n $cstacks_joined_soptions 2> $cstacks_log_outfile";
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

	# The standard out log file for the cstacks program.
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
	
	if(-d $input_dir){ # Check if $input_dir is a directory.
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
		die "Error $input_dir does not exist!\n";
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
	
	my ($fastq_filename, $fastq_dir) = fileparse($fastq_file, qr/\.fastq.gz/);

	# Split the fastq filename so that we can get the individual id.
	my @split_fastq_filename = split(/_/, $fastq_filename);
	my $individual_id = $split_fastq_filename[0];
	
	my $uncompressed_fastq_file = join('/', $fastq_dir, $individual_id . ".fastq");
	
	unless(-s $uncompressed_fastq_file){
		warn "Calling gunzip for $fastq_file....\n";
		my $gunzipCmd  = "$gunzip -c $fastq_file > $uncompressed_fastq_file";
		warn $gunzipCmd . "\n\n";
		system($gunzipCmd) == 0 or die "Error calling $gunzipCmd: $?";
	}
	return $uncompressed_fastq_file;
}
