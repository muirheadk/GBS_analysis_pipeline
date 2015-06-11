#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use File::Copy;

#### PROGRAM NAME ####
# fastq_quality_barcode_splitter.pl - A program that performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.

#### DESCRIPTION ####
# This program performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.

#### SAMPLE COMMAND ####
# ApeKI paired-end reads command using phred33 encoding
# perl fastq_quality_barcode_splitter.pl -1 ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/7_SBW_lane2_lane3_Pool_R1.fastq.gz -2 ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/7_SBW_lane2_lane3_Pool_R2.fastq.gz -e phred33 -r ApeKI -b ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/stacks_barcodes_ApeKI_paired-end_lane2_lane3-2015-01-10.txt -c inline_null -n 0 -o ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/ApeKI_paired-end_PROCESSED_RADTAGS

# ApeKI single-end reads commands using phred64 encoding
## ApeKI lane 1
# perl fastq_quality_barcode_splitter.pl -i ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/s_1_Lane1-1_sequence.txt.gz -e phred64 -r ApeKI -b ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/stacks_barcodes_ApeKI_lane1-2015-01-10.txt -c inline_null -n 0 -o ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/ApeKI_single-end_lane1_PROCESSED_RADTAGS

## ApeKI lane 2
# perl fastq_quality_barcode_splitter.pl -i ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/s_2_Lane2-1_sequence.txt.gz -e phred64 -r ApeKI -b ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/stacks_barcodes_ApeKI_single-end_lane2-2015-01-10.txt -c inline_null -n 0 -o ~/workspace/GBS_data-08-10-2013/ApeKI_GBS_Data/ApeKI_single-end_lane2_PROCESSED_RADTAGS

# PstI/MspI single-end reads command using phred33 encoding
# perl fastq_quality_barcode_splitter.pl -i ~/workspace/GBS_data-08-10-2013/PstI_MspI_GBS_Data/HI.1405.008.GQ03122013-5_R1.fastq.gz -e phred33 -r PstI/MspI -b ~/workspace/GBS_data-08-10-2013/PstI_MspI_GBS_Data/PstI_MspI_GBS_barcodes-2013-10-09.csv -n 0 -o ~/workspace/GBS_data-08-10-2013/PstI_MspI_GBS_Data/PstI_MspI_PROCESSED_RADTAGS
my ($single_end_fastq_infile, $first_paired_end_fastq_infile, $second_paired_end_fastq_infile, $restriction_enzymes, $encoded_phred_offset, $sliding_window_size, $quality_score_limit, $gbs_sequence_length, $barcode_infile, $barcode_option, $num_mismatches, $output_dir);
GetOptions(
	'i=s'    => \$single_end_fastq_infile, # The absolute path to the bulk fastq input file to split sequences based on barcode sequences if processing single-end seqeunces. Can either be *.fastq or *.fastq.gz extension.
	'1=s'    => \$first_paired_end_fastq_infile, # The absolute path to the first bulk fastq input file in a set of paired-end sequences to split sequences based on barcode sequences. Can either be *.fastq or *.fastq.gz extension.
	'2=s'    => \$second_paired_end_fastq_infile, # The absolute path to the second bulk fastq input file in a set of paired-end sequences to split sequences based on barcode sequences. Can either be *.fastq or *.fastq.gz extension.
	'r=s'    => \$restriction_enzymes, # The restriction enzyme(s) used to digest the genomic sequences. Default: PstI/MspI
	'e=s'    => \$encoded_phred_offset, # The fastq quality score encoding used in the Illumina sequencing run.  Use phred33 for Illumina 1.8+ and Sanger or phred64 for Illumina 1.3 to 1.5. Default: phred33
	'w=s'    => \$sliding_window_size, # The size of the sliding window as a fraction of the read length between 0 and 1. Default: 0.15
	's=s'    => \$quality_score_limit, # The quality score limit. If the average score within the sliding window drops below this value, the read is discarded. Default: 20
	't=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
	'b=s'    => \$barcode_infile, # The absolute path to the barcodes input file used to split sequences into individual fastq output files.
	'c=s'    => \$barcode_option, # The barcode option for whether or not single-end or paired-end barcodes are within the FASTQ header or inline with sequence. Default: inline_null
	'n=s'    => \$num_mismatches, # The number of mismatches allowed within the barcode sequences. Default: 0
	'o=s'    => \$output_dir, # The absolute path to the output directory to contain the split *.fastq output files.
);

# Display a usage message if the following parameters are not specified.
usage() unless (
      defined $barcode_infile
      and defined $output_dir
);

# The restriction enzyme(s) used to digest the genomic sequences. Default: PstI/MspI
$restriction_enzymes = 'PstI/MspI' unless defined $restriction_enzymes;

# The fastq quality score encoding used in the Illumina sequencing run.  Use phred33 for Illumina 1.8+ and Sanger or phred64 for Illumina 1.3 to 1.5. Default: phred33
$encoded_phred_offset = 'phred33' unless defined $encoded_phred_offset;

# The size of the sliding window as a fraction of the read length between 0 and 1. Default: 0.15
$sliding_window_size = 0.15 unless defined $sliding_window_size;

# The quality score limit. If the average score within the sliding window drops below this value, the read is discarded. Default: 20
$quality_score_limit = 20 unless defined $quality_score_limit;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

# The barcode option for whether or not single-end or paired-end barcodes are within the FASTQ header or inline with sequence. Default: inline_null
$barcode_option = 'inline_null' unless defined $barcode_option;

# The number of mismatches allowed within the barcode sequences. Default: 0
$num_mismatches = 0 unless defined $num_mismatches;

# Program dependencies - process_radtags program from the Stacks Software Suite.
my $process_radtags 				= '/usr/local/bin/process_radtags';

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i single_end_fastq_infile -1 first_paired_end_fastq_infile -2 second_paired_end_fastq_infile -r restriction_enzymes -e encoded_phred_offset -w sliding_window_size -s quality_score_limit -l gbs_sequence_length -b barcode_infile -c barcode_option -n num_mismatches -o output_dir
    
DESCRIPTION - This program performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.

OPTIONS:

-i single_end_fastq_infile - The absolute path to the bulk fastq input file to split sequences based on barcode sequences if processing single-end seqeunces. Can either be *.fastq or *.fastq.gz extension.

-1 first_paired_end_fastq_infile - The absolute path to the first bulk fastq input file in a set of paired-end sequences to split sequences based on barcode sequences. Can either be *.fastq or *.fastq.gz extension.

-2 second_paired_end_fastq_infile - The absolute path to the second bulk fastq input file in a set of paired-end sequences to split sequences based on barcode sequences. Can either be *.fastq or *.fastq.gz extension.

-r restriction_enzymes - The restriction enzyme(s) used to digest the genomic sequences. Default: PstI/MspI

	Can be one of the following;

	PstI/MspI
	ApeKI

-e encoded_phred_offset - The fastq quality score encoding used in the Illumina sequencing run.  Use phred33 for Illumina 1.8+ and Sanger or phred64 for Illumina 1.3 to 1.5. Default: phred33

-w sliding_window_size - The size of the sliding window as a fraction of the read length between 0 and 1. Default: 0.15

-s quality_score_limit - The quality score limit. If the average score within the sliding window drops below this value, the read is discarded. Default: 20

-t gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92

-b barcode_infile - The absolute path to the barcodes input file used to split sequences into individual fastq output files.

-c barcode_option - The barcode option for whether or not single-end or paired-end barcodes are within the FASTQ header or inline with sequence. Default: inline_null

	Can be one of the following;

	inline_null: barcode is inline with sequence, occurs only on single-end read (default).
	index_null: barcode is provded in FASTQ header, occurs only on single-end read.
	inline_inline: barcode is inline with sequence, occurs on single and paired-end read.
	index_index: barcode is provded in FASTQ header, occurs on single and paired-end read.
	inline_index: barcode is inline with sequence on single-end read, occurs in FASTQ header for paired-end read.
	index_inline: barcode occurs in FASTQ header for single-end read, is inline with sequence on paired-end read.

-n num_mismatches - The number of mismatches allowed within the barcode sequences. Default: 0

-o output_dir - The absolute path to the output directory to contain the split *.fastq output files.

USAGE
}

# Choose the restriction enzyme option based on the restriction enzyme(s) used.
my $renzyme_option = "";
if($restriction_enzymes eq 'ApeKI'){ # ApeKI
	$renzyme_option = 'apeKI';
}elsif($restriction_enzymes eq 'PstI/MspI'){ # PstI/MspI
	$renzyme_option = 'pstI';
}
else{ 
	die "Input $restriction_enzymes is not one of the recognized restriction enzyme(s)! Please specify either ApeKI or PstI/MspI on the command line."
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Parse the barcodes input file based on barcode sequence so that we can split the bulk fastq file using the process_radtags program variants.
my %project_leader_names = ();
open(INFILE, "<$barcode_infile") or die "Couldn't open file $barcode_infile for reading, $!";
my $i = 0;
while(<INFILE>){
	chomp $_;
	warn $_ . "\n";
	if($i ne 0){
		my @split_row_entry = split(/\t/, $_);
		
		my ($fastq_flowcell_name, $fastq_plate_num, $fastq_well_row, $fastq_well_column, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) = @split_row_entry;
		$project_leader_names{$fastq_run_id} = $fastq_project_leader;
	}
	$i++;
}
close(INFILE) or die "Couldn't close file $barcode_infile";

# Create project leader output directory if it doesn't already exist.
my $mismatch_dir = join("_", $num_mismatches, "MISMATCH");
$mismatch_dir = "NO_MISMATCHES" if($num_mismatches eq 0);
my $project_leader_dir = join('/', $output_dir, join("_", "PROJECT_LEADER_DIR", $mismatch_dir));
unless(-d $project_leader_dir){
	mkdir($project_leader_dir, 0777) or die "Can't make directory: $!";
}

my ($split_fastq_files, $split_fastq_file_counter);
if(defined($single_end_fastq_infile)){

	# Execute the process_radtags_single_end program for quality assessment and quality control filtering and demultiplexing the files based on barcode and project leader name.
	($split_fastq_files, $split_fastq_file_counter) = process_radtags_single_end($single_end_fastq_infile, $renzyme_option, $sliding_window_size, $quality_score_limit, $gbs_sequence_length, $barcode_infile, $barcode_option, $num_mismatches, $output_dir);

	# Iterate through each split fastq sample file from the process_radtags program and rename the files based on project leader.
	foreach my $file_name (sort keys %{$split_fastq_files}){
		warn "Processing " . $file_name . ".....\n";
		my $split_fastq_infile = $split_fastq_files->{$file_name};
		
		# Get the basename of the fastq filename without the .fastq extension.
		my $fastq_filename = fileparse($split_fastq_infile, qr/\.fq/);
		
		my $fastq_project_leader = $project_leader_names{$fastq_filename};
		
		# Create project output directory if it doesn't already exist.
		my $project_output_dir = join('/', $project_leader_dir, $fastq_project_leader);
		unless(-d $project_output_dir){
			mkdir($project_output_dir, 0777) or die "Can't make directory: $!";
		}
		
		my $fastq_project_leader_outfile = join('/', $project_output_dir, join("", $fastq_filename, ".fastq"));
		warn "Copying $file_name to $fastq_project_leader_outfile.....";
		copy($split_fastq_infile, $fastq_project_leader_outfile) or die "Copy failed: $!";
		unlink($split_fastq_infile) or die "Could not unlink $split_fastq_infile: $!";
	}
	
}elsif(defined($first_paired_end_fastq_infile) and defined($second_paired_end_fastq_infile)){

	# Execute the process_radtags_paired_end program for quality assessment and quality control filtering and demultiplexing the files based on barcode and project leader name.
	($split_fastq_files, $split_fastq_file_counter) = process_radtags_paired_end($first_paired_end_fastq_infile, $second_paired_end_fastq_infile, $renzyme_option, $sliding_window_size, $quality_score_limit, $gbs_sequence_length, $barcode_infile, $barcode_option, $num_mismatches, $output_dir);

	# Iterate through each split fastq sample file from the process_radtags program and rename the files based on project leader.
	foreach my $file_name (sort keys %{$split_fastq_files}){
		warn "Processing " . $file_name . ".....\n";
		my $split_fastq_infile = $split_fastq_files->{$file_name};
		
		# Get the basename of the fastq filename without the .fastq extension.
		my ($fastq_filename, $fastq_dirname, $fastq_suffix) = fileparse($split_fastq_infile, qr/\.\d+\.fq|\.rem\.\d+\.fq/);
		$fastq_suffix =~ s/\.fq//g;
		
		my $fastq_project_leader = $project_leader_names{$fastq_filename};
		
		# Create project output directory if it doesn't already exist.
		my $project_output_dir = join('/', $project_leader_dir, $fastq_project_leader);
		unless(-d $project_output_dir){
			mkdir($project_output_dir, 0777) or die "Can't make directory: $!";
		}
		
		my $fastq_project_leader_outfile = join('/', $project_output_dir, join("", $fastq_filename, "$fastq_suffix.fastq"));
		warn "Copying $file_name to $fastq_project_leader_outfile.....";
		copy($split_fastq_infile, $fastq_project_leader_outfile) or die "Copy failed: $!";
		unlink($split_fastq_infile) or die "Could not unlink $split_fastq_infile: $!";
	}
	
}else{
	die usage();
}

# ($split_fastq_files, $split_fastq_file_counter) = process_radtags_single_end($single_end_fastq_infile, $renzyme_option, $sliding_window_size, $quality_score_limit, $gbs_sequence_length, $barcode_infile, $barcode_option, $num_mismatches, $output_dir) - Executes the process_radtags program from the Stacks Software Suite using a single-end fastq input file for quality assessment and quality control filtering and demultiplexing the files based on barcode and project leader name.
#
# Input paramater(s):
# 
# $single_end_fastq_infile - The absolute path to the single-end fastq input file to split sequences based on barcode sequences.
#
# $renzyme_option - The restriction enzyme used (cut site occurs on single-end read).
# 
# $sliding_window_size - The size of the sliding window as a fraction of the read length.
#
# $quality_score_limit - The quality score limit. If the average score within the sliding window drops below this value, the read is discarded.
#
# $gbs_sequence_length - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
# 
# $barcode_infile - The absolute path to the barcodes input file used to split sequences into individual fastq output files.
# 
# $barcode_option - The barcode option for whether or not single-end barcodes are within the FASTQ header or inline with sequence.
#
# $num_mismatches - The number of mismatches allowed within the barcode sequences.
#
# $output_dir - The absolute path to the output directory to contain the split *.fq output files.
#
# Output paramater(s):
# 
# $split_fastq_files - A hash reference containing all the files with file extension *.fq in key/value pairs.
# 
# key => filename ( e.g. filename.suffix )
# value => absolue filepath ( e.g. /path/to/filename.suffix )
# 
# $split_fastq_file_counter - The number of split fastq files stored with file extension *.fq.
sub process_radtags_single_end{

	# The GBS single-end fastq input file to filter based on quality and demultiplex based on barcodes.
	my $single_end_fastq_infile = shift;
	die "Error lost fastq input file" unless defined $single_end_fastq_infile;
	
	# The restriction enzyme used (cut site occurs on single-end read).
	my $renzyme_option = shift;
	die "Error lost restriction enzyme option" unless defined $renzyme_option;
	
	# The size of the sliding window as a fraction of the read length.
	my $sliding_window_size = shift;
	die "Error lost size of the sliding window" unless defined $sliding_window_size;
	
	# The quality score limit. If the average score within the sliding window drops below this value, the read is discarded.
	my $quality_score_limit = shift;
	die "Error lost the quality score limit" unless defined $quality_score_limit;
	
	# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences.
	my $gbs_sequence_length = shift;
	die "Error lost the GBS fastq sequence length in base pairs (bps)" unless defined $gbs_sequence_length;
	
	# The absolute path to the barcodes input file used to split sequences into individual fastq output files.
	my $barcode_infile = shift;
	die "Error lost barcode sequence file" unless defined $barcode_infile;
	
	# The barcode option for whether or not single-end barcodes are within the FASTQ header or inline with sequence.
	my $barcode_option = shift;
	die "Error lost barcode option" unless defined $barcode_option;
	
	# The number of mismatches allowed within the barcode sequences.
	my $num_mismatches = shift;
	die "Error lost number of mismatches" unless defined $num_mismatches;
	
	# The absolute path to the output directory to contain the split *.fastq output files.
	my $output_dir = shift;
	die "Error lost output directory" unless defined $output_dir;
	
	# Get the flowcell name so that we can further demultiplex the fastq files.
	my $flowcell_name = fileparse($barcode_infile, qr/-barcodes\.\w+/);
	
	# Create output directory if it doesn't already exist.
	my $split_fastq_output_dir = join('/', $output_dir, join("_", "SPLIT", $flowcell_name, "FASTQ_FILES"));
	unless(-d $split_fastq_output_dir){
		mkdir($split_fastq_output_dir, 0777) or die "Can't make directory: $!";
	}
	
	# Check the $split_fastq_output_dir for files with extension *.fq and check to make sure the size of the resulting files are non-zero 
	# indicating that the $single_end_fastq_infile file hasn't been already filtered based on quality and demultiplexed into separate files based on barcodes.  
	my ($split_fastq_files, $split_fastq_file_counter) = find_files($split_fastq_output_dir, "fq");
	my $non_zero_split_fastq_files = 0;
	foreach my $file_name (sort keys %{$split_fastq_files}){
		warn $file_name . "\n";
		if(-s $split_fastq_files->{$file_name}){
			$non_zero_split_fastq_files++;
		}
	}
	
	# Check if the $single_end_fastq_infile file has been already filtered based on quality and demultiplexed into separate files based on barcodes.
	# If it hasn't, process the $single_end_fastq_infile file using the following process_radtags commands.
	unless(($non_zero_split_fastq_files eq $split_fastq_file_counter) and ($split_fastq_file_counter ne 0)){
	
		# Parse the barcodes input file based on barcode sequence and barcode sequence length so that we can split the bulk fastq file using the longest to shortest barcode so that we can use as input to the process_radtags program.
		open(INFILE, "<$barcode_infile") or die "Couldn't open file $barcode_infile for reading, $!";
		my @barcode_entries = ();
		my $i = 0;
		while(<INFILE>){
			chomp $_;
	# 		warn $_ . "\n";
			if($i ne 0){
				my @split_row_entry = split(/\t/, $_);
				my ($fastq_flowcell_name, $fastq_plate_num, $fastq_well_row, $fastq_well_column, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) = @split_row_entry;
				
				my $barcode_seq_length = length($fastq_barcode_seq);
				my $barcodes_entry = join("\t", $barcode_seq_length, $fastq_barcode_seq, $fastq_run_id);
	# 			warn $barcodes_entry . "\n";

				push(@barcode_entries, [split(/\t/, $barcodes_entry)]);
				
			}
			$i++;
		}
		close(INFILE) or die "Couldn't close file $barcode_infile";
		
		# Generate the barcodes input file for the process_radtags program.
		my $restriction_enzyme_names = $restriction_enzymes;
		$restriction_enzyme_names =~ s/\//_/g;
		my $barcode_outfile = join('/', $output_dir, join("_", $flowcell_name, $restriction_enzyme_names, "single-end", "process-radtags-barcodes.txt"));
		open(OUTFILE, ">$barcode_outfile") or die "Couldn't open file $barcode_outfile for writting, $!";
		foreach my $barcodes_entry (sort {$b->[0] <=> $a->[0] || $a->[1] cmp $b->[1]} @barcode_entries){
			my ($fastq_barcode_seq, $fastq_run_id) = (@$barcodes_entry[1], @$barcodes_entry[2]);
			print OUTFILE join("\t", $fastq_barcode_seq, $fastq_run_id) . "\n";
		}
		close(OUTFILE) or die "Couldn't close file $barcode_outfile";

		# Choose the input file type based on whether or not the bulk fastq file is compressed.
		my $infile_type = "";
		if($single_end_fastq_infile =~ m/\.gz$/){ # If the bulk fastq file is compressed with extension *.gz.
			$infile_type = 'gzfastq';
		}else{ # If the bulk fastq file is not compressed with extension *.fastq.
			$infile_type = 'fastq';
		}
		
		# Execute the process_radtags command using system.
		my $process_radtagsCmd = "$process_radtags -f $single_end_fastq_infile -i $infile_type -b $barcode_outfile -o $split_fastq_output_dir -y fastq -c -q -r -E $encoded_phred_offset -D -w $sliding_window_size -s $quality_score_limit -t $gbs_sequence_length --$barcode_option -e $renzyme_option --filter_illumina --barcode_dist $num_mismatches";
		warn $process_radtagsCmd . "\n\n";
		my $status = system($process_radtagsCmd) == 0 or die "Error calling $process_radtags: $?";
	}
	
	# Get all the split files with extension *.fq so that we can return them for further processing.
	($split_fastq_files, $split_fastq_file_counter) = find_files($split_fastq_output_dir, "fq");
	return ($split_fastq_files, $split_fastq_file_counter);
}

# ($split_fastq_files, $split_fastq_file_counter) = process_radtags_paired_end($first_paired_end_fastq_infile, $second_paired_end_fastq_infile, $renzyme_option, $sliding_window_size, $quality_score_limit, $gbs_sequence_length, $barcode_infile, $barcode_option, $num_mismatches, $output_dir) - Executes the process_radtags program from the Stacks Software Suite using paired-end fastq input files for quality assessment and quality control filtering and demultiplexing the files based on barcode and project leader name.
#
# Input paramater(s):
# 
# $first_paired_end_fastq_infile - The absolute path to the first GBS paired-end fastq fastq input file to split sequences based on barcode sequences.
#
# $second_paired_end_fastq_infile - The absolute path to the second GBS paired-end fastq fastq input file to split sequences based on barcode sequences.
#
# $renzyme_option - The restriction enzyme used (cut site occurs on single-end read).
# 
# $sliding_window_size - The size of the sliding window as a fraction of the read length.
#
# $quality_score_limit - The quality score limit. If the average score within the sliding window drops below this value, the read is discarded.
#
# $gbs_sequence_length - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
# 
# $barcode_infile - The absolute path to the barcodes input file used to split sequences into individual fastq output files.
# 
# $barcode_option - The barcode option for whether or not paired-end barcodes are within the FASTQ header or inline with sequence.
# 
# $num_mismatches - The number of mismatches allowed within the barcode sequences.
#
# $output_dir - The absolute path to the output directory to contain the split *.fq output files.
#
# Output paramater(s):
# 
# $split_fastq_files - A hash reference containing all the files with file extension *.fq in key/value pairs.
# 
# key => filename ( e.g. filename.suffix )
# value => absolue filepath ( e.g. /path/to/filename.suffix )
# 
# $split_fastq_file_counter - The number of split fastq files stored with file extension *.fq.
sub process_radtags_paired_end{

	# The first GBS paired-end fastq input file to filter based on quality and demultiplex based on barcodes.
	my $first_paired_end_fastq_infile = shift;
	die "Error lost first paired-end fastq input file" unless defined $first_paired_end_fastq_infile;
	
	# The second GBS paired-end fastq input file to filter based on quality and demultiplex based on barcodes.
	my $second_paired_end_fastq_infile = shift;
	die "Error lost second paired-end fastq input file" unless defined $second_paired_end_fastq_infile;
	
	# The restriction enzyme used (cut site occurs on single-end read).
	my $renzyme_option = shift;
	die "Error lost restriction enzyme option" unless defined $renzyme_option;
	
	# The size of the sliding window as a fraction of the read length.
	my $sliding_window_size = shift;
	die "Error lost size of the sliding window" unless defined $sliding_window_size;
	
	# The quality score limit. If the average score within the sliding window drops below this value, the read is discarded.
	my $quality_score_limit = shift;
	die "Error lost the quality score limit" unless defined $quality_score_limit;
	
	# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences.
	my $gbs_sequence_length = shift;
	die "Error lost the GBS fastq sequence length in base pairs (bps)" unless defined $gbs_sequence_length;
	
	# The absolute path to the barcodes input file used to split sequences into individual fastq output files.
	my $barcode_infile = shift;
	die "Error lost barcode sequence file" unless defined $barcode_infile;
	
	# The barcode option for whether or not paired-end barcodes are within the FASTQ header or inline with sequence.
	my $barcode_option = shift;
	die "Error lost barcode option" unless defined $barcode_option;
	
	# The number of mismatches allowed within the barcode sequences.
	my $num_mismatches = shift;
	die "Error lost number of mismatches" unless defined $num_mismatches;
	
	# The absolute path to the output directory to contain the split *.fastq output files.
	my $output_dir = shift;
	die "Error lost output directory" unless defined $output_dir;

	# Get the flowcell name so that we can further demultiplex the fastq files.
	my $flowcell_name = fileparse($barcode_infile, qr/-barcodes\.\w+/);
	
	# Create output directory if it doesn't already exist.
	my $split_fastq_output_dir = join('/', $output_dir, join("_", "SPLIT", $flowcell_name, "FASTQ_FILES"));
	unless(-d $split_fastq_output_dir){
		mkdir($split_fastq_output_dir, 0777) or die "Can't make directory: $!";
	}
	
	# Check the $split_fastq_output_dir for files with extension *.fq and check to make sure the size of the resulting files are non-zero 
	# indicating that the $single_end_fastq_infile file hasn't been already filtered based on quality and demultiplexed into separate files based on barcodes.  
	my ($split_fastq_files, $split_fastq_file_counter) = find_files($split_fastq_output_dir, "fq");
	my $non_zero_split_fastq_files = 0;
	foreach my $file_name (sort keys %{$split_fastq_files}){
		warn $file_name . "\n";
		if(-s $split_fastq_files->{$file_name}){
			$non_zero_split_fastq_files++;
		}
	}
	
	# Check if the $first_paired_end_fastq_infile and $second_paired_end_fastq_infile files have already been filtered based on quality and demultiplexed into separate files based on barcodes.
	# If it hasn't, process the $first_paired_end_fastq_infile and $second_paired_end_fastq_infile files using the following process_radtags commands.
	unless(($non_zero_split_fastq_files eq $split_fastq_file_counter) and ($split_fastq_file_counter ne 0)){
	
		# Parse the barcodes input file based on barcode sequence and barcode sequence length so that we can split the bulk fastq file using the longest to shortest barcode so that we can use as input to the process_radtags program.
		open(INFILE, "<$barcode_infile") or die "Couldn't open file $barcode_infile for reading, $!";
		my @barcode_entries = ();
		my $i = 0;
		while(<INFILE>){
			chomp $_;
	# 		warn $_ . "\n";
			if($i ne 0){
				my @split_row_entry = split(/\t/, $_);
				my ($fastq_flowcell_name, $fastq_plate_num, $fastq_well_row, $fastq_well_column, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) = @split_row_entry;

				my $barcode_seq_length = length($fastq_barcode_seq);
				my $barcodes_entry = join("\t", $barcode_seq_length, $fastq_barcode_seq, $fastq_run_id);
	# 			warn $barcodes_entry . "\n";

				push(@barcode_entries, [split(/\t/, $barcodes_entry)]);
				
			}
			$i++;
		}
		close(INFILE) or die "Couldn't close file $barcode_infile";
		
		# Generate the barcodes input file for the process_radtags program.
		my $restriction_enzyme_names = $restriction_enzymes;
		$restriction_enzyme_names =~ s/\//-/g;
        
		my $barcode_outfile = join('/', $output_dir, join("_", $flowcell_name, $restriction_enzyme_names, "paired-end", "process-radtags-barcodes.txt"));
		open(OUTFILE, ">$barcode_outfile") or die "Couldn't open file $barcode_outfile for writting, $!";
		foreach my $barcodes_entry (sort {$b->[0] <=> $a->[0] || $a->[1] cmp $b->[1]} @barcode_entries){
			my ($fastq_barcode_seq, $fastq_run_id) = (@$barcodes_entry[1], @$barcodes_entry[2]);
			print OUTFILE join("\t", $fastq_barcode_seq, $fastq_run_id) . "\n";
		}
		close(OUTFILE) or die "Couldn't close file $barcode_outfile";

		# Choose the input file type based on whether or not the bulk fastq file is compressed.
		my $infile_type = "";
		if(($first_paired_end_fastq_infile =~ m/\.gz$/) and ($second_paired_end_fastq_infile =~ m/\.gz$/)){ # If the bulk fastq file is compressed with extension *.gz.
			$infile_type = 'gzfastq';
		}elsif(($first_paired_end_fastq_infile =~ m/\.fastq$/) and ($second_paired_end_fastq_infile =~ m/\.fastq$/)){ # If the bulk fastq file is not compressed with extension *.fastq.
			$infile_type = 'fastq';
		}
		
		# Execute the process_radtags command using system.
		my $process_radtagsCmd = "$process_radtags -1 $first_paired_end_fastq_infile -2 $second_paired_end_fastq_infile -i $infile_type -b $barcode_outfile -o $split_fastq_output_dir -y fastq -c -q -r -E $encoded_phred_offset -D -w $sliding_window_size -s $quality_score_limit -t $gbs_sequence_length --$barcode_option -e $renzyme_option --filter_illumina --barcode_dist $num_mismatches";
		warn $process_radtagsCmd . "\n\n";
		my $status = system($process_radtagsCmd) == 0 or die "Error calling $process_radtags: $?";
	}
	
	# Get all the split files with extension *.fq so that we can return them for further processing.
	($split_fastq_files, $split_fastq_file_counter) = find_files($split_fastq_output_dir, "fq");
	return ($split_fastq_files, $split_fastq_file_counter);
}

# (\%files, $file_counter) = find_files($infile_dir) - Find all files in the specified input file directory with the file extension *.suffix.
# 
# Input paramater(s):
# 
# $infile_dir - The iinput file directory.
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