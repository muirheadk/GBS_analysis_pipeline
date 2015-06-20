#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;


#### PROGRAM NAME ####
# execute_multi_fastq_barcode_splitter.pl - A program that performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.

#### DESCRIPTION ####
# This program performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.

#### SAMPLE COMMAND ####
# PstI/MspI single-end reads command using phred33 encoding
my ($GBS_barcode_sample_index_file, $restriction_enzymes, $encoded_phred_offset, $sliding_window_size, $quality_score_limit, $gbs_sequence_length, $barcode_infile, $barcode_option, $num_mismatches,$num_threads, $output_dir);
GetOptions(
	'i=s'    => \$GBS_barcode_sample_index_file, # The GBS barcode sample index input file.
    'r=s'    => \$restriction_enzymes, # The restriction enzyme(s) used to digest the genomic sequences. Can be ApeKI, PstI/MspI, or SbfI/MspI. Default: PstI/MspI
    'e=s'    => \$encoded_phred_offset, # The fastq quality score encoding used in the Illumina sequencing run.  Use phred33 for Illumina 1.8+ and Sanger or phred64 for Illumina 1.3 to 1.5. Default: phred33
    'w=s'    => \$sliding_window_size, # The size of the sliding window as a fraction of the read length between 0 and 1. Default: 0.15
    's=s'    => \$quality_score_limit, # The quality score limit. If the average score within the sliding window drops below this value, the read is discarded. Default: 20
    't=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
    'b=s'    => \$barcode_infile, # The absolute path to the barcodes input file used to split sequences into individual fastq output files.
    'c=s'    => \$barcode_option, # The barcode option for whether or not single-end or paired-end barcodes are within the FASTQ header or inline with sequence. Default: inline_null
    'n=s'    => \$num_mismatches, # The number of mismatches allowed within the barcode sequences. Default: 0
    'h=s'    => \$num_threads, # The number of cpu threads to use for the stacks programs. Default: 1
    'o=s'    => \$output_dir, # The absolute path to the output directory to contain the split *.fastq output files.
);

# Display a usage message if the following parameters are not specified.
usage() unless (
	defined $GBS_barcode_sample_index_file
	and defined $output_dir
);

# The restriction enzyme(s) used to digest the genomic sequences. Can be ApeKI, PstI/MspI, or SbfI/MspI. Default: PstI/MspI
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

# The number of cpu threads to use for the stacks programs. Default: 1
$num_threads = 1 unless defined $num_threads;

# Program dependencies - fastq_quality_barcode_splitter.pl script that runs the process_radtags program from the Stacks Software Suite.
my $dirname = dirname($0);
my $fastq_quality_barcode_splitter 				= join('/', $dirname, 'fastq_quality_barcode_splitter.pl');

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i GBS_barcode_sample_index_file -r restriction_enzymes -e encoded_phred_offset -w sliding_window_size -s quality_score_limit -l gbs_sequence_length -b barcode_infile -c barcode_option -n num_mismatches -c num_threads -o output_dir

    
DESCRIPTION - This program performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.
        
-i GBS_barcode_sample_index_file - The GBS barcode sample index input file.
    
    
-r restriction_enzymes - The restriction enzyme(s) used to digest the genomic sequences. Can be ApeKI, PstI/MspI, or SbfI/MspI. Default: PstI/MspI

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

-h num_threads - The number of cpu threads to use for the stacks programs. Default: 1
        
-o output_dir - The absolute path to the output directory to contain the split *.fastq output files.

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}


# Parse the barcodes input file based on barcode sequence so that we can split the bulk fastq file using the process_radtags program variants.
my %GBS_fastq_index = ();
open(INFILE, "<$GBS_barcode_sample_index_file") or die "Couldn't open file $GBS_barcode_sample_index_file for reading, $!";

while(<INFILE>){
	chomp $_;
	warn $_ . "\n";
		my @split_row_entry = split(/\t/, $_);
		
		my ($flowcell_name, $raw_fastq_infile, $barcode_infile) = @split_row_entry;
		$GBS_fastq_index{$flowcell_name} = join("\t", $raw_fastq_infile, $barcode_infile);
}
close(INFILE) or die "Couldn't close file $GBS_barcode_sample_index_file";

# Create the process radtags directory if it doesn't already exist.
my $process_radtags_output_dir = join('/', $output_dir, "PROCESSED_RADTAGS");
unless(-d $process_radtags_output_dir){
    mkdir($process_radtags_output_dir, 0777) or die "Can't make directory: $!";
}

# Run fastq_quality_barcode_splitter.pl script in parallel if Parallel loops is available and if num_threads is defined. Otherwise,
if ((require Parallel::Loops)and($num_threads)){
    
    # Run these jobs in parallel.
    my $parallel = Parallel::Loops->new($num_threads);
	my @jobs = keys %GBS_fastq_index;
    $parallel->foreach(\@jobs,sub {
        
        # Get the flowcell name so that we can further demultiplex the fastq files.
		my $flowcell_name = $_;
        
        # Create output directory if it doesn't already exist.
        my $split_fastq_output_dir = join('/', $process_radtags_output_dir, join("_", "SPLIT", $flowcell_name, "FASTQ_FILES"));
        unless(-d $split_fastq_output_dir){
            mkdir($split_fastq_output_dir, 0777) or die "Can't make directory: $!";
        }
        
        # Split the fastq barcode index entry to obtain the fastq file and barcode files as input.
        my @split_index_entry = split(/\t/, $GBS_fastq_index{$flowcell_name});
        my $fastq_infile = $split_index_entry[0];
        my $barcode_infile = $split_index_entry[1];
        
        # Create a process_radtags log file of everything verbose.
        my $process_radtags_fastq_logfile = join('/', $split_fastq_output_dir, join("-", $flowcell_name, "process_radtags_log.txt"));
        my $fastq_quality_barcode_splitterCmd = "$fastq_quality_barcode_splitter -i $fastq_infile -b $barcode_infile -r PstI/MspI -e $encoded_phred_offset -w $sliding_window_size -s $quality_score_limit -t $gbs_sequence_length -c -n $num_mismatches -c inline_null -o $process_radtags_output_dir 2> $process_radtags_fastq_logfile";
        warn $fastq_quality_barcode_splitterCmd . "\n\n";
        my $status = system($fastq_quality_barcode_splitterCmd) == 0 or die "Error calling $fastq_quality_barcode_splitter: $?";
    });
    undef $parallel;
}else{
    # Iterate through each flowcell name.
    foreach my $flowcell_name (sort keys %GBS_fastq_index){

        # Create output directory if it doesn't already exist.
        my $split_fastq_output_dir = join('/', $process_radtags_output_dir, join("_", "SPLIT", $flowcell_name, "FASTQ_FILES"));
        unless(-d $split_fastq_output_dir){
            mkdir($split_fastq_output_dir, 0777) or die "Can't make directory: $!";
        }
        
        # Split the fastq barcode index entry to obtain the fastq file and barcode files as input.
        my @split_index_entry = split(/\t/, $GBS_fastq_index{$flowcell_name});
        my $fastq_infile = $split_index_entry[0];
        my $barcode_infile = $split_index_entry[1];
        
        # Create a process_radtags log file of everything verbose.
        my $process_radtags_fastq_logfile = join('/', $split_fastq_output_dir, join("-", $flowcell_name, "process_radtags.txt"));
        my $fastq_quality_barcode_splitterCmd = "$fastq_quality_barcode_splitter -i $fastq_infile -b $barcode_infile -r PstI/MspI -e $encoded_phred_offset -w $sliding_window_size -s $quality_score_limit -t $gbs_sequence_length -c -n $num_mismatches -c inline_null -o $process_radtags_output_dir 2> $process_radtags_fastq_logfile";
        warn $fastq_quality_barcode_splitterCmd . "\n\n";
        my $status = system($fastq_quality_barcode_splitterCmd) == 0 or die "Error calling $fastq_quality_barcode_splitter: $?";

    }
}
