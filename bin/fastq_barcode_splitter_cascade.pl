#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Copy;
use File::Basename;

# NEED TO ADDED COPY COMMANDS AFTER EACH ITERITION IN ORDER TO SAVE SPACE BECAUSE THE FILES GENERATED ARE TECHNICALLY DUPLICATED.

# perl fastq_barcode_splitter_cascade.pl -i ~/workspace/GBS_data-08-10-2013/HI.1405.008.GQ03122013-5_R1.fastq -b ~/workspace/GBS_data-08-10-2013/GBS_barcodes-2013-10-09.csv -n 1 -o ~/workspace/GBS_data-08-10-2013
my ($fastq_infile, $barcode_infile, $num_mismatches, $output_dir);
GetOptions(
      'i=s'    => \$fastq_infile, # The absolute path to the bulk fastq input file to split sequences based on barcode sequences. Can either be *.fastq or *.fastq.gz extension.
      'b=s'    => \$barcode_infile, # The absolute path to the barcodes input file used to split sequences into individual fastq output files.
      'n=s'    => \$num_mismatches, # The number of mismatches allowed within the barcode sequences.
      'o=s'    => \$output_dir, # The absolute path to the output directory to contain the split *.fastq output files.
);

# Display a usage message if the following parameters are not specified.
usage() unless (
      defined $fastq_infile
      and defined $barcode_infile
      and defined $output_dir
);

# The number of mismatches allowed within the barcode sequences. Default: 0
$num_mismatches = 0 unless defined $num_mismatches;

# Program dependencies - The absolute paths to gunzip to uncompress bulk compressed fastq.gz input file if present and the fastx barcode splitter from the fastx toolkit.
my ($gunzip, $fastx_barcode_splitter);
$gunzip				= '/bin/gunzip';
$fastx_barcode_splitter 	= '/usr/local/bin/fastx_barcode_splitter.pl';

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fastq_infile -b barcode_infile -n num_mismatches -o output_dir
    
DESCRIPTION - Split a bulk fastq input file into separate fastq files based on the longest to shortest barcode sequences for 
each individual specified by the barcodes input file and output the results based on project leader specified in the barcodes 
input file.
    
OPTIONS:

-i fastq_infile - The absolute path to the bulk fastq input file to split sequences based on barcode sequences. Can either be *.fastq or *.fastq.gz extension.

	e.g. /path/to/HI.1405.008.GQ03122013-5_R1.fastq    <---- fastq file format
	     /path/to/HI.1405.008.GQ03122013-5_R1.fastq.gz <---- compressed fastq file format
    
-b barcode_infile - The absolute path to the barcodes input file used to split sequences into individual fastq output files.

	e.g. /path/to/GBS_barcodes-2013-10-09.csv

-n num_mismatches - The number of mismatches allowed within the barcode sequences. Default: 0
      
-o output_dir - The absolute path to the output directory to contain the split *.fastq output files.

	e.g. /path/to/output_dir

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# If the bulk fastq file is compressed, uncompress the file and set the resulting fastq filename to be the fastq infile.
if($fastq_infile =~ m/\.gz$/){
	my $uncompressed_fastq_file = gunzip_fastq_file($fastq_infile);
	$fastq_infile = $uncompressed_fastq_file;
}

# Create output directory if it doesn't already exist.
my $split_fastq_output_dir = join('/', $output_dir, "SPLIT_FASTQ_OUTFILES");
unless(-d $split_fastq_output_dir){
      mkdir($split_fastq_output_dir, 0777) or die "Can't make directory: $!";
}

# Parse the barcodes input file based on barcode length so that we can split the bulk fastq file using the longest to shortest barcode.
my %barcode_data = ();
open(INFILE, "<$barcode_infile") or die "Couldn't open file $barcode_infile for reading, $!";
my $i = 0;
while(<INFILE>){
	chomp $_;
	warn $_ . "\n";
	if($i ne 0){
		my @split_row_entry = split(/\t/, $_);
		my ($fastq_plate_num, $fastq_well_row, $fastq_well_column, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) = ($split_row_entry[0], $split_row_entry[1], $split_row_entry[2], $split_row_entry[3], $split_row_entry[4], $split_row_entry[5]);

		my $barcode_seq_length = length($fastq_barcode_seq);
		my $fastq_barcode_name = join("_", $fastq_run_id, $fastq_barcode_seq, $fastq_plate_num, join("", $fastq_well_row, $fastq_well_column));
		push(@{$barcode_data{$barcode_seq_length}},  join("\t", $fastq_barcode_name, $fastq_barcode_seq, $fastq_project_leader));
		
	}
	$i++;
}
close(INFILE) or die "Couldn't close file $barcode_infile";

# Start from the longest to shortest barcode length and keep track of which fastq file belongs to which project leader.
my %barcode_filenames = ();
my $unmatched_file_counter = 0;
my $unmatched_fastq_infile = "";
foreach my $barcode_seq_length (sort {$b <=> $a} keys %barcode_data){
	
	my $barcodes_lengths_dir = join('/', $split_fastq_output_dir, "FASTQ_BARCODE_LENGTH_DIR");
	unless(-d $barcodes_lengths_dir){
		mkdir($barcodes_lengths_dir, 0777) or die "Can't make directory: $!";
	}
	my $barcodes_output_dir = join('/', $barcodes_lengths_dir, join("_", "FASTQ_BARCODES_LENGTH", $barcode_seq_length));
	unless(-d $barcodes_output_dir){
		mkdir($barcodes_output_dir, 0777) or die "Can't make directory: $!";
	}

	my $barcode_splitter_file = join('/', $split_fastq_output_dir, "fastq_barcodes_length_$barcode_seq_length" . ".txt");
	open(OUTFILE, ">$barcode_splitter_file") or die "Couldn't open file $barcode_splitter_file for writting, $!";
	foreach my $fastq_barcode_data (@{$barcode_data{$barcode_seq_length}}){
		my @split_fastq_barcode_data = split(/\t/, $fastq_barcode_data);
		my ($fastq_barcode_name, $fastq_barcode_seq, $fastq_project_leader) = ($split_fastq_barcode_data[0], $split_fastq_barcode_data[1], $split_fastq_barcode_data[2]);
		print OUTFILE join("\t", $fastq_barcode_name, $fastq_barcode_seq) . "\n";
		my $fastq_barcode_length_filename = join('/', $barcodes_output_dir, join("", $fastq_barcode_name, ".fastq"));
		
		my $project_leader_dir = join('/', $split_fastq_output_dir, "PROJECT_LEADER_DIR");
		unless(-d $project_leader_dir){
			mkdir($project_leader_dir, 0777) or die "Can't make directory: $!";
		}
		
		my $project_output_dir = join('/', $project_leader_dir, $fastq_project_leader);
		unless(-d $project_output_dir){
			mkdir($project_output_dir, 0777) or die "Can't make directory: $!";
		}
		
		my $fastq_barcode_project_filename = join('/', $project_output_dir, join("", $fastq_barcode_name, ".fastq"));
		
		push(@{$barcode_filenames{$fastq_project_leader}}, join("\t", $fastq_barcode_length_filename, $fastq_barcode_project_filename));
	}
	close(OUTFILE) or die "Couldn't close file $barcode_splitter_file";
	
	fastx_barcode_splitter($fastq_infile, $barcode_splitter_file, $num_mismatches, $barcodes_output_dir) if($unmatched_file_counter eq 0);
	
	fastx_barcode_splitter($unmatched_fastq_infile, $barcode_splitter_file, $num_mismatches, $barcodes_output_dir) if($unmatched_file_counter >= 1);
	
	$unmatched_fastq_infile = join('/', $barcodes_output_dir, "unmatched.fastq");
	$unmatched_file_counter++;
}

# Copy the fastq output files from the FASTQ_BARCODE_LENGTH_DIR/FASTQ_BARCODES_LENGTH_$barcode_seq_length directories based on project leader names.
foreach my $fastq_project_leader (sort keys %barcode_filenames){
	
	foreach my $fastq_barcode_filenames (@{$barcode_filenames{$fastq_project_leader}}){
		
		my @split_fastq_barcode_filenames = split(/\t/, $fastq_barcode_filenames);
		my ($fastq_barcode_length_filename, $fastq_barcode_project_filename) = ($split_fastq_barcode_filenames[0],$split_fastq_barcode_filenames[1]);
            	copy($fastq_barcode_length_filename, $fastq_barcode_project_filename) or die "Copy failed: $!";
	}
}

# execute the fastx_barcode_splitter.pl script from the fastx toolkit.
sub fastx_barcode_splitter{
	
	my $fastq_infile = shift;
	die "Error fastq input file" unless defined $fastq_infile;
	my $barcode_infile = shift;
	die "Error lost barcode sequence file" unless defined $barcode_infile;
	my $num_mismatches = shift;
	die "Error lost number of mismatches" unless defined $num_mismatches;
	my $fastq_output_dir = shift;
	die "Error lost fastq output directory" unless defined $fastq_output_dir;

	
# 	$fastx_barcode_splitter --bcfile mpb_barcodes_5nucls.txt --bol --mismatches 1 --prefix ~/MPB_GBS_Data-08-10-2013/gbs_ --suffix ".fastq" < HI.1405.008.GQ03122013-5_R1.fastq
	my $fastx_barcode_splitterCmd  = "$fastx_barcode_splitter --bcfile $barcode_infile --bol --mismatches $num_mismatches --prefix $fastq_output_dir/ --suffix \".fastq\" < $fastq_infile";
	warn $fastx_barcode_splitterCmd . "\n\n";

	my $status = system($fastx_barcode_splitterCmd) == 0 or die "Error calling $fastx_barcode_splitter: $?";
}

# execute the gunzip program to uncompress the compressed fastq file.
sub gunzip_fastq_file{
	
	my $fastq_file = shift;
	die "Error lost the fastq file to compress using gunzip" unless defined $fastq_file;
	
	my ($fastq_filename, $fastq_dir) = fileparse($fastq_file, ".gz");
	
	my $uncompressed_fastq_file = join('/', $fastq_dir, $fastq_filename);
	
	warn "Calling gunzip for $fastq_file....\n";
	my $gunzipCmd  = "$gunzip -c $fastq_file > $uncompressed_fastq_file";
	warn $gunzipCmd . "\n\n";
	system($gunzipCmd) == 0 or die "Error calling $gunzipCmd: $?";
	
	return $uncompressed_fastq_file;
}