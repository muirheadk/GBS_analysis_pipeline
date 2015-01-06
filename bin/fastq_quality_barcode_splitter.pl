#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Copy;
use File::Basename;
use Switch;

#### PROGRAM NAME ####
# fastq_quality_barcode_splitter.pl - A program that performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.

#### DESCRIPTION ####
# This program performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.

#### SAMPLE COMMAND ####
# perl fastq_quality_barcode_splitter.pl -i ~/workspace/GBS_data-08-10-2013/HI.1405.008.GQ03122013-5_R1.fastq.gz -b ~/workspace/GBS_data-08-10-2013/GBS_barcodes-2013-10-09.csv -n 0 -o ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS
my ($fastq_infile, $restriction_enzymes, $barcode_infile, $num_mismatches, $output_dir);
GetOptions(
	'i=s'    => \$fastq_infile, # The absolute path to the bulk fastq input file to split sequences based on barcode sequences. Can either be *.fastq or *.fastq.gz extension.
	'r=s'    => \$restriction_enzymes, # The restriction enzyme(s) used to digest the genomic sequences. Default: Pstl/MspI
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

# The restriction enzyme(s) used to digest the genomic sequences. Default: Pstl/MspI
$restriction_enzymes = 'Pstl/MspI' unless defined $restriction_enzymes;

# The number of mismatches allowed within the barcode sequences. Default: 0
$num_mismatches = 0 unless defined $num_mismatches;

# Program dependencies - process_radtags program from the Stacks Software Suite.
my $process_radtags 				= '/usr/local/bin/process_radtags';

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fastq_infile -r restriction_enzymes -b barcode_infile -n num_mismatches -o output_dir
    
DESCRIPTION - This program performs a quality assessment and quality control filtering step to filter out GBS fastq reads that do not meet the quality threshold given in the quality scores. It also demultiplexes the original raw bulk GBS fastq file based on barcode into separate samples that include the barcode in the filename. The quality filtering and quality threshold steps are performed using the process_radtags program in the STACKS software suite. Once the raw fastq file is demultiplexed by barcode the resulting files are renamed corresponding to the individual name, plate/well number, and barcode sequence and copied to the corresponding project leader directory.

OPTIONS:

-i fastq_infile - The absolute path to the bulk fastq input file to split sequences based on barcode sequences. Can either be *.fastq or *.fastq.gz extension.

-r restriction_enzymes - The restriction enzyme(s) used to digest the genomic sequences. Default: Pstl/MspI

-b barcode_infile - The absolute path to the barcodes input file used to split sequences into individual fastq output files.

-n num_mismatches - The number of mismatches allowed within the barcode sequences.

-o output_dir - The absolute path to the output directory to contain the split *.fastq output files.

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Execute the process_radtags program for quality assessment and quality control filtering and demultiplexing the files based on barcode and project leader name.
my ($split_fastq_files, $split_fastq_file_counter) = process_radtags($fastq_infile, $restriction_enzymes, $barcode_infile, $num_mismatches, $output_dir);
	
# Parse the barcodes input file based on barcode sequence so that we can split the bulk fastq file using the process_radtags program.
my %project_leader_barcodes = ();
open(INFILE, "<$barcode_infile") or die "Couldn't open file $barcode_infile for reading, $!";
my $i = 0;
while(<INFILE>){
	chomp $_;
	warn $_ . "\n";
	if($i ne 0){
		my @split_row_entry = split(/\t/, $_);
		my ($fastq_plate_num, $fastq_well_row, $fastq_well_column, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) = @split_row_entry;

		my $fastq_barcode_name = join("_", $fastq_run_id, $fastq_barcode_seq, $fastq_plate_num, join("", $fastq_well_row, $fastq_well_column));
		$project_leader_barcodes{$fastq_barcode_seq} = join("\t", $fastq_barcode_name, $fastq_barcode_seq, $fastq_project_leader);
		
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

# Iterate through each split fastq sample file from the process_radtags program and rename the files based on project leader.
foreach my $file_name (sort keys %{$split_fastq_files}){
	warn "Processing " . $file_name . ".....\n";
	my $split_fastq_infile = $split_fastq_files->{$file_name};
	
	# Get the basename of the fastq filename without the .fastq extension.
	my $fastq_filename = fileparse($split_fastq_infile, qr/\.fq/);
	
	my ($split_sample_name, $split_fastq_barcode_seq) = split(/_/, $fastq_filename);
	
	my @split_row_entry = split(/\t/, $project_leader_barcodes{$split_fastq_barcode_seq});
	my ($fastq_barcode_name, $fastq_barcode_seq, $fastq_project_leader) = @split_row_entry;
	
	# Create project output directory if it doesn't already exist.
	my $project_output_dir = join('/', $project_leader_dir, $fastq_project_leader);
	unless(-d $project_output_dir){
		mkdir($project_output_dir, 0777) or die "Can't make directory: $!";
	}
	
	my $fastq_barcode_project_filename = join('/', $project_output_dir, join("", $fastq_barcode_name, ".fastq"));
	warn "Copying $file_name to $fastq_barcode_project_filename.....";
	copy($split_fastq_infile, $fastq_barcode_project_filename) or die "Copy failed: $!";
# 	unlink($split_fastq_infile) or die "Could not unlink $split_fastq_infile: $!";
}

# ($split_fastq_files, $split_fastq_file_counter) = process_radtags($fastq_infile, $restriction_enzymes, $barcode_infile, $num_mismatches, $output_dir) - Executes the process_radtags program from the Stacks Software Suite for quality assessment and quality control filtering and demultiplexing the files based on barcode and project leader name.
#
# Input paramater(s):
# 
# $fastq_infile - The absolute path to the bulk fastq input file to split sequences based on barcode sequences.
#
# $restriction_enzymes - The restriction enzyme(s) used to digest the genomic sequences.
# 
# $barcode_infile - The absolute path to the barcodes input file used to split sequences into individual fastq output files.
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
sub process_radtags{

	# The raw bulk GBS fastq input file to filter based on quality and demultiplex based on barcodes.
	my $fastq_infile = shift;
	die "Error lost fastq input file" unless defined $fastq_infile;
	
	# The restriction enzyme(s) used to digest the genomic sequences.
	my $restriction_enzymes = shift;
	die "Error lost restriction enzyme(s) used for GBS restriction digest protocol" unless defined $restriction_enzymes;
	
	# The absolute path to the barcodes input file used to split sequences into individual fastq output files.
	my $barcode_infile = shift;
	die "Error lost barcode sequence file" unless defined $barcode_infile;
	
	# The number of mismatches allowed within the barcode sequences.
	my $num_mismatches = shift;
	die "Error lost number of mismatches" unless defined $num_mismatches;
	
	# The absolute path to the output directory to contain the split *.fastq output files.
	my $output_dir = shift;
	die "Error lost output directory" unless defined $output_dir;

	# Create output directory if it doesn't already exist.
	my $split_fastq_output_dir = join('/', $output_dir, "SPLIT_FASTQ_OUTFILES");
	unless(-d $split_fastq_output_dir){
		mkdir($split_fastq_output_dir, 0777) or die "Can't make directory: $!";
	}
	
	# Check the $split_fastq_output_dir for files with extension *.fq and check to make sure the size of the resulting files are non-zero 
	# indicating that the $fastq_infile file hasn't been already filtered based on quality and demultiplexed into separate files based on barcodes.  
	my ($split_fastq_files, $split_fastq_file_counter) = find_files($split_fastq_output_dir, "fq");
	my $non_zero_split_fastq_files = 0;
	foreach my $file_name (sort keys %{$split_fastq_files}){
		warn $file_name . "\n";
		if(-s $split_fastq_files->{$file_name}){
			$non_zero_split_fastq_files++;
		}
	}
	
	# Check if the $fastq_infile file has been already filtered based on quality and demultiplexed into separate files based on barcodes.
	# If it hasn't, process the $fastq_infile file using the following process_radtags commands.
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
				my ($fastq_plate_num, $fastq_well_row, $fastq_well_column, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) = ($split_row_entry[0], $split_row_entry[1], $split_row_entry[2], $split_row_entry[3], $split_row_entry[4], $split_row_entry[5]);

				my $barcode_seq_length = length($fastq_barcode_seq);
				my $barcodes_entry = join("\t", $barcode_seq_length, $fastq_barcode_seq);
	# 			warn $barcodes_entry . "\n";

				push(@barcode_entries, [split(/\t/, $barcodes_entry)]);
				
			}
			$i++;
		}
		close(INFILE) or die "Couldn't close file $barcode_infile";
		
		# Generate the barcodes input file for the process_radtags program.
		my $barcode_outfile = join('/', $output_dir, "process_radtags_barcodes.txt");
		open(OUTFILE, ">$barcode_outfile") or die "Couldn't open file $barcode_outfile for writting, $!";
		foreach my $barcodes_entry (sort {$b->[0] <=> $a->[0] || $a->[1] cmp $b->[1]} @barcode_entries){
			print OUTFILE @$barcodes_entry[1] . "\n";
		}
		close(OUTFILE) or die "Couldn't close file $barcode_outfile";

		# Choose the process_radtags command based on the restriction enzyme(s) used and whether or not the bulk fastq file is compressed.
		my $process_radtagsCmd = "";
		switch($restriction_enzymes){
			case 'ApeKI'{ # ApeKI: Not operational yet
				if($fastq_infile =~ m/\.gz$/){ # If the bulk fastq file is compressed with extension *.gz.
					# process_radtags -f ./HI.1405.008.GQ03122013-5_R1.fastq.gz -i gzfastq -b ./process_radtags_barcodes.txt -o ./GBS_filtered_fastq -y fastq -c -q -r -E phred33 -D -w 0.15 -s 10 --inline_null -e apeKI --filter_illumina --barcode_dist 0
					$process_radtagsCmd = "$process_radtags -f $fastq_infile -i gzfastq -b $barcode_outfile -o $split_fastq_output_dir -y fastq -c -q -r -E phred33 -D -w 0.15 -s 10 --inline_null -e apeKI --filter_illumina --barcode_dist $num_mismatches";
				}else{ # If the bulk fastq file is not compressed with extension *.fastq.
					# process_radtags -f ./HI.1405.008.GQ03122013-5_R1.fastq -i fastq -b ./process_radtags_barcodes.txt -o ./GBS_filtered_fastq -y fastq -c -q -r -E phred33 -D -w 0.15 -s 10 --inline_null -e apeKI --filter_illumina --barcode_dist 0
					$process_radtagsCmd = "$process_radtags -f $fastq_infile -i fastq -b $barcode_outfile -o $split_fastq_output_dir -y fastq -c -q -r -E phred33 -D -w 0.15 -s 10 --inline_null -e apeKI --filter_illumina --barcode_dist $num_mismatches";
				}
			}
			case 'Pstl/MspI'{ # Pstl/MspI
				if($fastq_infile =~ m/\.gz$/){ # If the bulk fastq file is compressed with extension *.gz.
					# process_radtags -f ./HI.1405.008.GQ03122013-5_R1.fastq.gz -i gzfastq -b ./process_radtags_barcodes.txt -o ./GBS_filtered_fastq -y fastq -c -q -r -E phred33 -D -w 0.15 -s 10 --inline_null --renz_1 pstI --renz_2 mspI --filter_illumina --barcode_dist 0
					$process_radtagsCmd = "$process_radtags -f $fastq_infile -i gzfastq -b $barcode_outfile -o $split_fastq_output_dir -y fastq -c -q -r -E phred33 -D -w 0.15 -s 10 --inline_null --renz_1 pstI --renz_2 mspI --filter_illumina --barcode_dist $num_mismatches";
				}else{ # If the bulk fastq file is not compressed with extension *.fastq.
					# process_radtags -f ./HI.1405.008.GQ03122013-5_R1.fastq -i fastq -b ./process_radtags_barcodes.txt -o ./GBS_filtered_fastq -y fastq -c -q -r -E phred33 -D -w 0.15 -s 10 --inline_null --renz_1 pstI --renz_2 mspI --filter_illumina --barcode_dist 0
					$process_radtagsCmd = "$process_radtags -f $fastq_infile -i fastq -b $barcode_outfile -o $split_fastq_output_dir -y fastq -c -q -r -E phred33 -D -w 0.15 -s 10 --inline_null --renz_1 pstI --renz_2 mspI --filter_illumina --barcode_dist $num_mismatches";
				}
			}
			else{ 
				die "Input $restriction_enzymes is not one of the recognized restriction enzyme(s)! Please specify either ApeKI or Pstl/MspI on the command line"
			}
		}
		warn $process_radtagsCmd . "\n\n";
		
		# Execute the process_radtags command using system.
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
	die "Error lost file extension suffix" unless defined $suffix;
	
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
}