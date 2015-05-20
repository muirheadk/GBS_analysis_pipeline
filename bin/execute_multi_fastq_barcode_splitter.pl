#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

my ($GBS_fastq_file_input_dir, $GBS_barcode_input_dir, $output_dir);
GetOptions(
	'i=s'    => \$GBS_fastq_file_input_dir, # The absolute path to the GBS raw fastq file directory.
	'b=s'    => \$GBS_barcode_input_dir, # The absolute path to the GBS barcode file directory.
	'o=s'    => \$output_dir,  # The absolute path to the output directory to contain demultiplexed files.
);

# Display a usage message if the following parameters are not specified.
usage() unless (
	defined $GBS_fastq_file_input_dir
	and defined $GBS_barcode_input_dir
	and defined $output_dir
);

# Program dependencies - process_radtags program from the Stacks Software Suite.
my $fastq_quality_barcode_splitter 				= '/home/muirheadk/fastq_quality_barcode_splitter.pl';

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i GBS_fastq_file_input_dir -b GBS_barcode_input_dir -o output_dir

-i GBS_fastq_file_input_dir - The absolute path to the GBS raw fastq file directory.

-b GBS_barcode_input_dir - The absolute path to the GBS barcode file directory.

-o output_dir - The absolute path to the output directory to contain demultiplexed files.

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my ($fastq_files, $fastq_file_counter) = find_raw_fastq_files($GBS_fastq_file_input_dir, "fastq.gz");

my ($barcode_files, $barcode_file_counter) = find_barcode_files($GBS_barcode_input_dir, "txt");

foreach my $flowcell_name (keys %{$barcode_files}){

	my $barcode_infile = $barcode_files->{$flowcell_name};
 	
 	my $fastq_infile = $fastq_files->{$flowcell_name};
	
	my $process_radtags_output_dir = join('/', $output_dir, "PROCESSED_RADTAGS");

	my $fastq_quality_barcode_splitterCmd = "$fastq_quality_barcode_splitter -i $fastq_infile -e phred33 -r PstI/MspI -b $barcode_infile -n 0 -o $process_radtags_output_dir";
	warn $fastq_quality_barcode_splitterCmd . "\n\n";
 	my $status = system($fastq_quality_barcode_splitterCmd) == 0 or die "Error calling $fastq_quality_barcode_splitter: $?";

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
sub find_raw_fastq_files{
    
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
			my $flowcell_name;
			if($file_name =~ m/.+(GQ\d+)(-|_)(L(\d+)|(\d+)).+\.fastq\.gz$/){
				my $flowcell = $1;
				my $plate = $3;
				$plate =~ s/L//g;
				$flowcell_name = join("_", $flowcell, $plate);
			}elsif($file_name =~ m/(GQ\d+)(-|_)(\d+)\.fastq\.gz$/){
				my $flowcell = $1;
				my $plate = $3;
				$flowcell_name = join("_", $flowcell, $plate);
			}
			my $infile_name = join('/', $infile_dir, $file_name) if ($file_name =~ m/\.$suffix$/);
			warn "$infile_name\n" if ($file_name =~ m/\.$suffix$/);
			$files{$flowcell_name} = $infile_name if ($file_name =~ m/\.$suffix$/);
			$file_counter++ if ($file_name =~ m/\.$suffix$/);
		}
		closedir(DIR);
		return (\%files, $file_counter);
	}else{
		die "Error $infile_dir does not exist!\n";
	}
}

sub find_barcode_files{
    
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
			my $flowcell_name;
			if($file_name =~ m/(\w+\_\d+)-barcodes.txt/){
				$flowcell_name = $1;
			}
			my $infile_name = join('/', $infile_dir, $file_name) if ($file_name =~ m/\.$suffix$/);
			warn "$infile_name\n" if ($file_name =~ m/\.$suffix$/);
			$files{$flowcell_name} = $infile_name if ($file_name =~ m/\.$suffix$/);
			$file_counter++ if ($file_name =~ m/\.$suffix$/);
		}
		closedir(DIR);
		
		return (\%files, $file_counter);
	}else{
		die "Error $infile_dir does not exist!\n";
	}
}