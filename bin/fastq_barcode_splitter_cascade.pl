#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Copy;

# perl fastq_barcode_splitter_cascade.pl -i ~/workspace/GBS_data-08-10-2013/HI.1405.008.GQ03122013-5_R1.fastq -b ~/workspace/GBS_data-08-10-2013/GBS_barcodes-2013-10-09.csv -o ~/workspace/GBS_data-08-10-2013
my ($fastq_infile, $barcode_infile, $num_mismatches, $output_dir);
GetOptions(
      'i=s'    => \$fastq_infile,
      'b=s'    => \$barcode_infile,
      'n=s'    => \$num_mismatches,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $fastq_infile
      and defined $barcode_infile
      and defined $output_dir
);
$num_mismatches = 1 unless defined $num_mismatches;

my ($gunzip, $fastx_barcode_splitter, $send_mail);
$gunzip				= '/bin/gunzip';
$fastx_barcode_splitter 	= '/usr/local/bin/fastx_barcode_splitter.pl';
$send_mail			= '/TRIA-NetUtils/bin/send_mail.pl';

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fastq_infile -b barcode_infile -m metadata_infile -n num_mismatches -o output_dir
    
Description - 
    
OPTIONS:
      -i fastq_infile - 
    
      -b barcode_infile -
      
      -n num_mismatches - 
      
      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}


# Create output directory if it doesn't already exist.
my $split_fastq_output_dir = join('/', $output_dir, "SPLIT_FASTQ_OUTFILES");
unless(-d $split_fastq_output_dir){
      mkdir($split_fastq_output_dir, 0777) or die "Can't make directory: $!";
}

my %barcode_data = ();
open(INFILE, "<$barcode_infile") or die "Couldn't open file $barcode_infile for reading, $!";
my $i = 0;
while(<INFILE>){
	chomp $_;
	warn $_ . "\n";
	if($i ne 0){
		my @split_row_entry = split(/\t/, $_);
		my ($fastq_plate_num, $fastq_well_num, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) = ($split_row_entry[0], $split_row_entry[1], $split_row_entry[2], $split_row_entry[3], $split_row_entry[4]);

		my $barcode_seq_length = length($fastq_barcode_seq);
		my $fastq_barcode_name = join("_", $fastq_run_id, $fastq_barcode_seq, $fastq_plate_num, $fastq_well_num);
		push(@{$barcode_data{$barcode_seq_length}},  join("\t", $fastq_barcode_name, $fastq_barcode_seq, $fastq_project_leader));
		
	}
	$i++;
}
close(INFILE) or die "Couldn't close file $barcode_infile";

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


foreach my $fastq_project_leader (sort keys %barcode_filenames){
	
	foreach my $fastq_barcode_filenames (@{$barcode_filenames{$fastq_project_leader}}){
		
		my @split_fastq_barcode_filenames = split(/\t/, $fastq_barcode_filenames);
		my ($fastq_barcode_length_filename, $fastq_barcode_project_filename) = ($split_fastq_barcode_filenames[0],$split_fastq_barcode_filenames[1]);
            	copy($fastq_barcode_length_filename, $fastq_barcode_project_filename) or die "Copy failed: $!";
	}
}

    my $subject = "$0 Process Complete";
    my $message = "$0 process finished successfully! You can find the fastq output files for each project leader in the $split_fastq_output_dir/PROJECT_LEADER_DIR directory.";
send_mail($subject, $message, 'email');


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
# 	warn "Generating blastp alignment file....\n";
	my $fastx_barcode_splitterCmd  = "$fastx_barcode_splitter --bcfile $barcode_infile --bol --mismatches $num_mismatches --prefix $fastq_output_dir/ --suffix \".fastq\" < $fastq_infile";
	warn $fastx_barcode_splitterCmd . "\n\n";

	my $status = system($fastx_barcode_splitterCmd) == 0 or die "Error calling $fastx_barcode_splitter: $?";
}

sub send_mail{
    my $subject = shift or die "lost email subject";
    my $message = shift or die "lost email message";
    my $email_type = shift or die "lost email type";
    
    
		warn "$send_mail -s $subject -m $message -t $email_type\n\n";
		system($send_mail, 
			'-s', "\"$subject\"", 
			'-m', "\"$message\"", 
			'-t', $email_type
		) == 0 or die "Error calling $send_mail -s $subject -m $message -t $email_type: $?";

}
