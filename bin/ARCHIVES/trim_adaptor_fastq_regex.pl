#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;
use List::Compare;

# perl trim_adaptor_seq_fastq-new.pl -i ~/workspace/GBS_data-08-10-2013/PROJECT_LEADER_DIR/CHRISTIANNE_MCDONALD -p CHRISTIANNE_MCDONALD -c 7 -o ~/workspace/GBS_data-08-10-2013/TRIMMED_ADAPTOR_FASTQ_DIR-2014-10-29
my ($fastq_file_dir, $project_name, $fastq_adaptor_sequence, $gbs_sequence_length, $adaptor_length_min_threshold, $adaptor_trim_offset, $min_trimmed_fastq_sequence_length, $blast_num_cpu, $output_dir);
GetOptions(
	'i=s'    => \$fastq_file_dir, # The *.fastq input file directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	'p=s'    => \$project_name, # The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	'a=s'    => \$fastq_adaptor_sequence, # The GBS common adaptor sequence to trim from the GBS fastq sequences. Default: CCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
	'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 100
	'm=s'    => \$adaptor_length_min_threshold, # The minimum GBS common adaptor sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adaptor regex searches. Default: 16
	't=s'    => \$adaptor_trim_offset, # The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adaptor sequence found in the adaptor regex searches. Default: 5
	'q=s'    => \$min_trimmed_fastq_sequence_length, # The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Default: 32
	'c=s'    => \$blast_num_cpu, # The number of cpu cores to use for the adaptor regex searches. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2
	'o=s'    => \$output_dir, # The absolute path to the output directory to contain the trimmed adaptor sequence fastq output files.
);

# Display a usage message if the following parameters are not specified.
usage() unless (
	defined $fastq_file_dir
	and $project_name
	and defined $output_dir
);

# The GBS common adaptor sequence to trim from the GBS fastq sequences. Default: CCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
$fastq_adaptor_sequence = 'CCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG' unless defined $fastq_adaptor_sequence;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 100
$gbs_sequence_length = 100 unless defined $gbs_sequence_length;

# The minimum GBS common adaptor sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adaptor regex searches. Default: 1
$adaptor_length_min_threshold = 1 unless defined $adaptor_length_min_threshold;

# The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adaptor sequence found in the adaptor regex searches. Default: 5
$adaptor_trim_offset = 5 unless defined $adaptor_trim_offset;

# The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Keep in mind that this is the minimum trimmed fastq sequence length used before we add the length of the barcode 
# used for splitting each individual fastq file. Tassel requires sequences at least 32 base pairs (bps) plus the length of a particular barcode that can be in the range of 4-8 base pairs (bps) in length. Default: 32
$min_trimmed_fastq_sequence_length = 32  unless defined $min_trimmed_fastq_sequence_length;

# The number of cpu cores to use for the adaptor regex searches. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2
$blast_num_cpu = 2 unless defined $blast_num_cpu;

# Program dependencies - The absolute paths to gzip to compress all project leader *.fastq input files.
my $gzip 				= '/bin/gzip';

sub usage {

die <<"USAGE";


Usage: $0 -i fastq_file_dir -p project_name -a fastq_adaptor_sequence -l gbs_sequence_length -m adaptor_length_min_threshold -t adaptor_trim_offset -q min_trimmed_fastq_sequence_length -c blast_num_cpu -o output_dir

DESCRIPTION - A program to trim the GBS common adaptor sequence from each GBS fastq file within a particular Genotyping by Sequencing (GBS) project. Fixes the misprimming issue where the GBS common adaptor is sequenced along with the DNA of an individual

OPTIONS:

-i fastq_file_dir - The *.fastq input file directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	e.g. /path/to/fastq_file_dir
	
-p project_name - The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	e.g. MPB_MALE_GBS
	
-a fastq_adaptor_sequence - The GBS common adaptor sequence to trim from the GBS fastq sequences. Default: CCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG

-l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 100

-m adaptor_length_min_threshold - The minimum GBS common adaptor sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adaptor regex searches. Default: 1

-t adaptor_trim_offset - The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adaptor sequence found in the adaptor regex searches. Default: 5

-q min_trimmed_fastq_sequence_length - The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Keep in mind that this is the minimum trimmed fastq sequence length used before we add the 
length of the barcode used for splitting each individual fastq file. Tassel requires sequences at least 32 base pairs (bps) plus the length of a particular barcode that can be in the range of 4-8 base pairs (bps) in 
length. Default: 32

-c blast_num_cpu - The number of cpu cores to use for the adaptor regex searches. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2

-o output_dir - The absolute path to the output directory to contain the trimmed adaptor sequence fastq output files.

	e.g. /path/to/output_dir
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Get the length of the adaptor so that we can verify adaptor blast results.
my $adaptor_sequence_length = length($fastq_adaptor_sequence);

# Find all files in the specified directory with the extension *.fastq.
my ($fastq_files, $file_count) = find_fastq_files($fastq_file_dir);

# Container for compressing fastq files using gzip.
my @fastq_files2gzip = ();

# Create output directory if it doesn't already exist.
my $project_dir = join('/', $output_dir, $project_name);
unless(-d $project_dir){
	mkdir($project_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $fasta_output_dir = join('/', $project_dir, "FASTA_FILES");
unless(-d $fasta_output_dir){
	mkdir($fasta_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $regex_output_dir = join('/', $project_dir, "ADAPTOR_REGEX_FILES");
unless(-d $regex_output_dir){
	mkdir($regex_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $trimmed_output_dir = join('/', $project_dir, "TRIMMED_OUTPUT_FILES");
unless(-d $trimmed_output_dir){# Need this to make the bulk fastq sequence file.
	mkdir($trimmed_output_dir, 0777) or die "Can't make directory: $!";
}

# Generate bulk trimmed fastq files for each Genotyping by Sequencing (GBS) project. 
# $trimmed_fastq_bulk_outfile contains the trimmed GBS sequences in fastq format.
my $trimmed_fastq_bulk_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adaptor_trim_offset) . ".fastq");
push(@fastq_files2gzip, $trimmed_fastq_bulk_outfile);
open(TRIMMED_BULK_OUTFILE, ">$trimmed_fastq_bulk_outfile") or die "Couldn't open file $trimmed_fastq_bulk_outfile for writting, $!";


# The $trimmed_seqs_layout_bulk_outfile contains alled the trimmed coordinates layout for each sequence trimmed of the GBS common adaptor sequence.
my $trimmed_seqs_layout_bulk_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adaptor_trim_offset, "all_trimmed_seqs_layout") . ".txt");
open(TRIMMED_BULK_LAYOUT_OUTFILE, ">$trimmed_seqs_layout_bulk_outfile") or die "Couldn't open file $trimmed_seqs_layout_bulk_outfile for writting, $!";
print TRIMMED_BULK_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adaptor_start", "trimmed_adaptor_end",
"trimmed_adaptor_length", "adaptor_seq_start", "adaptor_seq_end", "adaptor_seq_length") . "\n";

# The $trimmed_seqs_layout_outfile contains the trimmed coordinates layout for each sequence trimmed of the GBS common adaptor sequence that passed the retaining criteria.
my $trimmed_seqs_layout_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adaptor_trim_offset, "trimmed_seqs_layout") . ".txt");
open(TRIMMED_LAYOUT_OUTFILE, ">$trimmed_seqs_layout_outfile") or die "Couldn't open file $trimmed_seqs_layout_outfile for writting, $!";
print TRIMMED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adaptor_start", "trimmed_adaptor_end",
"trimmed_adaptor_length", "adaptor_seq_start", "adaptor_seq_end", "adaptor_seq_length") . "\n";

# The $removed_seqs_layout_outfile contains the trimmed coordinates layout for each sequence trimmed of the GBS common adaptor sequence that failed the retaining criteria and was therefore removed.
my $removed_seqs_layout_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adaptor_trim_offset, "removed_seqs_layout") . ".txt");
open(REMOVED_LAYOUT_OUTFILE, ">$removed_seqs_layout_outfile") or die "Couldn't open file $removed_seqs_layout_outfile for writting, $!";
print REMOVED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adaptor_start", "trimmed_adaptor_end",
"trimmed_adaptor_length", "adaptor_seq_start", "adaptor_seq_end", "adaptor_seq_length") . "\n";
my %fastq_seq_counter = ();
my %adaptor_length_counter = ();

# Iterate through the files with the extension *.fastq.
foreach my $fastq_filename (sort keys %{$fastq_files}){
	my $fastq_infile = $fastq_files->{$fastq_filename};
	warn "Processing " . $fastq_infile . ".....\n";
	
	my %fastq_sequences = ();
	
	# open the fastq input file for parsing.
	open(INFILE, "<$fastq_infile") or die "Couldn't open file $fastq_infile for reading, $!";
	
	# Get the basename of the fastq filename without the .fastq extension.
	my $fasta_filename = fileparse($fastq_infile, qr/\.fastq/);
	
	# Need to add the length of the barcode to the $min_trimmed_fastq_sequence_length so that we get the minimum trimmed fastq sequence length plus the length of the barcode sequence.
	my ($individual_id, $barcode, $plate_num, $well_num) = split(/_/, $fasta_filename);
	my $min_trimmed_fastq_sequence_length_plus_barcode = ($min_trimmed_fastq_sequence_length + length($barcode));
	
	# Generate files in fasta format from each fastq file for use in the adaptor regex searches
	my $fasta_target_outfile = join("/", $fasta_output_dir, $fasta_filename . ".fasta");
	open(OUTFILE, ">$fasta_target_outfile") or die "Couldn't open file $fasta_target_outfile for writting, $!";
	my ($fastq_header, $fastq_sequence, $fastq_plus, $fastq_quality_scores);
	my $i = 1;
	my $num_fastq_seqs = 0;
	
	# Parse the fastq files for the fastq header, sequence, plus, and quality scores and reformat to *.fasta format so that we can use adaptor regex.
	while(<INFILE>){
		chomp $_;
		#warn $_ . "\n";
		if(($_ =~ m/^\@[A-Za-z0-9-_]+:\d+:[A-Za-z0-9]+:\d+:\d+:\d+:\d+ \d:[A-Z]:\d:[ACGTRYKMSWBDHVN]*$/) 
			and ($i eq 1)){ # The fastq sequence header is on the first line. i.e. @HWI-ST767:215:C30VBACXX:8:1101:1801:1484 1:N:0:
			$fastq_header = $_;
#  				die $fastq_header;
		}elsif(($_ =~ m/^[ACGTRYKMSWBDHVN]+$/i) and ($i eq 2)){ # The fastq sequence is on the second line.
			$fastq_sequence = $_;
#  				die $fastq_sequence;
		}elsif(($_ =~ m/^\+$/) and ($i eq 3)){ # The fastq plus character is on the third line.
			$fastq_plus = $_;
#  				die $fastq_plus;
		}elsif(($_ =~ m/^.+$/) and (($i % 4) eq 0)){ # the fastq quality scores are on the fourth line.
			$fastq_quality_scores = $_;
#   				die $fastq_quality_scores;
		}
		
		if(($i % 4) eq 0){ # Once we are finished parsing a fastq sequence entry we check to make sure it was parsed correctly and that the lengths of the sequence and quality scores match. 
			
			die "Error: fastq_header is undefined" unless(defined($fastq_header));
			die "Error: fastq_sequence is undefined" unless(defined($fastq_sequence));
			die "Error: fastq_plus is undefined" unless(defined($fastq_plus));
			die "Error: fastq_quality_scores is undefined" unless(defined($fastq_quality_scores));
			
			my $fastq_sequence_length = length($fastq_sequence);
			my $fastq_quality_scores_length = length($fastq_quality_scores);
			
			# Store each fastq sequence entry in a hash variable indexed to the fastq header for identification and an attribute index for the fastq sequence, fastq plus character, and the fastq quality score.
			$fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $fastq_sequence;
			$fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
			$fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $fastq_quality_scores;
			die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($fastq_sequence_length ne $gbs_sequence_length);
			die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length ne fastq_quality_scores_length=$fastq_quality_scores_length" if($fastq_sequence_length ne $fastq_quality_scores_length);
			
			# Print the fastq header and sequence to the fasta file for the adaptor regex searches in fasta format.
			print OUTFILE join("\n", join("", ">", $fastq_header), $fastq_sequence) . "\n";
			
			$i = 1;
			$num_fastq_seqs++;
			
		}else{
			$i++;
		}
	}
	close(INFILE) or die "Couldn't close file $fastq_infile";
	close(OUTFILE) or die "Couldn't close file $fasta_target_outfile";
	
	# Execute the adaptor regex search using the GBS common adaptor sequence as the query and the fastq sequences in fasta format as the blast database.
	# Keep a maximum of target sequences equal to the number of sequences in the fastq file. Retain adaptor regex hits based on the minimum adaptor length threshold.
	# Keep adaptor regex hits that begin at position 1 of the full GBS common adaptor sequence up to the minimum adaptor length threshold.
	my $adaptor_regex_outfile = join('/', $regex_output_dir, join("_", $fasta_filename, "gbs_adaptor_regex.tsv"));
	open(OUTFILE, ">$adaptor_regex_outfile") or die "Couldn't open file $adaptor_regex_outfile for writting, $!";
	print OUTFILE join("\t", "query_name", "target_name", "align_length", "query_start", "query_end", "target_start", "target_end") . "\n";
	my %regex_align_hits = ();
	foreach my $fastq_header (sort keys %fastq_sequences){
		my $fastq_sequence = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
		my ($align_length, $query_start, $query_end, $target_start, $target_end);
		my $adaptor_length_count = 1;
		my $common_adaptor_sequence = $fastq_adaptor_sequence;
		my $align_found = "false";
		for(my $i = $adaptor_sequence_length; $i > $adaptor_length_min_threshold; $i--){
			my $common_adaptor_regex = qr($common_adaptor_sequence);
			if($fastq_sequence =~ /$common_adaptor_regex/g){
				print "$i eq $adaptor_sequence_length\n";
				if($i eq $adaptor_sequence_length){
					$target_start = ($-[0] + 1);
					$target_end = $+[0];
					$align_length = ($+[0] - $-[0]);
					print join("\t", $align_length, $target_start, $target_end) . "\n";
					$align_found = "true";
					$query_start = 1;
					$query_end = $adaptor_sequence_length;
					warn join("\t", $adaptor_sequence_length, $common_adaptor_sequence) . "\n";
					last;
				}elsif($i < $adaptor_sequence_length){
				
					if($+[0] eq $gbs_sequence_length){
						$target_start = ($-[0] + 1);
						$target_end = $+[0];
						$align_length = ($+[0] - $-[0]);
						print join("\t", $align_length, $target_start, $target_end) . "\n";
						$align_found = "true";
						last;
					}
				}
			}
			
			$query_start = 1;
			$query_end = ($adaptor_sequence_length - $adaptor_length_count);
			
			$common_adaptor_sequence = get_subseq($common_adaptor_sequence, $query_start, $query_end);
			warn join("\t", ($adaptor_sequence_length - $adaptor_length_count), $common_adaptor_sequence) . "\n";
			$adaptor_length_count++;
		}
		if($align_found eq "true"){
			my $regex_alignment = join("\t", join("_", "GBS_adaptor_sequence", $fastq_adaptor_sequence), $fastq_header, $align_length, $query_start, $query_end, $target_start, $target_end);
			print OUTFILE $regex_alignment . "\n";
		}
	}
	close(OUTFILE) or die "Couldn't close file $adaptor_regex_outfile";
  
  	# Create output directory if it doesn't already exist.
  	my $trimmed_regex_output_dir = join('/', $trimmed_output_dir, "TRIMMED_ADAPTOR_REGEX_FILES");
  	unless(-d $trimmed_regex_output_dir){
  		mkdir($trimmed_regex_output_dir, 0777) or die "Can't make directory: $!";
  	}
  	
  	# Generate the trimmed adaptor regex files filtered for further processing.
  	my $trimmed_adaptor_regex_outfile = join('/', $trimmed_regex_output_dir, join("_", $fasta_filename, "gbs_adaptor_regex.tsv"));
  	open(OUTFILE, ">$trimmed_adaptor_regex_outfile") or die "Couldn't open file $trimmed_adaptor_regex_outfile for writting, $!";
  	print OUTFILE join("\t", "query_name", "target_name", "align_length", "query_start", "query_end", "target_start", "target_end") . "\n";
  	my %trimmed_fastq_sequences = ();
  	my @trimmed_fastq_list = ();
  	my %removed_fastq_sequences = ();
  	my @removed_fastq_list = ();
  	# Parse the tab-delimited adaptor blast output so that we can trimm the fastq sequences that contain the GBS adaptor sequence.
  	# Create new adaptor regex files so that we can visualize where we trimmed the sequence. Store results in a hash array indexed to the target sequence so that 
  	# we can use the longest length of the GBS common adaptor sequence as the optimal alignment in the trimming step.
  	warn "Parsing $adaptor_regex_outfile to get most optimal adaptor regex hits.....\n";
  	open(INFILE, "<$adaptor_regex_outfile") or die "Couldn't open file $adaptor_regex_outfile for reading, $!";
  	my %adaptor_regex_hits = ();
  	$i = 0;
  	while(<INFILE>){
  		chomp $_;
  		if($i ne 0){
  			#warn $_ . "\n";
  			my @split_adaptor_regex_hit =  split(/\t/, $_);
  			my ($query_name, $target_name, $align_length, $query_start, $query_end, $target_start, $target_end) = @split_adaptor_regex_hit;
  		
			die "Error: query_name is undefined" unless(defined($query_name));
			die "Error: target_name is undefined" unless(defined($target_name));
			die "Error: align_length is undefined" unless(defined($align_length));
			die "Error: query_start is undefined" unless(defined($query_start));
			die "Error: query_end is undefined" unless(defined($query_end));
			die "Error: target_start is undefined" unless(defined($target_start));
			die "Error: target_end is undefined" unless(defined($target_end));
			
			my $fastq_header = $target_name;
			
			my $fastq_sequence = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
			my $fastq_plus = $fastq_sequences{$fastq_header}{'PLUS'}; 
			my $fastq_quality_scores = $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
			
			die "Error: fastq_sequence is undefined" unless(defined($fastq_sequence));
			die "Error: fastq_plus is undefined" unless(defined($fastq_plus));
			die "Error: fastq_quality_scores is undefined" unless(defined($fastq_quality_scores));
			
			my $fastq_sequence_length = length($fastq_sequence);
			my $predicted_adaptor_start = ($fastq_sequence_length - $adaptor_sequence_length);
			
			my ($trimmed_fastq_sequence, $trimmed_fastq_quality_scores, $trimmed_adaptor_sequence);
			# If the adaptor alignment length is equal to the length of the full GBS common adaptor then trim where that GBS common adaptor sequence is found subtracting the trimmed offset from the target start.
			# If the adaptor alignment length is less than the length of the full GBS common adaptor then trim only if the target end of the aligned adaptor sequence is equal the the common GBS sequence length.
			if(($align_length eq $adaptor_sequence_length) 
				or ((($align_length >= $adaptor_length_min_threshold) and ($align_length < $adaptor_sequence_length)) and ($query_start eq 1) and ($target_end eq $gbs_sequence_length))){
				
				# Get the trimmed fastq sequence trimmed of the GBS common adaptor sequence and trimmed offset.
				$trimmed_fastq_sequence = get_subseq($fastq_sequence, 1, (($target_start - $adaptor_trim_offset) - 1));
				# Get the trimmed fastq quality scores trimmed of the GBS common adaptor sequence and trimmed offset.
				$trimmed_fastq_quality_scores = get_subseq($fastq_quality_scores, 1, (($target_start - $adaptor_trim_offset) - 1));
				# Get the trimmed adaptor sequence and trimmed offset that was trimmed off the GBS fastq sequence.
				$trimmed_adaptor_sequence = get_subseq($fastq_sequence, ($target_start - $adaptor_trim_offset), $fastq_sequence_length);
				
				# Get the length of the trimmed GBS common adaptor sequence and trimmed offset.
				my $trimmed_adaptor_sequence_length = length($trimmed_adaptor_sequence);
				my $trimmed_fastq_sequence_length = length($trimmed_fastq_sequence);
				
				# If the trimmed sequence length is greater than or equal to the minimum trimmed fastq sequence length plus the length of the barcode then the trimmed and sequence count files.
				if($trimmed_fastq_sequence_length >= $min_trimmed_fastq_sequence_length_plus_barcode){
					# Store the trimmed fastq sequence entry in a hash variable indexed by the fastq sequence header.
					$trimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $trimmed_fastq_sequence;
					$trimmed_fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
					$trimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $trimmed_fastq_quality_scores;
					
					my $trimmed_fastq_regex = join("\t", join("_", "trimmed_fastq_sequence_offset", $adaptor_trim_offset), $target_name, $trimmed_fastq_sequence_length, 1, (($target_start - $adaptor_trim_offset) - 1), 1, (($target_start - $adaptor_trim_offset) - 1));
					my $trimmed_adaptor_regex = join("\t", join("_", "trimmed_adaptor_sequence_offset", $adaptor_trim_offset), $target_name, $trimmed_adaptor_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length);
					my $original_adaptor_regex = join("\t", $query_name, $target_name, $align_length, $query_start, $query_end, $target_start, $target_end);
						
					print OUTFILE $trimmed_fastq_regex . "\n";
					print OUTFILE $original_adaptor_regex . "\n";
					print OUTFILE $trimmed_adaptor_regex . "\n";
					
					my $trimmed_fastq_header = $fastq_header;
					my $trimmed_seqs_layout_header = join("_", $fasta_filename, $trimmed_fastq_header);
					print TRIMMED_LAYOUT_OUTFILE join("\t", $trimmed_seqs_layout_header, 1, (($target_start - $adaptor_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, $trimmed_adaptor_sequence_length, $target_start, $target_end, $align_length) . "\n";
					print TRIMMED_BULK_LAYOUT_OUTFILE  join("\t", $trimmed_seqs_layout_header, 1, (($target_start - $adaptor_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, $trimmed_adaptor_sequence_length, $target_start, $target_end, $align_length) . "\n";
					my $adaptor_concatenated_sequence;  
					if($query_name =~ m/([\w_]+)_([ACGT]+)/){
						my ($query_id, $adaptor_sequence) = ($1, $2);
						my $adaptor_end = length($adaptor_sequence);
						my $adaptor_sub_sequence = get_subseq($adaptor_sequence, $query_start, $query_end);
						my $adaptor_start_sequence = get_subseq($adaptor_sequence, 1, ($query_start - 1));
						my $adaptor_end_sequence = get_subseq($adaptor_sequence, ($query_end + 1), $adaptor_end);
						
						if($align_length eq $adaptor_sequence_length){
							$adaptor_concatenated_sequence  = join("-", $query_start, $adaptor_sub_sequence, $query_end);
						}elsif($adaptor_start_sequence eq ""){
							$adaptor_concatenated_sequence  = join("-", $query_start, $adaptor_sub_sequence, $query_end, $adaptor_end_sequence);
						}elsif($adaptor_end_sequence eq ""){
							$adaptor_concatenated_sequence  = join("-", $adaptor_start_sequence, $query_start, $adaptor_sub_sequence, $query_end);
						}else{
							$adaptor_concatenated_sequence  = join("-", $adaptor_start_sequence, $query_start, $adaptor_sub_sequence, $query_end, $adaptor_end_sequence);
						}
					}
					
					$adaptor_length_counter{$align_length}{$adaptor_concatenated_sequence}++;
					push(@trimmed_fastq_list, $fastq_header);
					
				}elsif($trimmed_fastq_sequence_length < $min_trimmed_fastq_sequence_length_plus_barcode){ #If this alignment doesn't meet our filtering criteria separate fastq sequence from trimmed fastq sequence data.
					$removed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
					$removed_fastq_sequences{$fastq_header}{'PLUS'} = $fastq_sequences{$fastq_header}{'PLUS'};
					$removed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
					
					my $removed_fastq_header = $fastq_header;
					my $removed_seqs_layout_header = join("_", $fasta_filename, $removed_fastq_header);
					
					print REMOVED_LAYOUT_OUTFILE join("\t", $removed_seqs_layout_header, 1, (($target_start - $adaptor_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, $trimmed_adaptor_sequence_length, $target_start, $target_end, $align_length) . "\n";
					print TRIMMED_BULK_LAYOUT_OUTFILE join("\t", $removed_seqs_layout_header, 1, (($target_start - $adaptor_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, $trimmed_adaptor_sequence_length, $target_start, $target_end, $align_length) . "\n";
					push(@removed_fastq_list, $fastq_header);
				}
			}
		}
		$i++;
  	}
  	close(OUTFILE) or die "Couldn't close file $trimmed_adaptor_regex_outfile";
  	close(INFILE) or die "Couldn't close file $adaptor_regex_outfile";
  	
  	# Grab the list of fastq sequence headers
  	my @fastq_sequence_list = keys %fastq_sequences;
  	
  	# Grab the sub list of untrimmed fastq sequences and put them into the trimmed_fastq_sequences hash.
  	my @fastq_list2remove = (@trimmed_fastq_list, @removed_fastq_list);
  	my $fastq_list_comparision = List::Compare->new(\@fastq_sequence_list, \@fastq_list2remove);
  	my @untrimmed_fastq_sequence_list = $fastq_list_comparision->get_unique;
  # 		die "unique fastq sequences: ", $fastq_list_comparision->get_unique, "\n";
  # 		die join("\t", scalar(@fastq_sequence_list), scalar(@trimmed_fastq_list), scalar(@untrimmed_fastq_sequence_list));
  	my %untrimmed_fastq_sequences = ();
  	foreach my $fastq_header (@untrimmed_fastq_sequence_list){
  		my $fastq_sequence = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
  		my $fastq_plus = $fastq_sequences{$fastq_header}{'PLUS'};
  		my $fastq_quality_scores = $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
  		
  		$untrimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $fastq_sequence;
  		$untrimmed_fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
  		$untrimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $fastq_quality_scores;
  	}
  
  	# Create output directory if it doesn't already exist.
  	my $trimmed_fastq_output_dir = join('/', $trimmed_output_dir, "TRIMMED_FASTQ_FILES");
  	unless(-d $trimmed_fastq_output_dir){
  		mkdir($trimmed_fastq_output_dir, 0777) or die "Can't make directory: $!";
  	}
  
  	# Print out trimmed fastq sequences first so that we can see what was trimmed.
  	my $trimmed_fastq_outfile = join("/", $trimmed_fastq_output_dir, join("_", $fasta_filename, "trimmed_offset", $adaptor_trim_offset) . ".fastq");
  	open(OUTFILE, ">$trimmed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for writting, $!";
  	foreach my $fastq_header (sort keys %trimmed_fastq_sequences){
  		
  		my $fastq_sequence = $trimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
  		my $fastq_plus = $trimmed_fastq_sequences{$fastq_header}{'PLUS'};
  		my $fastq_quality_scores = $trimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
  		my $fastq_sequence_length = length($fastq_sequence);
  		my ($fastq_header_prefix, $fastq_header_suffix) = split(" ", $fastq_header);
  		my $new_fastq_header = join(" ", join("_", $fastq_header_prefix, $fasta_filename, join("=", "length", $fastq_sequence_length)), $fastq_header_suffix);
  		print OUTFILE  $new_fastq_header . "\n";
  		print OUTFILE  $fastq_sequence . "\n";
  		print OUTFILE  $fastq_plus . "\n";
  		print OUTFILE  $fastq_quality_scores . "\n";
  
  		$fastq_seq_counter{$fasta_filename}{'TRIMMED'}++;
  		
  	}
  	
  	# Print out the rest of the sequences that were not trimmed because they either did not have an alignment or a significant alignment that passed the trimming threshold.
  	foreach my $fastq_header (sort keys %untrimmed_fastq_sequences){
  		
  		my $fastq_sequence = $untrimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
  		my $fastq_plus = $untrimmed_fastq_sequences{$fastq_header}{'PLUS'};
  		my $fastq_quality_scores = $untrimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
  		my $fastq_sequence_length = length($fastq_sequence);
  		my ($fastq_header_prefix, $fastq_header_suffix) = split(" ", $fastq_header);
  		my $new_fastq_header = join(" ", join("_", $fastq_header_prefix, $fasta_filename, join("=", "length", $fastq_sequence_length)), $fastq_header_suffix);
  		print OUTFILE  $new_fastq_header . "\n";
  		print OUTFILE  $fastq_sequence . "\n";
  		print OUTFILE  $fastq_plus . "\n";
  		print OUTFILE  $fastq_quality_scores . "\n";
  		
  		$fastq_seq_counter{$fasta_filename}{'UNTRIMMED'}++;
  
  	}
  	close(OUTFILE) or die "Couldn't close file $trimmed_fastq_outfile";
  	
  	# Parsed trimmed fastq file and concatenate to the bulk fastq output file for this project.
  	open(INFILE, "<$trimmed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for reading, $!";
  	while(<INFILE>){
  		chomp $_;
  		print TRIMMED_BULK_OUTFILE $_ . "\n";
  	}
  	close(INFILE) or die "Couldn't close file $trimmed_fastq_outfile";
  	
  	
  	# Create the removed fastq sequence output directory if it doesn't already exist.
  	my $removed_fastq_output_dir = join('/', $trimmed_output_dir, "REMOVED_FASTQ_SEQUENCES");
  	unless(-d $removed_fastq_output_dir){
  		mkdir($removed_fastq_output_dir, 0777) or die "Can't make directory: $!";
  	}
  	
  	# Print out trimmed fastq sequences that did not pass the filtering critera so that we can see what was trimmed and removed.
  	my $removed_fastq_outfile = join("/", $removed_fastq_output_dir, join("_", $fasta_filename, "trimmed_offset", $adaptor_trim_offset, "removed_sequences") . ".fastq");
  	open(OUTFILE, ">$removed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for writting, $!";
  	foreach my $fastq_header (sort keys %removed_fastq_sequences){
  	
  		my $fastq_sequence = $removed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
  		my $fastq_plus = $removed_fastq_sequences{$fastq_header}{'PLUS'};
  		my $fastq_quality_scores = $removed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
  		my $fastq_sequence_length = length($fastq_sequence);
  		my ($fastq_header_prefix, $fastq_header_suffix) = split(" ", $fastq_header);
  		my $new_fastq_header = join(" ", join("_", $fastq_header_prefix, $fasta_filename, join("=", "length", $fastq_sequence_length)), $fastq_header_suffix);
  		print OUTFILE  $new_fastq_header . "\n";
  		print OUTFILE  $fastq_sequence . "\n";
  		print OUTFILE  $fastq_plus . "\n";
  		print OUTFILE  $fastq_quality_scores . "\n";
  		
  		$fastq_seq_counter{$fasta_filename}{'REMOVED'}++;
  
  	}
  	close(OUTFILE) or die "Couldn't close file $removed_fastq_outfile";
  	push(@fastq_files2gzip, $trimmed_fastq_outfile);
  	
  	# Empty all hash and array containers so that we don't use fastq sequences or fastq headers from different files.
  	undef %trimmed_fastq_sequences;
  	undef %untrimmed_fastq_sequences;
  	undef %removed_fastq_sequences;
  	undef %fastq_sequences;
  	undef @trimmed_fastq_list;
  	undef @fastq_sequence_list;
  	undef @untrimmed_fastq_sequence_list;
  	undef @fastq_list2remove;
}
close(TRIMMED_BULK_OUTFILE) or die "Couldn't close file $trimmed_fastq_bulk_outfile";
close(TRIMMED_BULK_LAYOUT_OUTFILE) or die "Couldn't close file $trimmed_seqs_layout_bulk_outfile";
close(TRIMMED_LAYOUT_OUTFILE) or die "Couldn't close file $trimmed_seqs_layout_outfile";
close(REMOVED_LAYOUT_OUTFILE) or die "Couldn't close file $removed_seqs_layout_outfile";

# Generate a fastq sequence counts file so that we can see the percentages of untrimmed, trimmed, and removed fastq sequences each fastq input file.
my $fastq_seq_counts_outfile = join("/", $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adaptor_trim_offset, "fastq_seq_counts") . ".txt");
open(OUTFILE, ">$fastq_seq_counts_outfile") or die "Couldn't open file $fastq_seq_counts_outfile for writting, $!";
print OUTFILE join("\t", "fastq_file_name", "total_num_seqs", "num_untrimmed_seqs", "percent_untrimmed_seqs", "num_trimmed_seqs", "percent_trimmed_seqs", "num_removed_seqs", "percent_removed_seqs") . "\n";
foreach my $fasta_filename (sort keys %fastq_seq_counter){

    # Get the fastq untrimmed, trimmed, and removed counts.
    my $num_untrimmed_seqs = $fastq_seq_counter{$fasta_filename}{'UNTRIMMED'};
    my $num_trimmed_seqs = $fastq_seq_counter{$fasta_filename}{'TRIMMED'};
    my $num_removed_seqs = $fastq_seq_counter{$fasta_filename}{'REMOVED'};


    # If trimmed or removed sequence counters how zero count reflect in number of sequence variables.
    $num_untrimmed_seqs = 0 unless(defined($num_untrimmed_seqs));
    $num_trimmed_seqs = 0 unless(defined($num_trimmed_seqs));
    $num_removed_seqs = 0 unless(defined($num_removed_seqs));

    # Calculate the total number of sequences.
    my $total_num_seqs = ($num_untrimmed_seqs + $num_trimmed_seqs + $num_removed_seqs);

    my $percent_untrimmed_seqs = (($num_untrimmed_seqs/$total_num_seqs) * 100);
    my $percent_trimmed_seqs = (($num_trimmed_seqs/$total_num_seqs) * 100);
    my $percent_removed_seqs = (($num_removed_seqs/$total_num_seqs) * 100);

    # If trimmed or removed sequence counters how zero count reflect in either percent of sequence variables.
    $percent_untrimmed_seqs = 0.00 unless(defined($num_untrimmed_seqs));
    $percent_trimmed_seqs = 0.00 unless(defined($num_trimmed_seqs));
    $percent_removed_seqs = 0.00 unless(defined($num_removed_seqs));

    die "The number of sequences are not equal to 100.00 percent.....\n" . join("\t", $percent_untrimmed_seqs, $percent_trimmed_seqs, $percent_removed_seqs, ($percent_untrimmed_seqs + $percent_trimmed_seqs + $percent_removed_seqs)) if(($percent_untrimmed_seqs + $percent_trimmed_seqs + $percent_removed_seqs) ne 100.00);
    print OUTFILE join("\t", $fasta_filename, $total_num_seqs, $num_untrimmed_seqs, $percent_untrimmed_seqs, $num_trimmed_seqs, $percent_trimmed_seqs, $num_removed_seqs, $percent_removed_seqs) . "\n";
}
close(OUTFILE) or die "Couldn't close file $fastq_seq_counts_outfile";

# Generate the number of GBS common adaptor sequence found within the adaptor regex searches so that we can see the variation of the length of the adaptor sequences.
my $adaptor_length_counts_outfile = join("/", $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adaptor_trim_offset, "adaptor_length_counts") . ".txt");
open(OUTFILE, ">$adaptor_length_counts_outfile") or die "Couldn't open file $adaptor_length_counts_outfile for writting, $!";
print OUTFILE join("\t", "adapter_sequence_id", "adapter_length", "adaptor_sequence_count") . "\n";
foreach my $adaptor_length (sort {$b <=> $a} keys %adaptor_length_counter){
    foreach my $adaptor_concatenated_sequence (keys %{$adaptor_length_counter{$adaptor_length}}){
        my $adaptor_sequence_count = $adaptor_length_counter{$adaptor_length}{$adaptor_concatenated_sequence};
        print OUTFILE join("\t", $adaptor_concatenated_sequence, $adaptor_length, $adaptor_sequence_count) . "\n";
    }
}
close(OUTFILE) or die "Couldn't close file $adaptor_length_counts_outfile";

# Iterate through each trimmed fastq file and compress them using gzip.
foreach my $fastq_file (sort @fastq_files2gzip){
    gzip_fastq_file($fastq_file);
}
  
# # Grab all the blastdb and fasta formatted files to remove to save space because we don't need them anymore.
# my $blastdb_files = find_blastdb_files($fasta_output_dir);
# 
# # Iterate through all the blastdb and fasta formatted files and remove them.
# foreach my $blastdb_file (@{$blastdb_files}){
#     unlink($blastdb_file) or die "Could not unlink $blastdb_file: $!";
# }
# 
# # Remove the fasta file directory because we don't need it anymore.
# rmdir($fasta_output_dir) or die "Could not remove directory $fasta_output_dir: $!";

# Empty all hash and array containers for garbage collection purposes.
undef @fastq_files2gzip;
undef %fastq_seq_counter;
undef %adaptor_length_counter;

# find the fastq files found in the given fastq input file directory and return a hash list of the files and the number of files found.
# Will add a pod documentation for this sub routine.
sub find_fastq_files{
    
	my $fastq_file_dir = shift;
	die "Error lost fastq file directory" unless defined $fastq_file_dir;

	my (%fastq_files, $file_count);
	$file_count = 0;
	opendir(DIR, $fastq_file_dir) || die "Error in opening dir $fastq_file_dir\n";
	while( my $file_name = readdir(DIR)){
		my $fastq_file_name = join('/', $fastq_file_dir, $file_name) if ($file_name =~ m/\.fastq$/);
		warn "$fastq_file_name\n" if ($file_name =~ m/\.fastq$/);
		$fastq_files{$file_name} = $fastq_file_name if ($file_name =~ m/\.fastq$/);
		$file_count++ if ($file_name =~ m/\.fastq$/);
	}
	closedir(DIR);
	return (\%fastq_files, $file_count);
}

# find all the blast database related *.fasta, *.nin, *.nsq, and *.nhr files found in the given blastdb input file directory and return a hash list of the files
# Will add a pod documentation for this sub routine.
sub find_blastdb_files{
    
	my $blastdb_dir = shift;
	die "Error lost blastdb file directory" unless defined $blastdb_dir;
    
	my @blastdb_files = ();
	opendir(DIR, $blastdb_dir) || die "Error in opening dir $blastdb_dir\n";
	while( my $file_name = readdir(DIR)){
		my $blastdb_file_name = join('/', $blastdb_dir, $file_name) if ($file_name =~ m/\.fasta$/ or $file_name =~ m/\.nin$/ or $file_name =~ m/\.nsq$/ or $file_name =~ m/\.nhr$/);
		warn "$blastdb_file_name\n" if ($file_name =~ m/\.fasta$/ or $file_name =~ m/\.nin$/ or $file_name =~ m/\.nsq$/ or $file_name =~ m/\.nhr$/);
		push(@blastdb_files, $blastdb_file_name) if ($file_name =~ m/\.fasta$/ or $file_name =~ m/\.nin$/ or $file_name =~ m/\.nsq$/ or $file_name =~ m/\.nhr$/);
	}
	closedir(DIR);
	return \@blastdb_files;
}

# Compress a fastq file using gzip to save space.
# Will add a pod documentation for this sub routine.
sub gzip_fastq_file{
	
	my $fastq_file = shift;
	die "Error lost the fastq file to compress using gzip" unless defined $fastq_file ;

	warn "Calling gzip for $fastq_file....\n";
	my $gzipCmd  = "$gzip -9 $fastq_file";
	warn $gzipCmd . "\n\n";
	system($gzip, 
		'-9', $fastq_file,
	) == 0 or die "Error calling $gzipCmd: $?";
}

# Get the subsequence based on the input sequence, sequence start, and sequence end.
# Will add a pod documentation for this sub routine.
# my $seq = get_subseq("AGCTTGCGTT", 3, 8);
# warn $seq . "\n";
sub get_subseq{

        my $sequence = shift;
        die "Error lost sequence" unless defined $sequence;

        my $seq_start = shift;
        die "Error lost start of sequence" unless defined $seq_start;

        my $seq_end = shift;
        die "Error lost end of sequence" unless defined $seq_end;

        $seq_start = $seq_start - 1;
        $seq_end = $seq_end;

        my $length = ($seq_end - $seq_start);

        my $trimmed_seq = substr($sequence, $seq_start, $length);

        return uc($trimmed_seq);
}
