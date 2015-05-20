#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Parallel::Loops;
use File::Basename;
use IPC::Open2;

#### PROGRAM NAME ####
# trim_adapter_fastq_parallel_regex.pl - Program to trim the GBS common adapter sequence from each GBS fastq file within a particular Genotyping by Sequencing (GBS) project. Fixes the misprimming issue where the GBS common adapter is sequenced along with the DNA of an individual

#### DESCRIPTION ####
# This program trims the GBS common adapter sequence from each GBS fastq file within a particular Genotyping by Sequencing (GBS) project. Fixes the misprimming issue where the GBS common adapter is sequenced along with the DNA of an individual

#### SAMPLE COMMAND ####
# perl trim_adapter_fastq_parallel_regex.pl -i ~/workspace/GBS_data-08-10-2013/PstI_MspI_GBS_Data/PstI_MspI_PROCESSED_RADTAGS/PROJECT_LEADER_DIR_NO_MISMATCHES/CHRISTIANNE_MCDONALD -p CHRISTIANNE_MCDONALD -t 3 -q 32 -m 16 -l 92 -r PstI/MspI -c 16 -n true -o ~/workspace/GBS_data-08-10-2013/PstI_MspI_GBS_Data/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR
my ($fastq_file_dir, $project_name, $restriction_enzymes, $gbs_sequence_length, $adapter_length_min_threshold, $adapter_trim_offset, $min_trimmed_fastq_sequence_length, $regex_num_cpu, $pad_sequences, $output_dir);
GetOptions(
	'i=s'    => \$fastq_file_dir, # The *.fastq input file directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	'p=s'    => \$project_name, # The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	'r=s'    => \$restriction_enzymes, # The restriction enzyme(s) used to digest the genomic sequences. Default: PstI/MspI
	'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
	'm=s'    => \$adapter_length_min_threshold, # The minimum GBS common adapter sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adapter regex searches. Default: 16
	't=s'    => \$adapter_trim_offset, # The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adapter sequence found in the adapter regex searches. Default: 5
	'q=s'    => \$min_trimmed_fastq_sequence_length, # The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Default: 32
	'n=s'    => \$pad_sequences, # The padded sequence controller. Specify true for padded trimmed sequences or false for unpadded trimmed sequences. Default: false
	'c=s'    => \$regex_num_cpu, # The number of cpu cores to use for the adapter regex searches. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2
	'o=s'    => \$output_dir, # The absolute path to the output directory to contain the trimmed adapter sequence fastq output files.
);

# Display a usage message if the following parameters are not specified.
usage() unless (
	defined $fastq_file_dir
	and $project_name
	and defined $output_dir
);

# The restriction enzyme(s) used to digest the genomic sequences. Default: PstI/MspI
$restriction_enzymes = 'PstI/MspI' unless defined $restriction_enzymes;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

# The minimum GBS common adapter sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adapter regex searches. Default: 16
$adapter_length_min_threshold = 16 unless defined $adapter_length_min_threshold;

# The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adapter sequence found in the adapter regex searches. Default: 5
$adapter_trim_offset = 5 unless defined $adapter_trim_offset;

# The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Keep in mind that this is the minimum trimmed fastq sequence length used before we add the length of the barcode 
# used for splitting each individual fastq file. Tassel required sequences at least 32 base pairs (bps) plus the length of a particular barcode that can be in the range of 4-8 base pairs (bps) in length. Default: 32
$min_trimmed_fastq_sequence_length = 32  unless defined $min_trimmed_fastq_sequence_length;

# The padded sequence controller. Specify true for padded trimmed sequences or false for unpadded trimmed sequences. Default: false
$pad_sequences = 'false' unless defined $pad_sequences;

# The number of cpu cores to use for the adapter regex searches. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2
$regex_num_cpu = 2 unless defined $regex_num_cpu;

# Program dependencies - The absolute paths to gzip to compress all project leader *.fastq input files.
my $gzip 				= '/bin/gzip';

sub usage {

die <<"USAGE";

Usage: $0 -i fastq_file_dir -p project_name -r restriction_enzymes -l gbs_sequence_length -m adapter_length_min_threshold -t adapter_trim_offset -q min_trimmed_fastq_sequence_length -c regex_num_cpu -o output_dir

DESCRIPTION - This program trims the GBS common adapter sequence from each GBS fastq file within a particular Genotyping by Sequencing (GBS) project. Fixes the misprimming issue where the GBS common adapter is sequenced along with the DNA of an individual

OPTIONS:

-i fastq_file_dir - The *.fastq input file directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	e.g. /path/to/fastq_file_dir
	
-p project_name - The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	e.g. MPB-MALE-GBS
	
-r restriction_enzymes - The restriction enzyme(s) used to digest the genomic sequences. Default: PstI/MspI

-l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92

-m adapter_length_min_threshold - The minimum GBS common adapter sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adapter regex searches. Default: 16

-t adapter_trim_offset - The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adapter sequence found in the adapter regex searches. Default: 5

-q min_trimmed_fastq_sequence_length - The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Keep in mind that this is the minimum trimmed fastq sequence length used before we add the 
length of the barcode used for splitting each individual fastq file. Tassel required sequences at least 32 base pairs (bps) plus the length of a particular barcode that can be in the range of 4-8 base pairs (bps) in 
length. Default: 32

-n pad_sequences - The padded sequence controller. Specify true for padded trimmed sequences or false for unpadded trimmed sequences. Default: false

-c regex_num_cpu - The number of cpu cores to use for the adapter regex searches. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2

-o output_dir - The absolute path to the output directory to contain the trimmed adapter sequence fastq output files.

	e.g. /path/to/output_dir
USAGE
}

# Obtain the GBS common adapter sequence to trim from the GBS fastq sequences based on the restriction enzyme(s) used in the digest.
my $fastq_adapter_sequence;
if($restriction_enzymes eq 'ApeKI'){ # ApeKI
	$fastq_adapter_sequence = 'CWGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG';
}elsif($restriction_enzymes eq 'PstI/MspI'){ # PstI/MspI
	$fastq_adapter_sequence = 'CCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG';
}else{
	die "Input $restriction_enzymes is not one of the recognized restriction enzyme(s)! Please specify either ApeKI or PstI/MspI on the command line.";
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $project_dir = join('/', $output_dir, $project_name);
unless(-d $project_dir){
	mkdir($project_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $regex_output_dir = join('/', $project_dir, "ADAPTER_REGEX_FILES");
unless(-d $regex_output_dir){
	mkdir($regex_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $trimmed_output_dir = join('/', $project_dir, "TRIMMED_OUTPUT_FILES");
unless(-d $trimmed_output_dir){# Need this to make the bulk fastq sequence file.
	mkdir($trimmed_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $trimmed_regex_output_dir = join('/', $trimmed_output_dir, "TRIMMED_ADAPTER_REGEX_FILES");
unless(-d $trimmed_regex_output_dir){
	mkdir($trimmed_regex_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $trimmed_adapter_counts_output_dir = join('/', $trimmed_output_dir, "TRIMMED_ADAPTER_COUNTS_FILES");
unless(-d $trimmed_adapter_counts_output_dir){
	mkdir($trimmed_adapter_counts_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $trimmed_fastq_output_dir = join('/', $trimmed_output_dir, "TRIMMED_FASTQ_FILES");
unless(-d $trimmed_fastq_output_dir){
	mkdir($trimmed_fastq_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $trimmed_layout_output_dir = join('/', $trimmed_output_dir, "TRIMMED_LAYOUT_FILES");
unless(-d $trimmed_layout_output_dir){
	mkdir($trimmed_layout_output_dir, 0777) or die "Can't make directory: $!";
}

# Create the removed fastq sequence output directory if it doesn't already exist.
my $removed_fastq_output_dir = join('/', $trimmed_output_dir, "REMOVED_FASTQ_SEQUENCES");
unless(-d $removed_fastq_output_dir){
	mkdir($removed_fastq_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $removed_layout_output_dir = join('/', $trimmed_output_dir, "REMOVED_LAYOUT_FILES");
unless(-d $removed_layout_output_dir){
	mkdir($removed_layout_output_dir, 0777) or die "Can't make directory: $!";
}

# Get the length of the adapter so that we can verify adapter blast results.
my $adapter_sequence_length = length($fastq_adapter_sequence);

# Generate the adapter concatenated sequence positions string for annotation of the adapter counts file 
my %adapter_concatenated_sequences = ();
my $query_start = 1;
for(my $align_length = $adapter_sequence_length; $align_length >= $adapter_length_min_threshold; $align_length--){
	my $adapter_end = $adapter_sequence_length;
	my $query_end = $align_length;
	my $adapter_sub_sequence = get_subseq($fastq_adapter_sequence, $query_start, $query_end);
	my $adapter_start_sequence = get_subseq($fastq_adapter_sequence, 1, ($query_start - 1));
	my $adapter_end_sequence = get_subseq($fastq_adapter_sequence, ($query_end + 1), $adapter_end);

	my $adapter_concatenated_sequence;
	if($align_length eq $adapter_sequence_length){
		$adapter_concatenated_sequence  = join("-", $query_start, $adapter_sub_sequence, $query_end);
	}elsif($adapter_start_sequence eq ""){
		$adapter_concatenated_sequence  = join("-", $query_start, $adapter_sub_sequence, $query_end, $adapter_end_sequence);
	}elsif($adapter_end_sequence eq ""){
		$adapter_concatenated_sequence  = join("-", $adapter_start_sequence, $query_start, $adapter_sub_sequence, $query_end);
	}else{
		$adapter_concatenated_sequence  = join("-", $adapter_start_sequence, $query_start, $adapter_sub_sequence, $query_end, $adapter_end_sequence);
	}
	$adapter_concatenated_sequences{$align_length} = $adapter_concatenated_sequence;
}

my (%original_fastq_sequence_counter, %trimmed_fastq_sequence_counter) = ();
if ((require Parallel::Loops) and ($regex_num_cpu)){

        # Perform the adapter regex searches in parallel.
        my $parallel = Parallel::Loops->new($regex_num_cpu);
	my %original_fastq_seq_counter = ();
        my %trimmed_fastq_seq_counter = ();
	my %adapter_length_counter = ();
        $parallel->share(\%original_fastq_seq_counter, \%trimmed_fastq_seq_counter); # make sure that these are visible in the children.

        # Find all files in the specified directory with the extension *.fastq.
	my ($fastq_files, $fastq_file_count) = find_files($fastq_file_dir, "fastq");

	# Iterate through the files with the extension *.fastq.
	my @jobs = sort {$a cmp $b} values %{$fastq_files};
	$parallel->foreach(\@jobs,sub {

		# Get the full path to the GBS common adapter length counts file.
		my $fastq_infile = $_;
		
		warn "Processing " . $fastq_infile . ".....\n";
		
		# Get the basename of the fastq filename without the .fastq extension.
		my $fasta_filename = fileparse($fastq_infile, qr/\.fastq/);
		
		# Open the adapter regex output file.
		my $adapter_regex_outfile = join('/', $regex_output_dir, join("_", $fasta_filename, "gbs_adapter_regex.tsv"));
		open(ADAPTER_REGEX_OUTFILE, ">$adapter_regex_outfile") or die "Couldn't open file $adapter_regex_outfile for writting, $!";
		print ADAPTER_REGEX_OUTFILE join("\t", "query_name", "target_name", "align_length", "query_start", "query_end", "target_start", "target_end") . "\n";
		
		# Generate the trimmed adapter regex files filtered for further processing.
		my $trimmed_adapter_regex_outfile = join('/', $trimmed_regex_output_dir, join("_", $fasta_filename, "gbs_adapter_regex.tsv"));
		open(TRIMMED_ADAPTER_REGEX_OUTFILE, ">$trimmed_adapter_regex_outfile") or die "Couldn't open file $trimmed_adapter_regex_outfile for writting, $!";
		print TRIMMED_ADAPTER_REGEX_OUTFILE join("\t", "query_name", "target_name", "align_length", "query_start", "query_end", "target_start", "target_end") . "\n";
		
		# The $trimmed_seqs_layout_outfile contains the trimmed coordinates layout for each sequence trimmed of the GBS common adapter sequence that passed the retaining criteria.
		my $trimmed_seqs_layout_outfile = join('/', $trimmed_layout_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset, "trimmed_seqs_layout") . ".txt");
		open(TRIMMED_LAYOUT_OUTFILE, ">$trimmed_seqs_layout_outfile") or die "Couldn't open file $trimmed_seqs_layout_outfile for writting, $!";
		print TRIMMED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adapter_start", "trimmed_adapter_end",
		"trimmed_adapter_length", "adapter_seq_start", "adapter_seq_end", "adapter_seq_length") . "\n";

		# Print out trimmed fastq sequences first so that we can see what was trimmed.
		my $trimmed_fastq_outfile = join("/", $trimmed_fastq_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset) . ".fastq");
		open(TRIMMED_FASTQ_OUTFILE, ">$trimmed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for writting, $!";
		
		# The $removed_seqs_layout_outfile contains the trimmed coordinates layout for each sequence trimmed of the GBS common adapter sequence that failed the retaining criteria and was therefore removed.
		my $removed_seqs_layout_outfile = join('/', $removed_layout_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset, "removed_seqs_layout") . ".txt");
		open(REMOVED_LAYOUT_OUTFILE, ">$removed_seqs_layout_outfile") or die "Couldn't open file $removed_seqs_layout_outfile for writting, $!";
		print REMOVED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adapter_start", "trimmed_adapter_end",
		"trimmed_adapter_length", "adapter_seq_start", "adapter_seq_end", "adapter_seq_length") . "\n";
		
		# Print out trimmed fastq sequences that did not pass the filtering critera so that we can see what was trimmed and removed.
		my $removed_fastq_outfile = join("/", $removed_fastq_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset, "removed_sequences") . ".fastq");
		open(REMOVED_FASTQ_OUTFILE, ">$removed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for writting, $!";
		
		# Open the fastq input file for parsing.
		open(FASTQ_INFILE, "<$fastq_infile") or die "Couldn't open file $fastq_infile for reading, $!";
		# Parse the fastq files for the fastq header, sequence, plus, and quality scores and reformat to *.fasta format so that we can use adapter regex.
		my ($fastq_header, $fastq_sequence, $fastq_plus, $fastq_quality_scores);
		my $i = 1;
		my $num_fastq_seqs = 0;
		while(<FASTQ_INFILE>){
			chomp $_;
			#warn $_ . "\n";
			if(($_ =~ m/^\@[A-Za-z0-9-_]+:\d+:[A-Za-z0-9]+:\d+:\d+:\d+:\d+ \d:[A-Z]:\d:[ACGTRYKMSWBDHVN]*$/)
				or ($_ =~ m/^\@[A-Za-z0-9-_]+\_[A-Za-z0-9-_]+\_[A-Za-z0-9-_]+\_[A-Za-z0-9-_]+\_[A-Za-z0-9-_]+$/) 
				and ($i eq 1)){ # The fastq sequence header is on the first line. i.e. @HWI-ST767:215:C30VBACXX:8:1101:1801:1484 1:N:0:
				$fastq_header = $_;
				
			}elsif(($_ =~ m/^[ACGTRYKMSWBDHVN]+$/i) and ($i eq 2)){ # The fastq sequence is on the second line.
				$fastq_sequence = $_;
				
			}elsif(($_ =~ m/^\+$/) and ($i eq 3)){ # The fastq plus character is on the third line.
				$fastq_plus = $_;
				
			}elsif(($_ =~ m/^.+$/) and (($i % 4) eq 0)){ # the fastq quality scores are on the fourth line.
				$fastq_quality_scores = $_;
				
			}
			
			if(($i % 4) eq 0){ # Once we are finished parsing a fastq sequence entry we check to make sure it was parsed correctly and that the lengths of the sequence and quality scores match. 
				
				die "Error: fastq_header is undefined" unless(defined($fastq_header));
				die "Error: fastq_sequence is undefined" unless(defined($fastq_sequence));
				die "Error: fastq_plus is undefined" unless(defined($fastq_plus));
				die "Error: fastq_quality_scores is undefined" unless(defined($fastq_quality_scores));
				
				my $fastq_sequence_length = length($fastq_sequence);
				my $fastq_quality_scores_length = length($fastq_quality_scores);
				
				die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($fastq_sequence_length ne $gbs_sequence_length);
				die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length ne fastq_quality_scores_length=$fastq_quality_scores_length" if($fastq_sequence_length ne $fastq_quality_scores_length);
				
				# Execute the adapter regex search using the GBS common adapter sequence as the query and the fastq sequences as the target.
				# Keep a maximum of target sequences equal to the number of sequences in the fastq file. Retain adapter regex hits based on the minimum adapter length threshold.
				# Keep adapter regex hits that begin at position 1 of the full GBS common adapter sequence up to the minimum adapter length threshold.
				my ($align_length, $query_start, $query_end, $target_start, $target_end, $regex_alignment);
				my $adapter_length_count = 1;
				my $common_adapter_sequence = $fastq_adapter_sequence;
				my $alignment_found = "false";
				for(my $i = $adapter_sequence_length; $i >= $adapter_length_min_threshold; $i--){

# 					warn "$i eq $adapter_sequence_length\n";
					if($i eq $adapter_sequence_length){
						my $common_adapter_regex;
						if($restriction_enzymes eq 'ApeKI'){
							my $common_adapter_sequence_ApeKI = $common_adapter_sequence;
							$common_adapter_sequence_ApeKI =~ s/W/\[AT\]/g;
							$common_adapter_regex = qr($common_adapter_sequence_ApeKI);
						}elsif($restriction_enzymes eq 'PstI/MspI'){
							my $common_adapter_sequence_Pstl_MspI = $common_adapter_sequence;
							$common_adapter_regex = qr($common_adapter_sequence_Pstl_MspI);
						}
						
						if($fastq_sequence =~ /$common_adapter_regex/g){
							$target_start = ($-[0] + 1);
							$target_end = $+[0];
							$align_length = ($+[0] - $-[0]);
# 							warn join("\t", $align_length, $query_start, $query_end, $target_start, $target_end) . "\n";
							
							$query_start = 1;
							$query_end = $adapter_sequence_length;
# 							warn join("\t", $adapter_sequence_length, $common_adapter_sequence) . "\n";
							$regex_alignment = join("\t", join("_", "GBS_adapter_sequence", $fastq_adapter_sequence), $fastq_header, $align_length, $query_start, $query_end, $target_start, $target_end);
							print ADAPTER_REGEX_OUTFILE $regex_alignment . "\n";
							$alignment_found = "true";
							last;
						}
						
					}elsif($i < $adapter_sequence_length){
						my $common_adapter_regex;
						if($restriction_enzymes eq 'ApeKI'){
							my $common_adapter_sequence_ApeKI = $common_adapter_sequence;
							$common_adapter_sequence_ApeKI =~ s/W/\[AT\]/g;
							$common_adapter_regex = qr($common_adapter_sequence_ApeKI);
						}elsif($restriction_enzymes eq 'PstI/MspI'){
							my $common_adapter_sequence_Pstl_MspI = $common_adapter_sequence;
							$common_adapter_regex = qr($common_adapter_sequence_Pstl_MspI);
						}
						
						while($fastq_sequence =~ /$common_adapter_regex/g){
							if($+[0] eq $gbs_sequence_length){
								$target_start = ($-[0] + 1);
								$target_end = $+[0];
								$align_length = ($+[0] - $-[0]);
# 								warn join("\t", $align_length, $query_start, $query_end, $target_start, $target_end) . "\n";
								
								$regex_alignment = join("\t", join("_", "GBS_adapter_sequence", $fastq_adapter_sequence), $fastq_header, $align_length, $query_start, $query_end, $target_start, $target_end);
								print ADAPTER_REGEX_OUTFILE $regex_alignment . "\n";
								$alignment_found = "true";
								last;
							}
						}
						last if($alignment_found eq "true");
					}
					
					$query_start = 1;
					$query_end = ($adapter_sequence_length - $adapter_length_count);
					
					$common_adapter_sequence = get_subseq($common_adapter_sequence, $query_start, $query_end);
# 					warn join("\t", ($adapter_sequence_length - $adapter_length_count), $common_adapter_sequence) . "\n";
					$adapter_length_count++;
				}
				
				if($alignment_found eq "false"){ # If a regex alignment was not found.
				
					# Print out the rest of the sequences that were not trimmed because they either did not have an alignment or a significant alignment that passed the trimming threshold.
					my $new_fastq_header = join("\001", $fastq_header, $fasta_filename, join("=", "length", $fastq_sequence_length));
					print TRIMMED_FASTQ_OUTFILE  $new_fastq_header . "\n";
					print TRIMMED_FASTQ_OUTFILE  $fastq_sequence . "\n";
					print TRIMMED_FASTQ_OUTFILE  $fastq_plus . "\n";
					print TRIMMED_FASTQ_OUTFILE  $fastq_quality_scores . "\n";
					
					$trimmed_fastq_seq_counter{$fasta_filename}{'UNTRIMMED'}++;
					
				}elsif($alignment_found eq "true"){ # If a regex alignment was found.
				
					# Parse the tab-delimited adapter regex output so that we can trimm the fastq sequences that contain the GBS adapter sequence.
					# Create new adapter regex files so that we can visualize where we trimmed the sequence.
					my @split_adapter_regex_hit =  split(/\t/, $regex_alignment);
					my ($query_name, $target_name, $align_length, $query_start, $query_end, $target_start, $target_end) = @split_adapter_regex_hit;
				
					die "Error: query_name is undefined" unless(defined($query_name));
					die "Error: target_name is undefined" unless(defined($target_name));
					die "Error: align_length is undefined" unless(defined($align_length));
					die "Error: query_start is undefined" unless(defined($query_start));
					die "Error: query_end is undefined" unless(defined($query_end));
					die "Error: target_start is undefined" unless(defined($target_start));
					die "Error: target_end is undefined" unless(defined($target_end));
					
					die "Error: fastq_sequence is undefined" unless(defined($fastq_sequence));
					die "Error: fastq_plus is undefined" unless(defined($fastq_plus));
					die "Error: fastq_quality_scores is undefined" unless(defined($fastq_quality_scores));
					
					my $fastq_sequence_length = length($fastq_sequence);
					
					my ($trimmed_fastq_sequence, $trimmed_fastq_quality_scores, $trimmed_adapter_sequence);
					# If the adapter alignment length is equal to the length of the full GBS common adapter then trim where that GBS common adapter sequence is found subtracting the trimmed offset from the target start.
					# If the adapter alignment length is less than the length of the full GBS common adapter then trim only if the target end of the aligned adapter sequence is equal the the common GBS sequence length.
					if(($align_length eq $adapter_sequence_length) 
						or ((($align_length >= $adapter_length_min_threshold) and ($align_length < $adapter_sequence_length)) and ($query_start eq 1) and ($target_end eq $gbs_sequence_length))){
						
						# Get the trimmed fastq sequence trimmed of the GBS common adapter sequence and trimmed offset.
						$trimmed_fastq_sequence = get_subseq($fastq_sequence, 1, (($target_start - $adapter_trim_offset) - 1));
						# Get the trimmed fastq quality scores trimmed of the GBS common adapter sequence and trimmed offset.
						$trimmed_fastq_quality_scores = get_subseq($fastq_quality_scores, 1, (($target_start - $adapter_trim_offset) - 1));
						# Get the trimmed adapter sequence and trimmed offset that was trimmed off the GBS fastq sequence.
						$trimmed_adapter_sequence = get_subseq($fastq_sequence, ($target_start - $adapter_trim_offset), $fastq_sequence_length);
						
						# Get the length of the trimmed GBS common adapter sequence and trimmed offset.
						my $trimmed_adapter_sequence_length = length($trimmed_adapter_sequence);
						my $trimmed_fastq_sequence_length = length($trimmed_fastq_sequence);
						my $trimmed_fastq_quality_scores_length = length($trimmed_fastq_quality_scores);
						
						# If the trimmed sequence length is greater than or equal to the minimum trimmed fastq sequence length plus the length of the barcode then the trimmed and sequence count files.
						if($trimmed_fastq_sequence_length >= $min_trimmed_fastq_sequence_length){
							
							my $trimmed_fastq_regex = join("\t", join("_", "trimmed_fastq_sequence_offset", $adapter_trim_offset), $target_name, $trimmed_fastq_sequence_length, 1, (($target_start - $adapter_trim_offset) - 1), 1, (($target_start - $adapter_trim_offset) - 1));
							my $trimmed_adapter_regex = join("\t", join("_", "trimmed_adapter_sequence_offset", $adapter_trim_offset), $target_name, $trimmed_adapter_sequence_length, ($target_start - $adapter_trim_offset), $fastq_sequence_length, ($target_start - $adapter_trim_offset), $fastq_sequence_length);
							my $original_adapter_regex = join("\t", $query_name, $target_name, $align_length, $query_start, $query_end, $target_start, $target_end);
								
							print TRIMMED_ADAPTER_REGEX_OUTFILE $trimmed_fastq_regex . "\n";
							print TRIMMED_ADAPTER_REGEX_OUTFILE $original_adapter_regex . "\n";
							print TRIMMED_ADAPTER_REGEX_OUTFILE $trimmed_adapter_regex . "\n";
							
							print TRIMMED_LAYOUT_OUTFILE join("\t", join("_", $fasta_filename, $fastq_header), 1, (($target_start - $adapter_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adapter_trim_offset), $fastq_sequence_length, $trimmed_adapter_sequence_length, $target_start, $target_end, $align_length) . "\n";
							
							my $new_fastq_header = join("\001", $fastq_header, $fasta_filename, join("=", "length", $trimmed_fastq_sequence_length));
							if($pad_sequences eq "true"){ # Pad sequences with poly-Ns up to the common GBS sequence length.
								my $padded_N_length =  ($gbs_sequence_length - $trimmed_fastq_sequence_length);
								my $padded_N_seq = 'N' x $padded_N_length;
								my $padded_fastq_sequence = join("", $trimmed_fastq_sequence, $padded_N_seq);
								
								my $padded_score_length =  ($gbs_sequence_length - $trimmed_fastq_quality_scores_length);
								my $padded_score_seq = '#' x $padded_score_length;
								my $padded_fastq_quality_scores = join("", $trimmed_fastq_quality_scores, $padded_score_seq);
								
								my $padded_fastq_sequence_length = length($padded_fastq_sequence);
								my $padded_fastq_quality_scores_length = length($padded_fastq_quality_scores);
								
								die "Error: $new_fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($padded_fastq_sequence_length ne $gbs_sequence_length);
								die "Error: $new_fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length ne padded_fastq_quality_scores_length=$padded_fastq_quality_scores_length" if($padded_fastq_sequence_length ne $padded_fastq_quality_scores_length);
					
								print TRIMMED_FASTQ_OUTFILE  $new_fastq_header . "\n";
								print TRIMMED_FASTQ_OUTFILE  $padded_fastq_sequence . "\n";
								print TRIMMED_FASTQ_OUTFILE  $fastq_plus . "\n";
								print TRIMMED_FASTQ_OUTFILE  $padded_fastq_quality_scores . "\n";
							}else{
							
								print TRIMMED_FASTQ_OUTFILE  $new_fastq_header . "\n";
								print TRIMMED_FASTQ_OUTFILE  $trimmed_fastq_sequence . "\n";
								print TRIMMED_FASTQ_OUTFILE  $fastq_plus . "\n";
								print TRIMMED_FASTQ_OUTFILE  $trimmed_fastq_quality_scores . "\n";
							
							}
							
							$trimmed_fastq_seq_counter{$fasta_filename}{'TRIMMED'}++;
							
							$adapter_length_counter{$align_length}++;
							
						}elsif($trimmed_fastq_sequence_length < $min_trimmed_fastq_sequence_length){ #If this alignment doesn't meet our filtering criteria separate fastq sequence from trimmed fastq sequence data.
							
							print REMOVED_LAYOUT_OUTFILE join("\t", join("_", $fasta_filename, $fastq_header), 1, (($target_start - $adapter_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adapter_trim_offset), $fastq_sequence_length, $trimmed_adapter_sequence_length, $target_start, $target_end, $align_length) . "\n";

							my $new_fastq_header = join("\001", $fastq_header, $fasta_filename, join("=", "length", $trimmed_fastq_sequence_length));
							if($pad_sequences eq "true"){
								my $padded_N_length =  ($gbs_sequence_length - $trimmed_fastq_sequence_length);
								my $padded_N_seq = 'N' x $padded_N_length;
								my $padded_fastq_sequence = join("", $trimmed_fastq_sequence, $padded_N_seq);
								
								my $padded_score_length =  ($gbs_sequence_length - $trimmed_fastq_quality_scores_length);
								my $padded_score_seq = '#' x $padded_score_length;
								my $padded_fastq_quality_scores = join("", $trimmed_fastq_quality_scores, $padded_score_seq);
								
								my $padded_fastq_sequence_length = length($padded_fastq_sequence);
								my $padded_fastq_quality_scores_length = length($padded_fastq_quality_scores);
								
								die "Error: $new_fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($padded_fastq_sequence_length ne $gbs_sequence_length);
								die "Error: $new_fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length ne padded_fastq_quality_scores_length=$padded_fastq_quality_scores_length" if($padded_fastq_sequence_length ne $padded_fastq_quality_scores_length);
					
								print REMOVED_FASTQ_OUTFILE  $new_fastq_header . "\n";
								print REMOVED_FASTQ_OUTFILE  $padded_fastq_sequence . "\n";
								print REMOVED_FASTQ_OUTFILE  $fastq_plus . "\n";
								print REMOVED_FASTQ_OUTFILE  $padded_fastq_quality_scores . "\n";
							}else{
							
								print REMOVED_FASTQ_OUTFILE  $new_fastq_header . "\n";
								print REMOVED_FASTQ_OUTFILE  $trimmed_fastq_sequence . "\n";
								print REMOVED_FASTQ_OUTFILE  $fastq_plus . "\n";
								print REMOVED_FASTQ_OUTFILE  $trimmed_fastq_quality_scores . "\n";
							
							}
							
							$trimmed_fastq_seq_counter{$fasta_filename}{'REMOVED'}++;
						}
					}
				}
				$i = 1;
				$num_fastq_seqs++;
				
			}else{
				$i++;
			}
		}
		close(FASTQ_INFILE) or die "Couldn't close file $fastq_infile";
		close(ADAPTER_REGEX_OUTFILE) or die "Couldn't close file $adapter_regex_outfile";
		close(TRIMMED_ADAPTER_REGEX_OUTFILE) or die "Couldn't close file $trimmed_adapter_regex_outfile";
		close(TRIMMED_LAYOUT_OUTFILE) or die "Couldn't close file $trimmed_seqs_layout_outfile";
		close(TRIMMED_FASTQ_OUTFILE) or die "Couldn't close file $trimmed_fastq_outfile";
		close(REMOVED_LAYOUT_OUTFILE) or die "Couldn't close file $removed_seqs_layout_outfile";
		close(REMOVED_FASTQ_OUTFILE) or die "Couldn't close file $removed_fastq_outfile";
		
		# Get original fastq sequences counts for each file.
		$original_fastq_seq_counter{$fasta_filename} = $num_fastq_seqs;
		
		# Print the GBS common adapter length counts.
		my $trimmed_adapter_counts_outfile = join("/", $trimmed_adapter_counts_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset, "adapter_length_counts") . ".txt");
		open(TRIMMED_ADAPTER_COUNTS_OUTFILE, ">$trimmed_adapter_counts_outfile") or die "Couldn't open file $trimmed_adapter_counts_outfile for writting, $!";
		print TRIMMED_ADAPTER_COUNTS_OUTFILE join("\t", "adapter_sequence_id", "adapter_length", "adapter_sequence_count") . "\n";
		foreach my $adapter_length (sort {$b <=> $a} keys %adapter_length_counter){
			my $adapter_sequence_count = $adapter_length_counter{$adapter_length};
			my $adapter_concatenated_sequence = $adapter_concatenated_sequences{$adapter_length};
			print TRIMMED_ADAPTER_COUNTS_OUTFILE join("\t", $adapter_concatenated_sequence, $adapter_length, $adapter_sequence_count) . "\n";
		}
		close(TRIMMED_ADAPTER_COUNTS_OUTFILE) or die "Couldn't close file $trimmed_adapter_counts_outfile";

	});

	# Get the fastq untrimmed, trimmed, and removed counts.
	foreach my $fasta_filename (sort keys %trimmed_fastq_seq_counter){
		$trimmed_fastq_sequence_counter{$fasta_filename}{'UNTRIMMED'} = $trimmed_fastq_seq_counter{$fasta_filename}{'UNTRIMMED'};
		$trimmed_fastq_sequence_counter{$fasta_filename}{'TRIMMED'} = $trimmed_fastq_seq_counter{$fasta_filename}{'TRIMMED'};
		$trimmed_fastq_sequence_counter{$fasta_filename}{'REMOVED'} = $trimmed_fastq_seq_counter{$fasta_filename}{'REMOVED'};
	}
	
	# Get the original fastq sequence counts for each fastq file so that we can use them to make sure we have the correct proportions of untrimmed, trimmed, and removed counts.
	foreach my $fasta_filename (sort keys %original_fastq_seq_counter){
		$original_fastq_sequence_counter{$fasta_filename} = $original_fastq_seq_counter{$fasta_filename};
	}
	
	# close parallel loops to free the memory contained in the shared hash variables.
	undef $parallel;
}

# Generate bulk trimmed fastq file for this Genotyping by Sequencing (GBS) project. 
# $trimmed_fastq_bulk_outfile contains the trimmed GBS sequences in fastq format.
# Find all files in the specified directory with the extension *.fastq.
my ($trimmed_fastq_files, $trimmed_fastq_file_count) = find_files($trimmed_fastq_output_dir, "fastq");
# my $trimmed_fastq_bulk_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adapter_trim_offset) . ".fastq");
# open(TRIMMED_BULK_OUTFILE, ">$trimmed_fastq_bulk_outfile") or die "Couldn't open file $trimmed_fastq_bulk_outfile for writting, $!";
# Iterate through the files with extension *.fastq.
foreach my $trimmed_fastq_filename (sort keys %{$trimmed_fastq_files}){
	
	# Get the full path to the trimmed fastq file.
	my $trimmed_fastq_infile = $trimmed_fastq_files->{$trimmed_fastq_filename};
	
# 	# Parse trimmed fastq file and concatenate to the bulk fastq output file for this project.
# 	open(INFILE, "<$trimmed_fastq_infile") or die "Couldn't open file $trimmed_fastq_infile for reading, $!";
# 	while(<INFILE>){
# 		chomp $_;
# 		print TRIMMED_BULK_OUTFILE $_ . "\n";
# 	}
# 	close(INFILE) or die "Couldn't close file $trimmed_fastq_infile";
	
	# Compress the trimmed fastq file using gzip.
	gzip_file($trimmed_fastq_infile);
}
# close(TRIMMED_BULK_OUTFILE) or die "Couldn't close file $trimmed_fastq_bulk_outfile";

# Compress the bulk trimmed fastq file using gzip.
# gzip_file($trimmed_fastq_bulk_outfile);

# The $trimmed_seqs_layout_bulk_outfile contains all the trimmed coordinates layout for each sequence trimmed of the GBS common adapter sequence.
my $trimmed_seqs_layout_bulk_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adapter_trim_offset, "all_trimmed_seqs_layout") . ".txt");
open(TRIMMED_BULK_LAYOUT_OUTFILE, ">$trimmed_seqs_layout_bulk_outfile") or die "Couldn't open file $trimmed_seqs_layout_bulk_outfile for writting, $!";
print TRIMMED_BULK_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adapter_start", "trimmed_adapter_end",
"trimmed_adapter_length", "adapter_seq_start", "adapter_seq_end", "adapter_seq_length") . "\n";

# The $trimmed_seqs_layout_outfile contains the trimmed coordinates layout for each sequence trimmed of the GBS common adapter sequence that passed the retaining criteria.
# Find all files in the specified directory with the extension *.txt.
my ($trimmed_seqs_layout_files, $trimmed_seqs_layout_file_count) = find_files($trimmed_layout_output_dir, "txt");
my $trimmed_seqs_layout_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adapter_trim_offset, "trimmed_seqs_layout") . ".txt");
open(TRIMMED_LAYOUT_OUTFILE, ">$trimmed_seqs_layout_outfile") or die "Couldn't open file $trimmed_seqs_layout_outfile for writting, $!";
print TRIMMED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adapter_start", "trimmed_adapter_end",
"trimmed_adapter_length", "adapter_seq_start", "adapter_seq_end", "adapter_seq_length") . "\n";
foreach my $trimmed_seqs_layout_filename (sort keys %{$trimmed_seqs_layout_files}){

	# Get the full path to the trimmed sequence layout file.
	my $trimmed_seqs_layout_infile = $trimmed_seqs_layout_files->{$trimmed_seqs_layout_filename};
	
	# Parse trimmed sequence layout file and concatenate to the bulk sequence layout output file for this project.
	open(INFILE, "<$trimmed_seqs_layout_infile") or die "Couldn't open file $trimmed_seqs_layout_infile for reading, $!";
	my $i = 0;
	while(<INFILE>){
		chomp $_;
		if($i ne 0){
			print TRIMMED_LAYOUT_OUTFILE $_ . "\n";
			print TRIMMED_BULK_LAYOUT_OUTFILE $_ . "\n";
		}
		$i++;
	}
	close(INFILE) or die "Couldn't close file $trimmed_seqs_layout_infile";
	
	# Compress the trimmed sequence layout file using gzip.
	gzip_file($trimmed_seqs_layout_infile);
}
close(TRIMMED_LAYOUT_OUTFILE) or die "Couldn't close file $trimmed_seqs_layout_outfile";

# The $removed_seqs_layout_outfile contains the trimmed coordinates layout for each sequence trimmed of the GBS common adapter sequence that failed the retaining criteria and was therefore removed.
# Find all files in the specified directory with the extension *.txt.
my ($removed_seqs_layout_files, $removed_seqs_layout_file_count) = find_files($removed_layout_output_dir, "txt");
my $removed_seqs_layout_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adapter_trim_offset, "removed_seqs_layout") . ".txt");
open(REMOVED_LAYOUT_OUTFILE, ">$removed_seqs_layout_outfile") or die "Couldn't open file $removed_seqs_layout_outfile for writting, $!";
print REMOVED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adapter_start", "trimmed_adapter_end",
"trimmed_adapter_length", "adapter_seq_start", "adapter_seq_end", "adapter_seq_length") . "\n";
foreach my $removed_seqs_layout_filename (sort keys %{$removed_seqs_layout_files}){

	# Get the full path to the removed sequence layout file.
	my $removed_seqs_layout_infile = $removed_seqs_layout_files->{$removed_seqs_layout_filename};
	
	# Parse removed sequence layout file and concatenate to the bulk sequence layout output file for this project.
	open(INFILE, "<$removed_seqs_layout_infile") or die "Couldn't open file $removed_seqs_layout_infile for reading, $!";
	my $i = 0;
	while(<INFILE>){
		chomp $_;
		if($i ne 0){
			print REMOVED_LAYOUT_OUTFILE $_ . "\n";
			print TRIMMED_BULK_LAYOUT_OUTFILE $_ . "\n";
		}
		$i++;
	}
	close(INFILE) or die "Couldn't close file $removed_seqs_layout_infile";
	
	# Compress the trimmed sequence layout file using gzip.
	gzip_file($removed_seqs_layout_infile);
}
close(REMOVED_LAYOUT_OUTFILE) or die "Couldn't close file $removed_seqs_layout_outfile";
close(TRIMMED_BULK_LAYOUT_OUTFILE) or die "Couldn't close file $trimmed_seqs_layout_bulk_outfile";

# Compress the trimmed sequence layout file using gzip.
gzip_file($trimmed_seqs_layout_outfile);

# Compress the trimmed sequence layout file using gzip.
gzip_file($removed_seqs_layout_outfile);

# Compress the bulk trimmed sequence layout file using gzip.
gzip_file($trimmed_seqs_layout_bulk_outfile);

# Generate a fastq sequence counts file so that we can see the percentages of untrimmed, trimmed, and removed fastq sequences each fastq input file.
my $fastq_seq_counts_outfile = join("/", $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adapter_trim_offset, "fastq_seq_counts") . ".txt");
open(OUTFILE, ">$fastq_seq_counts_outfile") or die "Couldn't open file $fastq_seq_counts_outfile for writting, $!";
print OUTFILE join("\t", "fastq_file_name", "total_num_seqs", "num_untrimmed_seqs", "percent_untrimmed_seqs", "num_trimmed_seqs", "percent_trimmed_seqs", "num_removed_seqs", "percent_removed_seqs") . "\n";
foreach my $fasta_filename (sort keys %trimmed_fastq_sequence_counter){

    # Get the fastq untrimmed, trimmed, and removed counts.
    my $num_untrimmed_seqs = $trimmed_fastq_sequence_counter{$fasta_filename}{'UNTRIMMED'};
    my $num_trimmed_seqs = $trimmed_fastq_sequence_counter{$fasta_filename}{'TRIMMED'};
    my $num_removed_seqs = $trimmed_fastq_sequence_counter{$fasta_filename}{'REMOVED'};


    # If trimmed or removed sequence counters how zero count reflect in number of sequence variables.
    $num_untrimmed_seqs = 0 unless(defined($num_untrimmed_seqs));
    $num_trimmed_seqs = 0 unless(defined($num_trimmed_seqs));
    $num_removed_seqs = 0 unless(defined($num_removed_seqs));

    # Calculate the total number of sequences.
    my $total_num_seqs = $original_fastq_sequence_counter{$fasta_filename};
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

# Get the GBS common adapter length counts for each fastq file.
my ($trimmed_adapter_counts_files, $trimmed_adapter_counts_file_count) = find_files($trimmed_adapter_counts_output_dir, "txt");
my %adapter_length_counts = ();
foreach my $trimmed_adapter_counts_filename (sort keys %{$trimmed_adapter_counts_files}){

	# Get the full path to the GBS common adapter length counts file.
	my $trimmed_adapter_counts_infile = $trimmed_adapter_counts_files->{$trimmed_adapter_counts_filename};
	
	# Parse the GBS common adapter length counts file.
	open(INFILE, "<$trimmed_adapter_counts_infile") or die "Couldn't open file $trimmed_adapter_counts_infile for reading, $!";
	my $i = 0;
	while(<INFILE>){
		chomp $_;
		if($i ne 0){
			my @split_adapter_count_entries = split(/\t/, $_);
			
			my ($adapter_concatenated_sequence, $adapter_length, $adapter_sequence_count) = @split_adapter_count_entries;
			push(@{$adapter_length_counts{$adapter_length}}, $adapter_sequence_count);
		}
		$i++;
	}
	close(INFILE) or die "Couldn't close file $trimmed_adapter_counts_infile";
}

# Get the sum of adapter length counts obtained from each file.
my %adapter_regex_length_counter = ();
foreach my $adapter_length (sort {$b <=> $a} keys %adapter_length_counts){
	my $adaptor_length_counter = 0;
	foreach my $adapter_length_count (@{$adapter_length_counts{$adapter_length}}){
		$adaptor_length_counter += $adapter_length_count;
	}
	$adapter_regex_length_counter{$adapter_length} = $adaptor_length_counter;
}

# Generate the number of GBS common adapter sequence found within the adapter regex searches so that we can see the variation of the length of the adapter sequences.
my $adapter_length_counts_outfile = join("/", $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adapter_trim_offset, "adapter_length_counts") . ".txt");
open(OUTFILE, ">$adapter_length_counts_outfile") or die "Couldn't open file $adapter_length_counts_outfile for writting, $!";
print OUTFILE join("\t", "adapter_sequence_id", "adapter_length", "adapter_sequence_count") . "\n";
foreach my $adapter_length (sort {$b <=> $a} keys %adapter_regex_length_counter){
	my $adapter_sequence_count = $adapter_regex_length_counter{$adapter_length};
	my $adapter_concatenated_sequence = $adapter_concatenated_sequences{$adapter_length};
	print OUTFILE join("\t", $adapter_concatenated_sequence, $adapter_length, $adapter_sequence_count) . "\n";
}
close(OUTFILE) or die "Couldn't close file $adapter_length_counts_outfile";

# Empty all hash and array containers for garbage collection purposes.
undef %trimmed_fastq_sequence_counter;
undef %adapter_regex_length_counter;

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

# gzip_file($fastq_file) - Execute the gzip program to compress a fastq file to save space..
#
# Input paramater(s):
#
# $fastq_file - The fastq file to compress using gzip.
sub gzip_file{

	# The fastq file to compress using gzip.
	my $fastq_file = shift;
	die "Error lost the fastq file to compress using gzip" unless defined $fastq_file ;

	warn "Calling gzip for $fastq_file....\n";
	my $gzipCmd  = "$gzip -9 $fastq_file";
	warn $gzipCmd . "\n\n";
	system($gzip, 
		'-9', $fastq_file,
	) == 0 or die "Error calling $gzipCmd: $?";
}

# $trimmed_seq = get_subseq($sequence, $seq_start, $seq_end) - Get the subsequence based on the input sequence, sequence start, and sequence end.
#
# Input paramater(s):
#
# $sequence - The input sequence to obtain a subsequence.
#
# $seq_start - The start position of the sequence.
#
# $seq_end - The end position of the sequence.
#
# Output paramater(s):
#
# $trimmed_seq - The trimmed subsequence of the input sequence.
sub get_subseq{

	# The input sequence to obtain a subsequence.
        my $sequence = shift;
        die "Error lost input sequence to obtain a subsequence" unless defined $sequence;

        # The start position of the sequence.
        my $seq_start = shift;
        die "Error lost start position of the sequence" unless defined $seq_start;

        # The end position of the sequence.
        my $seq_end = shift;
        die "Error lost end position of the sequence" unless defined $seq_end;

        $seq_start = $seq_start - 1;
        $seq_end = $seq_end;

        my $length = ($seq_end - $seq_start);

        my $trimmed_seq = substr($sequence, $seq_start, $length);

        return uc($trimmed_seq);
}
