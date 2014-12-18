#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Parallel::Loops;
use Switch;
use File::Basename;
use IPC::Open2;
use List::Compare;

# perl trim_adapter_seq_fastq-new.pl -i ~/workspace/GBS_data-08-10-2013/PROJECT_LEADER_DIR/CHRISTIANNE_MCDONALD -p CHRISTIANNE_MCDONALD -c 7 -o ~/workspace/GBS_data-08-10-2013/TRIMMED_adapter_FASTQ_DIR-2014-10-29
my ($fastq_file_dir, $project_name, $restriction_enzymes, $gbs_sequence_length, $adapter_length_min_threshold, $adapter_trim_offset, $min_trimmed_fastq_sequence_length, $output_dir);
GetOptions(
	'i=s'    => \$fastq_file_dir, # The *.fastq input file directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	'p=s'    => \$project_name, # The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	'r=s'    => \$restriction_enzymes, # The restriction enzyme(s) used to digest the genomic sequences. Default: Pstl/MspI
	'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 100
	'm=s'    => \$adapter_length_min_threshold, # The minimum GBS common adapter sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adapter regex searches. Default: 16
	't=s'    => \$adapter_trim_offset, # The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adapter sequence found in the adapter regex searches. Default: 5
	'q=s'    => \$min_trimmed_fastq_sequence_length, # The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Default: 32
	'o=s'    => \$output_dir, # The absolute path to the output directory to contain the trimmed adapter sequence fastq output files.
);

# Display a usage message if the following parameters are not specified.
usage() unless (
	defined $fastq_file_dir
	and $project_name
	and defined $output_dir
);

# The restriction enzyme(s) used to digest the genomic sequences. Default: Pstl/MspI
$restriction_enzymes = 'Pstl/MspI' unless defined $restriction_enzymes;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 100
$gbs_sequence_length = 100 unless defined $gbs_sequence_length;

# The minimum GBS common adapter sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adapter regex searches. Default: 1
$adapter_length_min_threshold = 1 unless defined $adapter_length_min_threshold;

# The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adapter sequence found in the adapter regex searches. Default: 5
$adapter_trim_offset = 5 unless defined $adapter_trim_offset;

# The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Keep in mind that this is the minimum trimmed fastq sequence length used before we add the length of the barcode 
# used for splitting each individual fastq file. Tassel requires sequences at least 32 base pairs (bps) plus the length of a particular barcode that can be in the range of 4-8 base pairs (bps) in length. Default: 32
$min_trimmed_fastq_sequence_length = 32  unless defined $min_trimmed_fastq_sequence_length;

# Program dependencies - The absolute paths to gzip to compress all project leader *.fastq input files.
my $gzip 				= '/bin/gzip';

sub usage {

die <<"USAGE";


Usage: $0 -i fastq_file_dir -p project_name -r restriction_enzymes -l gbs_sequence_length -m adapter_length_min_threshold -t adapter_trim_offset -q min_trimmed_fastq_sequence_length -c blast_num_cpu -o output_dir

DESCRIPTION - A program to trim the GBS common adapter sequence from each GBS fastq file within a particular Genotyping by Sequencing (GBS) project. Fixes the misprimming issue where the GBS common adapter is sequenced along with the DNA of an individual

OPTIONS:

-i fastq_file_dir - The *.fastq input file directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	e.g. /path/to/fastq_file_dir
	
-p project_name - The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	e.g. MPB_MALE_GBS
	
-r restriction_enzymes - The restriction enzyme(s) used to digest the genomic sequences. Default: Pstl/MspI

-l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 100

-m adapter_length_min_threshold - The minimum GBS common adapter sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adapter regex searches. Default: 1

-t adapter_trim_offset - The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adapter sequence found in the adapter regex searches. Default: 5

-q min_trimmed_fastq_sequence_length - The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Keep in mind that this is the minimum trimmed fastq sequence length used before we add the 
length of the barcode used for splitting each individual fastq file. Tassel requires sequences at least 32 base pairs (bps) plus the length of a particular barcode that can be in the range of 4-8 base pairs (bps) in 
length. Default: 32

-o output_dir - The absolute path to the output directory to contain the trimmed adapter sequence fastq output files.

	e.g. /path/to/output_dir
USAGE
}

# Obtain the GBS common adapter sequence to trim from the GBS fastq sequences based on the restriction enzyme(s) used in the digest.
my $fastq_adapter_sequence = '';
switch($restriction_enzymes){
   case 'ApeKI'       { $fastq_adapter_sequence = 'CWGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG' } # Not operational yet
   case 'Pstl/MspI'   { $fastq_adapter_sequence = 'CCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG' }
   else               { die "Input $restriction_enzymes is not one of the recognized restriction enzyme(s)! Please specify either ApeKI or Pstl/MspI on the command line" }
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

# Initialize all the sequence hash counters.
my %original_fastq_seq_counter = ();
my %trimmed_fastq_seq_counter = ();
my %adapter_length_counter = ();

# Find all files in the specified directory with the extension *.fastq.
my ($fastq_files, $fastq_file_count) = find_files($fastq_file_dir, "fastq");

# Iterate through the files with the extension *.fastq.
foreach my $fastq_filename (sort keys %{$fastq_files}){

	# Get the full path to the GBS common adapter length counts file.
	my $fastq_infile = $fastq_files->{$fastq_filename};
	
	warn "Processing " . $fastq_infile . ".....\n";
	
	# Get the basename of the fastq filename without the .fastq extension.
	my $fasta_filename = fileparse($fastq_infile, qr/\.fastq/);
	
	# Need to add the length of the barcode to the $min_trimmed_fastq_sequence_length so that we get the minimum trimmed fastq sequence length plus the length of the barcode sequence.
	my ($individual_id, $barcode, $plate_num, $well_num) = split(/_/, $fasta_filename);
	my $min_trimmed_fastq_sequence_length_plus_barcode = ($min_trimmed_fastq_sequence_length + length($barcode));
	
	# open the fastq input file for parsing.
	open(INFILE, "<$fastq_infile") or die "Couldn't open file $fastq_infile for reading, $!";
	# Parse the fastq files for the fastq header, sequence, plus, and quality scores and reformat to *.fasta format so that we can use adapter regex.
	my ($fastq_header, $fastq_sequence, $fastq_plus, $fastq_quality_scores);
	my $i = 1;
	my $num_fastq_seqs = 0;
	my %fastq_sequences = ();
	while(<INFILE>){
		chomp $_;
		#warn $_ . "\n";
		if(($_ =~ m/^\@[A-Za-z0-9-_]+:\d+:[A-Za-z0-9]+:\d+:\d+:\d+:\d+ \d:[A-Z]:\d:[ACGTRYKMSWBDHVN]*$/) 
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
			
			# Store each fastq sequence entry in a hash variable indexed to the fastq header for identification and an attribute index for the fastq sequence, fastq plus character, and the fastq quality score.
			$fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $fastq_sequence;
			$fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
			$fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $fastq_quality_scores;
			die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($fastq_sequence_length ne $gbs_sequence_length);
			die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length ne fastq_quality_scores_length=$fastq_quality_scores_length" if($fastq_sequence_length ne $fastq_quality_scores_length);
			
			$i = 1;
			$num_fastq_seqs++;
			
		}else{
			$i++;
		}
	}
	close(INFILE) or die "Couldn't close file $fastq_infile";
	
	# Get original fastq sequences counts for each file.
	$original_fastq_seq_counter{$fasta_filename} = $num_fastq_seqs;

	# Execute the adapter regex search using the GBS common adapter sequence as the query and the fastq sequences in fasta format as the blast database.
	# Keep a maximum of target sequences equal to the number of sequences in the fastq file. Retain adapter regex hits based on the minimum adapter length threshold.
	# Keep adapter regex hits that begin at position 1 of the full GBS common adapter sequence up to the minimum adapter length threshold.
	my $adapter_regex_outfile = join('/', $regex_output_dir, join("_", $fasta_filename, "gbs_adapter_regex.tsv"));
	open(OUTFILE, ">$adapter_regex_outfile") or die "Couldn't open file $adapter_regex_outfile for writting, $!";
	print OUTFILE join("\t", "query_name", "target_name", "align_length", "query_start", "query_end", "target_start", "target_end") . "\n";
	foreach my $fastq_header (sort keys %fastq_sequences){
		my $fastq_sequence = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
		my ($align_length, $query_start, $query_end, $target_start, $target_end);
		my $adapter_length_count = 1;
		my $common_adapter_sequence = $fastq_adapter_sequence;
		for(my $i = $adapter_sequence_length; $i >= $adapter_length_min_threshold; $i--){

# 				warn "$i eq $adapter_sequence_length\n";
			if($i eq $adapter_sequence_length){
				my $common_adapter_regex = qr($common_adapter_sequence);
				if($fastq_sequence =~ /$common_adapter_regex/g){
					$target_start = ($-[0] + 1);
					$target_end = $+[0];
					$align_length = ($+[0] - $-[0]);
					warn join("\t", $align_length, $query_start, $query_end, $target_start, $target_end) . "\n";
					
					$query_start = 1;
					$query_end = $adapter_sequence_length;
					warn join("\t", $adapter_sequence_length, $common_adapter_sequence) . "\n";
					my $regex_alignment = join("\t", join("_", "GBS_adapter_sequence", $fastq_adapter_sequence), $fastq_header, $align_length, $query_start, $query_end, $target_start, $target_end);
					print OUTFILE $regex_alignment . "\n";
					last;
				}
			}elsif($i < $adapter_sequence_length){
				my $common_adapter_regex = qr($common_adapter_sequence);
				my $alignment_found = "false";
				while($fastq_sequence =~ /$common_adapter_regex/g){
					if($+[0] eq $gbs_sequence_length){
						$target_start = ($-[0] + 1);
						$target_end = $+[0];
						$align_length = ($+[0] - $-[0]);
						warn join("\t", $align_length, $query_start, $query_end, $target_start, $target_end) . "\n";
						
						my $regex_alignment = join("\t", join("_", "GBS_adapter_sequence", $fastq_adapter_sequence), $fastq_header, $align_length, $query_start, $query_end, $target_start, $target_end);
						print OUTFILE $regex_alignment . "\n";
						$alignment_found = "true";
						last;
					}
				}
				last if($alignment_found eq "true");
			}
			
			$query_start = 1;
			$query_end = ($adapter_sequence_length - $adapter_length_count);
			
			$common_adapter_sequence = get_subseq($common_adapter_sequence, $query_start, $query_end);
			warn join("\t", ($adapter_sequence_length - $adapter_length_count), $common_adapter_sequence) . "\n";
			$adapter_length_count++;
		}
	}
	close(OUTFILE) or die "Couldn't close file $adapter_regex_outfile";
	
	# Generate the trimmed adapter regex files filtered for further processing.
	my $trimmed_adapter_regex_outfile = join('/', $trimmed_regex_output_dir, join("_", $fasta_filename, "gbs_adapter_regex.tsv"));
	open(ADAPTER_REGEX_OUTFILE, ">$trimmed_adapter_regex_outfile") or die "Couldn't open file $trimmed_adapter_regex_outfile for writting, $!";
	print ADAPTER_REGEX_OUTFILE join("\t", "query_name", "target_name", "align_length", "query_start", "query_end", "target_start", "target_end") . "\n";
	
	# The $trimmed_seqs_layout_outfile contains the trimmed coordinates layout for each sequence trimmed of the GBS common adapter sequence that passed the retaining criteria.
	my $trimmed_seqs_layout_outfile = join('/', $trimmed_layout_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset, "trimmed_seqs_layout") . ".txt");
	open(TRIMMED_LAYOUT_OUTFILE, ">$trimmed_seqs_layout_outfile") or die "Couldn't open file $trimmed_seqs_layout_outfile for writting, $!";
	print TRIMMED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adapter_start", "trimmed_adapter_end",
	"trimmed_adapter_length", "adapter_seq_start", "adapter_seq_end", "adapter_seq_length") . "\n";

	# The $removed_seqs_layout_outfile contains the trimmed coordinates layout for each sequence trimmed of the GBS common adapter sequence that failed the retaining criteria and was therefore removed.
	my $removed_seqs_layout_outfile = join('/', $removed_layout_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset, "removed_seqs_layout") . ".txt");
	open(REMOVED_LAYOUT_OUTFILE, ">$removed_seqs_layout_outfile") or die "Couldn't open file $removed_seqs_layout_outfile for writting, $!";
	print REMOVED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adapter_start", "trimmed_adapter_end",
	"trimmed_adapter_length", "adapter_seq_start", "adapter_seq_end", "adapter_seq_length") . "\n";
	
	my %trimmed_fastq_sequences = ();
	my @trimmed_fastq_list = ();
	my %removed_fastq_sequences = ();
	my @removed_fastq_list = ();
	
	# Parse the tab-delimited adapter blast output so that we can trimm the fastq sequences that contain the GBS adapter sequence.
	# Create new adapter regex files so that we can visualize where we trimmed the sequence. Store results in a hash array indexed to the target sequence so that 
	# we can use the longest length of the GBS common adapter sequence as the optimal alignment in the trimming step.
	warn "Parsing $adapter_regex_outfile to get most optimal adapter regex hits.....\n";
	open(INFILE, "<$adapter_regex_outfile") or die "Couldn't open file $adapter_regex_outfile for reading, $!";
	my %adapter_regex_hits = ();
	$i = 0;
	while(<INFILE>){
		chomp $_;
		if($i ne 0){
			#warn $_ . "\n";
			my @split_adapter_regex_hit =  split(/\t/, $_);
			my ($query_name, $target_name, $align_length, $query_start, $query_end, $target_start, $target_end) = @split_adapter_regex_hit;
		
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
				
				# If the trimmed sequence length is greater than or equal to the minimum trimmed fastq sequence length plus the length of the barcode then the trimmed and sequence count files.
				if($trimmed_fastq_sequence_length >= $min_trimmed_fastq_sequence_length_plus_barcode){
					# Store the trimmed fastq sequence entry in a hash variable indexed by the fastq sequence header.
					$trimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $trimmed_fastq_sequence;
					$trimmed_fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
					$trimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $trimmed_fastq_quality_scores;
					
					my $trimmed_fastq_regex = join("\t", join("_", "trimmed_fastq_sequence_offset", $adapter_trim_offset), $target_name, $trimmed_fastq_sequence_length, 1, (($target_start - $adapter_trim_offset) - 1), 1, (($target_start - $adapter_trim_offset) - 1));
					my $trimmed_adapter_regex = join("\t", join("_", "trimmed_adapter_sequence_offset", $adapter_trim_offset), $target_name, $trimmed_adapter_sequence_length, ($target_start - $adapter_trim_offset), $fastq_sequence_length, ($target_start - $adapter_trim_offset), $fastq_sequence_length);
					my $original_adapter_regex = join("\t", $query_name, $target_name, $align_length, $query_start, $query_end, $target_start, $target_end);
						
					print ADAPTER_REGEX_OUTFILE $trimmed_fastq_regex . "\n";
					print ADAPTER_REGEX_OUTFILE $original_adapter_regex . "\n";
					print ADAPTER_REGEX_OUTFILE $trimmed_adapter_regex . "\n";
					
					my $trimmed_fastq_header = $fastq_header;
					my $trimmed_seqs_layout_header = join("_", $fasta_filename, $trimmed_fastq_header);
					print TRIMMED_LAYOUT_OUTFILE join("\t", $trimmed_seqs_layout_header, 1, (($target_start - $adapter_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adapter_trim_offset), $fastq_sequence_length, $trimmed_adapter_sequence_length, $target_start, $target_end, $align_length) . "\n";
					
					$adapter_length_counter{$align_length}++;
					push(@trimmed_fastq_list, $fastq_header);
					
				}elsif($trimmed_fastq_sequence_length < $min_trimmed_fastq_sequence_length_plus_barcode){ #If this alignment doesn't meet our filtering criteria separate fastq sequence from trimmed fastq sequence data.
					$removed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $trimmed_fastq_sequence;
					$removed_fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
					$removed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $trimmed_fastq_quality_scores;
					
					my $removed_fastq_header = $fastq_header;
					my $removed_seqs_layout_header = join("_", $fasta_filename, $removed_fastq_header);
					
					print REMOVED_LAYOUT_OUTFILE join("\t", $removed_seqs_layout_header, 1, (($target_start - $adapter_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adapter_trim_offset), $fastq_sequence_length, $trimmed_adapter_sequence_length, $target_start, $target_end, $align_length) . "\n";

					push(@removed_fastq_list, $fastq_header);
				}
			}
		}
		$i++;
	}
	close(ADAPTER_REGEX_OUTFILE) or die "Couldn't close file $trimmed_adapter_regex_outfile";
	close(TRIMMED_LAYOUT_OUTFILE) or die "Couldn't close file $trimmed_seqs_layout_outfile";
	close(REMOVED_LAYOUT_OUTFILE) or die "Couldn't close file $removed_seqs_layout_outfile";
	close(INFILE) or die "Couldn't close file $adapter_regex_outfile";
	
	# Grab the list of fastq sequence headers
	my @fastq_sequence_list = keys %fastq_sequences;
	
	# Grab the sub list of untrimmed fastq sequences and put them into the trimmed_fastq_sequences hash.
	my @fastq_list2remove = ();
	push(@fastq_list2remove, @trimmed_fastq_list);
	push(@fastq_list2remove, @removed_fastq_list);
	
	my $fastq_list_comparision = List::Compare->new(\@fastq_sequence_list, \@fastq_list2remove);
	my @untrimmed_fastq_sequence_list = $fastq_list_comparision->get_unique;

	my %untrimmed_fastq_sequences = ();
	foreach my $fastq_header (@untrimmed_fastq_sequence_list){
		my $fastq_sequence = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
		my $fastq_plus = $fastq_sequences{$fastq_header}{'PLUS'};
		my $fastq_quality_scores = $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
		
		$untrimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $fastq_sequence;
		$untrimmed_fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
		$untrimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $fastq_quality_scores;
	}
	
	# Print out trimmed fastq sequences first so that we can see what was trimmed.
	my $trimmed_fastq_outfile = join("/", $trimmed_fastq_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset) . ".fastq");
	open(OUTFILE, ">$trimmed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for writting, $!";
	foreach my $fastq_header (sort keys %trimmed_fastq_sequences){
		
		my $fastq_sequence = $trimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
		my $fastq_plus = $trimmed_fastq_sequences{$fastq_header}{'PLUS'};
		my $fastq_quality_scores = $trimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
		my $fastq_sequence_length = length($fastq_sequence);
		my $new_fastq_header = join("\001", $fastq_header, $fasta_filename, join("=", "length", $fastq_sequence_length));
		print OUTFILE  $new_fastq_header . "\n";
		print OUTFILE  $fastq_sequence . "\n";
		print OUTFILE  $fastq_plus . "\n";
		print OUTFILE  $fastq_quality_scores . "\n";

		$trimmed_fastq_seq_counter{$fasta_filename}{'TRIMMED'}++;
		
	}
	
	# Print out the rest of the sequences that were not trimmed because they either did not have an alignment or a significant alignment that passed the trimming threshold.
	foreach my $fastq_header (sort keys %untrimmed_fastq_sequences){
		
		my $fastq_sequence = $untrimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
		my $fastq_plus = $untrimmed_fastq_sequences{$fastq_header}{'PLUS'};
		my $fastq_quality_scores = $untrimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
		my $fastq_sequence_length = length($fastq_sequence);
		my $new_fastq_header = join("\001", $fastq_header, $fasta_filename, join("=", "length", $fastq_sequence_length));
		print OUTFILE  $new_fastq_header . "\n";
		print OUTFILE  $fastq_sequence . "\n";
		print OUTFILE  $fastq_plus . "\n";
		print OUTFILE  $fastq_quality_scores . "\n";
		
		$trimmed_fastq_seq_counter{$fasta_filename}{'UNTRIMMED'}++;

	}
	close(OUTFILE) or die "Couldn't close file $trimmed_fastq_outfile";
	
	# Print out trimmed fastq sequences that did not pass the filtering critera so that we can see what was trimmed and removed.
	my $removed_fastq_outfile = join("/", $removed_fastq_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset, "removed_sequences") . ".fastq");
	open(OUTFILE, ">$removed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for writting, $!";
	foreach my $fastq_header (sort keys %removed_fastq_sequences){
	
		my $fastq_sequence = $removed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
		my $fastq_plus = $removed_fastq_sequences{$fastq_header}{'PLUS'};
		my $fastq_quality_scores = $removed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
		my $fastq_sequence_length = length($fastq_sequence);
		my $new_fastq_header = join("\001", $fastq_header, $fasta_filename, join("=", "length", $fastq_sequence_length));
		print OUTFILE  $new_fastq_header . "\n";
		print OUTFILE  $fastq_sequence . "\n";
		print OUTFILE  $fastq_plus . "\n";
		print OUTFILE  $fastq_quality_scores . "\n";
		
		$trimmed_fastq_seq_counter{$fasta_filename}{'REMOVED'}++;

	}
	close(OUTFILE) or die "Couldn't close file $removed_fastq_outfile";
	
	# Print the GBS common adapter length counts.
	my $trimmed_adapter_counts_outfile = join("/", $trimmed_adapter_counts_output_dir, join("_", $fasta_filename, "trimmed_offset", $adapter_trim_offset, "adapter_length_counts") . ".txt");
	open(OUTFILE, ">$trimmed_adapter_counts_outfile") or die "Couldn't open file $trimmed_adapter_counts_outfile for writting, $!";
	print OUTFILE join("\t", "adapter_sequence_id", "adapter_length", "adapter_sequence_count") . "\n";
	foreach my $adapter_length (sort {$b <=> $a} keys %adapter_length_counter){
		my $adapter_sequence_count = $adapter_length_counter{$adapter_length};
		my $adapter_concatenated_sequence = $adapter_concatenated_sequences{$adapter_length};
		print OUTFILE join("\t", $adapter_concatenated_sequence, $adapter_length, $adapter_sequence_count) . "\n";
	}
	close(OUTFILE) or die "Couldn't close file $trimmed_adapter_counts_outfile";
	
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

# Generate bulk trimmed fastq file for this Genotyping by Sequencing (GBS) project. 
# $trimmed_fastq_bulk_outfile contains the trimmed GBS sequences in fastq format.
# Find all files in the specified directory with the extension *.fastq.
my ($trimmed_fastq_files, $trimmed_fastq_file_count) = find_files($trimmed_fastq_output_dir, "fastq");
my $trimmed_fastq_bulk_outfile = join('/', $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adapter_trim_offset) . ".fastq");
open(TRIMMED_BULK_OUTFILE, ">$trimmed_fastq_bulk_outfile") or die "Couldn't open file $trimmed_fastq_bulk_outfile for writting, $!";
# Iterate through the files with extension *.fastq.
foreach my $trimmed_fastq_filename (sort keys %{$trimmed_fastq_files}){
	
	# Get the full path to the trimmed fastq file.
	my $trimmed_fastq_infile = $trimmed_fastq_files->{$trimmed_fastq_filename};
	
	# Parse trimmed fastq file and concatenate to the bulk fastq output file for this project.
	open(INFILE, "<$trimmed_fastq_infile") or die "Couldn't open file $trimmed_fastq_infile for reading, $!";
	while(<INFILE>){
		chomp $_;
		print TRIMMED_BULK_OUTFILE $_ . "\n";
	}
	close(INFILE) or die "Couldn't close file $trimmed_fastq_infile";
	
	# Compress the trimmed fastq file using gzip.
	gzip_file($trimmed_fastq_infile);
}
close(TRIMMED_BULK_OUTFILE) or die "Couldn't close file $trimmed_fastq_bulk_outfile";

# Compress the bulk trimmed fastq file using gzip.
gzip_file($trimmed_fastq_bulk_outfile);

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
foreach my $fasta_filename (sort keys %trimmed_fastq_seq_counter){

    # Get the fastq untrimmed, trimmed, and removed counts.
    my $num_untrimmed_seqs = $trimmed_fastq_seq_counter{$fasta_filename}{'UNTRIMMED'};
    my $num_trimmed_seqs = $trimmed_fastq_seq_counter{$fasta_filename}{'TRIMMED'};
    my $num_removed_seqs = $trimmed_fastq_seq_counter{$fasta_filename}{'REMOVED'};


    # If trimmed or removed sequence counters how zero count reflect in number of sequence variables.
    $num_untrimmed_seqs = 0 unless(defined($num_untrimmed_seqs));
    $num_trimmed_seqs = 0 unless(defined($num_trimmed_seqs));
    $num_removed_seqs = 0 unless(defined($num_removed_seqs));

    # Calculate the total number of sequences.
#    my $total_num_seqs = ($num_untrimmed_seqs + $num_trimmed_seqs + $num_removed_seqs);
    my $total_num_seqs = $original_fastq_seq_counter{$fasta_filename};
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
undef %trimmed_fastq_seq_counter;
undef %adapter_regex_length_counter;

# find the fastq files found in the given fastq input file directory and return a hash list of the files and the number of files found.
# Will add a pod documentation for this sub routine.
sub find_files{
    
	my $infile_dir = shift;
	die "Error lost input file directory" unless defined $infile_dir;
	
	my $suffix = shift;
	die "Error lost file extension suffix directory" unless defined $suffix;
	
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

# Compress a file using gzip to save space.
# Will add a pod documentation for this sub routine.
sub gzip_file{
	
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
