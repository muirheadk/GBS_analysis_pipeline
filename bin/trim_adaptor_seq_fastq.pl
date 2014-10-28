#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;
use Math::Round;
use List::Compare;

# perl trim_adaptor_seq_fastq.pl -i ~/workspace/GBS_data-08-10-2013/PROJECT_LEADER_DIR -c 7 -o ~/workspace/GBS_data-08-10-2013/TRIMMED_ADAPTOR_FASTQ_DIR
my ($project_leader_dir, $fastq_adaptor_sequence, $gbs_sequence_length, $adaptor_length_min_threshold, $adaptor_length_max_threshold, $adaptor_trim_offset, $blast_num_cpu, $output_dir);
GetOptions(
      'i=s'    => \$project_leader_dir,
      'a=s'    => \$fastq_adaptor_sequence,
      'l=s'    => \$gbs_sequence_length,
      'n=s'    => \$adaptor_length_min_threshold,
#       'm=s'    => \$adaptor_length_max_threshold,
      't=s'    => \$adaptor_trim_offset,
      'c=s'    => \$blast_num_cpu,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $project_leader_dir
      and defined $output_dir
);

$fastq_adaptor_sequence = 'CCGAGATCGGAAGAGCGGGGACTTTAAGC' unless defined $fastq_adaptor_sequence;
$gbs_sequence_length = 100 unless defined $gbs_sequence_length;
$adaptor_length_min_threshold = 18 unless defined $adaptor_length_min_threshold;
$adaptor_trim_offset = 5 unless defined $adaptor_trim_offset;
$blast_num_cpu = 2 unless defined $blast_num_cpu;

my ($makeblastdb, $blastn, $gbs_adaptor_align_graphics, $send_mail);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$blastn				= '/usr/local/bin/blastn';
$gbs_adaptor_align_graphics	= '/GBS_analysis_pipeline/bin/gbs_adaptor_align_graphics.pl';
$send_mail			= '/GBS_analysis_pipeline/bin/send_mail.pl';

sub usage {

die <<"USAGE";


Usage: $0 -i project_leader_dir -a fastq_adaptor_sequence -l gbs_sequence_length -n adaptor_length_min_threshold -m adaptor_length_max_threshold -t adaptor_trim_offset -c blast_num_cpu -o output_dir

Description - 

OPTIONS:

	-i project_leader_dir - 

	-a fastq_adaptor_sequence -

	-l gbs_sequence_length - 

	-n adaptor_length_min_threshold -

# 	-m adaptor_length_max_threshold - 

	-t adaptor_trim_offset - 

	-c blast_num_cpu -

	-o output_dir -



USAGE
}


# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Find all files in the project leader directory with the extension *.fastq.
my ($project_fastq_files, $file_count) = find_fastq_files($project_leader_dir);

# Iterate through the files with the extension *.fastq, rename them to *.fasta so that we can use adaptor blastn.
# my $num_of_files = 0;
foreach my $project_leader (sort keys %{$project_fastq_files}){

	warn "Processing $project_leader\'s fastq files.....\n";
	# Create output directory if it doesn't already exist.
	my $project_leader_dir = join('/', $output_dir, $project_leader);
	unless(-d $project_leader_dir){
		mkdir($project_leader_dir, 0777) or die "Can't make directory: $!";
	}

	# Create output directory if it doesn't already exist.
	my $fasta_output_dir = join('/', $project_leader_dir, "FASTA_FILES");
	unless(-d $fasta_output_dir){
		mkdir($fasta_output_dir, 0777) or die "Can't make directory: $!";
	}

	# Create output directory if it doesn't already exist.
	my $blastn_output_dir = join('/', $project_leader_dir, "ADAPTOR_BLASTN_FILES");
	unless(-d $blastn_output_dir){
		mkdir($blastn_output_dir, 0777) or die "Can't make directory: $!";
	}

	# Create output directory if it doesn't already exist.
	my $trimmed_output_dir = join('/', $project_leader_dir, "TRIMMED_OUTPUT_FILES");
	unless(-d $trimmed_output_dir){# Need this to make the bulk fastq sequence file.
		mkdir($trimmed_output_dir, 0777) or die "Can't make directory: $!";
	}

	my $trimmed_fastq_bulk_outfile = join('/', $trimmed_output_dir, join("_", $project_leader, "trimmed") . ".fastq");
 	open(TRIMMED_BULK_OUTFILE, ">$trimmed_fastq_bulk_outfile") or die "Couldn't open file $trimmed_fastq_bulk_outfile for writting, $!";
 	my $trimmed_seqs_layout_outfile = join('/', $trimmed_output_dir, join("_", $project_leader, "trimmed_seqs_layout") . ".txt");
 	open(TRIMMED_LAYOUT_OUTFILE, ">$trimmed_seqs_layout_outfile") or die "Couldn't open file $trimmed_seqs_layout_outfile for writting, $!";
 	print TRIMMED_LAYOUT_OUTFILE join("\t", "sequence_id", "trimmed_fastq_start", "trimmed_fastq_end", "trimmed_fastq_length", "trimmed_adaptor_start", "trimmed_adaptor_end", 
 		"trimmed_adaptor_length", "adaptor_seq_start", "adaptor_seq_end", "adaptor_seq_length") . "\n";
	my @trimmed_fastq_files = ();
	my %fastq_seq_counter = ();
	my %adaptor_length_counter = ();
	foreach my $fastq_infile (sort @{$project_fastq_files->{$project_leader}}){

		warn $fastq_infile . "\n";
		
		my %fastq_sequences = ();
		open(INFILE, "<$fastq_infile") or die "Couldn't open file $fastq_infile for reading, $!";
		my $fasta_filename = fileparse($fastq_infile, qr/\.fastq/);
		my $fasta_target_outfile = join("/", $fasta_output_dir, $fasta_filename . ".fasta");
		open(OUTFILE, ">$fasta_target_outfile") or die "Couldn't open file $fasta_target_outfile for writting, $!";
		#my $fastq_test_outfile = join("/", $fasta_output_dir, $fasta_filename . ".test.fastq");
		#open(OUTFILE2, ">$fastq_test_outfile") or die "Couldn't open file $fastq_test_outfile for writting, $!";
		my ($fastq_header, $fastq_sequence, $fastq_plus, $fastq_quality_scores);
		my $i = 1;
		my $fastq_counter = 0;
		while(<INFILE>){
			chomp $_;
			#warn $_ . "\n";
			if($_ =~ m/^\@[A-Za-z0-9-_]+:\d+:[A-Za-z0-9]+:\d+:\d+:\d+:\d+ \d:[A-Z]:\d:[ACGTRYKMSWBDHVN]*$/){ # @HWI-ST767:215:C30VBACXX:8:1101:1801:1484 1:N:0:
				$fastq_header = $_;
#  				die $fastq_header;
			}elsif($_ =~ m/^[ACGTRYKMSWBDHVN]+$/i){
				$fastq_sequence = $_;
#  				die $fastq_sequence;
			}elsif($_ =~ m/^\+$/){
				$fastq_plus = $_;
#  				die $fastq_plus;
			}elsif($_ =~ m/^.+$/){
				$fastq_quality_scores = $_;
#   				die $fastq_quality_scores;
			}
			
			if(($i % 4) eq 0){
				
				die "Error: fastq_header is undefined" unless(defined($fastq_header));
				die "Error: fastq_sequence is undefined" unless(defined($fastq_sequence));
				die "Error: fastq_plus is undefined" unless(defined($fastq_plus));
				die "Error: fastq_quality_scores is undefined" unless(defined($fastq_quality_scores));
				
				my $fastq_sequence_length = length($fastq_sequence);
				my $fastq_quality_scores_length = length($fastq_quality_scores);
				
				$fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $fastq_sequence;
				$fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
				$fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $fastq_quality_scores;
				die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($fastq_sequence_length ne $gbs_sequence_length);
				die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length ne fastq_quality_scores_length=$fastq_quality_scores_length" if($fastq_sequence_length ne $fastq_quality_scores_length);
				
				my $fasta_header = $fastq_header;
				$fasta_header =~ s/\s/_/g;
				$fasta_header = join("_", $fasta_header, "length=$fastq_sequence_length");
				print OUTFILE join("\n", join("", ">", $fasta_header), $fastq_sequence) . "\n";
				
				#print OUTFILE2  $fastq_header . "\n";
				#print OUTFILE2  $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} . "\n";
				#print OUTFILE2  $fastq_sequences{$fastq_header}{'PLUS'} . "\n";
				#print OUTFILE2  $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} . "\n";
				$i = 1;
				$fastq_counter++;
			}else{
				$i++;
			}
		}
		close(INFILE) or die "Couldn't close file $fastq_infile";
		close(OUTFILE) or die "Couldn't close file $fasta_target_outfile";
		#close(OUTFILE2) or die "Couldn't close file $fastq_test_outfile";

		my ($adaptor_blastn_tsv_outfile, $adaptor_blastn_aln_outfile, $adaptor_sequence_length) = generate_adaptor_blastn($fastq_adaptor_sequence, $fasta_target_outfile, $blast_num_cpu, $blastn_output_dir);

		# Create output directory if it doesn't already exist.
# 		my $adaptor_graphics_output_dir = join('/', $blastn_output_dir, "ADAPTOR_BLASTN_GRAPHICS_FILES");
# 		unless(-d $adaptor_graphics_output_dir){
# 			mkdir($adaptor_graphics_output_dir, 0777) or die "Can't make directory: $!";
# 		}

# 		generate_adaptor_blast_graphics($adaptor_blastn_tsv_outfile, $adaptor_graphics_output_dir);

		# Create output directory if it doesn't already exist.
		my $trimmed_blastn_output_dir = join('/', $trimmed_output_dir, "TRIMMED_ADAPTOR_BLASTN_FILES");
		unless(-d $trimmed_blastn_output_dir){
			mkdir($trimmed_blastn_output_dir, 0777) or die "Can't make directory: $!";
		}
		
		# Parse the tab-delimited adaptor blast output so that we can trimm the fastq sequences that contain the GBS adaptor sequence.
		# Create new adaptor blastn files so that we can visualize where we trimmed the sequence.
		my $trimmed_adaptor_blastn_outfile = join('/', $trimmed_blastn_output_dir, $fasta_filename . ".gbs_adaptor_blastn.tsv");
		open(INFILE, "<$adaptor_blastn_tsv_outfile") or die "Couldn't open file $adaptor_blastn_tsv_outfile for reading, $!";
		open(OUTFILE, ">$trimmed_adaptor_blastn_outfile") or die "Couldn't open file $trimmed_adaptor_blastn_outfile for writting, $!";
		print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch", 
		"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score", "graphics_colour") . "\n";
		my %trimmed_fastq_sequences = ();
		my @trimmed_fastq_list = ();
		$i = 0;
		while(<INFILE>){
			chomp $_;
			if($i ne 0){
		 		#warn $_ . "\n";
				my @adaptor_blastn_hit =  split(/\t/, $_);
				my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
					$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, $gylph_colour) = @adaptor_blastn_hit;
					
    				if($query_start > $query_end){
        				warn "Found hit in antisense direction for query_id, so skipping because we aren't interested in this alignment....\n";
        				next;
    				}
    				if($target_start > $target_end){
        				warn "Found hit in antisense direction for target_id, so skipping because we aren't interested in this alignment....\n";
        				next;
    				}
    				if($num_mismatch > 0){
        				warn "Found hit with number of mismatches greater than 0, so skipping because we aren't interested in this alignment....\n";
        				next;
    				}
    				if($num_gaps > 0){
        				warn "Found hit with number of gaps greater than 0, so skipping because we aren't interested in this alignment....\n";
        				next;
				}

				my $fasta_header = $target_name;
				my ($fastq_header_part1, $fastq_header_part2, $fasta_header_length) = split(/_/, $fasta_header);
				my $fastq_header = join(" ", $fastq_header_part1, $fastq_header_part2);
				
				my $fastq_sequence = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
				my $trimmed_fastq_plus = $fastq_sequences{$fastq_header}{'PLUS'}; 
				my $fastq_quality_scores = $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
				
				#die $fastq_sequence;
				my $fastq_sequence_length = length($fastq_sequence);
				my $predicted_adaptor_start = ($fastq_sequence_length - $adaptor_sequence_length);
				
# 				my ($min_threshold_coverage, $max_threshold_coverage, $trimmed_fastq_sequence, $trimmed_fastq_quality_scores, $trimmed_adaptor_sequence);
				my ($trimmed_fastq_sequence, $trimmed_fastq_quality_scores, $trimmed_adaptor_sequence);
				# We are using ($adaptor_sequence_length + 1) because blast thinks that the adaptor sequence is 30bp when it is really 29. 
				# So any correct threshold value would give incorrect results or nothing going through the if statement.
# 				$min_threshold_coverage = round((($adaptor_length_min_threshold/($adaptor_sequence_length + 1)) * 100));
# 				$max_threshold_coverage = round((($adaptor_length_max_threshold/($adaptor_sequence_length + 1)) * 100));
# 				die join("\t", "query_coverage=$query_coverage", "threshold_coverage=$max_threshold_coverage");
				if($align_length >= $adaptor_length_min_threshold){
					$trimmed_fastq_sequence = get_subseq($fastq_sequence, 1, (($target_start - $adaptor_trim_offset) - 1));
					$trimmed_fastq_quality_scores = get_subseq($fastq_quality_scores, 1, (($target_start - $adaptor_trim_offset) - 1));
					
					$trimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $trimmed_fastq_sequence;
					$trimmed_fastq_sequences{$fastq_header}{'PLUS'} = $trimmed_fastq_plus;
					$trimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $trimmed_fastq_quality_scores;
					
					$trimmed_adaptor_sequence = get_subseq($fastq_sequence, ($target_start - $adaptor_trim_offset), $fastq_sequence_length);
					
					my $trimmed_fastq_sequence_length = length($trimmed_fastq_sequence);
					my $trimmed_adaptor_sequence_length = length($trimmed_adaptor_sequence);
					my $trimmed_fastq_blastn = join("\t", join("_", "trimmed_fastq_sequence_offset", $adaptor_trim_offset), $target_name, $query_coverage, $percent_identity, $trimmed_fastq_sequence_length, $num_mismatch, $num_gaps, 1, (($target_start - $adaptor_trim_offset) - 1), 1, (($target_start - $adaptor_trim_offset) - 1), $e_value, $bit_score, "blue");
					my $trimmed_adaptor_blastn = join("\t", join("_", "trimmed_adaptor_sequence_offset", $adaptor_trim_offset), $target_name, $query_coverage, $percent_identity, $trimmed_adaptor_sequence_length, $num_mismatch, $num_gaps, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, $e_value, $bit_score, "red");
					my $original_adaptor_blastn = join("\t", $query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
						$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, "orange");
					print OUTFILE $trimmed_fastq_blastn . "\n";
					print OUTFILE $original_adaptor_blastn . "\n";
					print OUTFILE $trimmed_adaptor_blastn . "\n";
					
					my $trimmed_fastq_header = $fastq_header;
					my $trimmed_seqs_layout_header = join("_", $fasta_filename, $trimmed_fastq_header);
					print TRIMMED_LAYOUT_OUTFILE join("\t", $trimmed_seqs_layout_header, 1, (($target_start - $adaptor_trim_offset) - 1), $trimmed_fastq_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, $trimmed_adaptor_sequence_length, $target_start, $target_end, $align_length) . "\n";
					
					my $adaptor_concatenated_sequence;  
					if($query_name =~ m/([\w_]+)_([ACGT]+)/){
						my ($query_id, $adaptor_sequence) = ($1, $2);
						my $adaptor_end = length($adaptor_sequence);
						my $adaptor_sub_sequence = get_subseq($adaptor_sequence, $query_start, $query_end);
						my $adaptor_start_sequence = get_subseq($adaptor_sequence, 1, ($query_start - 1));
						my $adaptor_end_sequence = get_subseq($adaptor_sequence, ($query_end + 1), $adaptor_end);

						if($adaptor_start_sequence eq ""){
							$adaptor_concatenated_sequence  = join("-", $query_start, $adaptor_sub_sequence, $query_end, $adaptor_end_sequence);
						}elsif($adaptor_end_sequence eq ""){
							$adaptor_concatenated_sequence  = join("-", $adaptor_start_sequence, $query_start, $adaptor_sub_sequence, $query_end);
						}else{
							$adaptor_concatenated_sequence  = join("-", $adaptor_start_sequence, $query_start, $adaptor_sub_sequence, $query_end, $adaptor_end_sequence);
						}
					}
					
					$adaptor_length_counter{$align_length}{$adaptor_concatenated_sequence}++;
					push(@trimmed_fastq_list, $fastq_header);
				}
			}
			$i++;

		}
		close(INFILE) or die "Couldn't close file $adaptor_blastn_tsv_outfile";
		close(OUTFILE) or die "Couldn't close file $trimmed_adaptor_blastn_outfile";
		
		# Grab the list of fastq sequence headers
		my @fastq_sequence_list = keys %fastq_sequences;
		
		# Grab the sub list of untrimmed fastq sequences and put them into the trimmed_fastq_sequences hash.
		my $fastq_list_comparision = List::Compare->new(\@fastq_sequence_list, \@trimmed_fastq_list);
# 		die "unique fastq sequences: ", $fastq_list_comparision->get_unique, "\n";		
		my @untrimmed_fastq_sequence_list = $fastq_list_comparision->get_unique;
		
# 		die join("\t", scalar(@fastq_sequence_list), scalar(@trimmed_fastq_list), scalar(@untrimmed_fastq_sequence_list));
		my %untrimmed_fastq_sequences = ();
		foreach my $fastq_header (@untrimmed_fastq_sequence_list){
			my $fastq_sequence = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
			my $fastq_plus = $fastq_sequences{$fastq_header}{'PLUS'};
			my $fastq_quality_scores = $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
			
			$untrimmed_fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $fastq_sequence;
			$untrimmed_fastq_sequences{$fastq_header}{'PLUS'} = $fastq_quality_scores;
			$untrimmed_fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $fastq_plus;		
		}
		
		# Create output directory if it doesn't already exist.
# 		my $trimmed_adaptor_graphics_output_dir = join('/', $trimmed_output_dir, "TRIMMED_ADAPTOR_BLASTN_GRAPHICS_FILES");
# 		unless(-d $trimmed_adaptor_graphics_output_dir){
# 			mkdir($trimmed_adaptor_graphics_output_dir, 0777) or die "Can't make directory: $!";
# 		}
		
# 		generate_adaptor_blast_graphics($trimmed_adaptor_blastn_outfile, $trimmed_adaptor_graphics_output_dir);

		# Create output directory if it doesn't already exist.
		my $trimmed_fastq_output_dir = join('/', $trimmed_output_dir, "TRIMMED_FASTQ_FILES");
		unless(-d $trimmed_fastq_output_dir){
			mkdir($trimmed_fastq_output_dir, 0777) or die "Can't make directory: $!";
		}

		my $trimmed_fastq_outfile = join("/", $trimmed_fastq_output_dir, join("_", $fasta_filename, "trimmed") . ".fastq");
		open(OUTFILE, ">$trimmed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for writting, $!";
		# Print out trimmed fastq sequences first so that we can see what was trimmed.
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

		open(INFILE, "<$trimmed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for reading, $!";
		while(<INFILE>){
			chomp $_;
			print TRIMMED_BULK_OUTFILE $_ . "\n";
		}
		close(INFILE) or die "Couldn't close file $trimmed_fastq_outfile";

		# Empty all hash and array containers so that we don't use fastq sequences or fastq headers from different files.
		%trimmed_fastq_sequences = ();
		%untrimmed_fastq_sequences = ();
		%fastq_sequences = ();
		@trimmed_fastq_list = ();
		@fastq_sequence_list = ();
		@untrimmed_fastq_sequence_list = ();
 	}
 	close(TRIMMED_BULK_OUTFILE) or die "Couldn't close file $trimmed_fastq_bulk_outfile";
 	close(TRIMMED_LAYOUT_OUTFILE) or die "Couldn't close file $trimmed_seqs_layout_outfile";
 	
 	my $fastq_seq_counts_outfile = join("/", $trimmed_output_dir, join("_", $project_leader, "fastq_seq_counts") . ".txt");
	open(OUTFILE, ">$fastq_seq_counts_outfile") or die "Couldn't open file $fastq_seq_counts_outfile for writting, $!";
	print OUTFILE join("\t", "fastq_file_name", "total_num_seqs", "num_untrimmed_seqs", "percent_untrimmed_seqs", "num_trimmed_seqs", "percent_trimmed_seqs") . "\n";
	foreach my $fasta_filename (sort keys %fastq_seq_counter){
	
		my $num_untrimmed_seqs = $fastq_seq_counter{$fasta_filename}{'UNTRIMMED'};
		my $num_trimmed_seqs = $fastq_seq_counter{$fasta_filename}{'TRIMMED'};
		
		my $total_num_seqs = ($num_untrimmed_seqs + $num_trimmed_seqs);
		
		my $percent_untrimmed_seqs = (($num_untrimmed_seqs/$total_num_seqs) * 100);
		my $percent_trimmed_seqs = (($num_trimmed_seqs/$total_num_seqs) * 100);
		
		print OUTFILE join("\t", $fasta_filename, $total_num_seqs, $num_untrimmed_seqs, $percent_untrimmed_seqs, $num_trimmed_seqs, $percent_trimmed_seqs) . "\n";
	}
 	close(OUTFILE) or die "Couldn't close file $fastq_seq_counts_outfile";
 	
	my $adaptor_length_counts_outfile = join("/", $trimmed_output_dir, join("_", $project_leader, "adaptor_length_counts") . ".txt");
	open(OUTFILE, ">$adaptor_length_counts_outfile") or die "Couldn't open file $adaptor_length_counts_outfile for writting, $!";
	print OUTFILE join("\t", "adapter_sequence_id", "adapter_length", "adaptor_sequence_count") . "\n";
	foreach my $adaptor_length (sort {$b <=> $a} keys %adaptor_length_counter){
		foreach my $adaptor_concatenated_sequence (keys %{$adaptor_length_counter{$adaptor_length}}){
			my $adaptor_sequence_count = $adaptor_length_counter{$adaptor_length}{$adaptor_concatenated_sequence};
	 		print OUTFILE join("\t", $adaptor_concatenated_sequence, $adaptor_length, $adaptor_sequence_count) . "\n";
		}
		
	}
 	close(OUTFILE) or die "Couldn't close file $adaptor_length_counts_outfile";

 	%fastq_seq_counter = ();
 	%adaptor_length_counter = ();
}

my $subject = "\"$0 Process Complete\"";
my $message = "\"$0 process finished successfully! You can find all the trimmed adaptor fastq output files for each project leader in the $output_dir directory.\"";
send_mail($subject, $message, 'both');

# warn "$num_of_files out of $file_count .fastq files copied successfully....\n";

# print $fastq_header_counter . "\n";
sub find_fastq_files{
    
	my $project_leader_dir = shift;
	die "Error lost project leader directory" unless defined $project_leader_dir;

	my (%project_fastq_files, $file_count);
	$file_count = 0;
	opendir(DIR1, $project_leader_dir) || die "Error in opening dir $project_leader_dir\n";
	while( my $dir_name = readdir(DIR1)){
		if($dir_name =~ m/[\w_]+/){
			my $fastq_file_dir = join('/', $project_leader_dir, $dir_name);
			opendir(DIR2, $fastq_file_dir) || die "Error in opening dir $fastq_file_dir\n";
			while( my $file_name = readdir(DIR2)){
				my $fastq_file_name = join('/', $fastq_file_dir, $file_name) if ($file_name =~ m/.fastq$/);
				warn "$fastq_file_name\n" if ($file_name =~ m/.fastq$/);
				push(@{$project_fastq_files{$dir_name}}, $fastq_file_name) if ($file_name =~ m/.fastq$/);
				$file_count++ if ($file_name =~ m/.fastq$/);
			}
		closedir(DIR2);
		}
	}
	closedir(DIR1);
	return (\%project_fastq_files, $file_count);
}

# makeblastdb -in ncbi_nr_db_2014-05-30_organisms.fasta -dbtype 'nucl' -out ncbi_nr_db_2014-05-30_organisms.fasta
sub makeblastdb_nuc{
	my $fastadb = shift;
	die "Error lost fastadb to makeblastdb" unless defined $fastadb;

	# format the database file into .nin .nsq .nhr files.
	my ($fastadbNIN, $fastadbNSQ, $fastadbNHR);
	$fastadbNIN = $fastadb . '.nin';
	$fastadbNSQ = $fastadb . '.nsq';
	$fastadbNHR = $fastadb . '.nhr';
	unless(-s $fastadbNIN and -s $fastadbNSQ and -s $fastadbNHR){
		warn "Calling makeblastdb for $fastadb....\n";
		my $makeblastdbCmd = "$makeblastdb -in $fastadb -dbtype nucl";
		warn $makeblastdbCmd . "\n\n";
		system($makeblastdb, 
			'-in', $fastadb, 
			'-dbtype', 'nucl'
		) == 0 or die "Error calling $makeblastdbCmd: $?";
	}

}

sub generate_adaptor_blastn{
	my $fastq_adaptor_sequence = shift;
	die "Error lost fastq adaptor sequence" unless defined $fastq_adaptor_sequence;
	my $fasta_target = shift;
	die "Error lost fasta database target file" unless defined $fasta_target;
	my $blast_num_cpu = shift;
	die "Error lost number of cpus to allocate" unless defined $blast_num_cpu;
	my $blastn_output_dir = shift;
	die "Error lost adaptor blastn output directory" unless defined $blastn_output_dir;

	# Get the length of the adaptor so that we can verify adaptor blast results.
	my $adaptor_sequence_length = length($fastq_adaptor_sequence);
	
	my $fasta_target_name = fileparse($fasta_target);
	
	makeblastdb_nuc($fasta_target);
	my $blastn_tab_output_dir = join('/', $blastn_output_dir, "ADAPTOR_BLASTN_TAB_FILES");
	unless(-d $blastn_tab_output_dir){
		mkdir($blastn_tab_output_dir, 0777) or die "Can't make directory: $!";
	}
	my $adaptor_blastn_tsv_outfile = join('/', $blastn_tab_output_dir, $fasta_target_name . ".gbs_adaptor_blastn.tsv");
	unless(-s $adaptor_blastn_tsv_outfile){
		warn join(" ", "Generating blast tab-delimited file", join("", $fasta_target_name, ".gbs_adaptor_blastn", ".tsv"), "using adaptor blastn") . "\n";
		my $adaptorBlastnCmd  = "echo -e \"$fastq_adaptor_sequence\" | $blastn -query - -db $fasta_target -task blastn -word_size 7 -dust no -evalue 1000 -outfmt '6 qseqid salltitles qcovhsp pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
		warn $adaptorBlastnCmd . "\n\n";


		open(OUTFILE, ">$adaptor_blastn_tsv_outfile") or die "Couldn't open file $adaptor_blastn_tsv_outfile for writting, $!";
                print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch",
                "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score", "graphics_colour") . "\n";
		local (*ADAPTOR_BLASTN_OUT, *ADAPTOR_BLASTN_IN);
		my $pid = open2(\*ADAPTOR_BLASTN_OUT,\*ADAPTOR_BLASTN_IN, $adaptorBlastnCmd) or die "Error calling open2: $!";
		close ADAPTOR_BLASTN_IN or die "Error closing STDIN to adaptor blastn process: $!";
		while(<ADAPTOR_BLASTN_OUT>){
			chomp $_;
			my @adaptor_blastn_hit =  split(/\t/, $_);
			my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
			$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @adaptor_blastn_hit;
			my $glyph_colour = "blue";

	    		print OUTFILE join("\t", join("_", "GBS_adaptor_sequence", $fastq_adaptor_sequence), $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch, 
				$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, $glyph_colour) . "\n";

		
		}
		close ADAPTOR_BLASTN_OUT or die "Error closing STDOUT from adaptor blastn process: $!";
		wait;
		close(OUTFILE) or die "Couldn't close file $adaptor_blastn_tsv_outfile";
	}
	
	my $blastn_aln_output_dir = join('/', $blastn_output_dir, "ADAPTOR_BLASTN_ALIGN_FILES");
	unless(-d $blastn_aln_output_dir){
		mkdir($blastn_aln_output_dir, 0777) or die "Can't make directory: $!";
	}
	
	my $adaptor_blastn_aln_outfile = join('/', $blastn_aln_output_dir, $fasta_target_name . ".gbs_adaptor_blastn.aln.txt");
	unless(-s $adaptor_blastn_aln_outfile){
		warn "Generating primer blastn alignment file....\n";
		my $adaptorBlastnCmd  = "echo -e \"$fastq_adaptor_sequence\" | $blastn -query - -db $fasta_target -task blastn -word_size 7 -dust no -evalue 1000 -out $adaptor_blastn_aln_outfile -num_threads $blast_num_cpu";
		warn $adaptorBlastnCmd . "\n\n";

		my $status = system($adaptorBlastnCmd) == 0 or die "Error calling $blastn: $?";

	}
	return ($adaptor_blastn_tsv_outfile, $adaptor_blastn_aln_outfile, $adaptor_sequence_length);
}

sub generate_adaptor_blast_graphics{
	
	my $adaptor_blastn_tsv_outfile = shift;
	die "Error lost the tab-delimited adaptor blastn file" unless defined $adaptor_blastn_tsv_outfile ;
	my $adaptor_graphics_output_dir = shift;
	die "Error lost adaptor blastn graphics output directory" unless defined $adaptor_graphics_output_dir;

	
	warn "Calling gbs_adaptor_align_graphics.pl for $adaptor_blastn_tsv_outfile....\n";
	my $gbsAdaptorAlignGraphicsCmd  = "$gbs_adaptor_align_graphics -i $adaptor_blastn_tsv_outfile -o $adaptor_graphics_output_dir";
	warn $gbsAdaptorAlignGraphicsCmd . "\n\n";
	system($gbs_adaptor_align_graphics, 
		'-i', $adaptor_blastn_tsv_outfile, 
		'-o', $adaptor_graphics_output_dir
	) == 0 or die "Error calling $gbsAdaptorAlignGraphicsCmd: $?";
}

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

sub send_mail{

	my $subject = shift or die "lost email subject";
	my $message = shift or die "lost email message";
	my $email_type = shift or die "lost email type";
	
	my $sendMailCmd  = "$send_mail -s $subject -m $message -t $email_type";
	warn $sendMailCmd . "\n\n";
	system($sendMailCmd) == 0 or die "Error calling $sendMailCmd: $?";

}
