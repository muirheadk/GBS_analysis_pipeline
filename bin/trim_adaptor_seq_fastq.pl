#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;
# perl trim_adaptor_seq_fastq.pl -i ~/workspace/GBS_data-08-10-2013/PROJECT_LEADER_DIR -c 7 -o ~/workspace/GBS_data-08-10-2013/TRIM_ADAPTOR_SEQ_FASTQ_DIR

my ($project_leader_dir, $fastq_adaptor_sequence, $adaptor_length_threshold, $adaptor_trim_offset, $blast_num_cpu, $output_dir);
GetOptions(
      'i=s'    => \$project_leader_dir,
      'a=s'    => \$fastq_adaptor_sequence,
      't=s'    => \$adaptor_length_threshold,
      'm=s'    => \$adaptor_trim_offset,
      'c=s'    => \$blast_num_cpu,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $project_leader_dir
      and defined $output_dir
);

$blast_num_cpu = 2 unless defined $blast_num_cpu;
$fastq_adaptor_sequence = 'CCGAGATCGGAAGAGCGGGGACTTTAAGC' unless defined $fastq_adaptor_sequence;
$adaptor_length_threshold = 16 unless defined $adaptor_length_threshold;
$adaptor_trim_offset = 5 unless defined $adaptor_trim_offset;

my ($makeblastdb, $blastn, $gbs_adaptor_align_graphics, $send_mail);
$makeblastdb 			= '/usr/bin/makeblastdb';
$blastn				= '/usr/bin/blastn';
$gbs_adaptor_align_graphics	= '/TRIA-NetUtils/bin/gbs_adaptor_align_graphics.pl';
$send_mail			= '/TRIA-NetUtils/bin/send_mail.pl';

sub usage {

die <<"USAGE";

Usage: $0 -i project_leader_dir -a fastq_adaptor_sequence -t adaptor_length_threshold -m adaptor_trim_offset -c blast_num_cpu -o output_dir

Description - 

OPTIONS:

      -i project_leader_dir - 

      -a fastq_adaptor_sequence -

      -t adaptor_length_threshold - 

      -m adaptor_trim_offset - 

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
	
	
	foreach my $fastq_infile (sort @{$project_fastq_files->{$project_leader}}){
		my %fastq_sequences = ();

		warn $fastq_infile . "\n";
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
		
		open(INFILE, "<$fastq_infile") or die "Couldn't open file $fastq_infile for reading, $!";
		my $fasta_filename = fileparse($fastq_infile, qr/\.fastq/);
		my $fasta_target_outfile = join("/", $fasta_output_dir, $fasta_filename . ".fasta");
		open(OUTFILE, ">$fasta_target_outfile") or die "Couldn't open file $fasta_target_outfile for writting, $!" unless(-s $fasta_target_outfile);
		my ($fastq_header, $fastq_sequence, $fastq_plus, $fastq_quality_scores);
		my $i = 0;
		# my $fastq_header_counter = 0;
		while(<INFILE>){
			chomp $_;
		# 	warn $_ . "\n";
			$fastq_header = $_ if($i eq 0);
		# 	$fastq_header_counter++ if($i eq 0);
			$fastq_sequence = $_ if($i eq 1);
			$fastq_plus = $_ if($i eq 2);
			$fastq_quality_scores = $_ if($i eq 3);
			
			if($i eq 3){

				my $fasta_header = $fastq_header;
				$fasta_header =~ s/\s/_/g;
				my $fastq_sequence_length = length($fastq_sequence);
				$fasta_header = join("_", $fasta_header, "length=$fastq_sequence_length");

				$fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $fastq_sequence;
				$fastq_sequences{$fastq_header}{'PLUS'} = $fastq_plus;
				$fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $fastq_quality_scores;
				die if($fastq_sequence_length ne 100);

				print OUTFILE join("\n", join("", ">", $fasta_header), $fastq_sequence) . "\n" unless(-s $fasta_target_outfile);
			}
			
			if($i eq 3){
				$i = 0;
			}elsif(($i >= 0) and ($i < 3)){
				$i++;
			}
		}
		close(INFILE) or die "Couldn't close file $fastq_infile";
		close(OUTFILE) or die "Couldn't close file $fasta_target_outfile" unless(-s $fasta_target_outfile);

		my ($adaptor_blastn_tsv_outfile, $adaptor_blastn_aln_outfile, $adaptor_sequence_length) = generate_adaptor_blastn($fastq_adaptor_sequence, $fasta_target_outfile, $blast_num_cpu, $blastn_output_dir);

		# Create output directory if it doesn't already exist.
		my $adaptor_graphics_output_dir = join('/', $blastn_output_dir, "ADAPTOR_BLASTN_GRAPHICS_FILES");
		unless(-d $adaptor_graphics_output_dir){
			mkdir($adaptor_graphics_output_dir, 0777) or die "Can't make directory: $!";
		}

		generate_adaptor_blast_graphics($adaptor_blastn_tsv_outfile, $adaptor_graphics_output_dir);

		# Create output directory if it doesn't already exist.
		my $trimmed_output_dir = join('/',$project_leader_dir, "TRIMMED_OUTPUT_FILES");
		unless(-d $trimmed_output_dir){
			mkdir($trimmed_output_dir, 0777) or die "Can't make directory: $!";
		}

		# Create output directory if it doesn't already exist.
		my $trimmed_blastn_output_dir = join('/', $trimmed_output_dir, "TRIMMED_ADAPTOR_BLASTN_FILES");
		unless(-d $trimmed_blastn_output_dir){
			mkdir($trimmed_blastn_output_dir, 0777) or die "Can't make directory: $!";
		}

		
		# Parse the tab-delimited adaptor blast output so that we can trimm the fastq sequences that contain the GBS adaptor sequence.
		# Create new adaptor blastn files so that we can visualize where we trimmed the sequence.
		my $trimmed_adaptor_blastn_outfile = join('/', $trimmed_blastn_output_dir, $fasta_filename . ".gbs_adaptor_blastn.tsv");
		open(INFILE, "<$adaptor_blastn_tsv_outfile") or die "Couldn't open file $adaptor_blastn_tsv_outfile for reading, $!";
		open(OUTFILE, ">$trimmed_adaptor_blastn_outfile") or die "Couldn't open file $trimmed_adaptor_blastn_outfile for writting, $!" unless(-s $trimmed_adaptor_blastn_outfile);

		print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch", 
		"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score", "graphics_colour") . "\n" unless(-s $trimmed_adaptor_blastn_outfile); 
		$i = 0;
		while(<INFILE>){
			chomp $_;
			if($i ne 0){

		 		warn $_ . "\n";
				
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
				my $fastq_quality_scores = $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
				#die $fastq_sequence;
				my $fastq_sequence_length = length($fastq_sequence);
				my $predicted_adaptor_start = ($fastq_sequence_length - $adaptor_sequence_length);
				my ($threshold_coverage, $trimmed_fastq_sequence, $trimmed_fastq_quality_scores, $trimmed_adaptor_sequence);
				$threshold_coverage = (($adaptor_length_threshold/$adaptor_sequence_length) * 100);
				#die join("\t", "query_coverage=$query_coverage", "threshold_coverage=$threshold_coverage");
				if(($target_start >= $predicted_adaptor_start) and ($target_end <= $fastq_sequence_length) 
					and ($query_coverage >= $threshold_coverage)){
					$trimmed_fastq_sequence = get_subseq($fastq_sequence, 1, (($target_start - $adaptor_trim_offset) - 1));
					$trimmed_fastq_quality_scores = get_subseq($fastq_quality_scores, 1, (($target_start - $adaptor_trim_offset) - 1));
					my $trimmed_fastq_sequence_length = length($trimmed_fastq_sequence);
					$fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'} = $trimmed_fastq_sequence;
					$fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'} = $trimmed_fastq_quality_scores;
					$trimmed_adaptor_sequence = get_subseq($fastq_sequence, ($target_start - $adaptor_trim_offset), $fastq_sequence_length);
					my $trimmed_fastq_sequence_length = length($trimmed_fastq_sequence);
					my $trimmed_adaptor_sequence_length = length($trimmed_adaptor_sequence);
					my $trimmed_fastq_blastn = join("\t", join("_", "trimmed_fastq_sequence_offset", $adaptor_trim_offset), $target_name, $query_coverage, $percent_identity, $trimmed_fastq_sequence_length, $num_mismatch, $num_gaps, 1, (($target_start - $adaptor_trim_offset) - 1), 1, (($target_start - $adaptor_trim_offset) - 1), $e_value, $bit_score, "blue");
					my $trimmed_adaptor_blastn = join("\t", join("_", "trimmed_adaptor_sequence_offset", $adaptor_trim_offset), $target_name, $query_coverage, $percent_identity, $trimmed_adaptor_sequence_length, $num_mismatch, $num_gaps, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, ($target_start - $adaptor_trim_offset), $fastq_sequence_length, $e_value, $bit_score, "red");
					my $original_adaptor_sequence = join("\t", $query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
				$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, "orange");
					print OUTFILE $trimmed_fastq_blastn . "\n" unless(-s $trimmed_adaptor_blastn_outfile);
					print OUTFILE $original_adaptor_sequence . "\n" unless(-s $trimmed_adaptor_blastn_outfile);
					print OUTFILE $trimmed_adaptor_blastn . "\n" unless(-s $trimmed_adaptor_blastn_outfile);

				}
			}
			$i++;

		}
		close(INFILE) or die "Couldn't close file $adaptor_blastn_tsv_outfile";
		close(OUTFILE) or die "Couldn't close file $trimmed_adaptor_blastn_outfile" unless(-s $trimmed_adaptor_blastn_outfile);

		# Create output directory if it doesn't already exist.
		my $trimmed_adaptor_graphics_output_dir = join('/', $trimmed_output_dir, "TRIMMED_ADAPTOR_BLASTN_GRAPHICS_FILES");
		unless(-d $trimmed_adaptor_graphics_output_dir){
			mkdir($trimmed_adaptor_graphics_output_dir, 0777) or die "Can't make directory: $!";
		}
		generate_adaptor_blast_graphics($trimmed_adaptor_blastn_outfile, $trimmed_adaptor_graphics_output_dir);

		# Create output directory if it doesn't already exist.
		my $trimmed_fastq_output_dir = join('/', $trimmed_output_dir, "TRIMMED_FASTQ_FILES");
		unless(-d $trimmed_fastq_output_dir){
			mkdir($trimmed_fastq_output_dir, 0777) or die "Can't make directory: $!";
		}

		my $trimmed_fastq_outfile = join("/", $trimmed_fastq_output_dir, $fasta_filename . ".fastq");
		open(OUTFILE, ">$trimmed_fastq_outfile") or die "Couldn't open file $trimmed_fastq_outfile for writting, $!" unless(-s $trimmed_fastq_outfile);
		foreach my $fastq_header (keys %fastq_sequences){
			
			my $fastq_sequence = $fastq_sequences{$fastq_header}{'FASTQ_SEQUENCE'};
			my $fastq_plus = $fastq_sequences{$fastq_header}{'PLUS'};
			my $fastq_quality_scores = $fastq_sequences{$fastq_header}{'FASTQ_QUALITY_SCORES'};
			my $fastq_sequence_length = length($fastq_sequence);

			my $new_fastq_header = join("_", $fastq_header, join("=", "length", $fastq_sequence_length));
			print OUTFILE  $new_fastq_header . "\n" unless(-s $trimmed_fastq_outfile);
			print OUTFILE  $fastq_sequence . "\n" unless(-s $trimmed_fastq_outfile);
			print OUTFILE  $fastq_plus . "\n" unless(-s $trimmed_fastq_outfile);
			print OUTFILE  $fastq_quality_scores . "\n" unless(-s $trimmed_fastq_outfile);

		}
		close(OUTFILE) or die "Couldn't close file $trimmed_fastq_outfile" unless(-s $trimmed_fastq_outfile);
		%fastq_sequences = ();
	}
}

my $subject = "$0 Process Complete";
my $message = "$0 process finished successfully! You can find all the trimmed adaptor fastq output files for each project leader in the $output_dir directory.";
send_mail($subject, $message, 'email');

# warn "$num_of_files out of $file_count .fastq files copied successfully....\n";

# print $fastq_header_counter . "\n";

=head1 SUBROUTINES
 
(\%music_files, $file_count) = find_music_files($ipod_music_dir) - Find all music files in the specified iPod music folder with the extension *.m4a.

Input paramater(s):
 
$ipod_music_dir - iPod music input directory.
 
Output paramater(s):

\%music_files - A hash reference containing all the files with file extension *.m4a in key/value pairs.
 
 key => filename [A-Z]{4}.m4a ( e.g. ABCD.m4a )
 value => absolue filepath /path/to/[A-Z]{4} ( e.g. /path/to/ABCD.m4a )
 
$file_count - The number of files stored with file extension *.m4a.
 
=cut
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
		warn "$makeblastdb -in $fastadb -dbtype nucl\n\n";
		system($makeblastdb, 
			'-in', $fastadb, 
			'-dbtype', 'nucl'
		) == 0 or die "Error calling $makeblastdb -in $fastadb -dbtype nucl: $?";
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

		my @adaptor_blastn_results = ();
		local (*ADAPTOR_BLASTN_OUT, *ADAPTOR_BLASTN_IN);
		my $pid = open2(\*ADAPTOR_BLASTN_OUT,\*ADAPTOR_BLASTN_IN, $adaptorBlastnCmd) or die "Error calling open2: $!";
		close ADAPTOR_BLASTN_IN or die "Error closing STDIN to adaptor blastn process: $!";
		while(<ADAPTOR_BLASTN_OUT>){
			chomp $_;
			my @adaptor_blastn_hit =  split(/\t/, $_);
			my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
			$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @adaptor_blastn_hit;
			my $glyph_colour = "blue";

	    		my $adaptor_blastn_entry = join("\t", join("_", "GBS_adaptor_sequence", $fastq_adaptor_sequence), $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch, 
				$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, $glyph_colour);
# 			warn $adaptor_blastn_entry . "\n";

	    		push(@adaptor_blastn_results, [split(/\t/, $adaptor_blastn_entry)]);
		
		}
		close ADAPTOR_BLASTN_OUT or die "Error closing STDOUT from adaptor blastn process: $!";
		wait;
		
		open(OUTFILE, ">$adaptor_blastn_tsv_outfile") or die "Couldn't open file $adaptor_blastn_tsv_outfile for writting, $!";
		print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch", 
		"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score", "graphics_colour") . "\n"; 
		foreach my $adaptor_blastn_entry (sort {$b->[12] <=> $a->[12]} @adaptor_blastn_results){
			print OUTFILE join("\t", @$adaptor_blastn_entry) . "\n";
		}
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
	warn "$gbs_adaptor_align_graphics -i $adaptor_blastn_tsv_outfile -o $adaptor_graphics_output_dir\n\n";
	system($gbs_adaptor_align_graphics, 
		'-i', $adaptor_blastn_tsv_outfile, 
		'-o', $adaptor_graphics_output_dir
	) == 0 or die "Error calling $gbs_adaptor_align_graphics -i $adaptor_blastn_tsv_outfile -o $adaptor_graphics_output_dir: $?";
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
    
    
		warn "$send_mail -s $subject -m $message -t $email_type\n\n";
		system($send_mail, 
			'-s', "\"$subject\"", 
			'-m', "\"$message\"", 
			'-t', $email_type
		) == 0 or die "Error calling $send_mail -s $subject -m $message -t $email_type: $?";

}
