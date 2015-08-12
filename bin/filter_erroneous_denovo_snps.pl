#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

#### PROGRAM NAME ####
# trim_adapter_fastq_parallel_regex.pl - Program to trim the GBS common adapter sequence from each GBS fastq file within a particular Genotyping by Sequencing (GBS) project. Fixes the misprimming issue where the GBS common adapter is sequenced along with the DNA of an individual

#### DESCRIPTION ####
# This program trims the GBS common adapter sequence from each GBS fastq file within a particular Genotyping by Sequencing (GBS) project. Fixes the misprimming issue where the GBS common adapter is sequenced along with the DNA of an individual

#### SAMPLE COMMAND ####
# perl trim_adapter_fastq_parallel_regex.pl -i ~/workspace/GBS_data-08-10-2013/PstI_MspI_GBS_Data/PstI_MspI_PROCESSED_RADTAGS/PROJECT_LEADER_DIR_NO_MISMATCHES/CHRISTIANNE_MCDONALD -p CHRISTIANNE_MCDONALD -t 3 -q 32 -m 16 -l 92 -r PstI/MspI -c 16 -n true -o ~/workspace/GBS_data-08-10-2013/PstI_MspI_GBS_Data/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR
my ($fastq_file_dir, $stacks_fasta_infile, $project_name, $restriction_enzymes, $gbs_sequence_length, $adapter_length_min_threshold, $adapter_trim_offset, $min_trimmed_fastq_sequence_length, $regex_num_cpu, $pad_sequences, $gzip_files_switch, $output_dir);
GetOptions(
	'i=s'    => \$fastq_file_dir, # The *.fastq input file directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.

    'f=s'    => \$stacks_fasta_infile, # The Population Stacks write_single_snp fasta formatted input file.
	'p=s'    => \$project_name, # The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
#	'r=s'    => \$restriction_enzymes, # The restriction enzyme(s) used to digest the genomic sequences. Can be ApeKI, PstI/MspI, or SbfI/MspI. Default: PstI/MspI
	'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
#'m=s'    => \$adapter_length_min_threshold, # The minimum GBS common adapter sequence length cut-off in base pairs (bps) to retain for trimming if found in a given GBS fastq sequence hit found in the adapter regex searches. Default: 16
#'t=s'    => \$adapter_trim_offset, # The trimming offset length in base pairs (bps) to trim upstream of the start of the GBS common adapter sequence found in the adapter regex searches. Default: 5
	'q=s'    => \$min_trimmed_fastq_sequence_length, # The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Default: 32
	'n=s'    => \$pad_sequences, # The padded sequence controller. Specify true for padded trimmed sequences or false for unpadded trimmed sequences. Default: false
	'c=s'    => \$regex_num_cpu, # The number of cpu cores to use for the adapter regex searches. Default: 1
	's=s'    => \$gzip_files_switch, # The gzip compression switch. Specify true for compressing file or false for not compressing files (saves time). Default: false
	'o=s'    => \$output_dir, # The absolute path to the output directory to contain the trimmed adapter sequence fastq output files.
);

# Display a usage message if the following parameters are not specified.
usage() unless (
	defined $fastq_file_dir
    and defined $stacks_fasta_infile
	and defined $project_name
	and defined $output_dir
);


# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

# The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Keep in mind that this is the minimum trimmed fastq sequence length used before we add the length of the barcode 
# used for splitting each individual fastq file. Tassel required sequences at least 32 base pairs (bps) plus the length of a particular barcode that can be in the range of 4-8 base pairs (bps) in length. Default: 32
$min_trimmed_fastq_sequence_length = 32  unless defined $min_trimmed_fastq_sequence_length;

# The padded sequence controller. Specify true for padded trimmed sequences or false for unpadded trimmed sequences. Default: false
$pad_sequences = 'false' unless defined $pad_sequences;

# The gzip compression switch. Specify true for compressing file or false for not compressing files (saves time). Default: false
$gzip_files_switch = 'false' unless defined $gzip_files_switch;

# The number of cpu cores to use for the adapter regex searches. Default: 1
$regex_num_cpu = 1 unless defined $regex_num_cpu;

# Program dependencies - The absolute paths to gzip to compress all project leader *.fastq input files.
my $gzip 				= '/bin/gzip';

sub usage {

die <<"USAGE";

Usage: $0 -i fastq_file_dir -p project_name -r restriction_enzymes -l gbs_sequence_length -m adapter_length_min_threshold -t adapter_trim_offset -q min_trimmed_fastq_sequence_length -s gzip_files_switch -c regex_num_cpu -o output_dir

DESCRIPTION - This program trims the GBS common adapter sequence from each GBS fastq file within a particular Genotyping by Sequencing (GBS) project. Fixes the misprimming issue where the GBS common adapter is sequenced along with the DNA of an individual

OPTIONS:

-i fastq_file_dir - The *.fastq input file directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	e.g. /path/to/fastq_file_dir
	
-p project_name - The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	e.g. SBW_TUTORIAL
	
-l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92

-q min_trimmed_fastq_sequence_length - The minimum trimmed fastq sequence length in base pairs (bps) to retain after trimming. Keep in mind that this is the minimum trimmed fastq sequence length used before we add the 
length of the barcode used for splitting each individual fastq file. Tassel required sequences at least 32 base pairs (bps) plus the length of a particular barcode that can be in the range of 4-8 base pairs (bps) in 
length. Default: 32

-n pad_sequences - The padded sequence controller. Specify true for padded trimmed sequences or false for unpadded trimmed sequences. Default: false
    
-s gzip_files_switch - The gzip compression switch. Specify true for compressing file or false for not compressing files (saves time). Default: false
    
-c regex_num_cpu - The number of cpu cores to use for the adapter regex searches. Default: 1

-o output_dir - The absolute path to the output directory to contain the trimmed adapter sequence fastq output files.

	e.g. /path/to/output_dir
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my (%original_fastq_sequence_counter, %trimmed_fastq_sequence_counter) = ();
if($regex_num_cpu >= 2){

    # Perform the adapter regex searches in parallel.
    require Parallel::Loops;
    
    my $parallel = Parallel::Loops->new($regex_num_cpu);
    
    my %flagged_trimmed_fastq_seqs = ();
    my %flagged_trimmed_fastq_lengths = ();
    
    $parallel->share(\%flagged_trimmed_fastq_seqs, \%flagged_trimmed_fastq_lengths); # make sure that these are visible in the children.

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
		
		# Open the fastq input file for parsing.
		open(FASTQ_INFILE, "<$fastq_infile") or die "Couldn't open file $fastq_infile for reading, $!";
		# Parse the fastq files for the fastq header, sequence, plus, and quality scores and reformat to *.fasta format so that we can use adapter regex.
		my ($fastq_header, $fastq_sequence, $fastq_plus, $fastq_quality_scores);
		my $i = 1;
		my $num_fastq_seqs = 0;
		while(<FASTQ_INFILE>){
			chomp $_;
			#warn $_ . "\n";
			if(($_ =~ m/^@.+length=\d+$/)
				and ($i eq 1)){ # The fastq sequence header is on the first line.
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
				
                if($fastq_sequence =~ m/N{1,}$/){
                    $fastq_sequence =~ s/N/A/g;
                    push(@{$flagged_trimmed_fastq_seqs{$fastq_sequence}}, $fastq_header);
                    
                    my @split_fastq_header = split(/\001/, $fastq_header);
                    #die join("\t", @split_fastq_header);
                    my $length = $split_fastq_header[2];
                    $length =~ s/length=//g;
                    $flagged_trimmed_fastq_lengths{$fastq_sequence} = $length;
                }
				
				$i = 1;
				
			}else{
				$i++;
			}
		}
		close(FASTQ_INFILE) or die "Couldn't close file $fastq_infile";

	});

#	# Get the fastq untrimmed, trimmed, and removed counts.
#	foreach my $fasta_filename (sort keys %trimmed_fastq_seq_counter){
#		$trimmed_fastq_sequence_counter{$fasta_filename}{'UNTRIMMED'} = $trimmed_fastq_seq_counter{$fasta_filename}{'UNTRIMMED'};
#		$trimmed_fastq_sequence_counter{$fasta_filename}{'TRIMMED'} = $trimmed_fastq_seq_counter{$fasta_filename}{'TRIMMED'};
#		$trimmed_fastq_sequence_counter{$fasta_filename}{'REMOVED'} = $trimmed_fastq_seq_counter{$fasta_filename}{'REMOVED'};
#	}
	
#    foreach my $fastq_sequence (sort keys %flagged_trimmed_fastq_lengths){
#        die join("\t", @{$flagged_trimmed_fastq_lengths{$fastq_sequence}});
#    }
    
    # Parse the contents of the Stacks fasta file to obtain the reference genome sequence ids.
    open(INFILE, "<$stacks_fasta_infile") or die "Couldn't open file $stacks_fasta_infile for reading, $!";
    my %all_loci_sequences = ();
    my %flagged_loci_sequences = ();
    my %flagged_loci_seq_lengths = ();
    my ($locus_id, $allele_id, $individual_id, $locus_seq_header) = "";
    while(<INFILE>){
        chomp $_;
        
        if($_ =~ /^>/){
            my $fasta_header = $_;
            # Parse stacks fasta output file for the locus, allele number, individual name id, reference genome sequence id, SNP position, and reference genome sequence strand orientation.
            if($fasta_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+)\]/){
                #warn $fasta_header . "\n";
                ($locus_id, $allele_id, $individual_id, $locus_seq_header) = ($1, $2, $3, $fasta_header);
                
            }else{
                die "Error: $fasta_header is not in the correct format!";
            }
        }elsif($_ =~ m/[ACGTN]+/){ # The fasta sequence entry.
            my $sequence = $_;
            if(defined($flagged_trimmed_fastq_seqs{$sequence})){
                push(@{$flagged_loci_sequences{$locus_id}}, $sequence);
                $flagged_loci_seq_lengths{$locus_id} = $flagged_trimmed_fastq_lengths{$sequence};
            }
            $all_loci_sequences{$locus_seq_header} = $sequence;
        }
    }
    close(INFILE) or die "Couldn't close file $stacks_fasta_infile";
    
    my %filtered_locus_ids = ();
    foreach my $locus_id (sort keys %flagged_loci_sequences){
        if(defined($flagged_loci_sequences{$locus_id})){
            my $sequence1 = @{$flagged_loci_sequences{$locus_id}}[0];
            #            die $sequence1;
            my $flagged_loci_sequence_1_length = length($sequence1);
            my @flagged_loci_sequence_1 = split('', $sequence1);
            
            my $flagged_loci_size1 = scalar(@flagged_loci_sequence_1);
            for(my $i = 1; $i < $flagged_loci_size1; $i++){
                
                if(defined(@{$flagged_loci_sequences{$locus_id}}[$i])){
                    my $sequencei = @{$flagged_loci_sequences{$locus_id}}[$i];
                    my @flagged_loci_sequence_i = split('', $sequencei);
                
                    my $flagged_loci_sequence_i_length = length($sequencei);
                    
                    
#                    die $flagged_loci_sequence_i_length;
                    die "Error: $locus_id: flagged_loci_sequence_1_length=$flagged_loci_sequence_1_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($flagged_loci_sequence_1_length ne $gbs_sequence_length);
                    
                    die "Error: $locus_id: flagged_loci_sequence_i_length=$flagged_loci_sequence_i_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($flagged_loci_sequence_i_length ne $gbs_sequence_length);
                    
                    for(my $j = $flagged_loci_seq_lengths{$locus_id}; $j < $flagged_loci_sequence_i_length; $j++){
                        if($flagged_loci_sequence_1[$j] ne $flagged_loci_sequence_i[$j]){
                            
                            $filtered_locus_ids{$locus_id} = $locus_id;
                        }
                    }
                }
            }
        }
    }
    
    # Generate a fastq sequence counts file so that we can see the percentages of untrimmed, trimmed, and removed fastq sequences each fastq input file.
    my $filtered_snps_outfile = join("/", $output_dir, join("_", $project_name, "filtered_snps") . ".fasta");
    open(OUTFILE, ">$filtered_snps_outfile") or die "Couldn't open file $filtered_snps_outfile for writting, $!";
    foreach my $locus_seq_header (sort keys %all_loci_sequences){
        
        if($locus_seq_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+)\]/){
            #warn $fasta_header . "\n";
            ($locus_id, $allele_id, $individual_id) = ($1, $2, $3);
            if(!defined($filtered_locus_ids{$locus_id})){
                print OUTFILE join("\n", "$locus_seq_header", $all_loci_sequences{$locus_seq_header}) . "\n";
            }
        }
        
    }
    close(OUTFILE) or die "Couldn't close file $filtered_snps_outfile";
    
	# close parallel loops to free the memory contained in the shared hash variables.
	undef $parallel;
}elsif($regex_num_cpu eq 1){

    
}



# Compress the trimmed fastq files using gzip.
# Find all files in the specified directory with the extension *.fastq.
#if($gzip_files_switch eq "true"){
#    my ($trimmed_fastq_files, $trimmed_fastq_file_count) = find_files($trimmed_fastq_output_dir, "fastq");
#    foreach my $trimmed_fastq_filename (sort keys %{$trimmed_fastq_files}){
#        
#        # Get the full path to the trimmed fastq file.
#        my $trimmed_fastq_infile = $trimmed_fastq_files->{$trimmed_fastq_filename};
#
#        # Compress the trimmed fastq file using gzip.
#        gzip_file($trimmed_fastq_infile);
#    }
#}

#
## Compress the bulk trimmed sequence layout file using gzip if $gzip_files_switch eq "true".
#gzip_file($trimmed_seqs_layout_bulk_outfile) if($gzip_files_switch eq "true");

# Generate a fastq sequence counts file so that we can see the percentages of untrimmed, trimmed, and removed fastq sequences each fastq input file.
#my $fastq_seq_counts_outfile = join("/", $trimmed_output_dir, join("_", $project_name, "trimmed_offset", $adapter_trim_offset, "fastq_seq_counts") . ".txt");
#open(OUTFILE, ">$fastq_seq_counts_outfile") or die "Couldn't open file $fastq_seq_counts_outfile for writting, $!";
#print OUTFILE join("\t", "fastq_file_name", "total_num_seqs", "num_untrimmed_seqs", "percent_untrimmed_seqs", "num_trimmed_seqs", "percent_trimmed_seqs", "num_removed_seqs", "percent_removed_seqs") . "\n";
#foreach my $fasta_filename (sort keys %trimmed_fastq_sequence_counter){
#
#    # Get the fastq untrimmed, trimmed, and removed counts.
#    my $num_untrimmed_seqs = $trimmed_fastq_sequence_counter{$fasta_filename}{'UNTRIMMED'};
#    my $num_trimmed_seqs = $trimmed_fastq_sequence_counter{$fasta_filename}{'TRIMMED'};
#    my $num_removed_seqs = $trimmed_fastq_sequence_counter{$fasta_filename}{'REMOVED'};
#
#
#    # If trimmed or removed sequence counters how zero count reflect in number of sequence variables.
#    $num_untrimmed_seqs = 0 unless(defined($num_untrimmed_seqs));
#    $num_trimmed_seqs = 0 unless(defined($num_trimmed_seqs));
#    $num_removed_seqs = 0 unless(defined($num_removed_seqs));
#
#    # Calculate the total number of sequences.
#    my $total_num_seqs = $original_fastq_sequence_counter{$fasta_filename};
#    my $percent_untrimmed_seqs = (($num_untrimmed_seqs/$total_num_seqs) * 100);
#    my $percent_trimmed_seqs = (($num_trimmed_seqs/$total_num_seqs) * 100);
#    my $percent_removed_seqs = (($num_removed_seqs/$total_num_seqs) * 100);
#
#    # If trimmed or removed sequence counters how zero count reflect in either percent of sequence variables.
#    $percent_untrimmed_seqs = 0.00 unless(defined($num_untrimmed_seqs));
#    $percent_trimmed_seqs = 0.00 unless(defined($num_trimmed_seqs));
#    $percent_removed_seqs = 0.00 unless(defined($num_removed_seqs));
#
#    die "The number of sequences are not equal to 100.00 percent.....\n" . join("\t", $percent_untrimmed_seqs, $percent_trimmed_seqs, $percent_removed_seqs, ($percent_untrimmed_seqs + $percent_trimmed_seqs + $percent_removed_seqs)) if(($percent_untrimmed_seqs + $percent_trimmed_seqs + $percent_removed_seqs) ne 100.00);
#    print OUTFILE join("\t", $fasta_filename, $total_num_seqs, $num_untrimmed_seqs, $percent_untrimmed_seqs, $num_trimmed_seqs, $percent_trimmed_seqs, $num_removed_seqs, $percent_removed_seqs) . "\n";
#}
#close(OUTFILE) or die "Couldn't close file $fastq_seq_counts_outfile";
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
