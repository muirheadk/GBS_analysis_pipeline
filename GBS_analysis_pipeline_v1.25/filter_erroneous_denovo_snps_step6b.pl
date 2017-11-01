#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

# VERSION 1.22

#### PROGRAM NAME ####
# filter_erroneous_denovo_snps.pl - A program that filters erroneous SNPs called within the denovo stacks analysis pipeline. SNPs called within the artificial poly-A sequence are flagged and filtered from the fasta formatted file that was generated from population stacks program. The program generates a filtered fasta formatted file and locus id list as output. It also generates a list of locus ids that were removed from the original fasta file and a fasta formatted file containing all locus sequences that were flagged and removed.

#### VERSION 1.20 ####

#### DESCRIPTION ####
# This program filters erroneous SNPs called within the denovo stacks analysis pipeline. SNPs called within the artificial poly-A sequence are flagged and filtered from the fasta formatted file that was generated from population stacks program. The program generates a filtered fasta formatted file and locus id list as output. It also generates a list of locus ids that were removed from the original fasta file and a fasta formatted file containing all locus sequences that were flagged and removed.

#### SAMPLE COMMAND ####
# perl ~/filter_erroneous_denovo_snps.pl -i ~/scratch/SBW20INDIVS -f ~/scratch/SBW20INDIVS/POPULATIONS_FASTA_DIR/batch_1.fa -p SBW20IndvsFiltered -l 92 -c 24 -o ~/output_dir
my ($fastq_file_dir, $stacks_fasta_infile, $project_name, $gbs_sequence_length, $trim_polyA, $num_cpu_cores, $output_dir);
GetOptions(
	'i=s'    => \$fastq_file_dir, # The absolute path to the trimmed *.fastq input file directory that contains files trimmed of the common adapter sequence with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	'f=s'    => \$stacks_fasta_infile, # The absolute path to the fasta formatted input file generated from the Stacks populations program.
	'p=s'    => \$project_name, # The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
	't=s'    => \$trim_polyA, # The trim poly-A sequence flag. If true, trim the poly-A sequence off the 3' end of the locus sequence of each individual. Otherwise, do not trim the poly-A sequence at a particular locus sequence. Default: true
	'c=s'    => \$num_cpu_cores, # The number of cpu cores to use for filtering the erroneous SNPs. Default: 1
	'o=s'    => \$output_dir, # The absolute path to the output directory to write the filtered and removed SNPs locus id lists and fasta files.
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

# The trim poly-A sequence flag. If true, trim the poly-A sequence off the 3' end of the locus sequence of each individual. Otherwise, do not trim the poly-A sequence at a particular locus sequence. Default: true
$trim_polyA = 'true' unless defined $trim_polyA;

# The number of cpu cores to use for filtering erroneous SNPs. Default: 1
$num_cpu_cores = 1 unless defined $num_cpu_cores;

sub usage {

die <<"USAGE";

Usage: $0 -i fastq_file_dir -f stacks_fasta_infile -p project_name -l gbs_sequence_length -t trim_polyA -c num_cpu_cores -o output_dir

VERSION 1.22

DESCRIPTION - This program filters erroneous SNPs called within the denovo stacks analysis pipeline. SNPs called within the artifical poly-A sequence are flagged and filtered from the fasta formatted file generated from population stacks program. The program generates a filtered fasta formatted file and locus id list as output. It also generates a locus id list of locus ids that were removed from the original fasta file and a fasta formatted file containing all locus sequences that were flagged and removed.

OPTIONS:

-i fastq_file_dir - The absolute path to the trimmed *.fastq input file directory that contains files trimmed of the common adapter sequence with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	e.g. /path/to/fastq_file_dir

-f stacks_fasta_infile - The absolute path to the fasta formatted input file generated from the Stacks populations program.

-p project_name - The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output directories and files with the specifed output directory.
	e.g. SBW_TUTORIAL
	
-l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92

-t trim_polyA - The trim poly-A sequence flag. If true, trim the poly-A sequence off the 3' end of the locus sequence of each individual. Otherwise, do not trim the poly-A sequence at a particular locus sequence. Default: true

-c num_cpu_cores - The number of cpu cores to use for filtering erroneous SNPs. Default: 1

-o output_dir - The absolute path to the output directory to write the filtered and removed SNPs locus id lists and fasta files.

	e.g. /path/to/output_dir
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my (%original_fastq_sequence_counter, %trimmed_fastq_sequence_counter) = ();
if($num_cpu_cores >= 1){

    # Perform the adapter regex searches in parallel.
    require Parallel::Loops;
    
    my $parallel = Parallel::Loops->new($num_cpu_cores);
    
    my %unique_trimmed_fastq_seqs = ();
    my %flagged_trimmed_fastq_lengths = ();
    
    $parallel->share(\%unique_trimmed_fastq_seqs, \%flagged_trimmed_fastq_lengths); # make sure that these are visible in the children.

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
                    push(@{$unique_trimmed_fastq_seqs{$fastq_sequence}}, $fastq_header);
                    
                    my @split_fastq_header = split(/\001/, $fastq_header);
                    #die join("\t", @split_fastq_header);
                    my $length = $split_fastq_header[2];
                    $length =~ s/length=//g;
                    push(@{$flagged_trimmed_fastq_lengths{$fastq_sequence}}, $length);
                    
                }
				
				$i = 1;
				
			}else{
				$i++;
			}
		}
		close(FASTQ_INFILE) or die "Couldn't close file $fastq_infile";

	});

    my %unique_trimmed_fasta_lengths = ();
    foreach my $fastq_sequence (sort keys %flagged_trimmed_fastq_lengths){
        
        my @flagged_sequence_lengths = ();
        push(@flagged_sequence_lengths, @{$flagged_trimmed_fastq_lengths{$fastq_sequence}});
        my @unique_sequence_lengths = do { my %seen; grep { !$seen{$_}++ } @flagged_sequence_lengths };
        
        my $uniq_seq_length_size = scalar(@unique_sequence_lengths);
        my $min_fastq_length = 0;
        if($uniq_seq_length_size > 1){
            foreach my $length (sort {$a <=> $b} @unique_sequence_lengths){
                $min_fastq_length = $length;
                last;
            }
        }elsif($uniq_seq_length_size eq 1){
            $min_fastq_length = $unique_sequence_lengths[0];
        }
        $unique_trimmed_fasta_lengths{$fastq_sequence} = $min_fastq_length;
    }
    
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
            if(defined($unique_trimmed_fastq_seqs{$sequence})){
                push(@{$flagged_loci_sequences{$locus_id}}, $sequence);
                push(@{$flagged_loci_seq_lengths{$locus_id}}, $unique_trimmed_fasta_lengths{$sequence});
            }
            $all_loci_sequences{$locus_seq_header} = $sequence;
        }
    }
    close(INFILE) or die "Couldn't close file $stacks_fasta_infile";
    
    my %unique_fastq_read_lengths = ();
    foreach my $locus_id (sort keys %flagged_loci_seq_lengths) {
        my $min_fastq_length = 0;
        foreach my $length (sort {$a <=> $b} @{$flagged_loci_seq_lengths{$locus_id}}){
            $min_fastq_length = $length;
            $unique_fastq_read_lengths{$locus_id} = $min_fastq_length;
            last;
        }
    }
    
    my %removed_locus_ids = ();
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
                    die "Error: $locus_id: flagged_loci_sequence_1_length=$flagged_loci_sequence_1_length bp ne flagged_loci_sequence_i_length=$flagged_loci_sequence_i_length bp" if($flagged_loci_sequence_1_length ne $flagged_loci_sequence_i_length);
                    
                    for(my $j = $unique_fastq_read_lengths{$locus_id}; $j < $flagged_loci_sequence_i_length; $j++){
                        if($flagged_loci_sequence_1[$j] ne $flagged_loci_sequence_i[$j]){
                            
                            $removed_locus_ids{$locus_id} = $locus_id;
                        }
                    }
                }
            }
        }
    }
    
	my (%filtered_locus_sequences, %removed_locus_sequences) = ();
    foreach my $locus_seq_header (sort keys %all_loci_sequences){
        
        if($locus_seq_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+)\]/){
            #warn $fasta_header . "\n";
            ($locus_id, $allele_id, $individual_id) = ($1, $2, $3);
            if(!defined($removed_locus_ids{$locus_id})){
				if(($trim_polyA eq "true") and defined($unique_fastq_read_lengths{$locus_id})){
					my $trimmed_sequence = get_subseq($all_loci_sequences{$locus_seq_header}, 1, ($unique_fastq_read_lengths{$locus_id} - 1));
					$filtered_locus_sequences{$locus_id}{$individual_id}{$allele_id} = join("\n", "$locus_seq_header", $trimmed_sequence);
				}else{
					$filtered_locus_sequences{$locus_id}{$individual_id}{$allele_id} = join("\n", "$locus_seq_header", $all_loci_sequences{$locus_seq_header});
				}
			}else{
                $removed_locus_sequences{$locus_id}{$individual_id}{$allele_id} = join("\n", "$locus_seq_header", $all_loci_sequences{$locus_seq_header});
           	}
        }
        
    }
	
	# Generate a fastq sequence counts file so that we can see the percentages of untrimmed, trimmed, and removed fastq sequences each fastq input file.
    my $filtered_snps_outfile = join("/", $output_dir, join("_", $project_name, "filtered_snps") . ".fasta");
    open(FILTERED_OUTFILE, ">$filtered_snps_outfile") or die "Couldn't open file $filtered_snps_outfile for writting, $!";
    foreach my $locus_id (sort keys %filtered_locus_sequences){
		foreach my $individual_id (sort keys %{$filtered_locus_sequences{$locus_id}}){
			foreach my $allele_id (sort keys %{$filtered_locus_sequences{$locus_id}{$individual_id}}){
					print FILTERED_OUTFILE $filtered_locus_sequences{$locus_id}{$individual_id}{$allele_id} . "\n";
			}
		}
	}
    close(FILTERED_OUTFILE) or die "Couldn't close file $filtered_snps_outfile";
	
	my $filtered_locus_id_outfile = join("/", $output_dir, join("_", $project_name, "filtered_locus_id_list") . ".txt");
    open(OUTFILE, ">$filtered_locus_id_outfile") or die "Couldn't open file $filtered_locus_id_outfile for writting, $!";
    foreach my $locus_id (sort keys %filtered_locus_sequences){
        print OUTFILE $locus_id . "\n";
    }
    close(OUTFILE) or die "Couldn't close file $filtered_locus_id_outfile";
	
    my $removed_snps_outfile = join("/", $output_dir, join("_", $project_name, "removed_snps") . ".fasta");
    open(REMOVED_OUTFILE, ">$removed_snps_outfile") or die "Couldn't open file $removed_snps_outfile for writting, $!";
    foreach my $locus_id (sort keys %removed_locus_sequences){
		foreach my $individual_id (sort keys %{$removed_locus_sequences{$locus_id}}){
			foreach my $allele_id (sort keys %{$removed_locus_sequences{$locus_id}{$individual_id}}){
					print REMOVED_OUTFILE $removed_locus_sequences{$locus_id}{$individual_id}{$allele_id} . "\n";
			}
		}
	}
    close(REMOVED_OUTFILE) or die "Couldn't close file $removed_snps_outfile";
    
	my $removed_locus_id_outfile = join("/", $output_dir, join("_", $project_name, "removed_locus_id_list") . ".txt");
    open(OUTFILE, ">$removed_locus_id_outfile") or die "Couldn't open file $removed_locus_id_outfile for writting, $!";
    foreach my $locus_id (sort keys %removed_locus_sequences){
        print OUTFILE $locus_id . "\n";
    }
    close(OUTFILE) or die "Couldn't close file $removed_locus_id_outfile";
    
	# Print all the unique flagged sequences that contain an erroneous SNP from each individual and the  
    my $unique_trimmed_seq_outfile = join("/", $output_dir, join("_", $project_name, "unique_sequence_list") . ".txt");
    open(OUTFILE, ">$unique_trimmed_seq_outfile") or die "Couldn't open file $unique_trimmed_seq_outfile for writting, $!";
    print OUTFILE join("\t", "fastq_sequence", "flagged_seq_length", "fastq_header_list") . "\n";
    foreach my $fastq_sequence (sort keys %unique_trimmed_fastq_seqs){
        print OUTFILE join("\t", $fastq_sequence, $unique_trimmed_fasta_lengths{$fastq_sequence}, join(", ", @{$unique_trimmed_fastq_seqs{$fastq_sequence}})) . "\n";
    }
    close(OUTFILE) or die "Couldn't close file $unique_trimmed_seq_outfile";
    
	# close parallel loops to free the memory contained in the shared hash variables.
	undef $parallel;
}

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
