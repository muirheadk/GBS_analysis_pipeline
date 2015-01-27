#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

#### PROGRAM NAME ####
# convert_stacks_fasta2nexus.pl - Program that generates files in NEXUS format from a Stacks fasta formatted file.

#### DESCRIPTION ####
# This program generates files in NEXUS format from a Stacks fasta formatted file.

#### SAMPLE COMMAND ####
# perl convert_stacks_fasta2nexus.pl -i ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -f ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/batch_1.fa -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/NEXUS_STACKS_OUTFILES
my ($stacks_input_dir, $stacks_fasta_infile, $gbs_sequence_length, $output_dir);
GetOptions(
	'i=s'    => \$stacks_input_dir, # The directory that contains all the Stacks generated output files.
	'f=s'    => \$stacks_fasta_infile, # The Stacks fasta formatted input file.
	'l=s'    => \$gbs_sequence_length, # The GBS fasta sequence length in base pairs (bps) common to all GBS fasta sequences. Default: 92
	'o=s'    => \$output_dir, # The directory to contain all the nexus output files.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
	defined $stacks_input_dir
	and defined $stacks_fasta_infile
	and defined $output_dir
);

# The GBS fasta sequence length in base pairs (bps) common to all GBS fasta sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

sub usage {

die <<"USAGE";


Usage: $0 -i stacks_input_dir -f stacks_fasta_infile -l gbs_sequence_length -o output_dir

DESCRIPTION - This program generates files in NEXUS format from a Stacks fasta formatted file.

OPTIONS:

-i stacks_input_dir - The directory that contains all the Stacks generated output files.

-f stacks_fasta_infile - The Stacks fasta formatted input file.

-l gbs_sequence_length - The GBS fasta sequence length in base pairs (bps) common to all GBS fasta sequences. Default: 92

-o output_dir - The directory to contain all the nexus output files.

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Find all stacks alleles output files from the stacks output directory with the extension *.alleles.tsv.
my ($stacks_alleles_files, $stacks_alleles_file_count) = find_files($stacks_input_dir, "alleles.tsv");

# Iterate through each alleles output file with extension *.alleles.tsv and obtain the sample id number.r to reference the individual name.
my %sample_names = ();
foreach my $file_name (sort keys %{$stacks_alleles_files}){
	if($file_name !~ m/batch_\d+\.catalog\.alleles\.tsv/){ # If *.alleles.tsv file does not match batch_*.alleles.tags.tsv.
		my $stacks_alleles_infile = $stacks_alleles_files->{$file_name};

		# Get the basename of the alleles filename without the .alleles.tsv extension.
		my $stacks_alleles_filename = fileparse($stacks_alleles_infile, qr/\.alleles\.tsv/);

		warn "Processing " . $stacks_alleles_filename . ".....\n";
		
		# Grab the first line of the stacks alleles file and parse for the sample id number.
		open(INFILE, "<$stacks_alleles_infile") or die "Couldn't open file $stacks_alleles_infile for reading, $!";
		my $line = <INFILE>;
		chomp($line);
		my @split_alleles_entry = split(/\t/, $line);
		my $sample_id = $split_alleles_entry[1];

		my ($sample_name, $project_id) = split(/_/, $stacks_alleles_filename);
		$sample_names{$sample_id} = $sample_name;
		close(INFILE) or die "Couldn't close file $stacks_alleles_infile";
	}
}

# Parse the contents of the Stacks fasta file to generate the nexus files.
my $seqio = Bio::SeqIO->new(-file => $stacks_fasta_infile, '-format' => 'Fasta');
my (%fasta_locus_nexus, %fasta_locus_seq_counts, %fasta_locus_refgen) = ();
while(my $seq_entry = $seqio->next_seq) {
    
    my $seq_id = $seq_entry->id;
    my $sequence = $seq_entry->seq;
    my $seq_desc = $seq_entry->desc;
    
    my $fasta_header = join(" ", $seq_id, $seq_desc);
    
    # Parse stacks fasta output file for the locus, sample, and allele number and the reference genome sequence id, SNP position, and reference genome sequence strand orientation.
    if($fasta_header =~ m/CLocus_(\d+)_Sample_(\d+)_Locus_\d+_Allele_(\d+) \[(\d+), (\d+), (\+|\-)\]/){
        warn $fasta_header . "\n";
        my ($locus_id, $sample_id, $allele_id, $refgen_seq_id, $snp_position, $refgen_seq_strand) = ($1, $2, $3, $4, $5, $6);
        
        $fasta_locus_nexus{$locus_id}{$sample_id}{$allele_id} = $sequence;
        $fasta_locus_refgen{$locus_id} = join("\t", $refgen_seq_id, $snp_position, $refgen_seq_strand);
        $fasta_locus_seq_counts{$locus_id}++;
        #die join("\t", $locus_id, $sample_id, $allele_id, $refgen_seq_id, $snp_position, $refgen_seq_strand) . "\n";
    }else{
        
        die "Error: $fasta_header is mot in the correct format!";
    }
}

# Clean out the sequence I/O object.
$seqio = ();

# Iterate through each locus id to obtain the metadata for that locus and print out the contents in nexus file format.
foreach my $locus_id (sort {$a <=> $b} keys %fasta_locus_nexus){
	# warn $locus_id . "\n";
	if($fasta_locus_seq_counts{$locus_id} > 1){ # Print nexus file if there is more than one sequence present.
		my @nexus_sequence_list = ();
		# Iterate through each sample id to grab the allele id for naming the file and printing out the metadata associated with that particular locus.
		foreach my $sample_id (sort {$a <=> $b} keys %{$fasta_locus_nexus{$locus_id}}){
			foreach my $allele_id (sort {$a <=> $b} keys %{$fasta_locus_nexus{$locus_id}{$sample_id}}){
				push(@nexus_sequence_list, join("  ", join("_", $sample_names{$sample_id}, "Allele", $allele_id), $fasta_locus_nexus{$locus_id}{$sample_id}{$allele_id}));
			}
		}

		# Grab the reference genome sequence id and strand orientation and the SNP position in the sequence.
		my ($refgen_seq_id, $snp_position, $refgen_seq_strand) = split(/\t/, $fasta_locus_refgen{$locus_id});

		my $refgen_strand = "";
		$refgen_strand = "Plus" if($refgen_seq_strand eq "+"); # If the sequence is in the sense strand orientation substitute "+" for "plus".
		$refgen_strand = "Minus" if($refgen_seq_strand eq "-"); # If the sequence is in the antisense strand orientation substitute "-" for "minus".

		# Print out the locus metadata contents in nexus file format.
		my $nexus_outfile = join('/', $output_dir, join("_", "locus$locus_id", "refSeq$refgen_seq_id", "pos$snp_position", "strand$refgen_strand") . ".nex");
		open(OUTFILE, ">$nexus_outfile") or die "Couldn't open file $nexus_outfile for writting, $!";
		print OUTFILE "#NEXUS" . "\n";
		print OUTFILE "Begin data;" . "\n";
		print OUTFILE "Dimensions ntax=$fasta_locus_seq_counts{$locus_id} nchar=$gbs_sequence_length;" . "\n";
		print OUTFILE "Format datatype=dna symbols=\"ACTG\" missing=N gap=-;" . "\n";
		print OUTFILE "Matrix" . "\n";
		print OUTFILE join("\n", @nexus_sequence_list) . "\n";
		print OUTFILE ";" . "\n";
		print OUTFILE "End;" . "\n";
		close(OUTFILE) or die "Couldn't close file $nexus_outfile";
	}
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

