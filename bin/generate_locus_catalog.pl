#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

#### PROGRAM NAME ####
# generate_locus_catalog.pl - Program that generates locus catalog files from a Stacks fasta formatted file from population stacks.

#### DESCRIPTION ####
# This program generates generates locus catalog files from a Stacks fasta formatted file.

#### SAMPLE COMMAND ####
# perl generate_locus_catalog.pl -i ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -f ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/batch_1.fa -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/NEXUS_STACKS_OUTFILES
my ($stacks_input_dir, $stacks_fasta_infile, $gbs_sequence_length, $output_dir);
GetOptions(
	'i=s'    => \$stacks_input_dir, # The directory that contains all the Stacks generated output files.
	'f=s'    => \$stacks_fasta_infile, # The Stacks fasta formatted input file.
	'l=s'    => \$gbs_sequence_length, # The GBS fasta sequence length in base pairs (bps) common to all GBS fasta sequences. Default: 92
	'o=s'    => \$output_dir, # The directory to contain all locus catalog output files.
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

# Find all stacks tags output files from the stacks output directory with the extension *.tags.tsv.
my ($stacks_tags_files, $stacks_tags_file_count) = find_files($stacks_input_dir, "tags.tsv");

# Iterate through each tags output file with extension *.tags.tsv and obtain each individual name.
my %sample_names = ();
foreach my $file_name (sort keys %{$stacks_tags_files}){
	if($file_name !~ m/batch_\d+\.catalog\.tags\.tsv/){ # If *.tags.tsv file does not match batch_*.tags.tags.tsv.
		my $stacks_tags_infile = $stacks_tags_files->{$file_name};

		# Get the basename of the tags filename without the .tags.tsv extension.
		my $stacks_tags_filename = fileparse($stacks_tags_infile, qr/\.tags\.tsv/);

		warn "Processing " . $stacks_tags_filename . ".....\n";
        $sample_names{$stacks_tags_filename} = $stacks_tags_filename;
	}
}

# Parse the contents of the Stacks fasta file to generate the nexus files.
my $seqio = Bio::SeqIO->new(-file => $stacks_fasta_infile, '-format' => 'Fasta');
my (%locus_metadata, %fasta_locus_seq_counts, %fasta_locus_refgen) = ();
while(my $seq_entry = $seqio->next_seq) {
    
    my $seq_id = $seq_entry->id;
    my $sequence = $seq_entry->seq;
    my $seq_desc = $seq_entry->desc;
    
    my $fasta_header = join(" ", $seq_id, $seq_desc);
    
    # Parse stacks fasta output file for the locus, sample, and allele number and the reference genome sequence id, SNP position, and reference genome sequence strand orientation.
    if($fasta_header =~ m/CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+); (\d+), (\d+), (\+|\-)\]/){
        warn $fasta_header . "\n";
        my ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = ($1, $2, $3, $4, $5, $6);
        
        $locus_metadata{$locus_id}{$individual_id}{$allele_id} = $sequence;
        $fasta_locus_refgen{$locus_id} = join("\t", $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand);
        $fasta_locus_seq_counts{$locus_id}++;
        #die join("\t", $locus_id, $individual_id, $allele_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) . "\n";
    }else{
        die "Error: $fasta_header is not in the correct format!";
    }
}

# Clean out the sequence I/O object.
$seqio = ();

# Iterate through each locus id to obtain the metadata for that locus and print out the contents in nexus file format.
my (%processed_loci, %snp_loci_positions) = ();
foreach my $locus_id (sort {$a <=> $b} keys %locus_metadata){
	# warn $locus_id . "\n";
		my @nexus_sequence_list = ();
		# Iterate through each sample id to grab the allele id for naming the file and printing out the metadata associated with that particular locus.
		foreach my $individual_id (sort {$a cmp $b} keys %{$locus_metadata{$locus_id}}){
            my @allele_sequences = ();
            my $allele_counts = 0;
			foreach my $allele_id (sort {$a <=> $b} keys %{$locus_metadata{$locus_id}{$individual_id}}){
#
				push(@nexus_sequence_list, join("  ", join("_", $sample_names{$individual_id}, "Allele", $allele_id), $locus_metadata{$locus_id}{$individual_id}{$allele_id}));
                
                
                push(@allele_sequences, $locus_metadata{$locus_id}{$individual_id}{$allele_id});
                $allele_counts++;
                
			}
            
            
            my @allele_0 = split('', $allele_sequences[0]);
            if($allele_counts > 1){
                die "Error: Allele_0 is not the same length as gbs_sequence_length=$gbs_sequence_length" if(length($allele_sequences[0]) ne $gbs_sequence_length);
                for(my $i = 1; $i < scalar(@allele_sequences); $i++){
                    
                    print join("\t", $locus_id, $individual_id, $i, $allele_sequences[$i]) . "\n";
                    
                    my @allele_i = split('', $allele_sequences[$i]);
                    die "Error: Allele_0 is not the same length as Allele_$i" if(length($allele_sequences[0]) ne length($allele_sequences[$i]));
                    die "Error: Allele_0 is not the same length as gbs_sequence_length=$gbs_sequence_length" if(length($allele_sequences[0]) ne $gbs_sequence_length);
                    
                    my %nucleotide_counts = ();
                    for(my $nuc_position = 0; $nuc_position < length($allele_sequences[$i]); $nuc_position++){
                        if($allele_0[$nuc_position] ne $allele_i[$nuc_position]){
                            
                            print "$allele_0[$nuc_position] ne $allele_i[$nuc_position] $nuc_position" . "\n";
                            my $concat_snps = join(",", $allele_0[$nuc_position], $allele_i[$nuc_position]);
                            push(@{$nucleotide_counts{$nuc_position}}, $concat_snps);
                        }
                        
                    }
                    
                    my @iupac_locus_sequence = ();
                    @iupac_locus_sequence = @allele_0;
                    my @snp_positions = ();
                    foreach my $nuc_position (sort keys %nucleotide_counts){
                        my @snp_nuc_chars = ();
                        foreach my $concat_snps (@{$nucleotide_counts{$nuc_position}}){
                            warn join("\t", $nuc_position, $concat_snps) . "\n";
                            my @split_concat_snps = split(/,/, $concat_snps);
                            push(@snp_nuc_chars, @split_concat_snps);
                        }
                        my @unique_snp_nuc_chars = do { my %seen; grep { !$seen{$_}++ } sort {$a cmp $b} @snp_nuc_chars };
                        my $unique_concat_snps = join(",", @unique_snp_nuc_chars);
                        
                        my $iupac_code = get_iupac_code($unique_concat_snps);
                        die "Error: iupac_code is undefined" if(!defined($iupac_code));
                        $iupac_locus_sequence[$nuc_position] = $iupac_code;
                        
                        push(@snp_positions, ($nuc_position + 1));
                    }
                    $snp_loci_positions{$locus_id}{$individual_id} = join(",", @snp_positions);
                    $processed_loci{$locus_id}{$individual_id} = join('', @iupac_locus_sequence);
    
                }
            }elsif($allele_counts eq 1){
                $processed_loci{$locus_id}{$individual_id} = join('', @allele_0);
                $snp_loci_positions{$locus_id}{$individual_id} = "N/A";
            }
		}
	
}

my $locus_catalog_outfile = join('/', $output_dir, "sample_locus_catalog.txt");
# Print out the locus metadata contents in nexus file format.
open(OUTFILE, ">$locus_catalog_outfile") or die "Couldn't open file $locus_catalog_outfile for writting, $!";
print OUTFILE join("\t", "locus_id", "refgen_seq_id", "refgen_seq_start", "strand", "individual_id", "snp_position", "sequence") . "\n";
foreach my $locus_id (sort {$a <=> $b} keys %locus_metadata){
    my @locus_catalog_sequence_list = ();
    foreach my $individual_id (sort {$a cmp $b} keys %sample_names){
        warn $individual_id . "\n";
        if(!defined($processed_loci{$locus_id}{$individual_id})){
            
            $processed_loci{$locus_id}{$individual_id} = '-' x $gbs_sequence_length;
            $snp_loci_positions{$locus_id}{$individual_id} = "N/A";
        }
        
		# Grab the reference genome sequence id and strand orientation and the SNP position in the sequence.
        push(@locus_catalog_sequence_list, join("\t", $locus_id, $fasta_locus_refgen{$locus_id}, $sample_names{$individual_id}, $snp_loci_positions{$locus_id}{$individual_id}, $processed_loci{$locus_id}{$individual_id}));
    }
    
    
    print OUTFILE join("\n", @locus_catalog_sequence_list) . "\n";
}
close(OUTFILE) or die "Couldn't close file $locus_catalog_outfile";

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

sub get_iupac_code{
    
    # The input nucleotide list.
	my $nuc_list = shift;
	die "Error lost input nucleotide list" unless defined $nuc_list;
    
    my %iupac_codes = (
        "A,G" => "R",
        "C,T" => "Y",
        "C,G" => "S",
        "A,T" => "W",
        "G,T" => "K",
        "A,C" => "M",
        "C,G,T" => "B",
        "A,G,T" => "D",
        "A,C,T" => "H",
        "A,C,G" => "V",
    );
 
    return $iupac_codes{$nuc_list};
}
