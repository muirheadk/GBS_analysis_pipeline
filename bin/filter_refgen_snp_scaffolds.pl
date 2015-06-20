#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

#### PROGRAM NAME ####
# filter_refgen_snp_scaffolds.pl - Program that generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file.

#### DESCRIPTION ####
# This program generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file.

#### SAMPLE COMMAND ####
#
my ($stacks_fasta_infile, $stacks_refgen_infile, $output_dir);
GetOptions(
'i=s'    => \$stacks_fasta_infile, # The Population Stacks write_single_snp fasta formatted input file.
'r=s'    => \$stacks_refgen_infile, # The Stacks population map formatted input file.
'o=s'    => \$output_dir, # The directory to contain all the output files.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
defined $stacks_fasta_infile
and defined $stacks_refgen_infile
and defined $output_dir
);

sub usage {
    
    die <<"USAGE";
    
    
Usage: $0 -i stacks_fasta_infile -f stacks_refgen_infile -o output_dir
    
    DESCRIPTION - A Program that generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file.
    
OPTIONS:
    
    -i stacks_fasta_infile - The Population Stacks write_single_snp fasta formatted input file.
    
    -r stacks_refgen_infile - The Stacks reference genome fasta formatted input file.
    
    -o output_dir - The directory to contain all the output files.
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
    mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Parse the contents of the Stacks fasta file to obtain the reference genome sequence ids.
open(INFILE, "<$stacks_fasta_infile") or die "Couldn't open file $stacks_fasta_infile for reading, $!";
my %fasta_refgen_ids = ();
my ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
while(<INFILE>){
    chomp $_;
    
    if($_ =~ /^>/){
        my $fasta_header = $_;
        # Parse stacks fasta output file for the locus, allele number, individual name id, reference genome sequence id, SNP position, and reference genome sequence strand orientation.
        if($fasta_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+); (\d+), (\d+), (\+|\-)\]/){
            warn $fasta_header . "\n";
            ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = ($1, $2, $3, $4, $5, $6);
            
            $fasta_refgen_ids{$refgen_seq_id} = $refgen_seq_id;
            #die join("\t", $locus_id, $individual_id, $allele_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) . "\n";
        }else{
            die "Error: $fasta_header is not in the correct format!";
        }
    }elsif($_ =~ m/[ACGTN]+/){ # The fasta sequence entry.
        my $sequence = $_;
        #$locus_metadata{$locus_id}{$individual_id}{$allele_id} = $sequence;
        ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
    }
}
close(INFILE) or die "Couldn't close file $stacks_fasta_infile";

# Parse the contents of the Stacks reference genome formatted fasta file to filter sequences containing SNPs the from refgen_stacks_analysis_pipeline.pl script.
open(INFILE, "<$stacks_refgen_infile") or die "Couldn't open file $stacks_refgen_infile for reading, $!";
my $filtered_refgen_filename = fileparse($stacks_refgen_infile, qr/\.fna|\.fasta|\.fna\.fasta/);
my $filtered_refgen_outfile = join('/', $output_dir, "$filtered_refgen_filename-filtered.fasta");
# Print out the filtered sequence to the new reference genome filtered outfile.
open(OUTFILE, ">$filtered_refgen_outfile") or die "Couldn't open file $filtered_refgen_outfile for writting, $!";
$refgen_seq_id = "";
while(<INFILE>){
    chomp $_;
    
    if($_ =~ /^>/){
        my $fasta_header = $_;
        # Parse stacks reference genome fasta input file for the reference genome.
        if($fasta_header =~ m/^>(\d+)/){
            warn $fasta_header . "\n";
            $refgen_seq_id = $1;
            
            #die $refgen_seq_id . "\n";
        }else{
            die "Error: $fasta_header is not in the correct format!";
        }
    }elsif($_ =~ m/[ACGTN]+/){
        my $sequence = $_;
        print OUTFILE join("\n", ">$refgen_seq_id", $sequence) . "\n" if(defined($fasta_refgen_ids{$refgen_seq_id}));
        $refgen_seq_id = "";
    }else{
        die "Error: $stacks_refgen_infile is not in the correct format!";
    }
}
close(INFILE) or die "Couldn't close file $stacks_refgen_infile";
close(OUTFILE) or die "Couldn't close file $filtered_refgen_outfile";


