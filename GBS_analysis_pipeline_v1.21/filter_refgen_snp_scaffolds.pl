#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

# VERSION 1.21

#### PROGRAM NAME ####
# filter_refgen_snp_scaffolds.pl - Program that generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file. You can also use a whitelist of locus ids to get only the reference genome scaffolds that pass the present data filter in convert_fasta2nexus.pl. Also prints out the populations fasta sequences that passed the present data filter to an output file.

#### DESCRIPTION ####
# This program generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file. You can also use a whitelist of locus ids to get only the reference genome scaffolds that pass the present data filter in convert_fasta2nexus.pl. Also prints out the populations fasta sequences that passed the present data filter to an output file.

#### SAMPLE COMMAND ####
#
my ($stacks_fasta_infile, $stacks_refgen_infile, $whitelist_infile, $output_dir);
GetOptions(
    'i=s'    => \$stacks_fasta_infile, # The Population Stacks write_single_snp fasta formatted input file.
    'r=s'    => \$stacks_refgen_infile, # The Stacks population map formatted input file.

    'w=s'    => \$whitelist_infile, # The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the filter reference genome fasta file.
    'o=s'    => \$output_dir, # The directory to contain all the output files.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
    defined $stacks_fasta_infile
    and defined $stacks_refgen_infile
    and defined $output_dir
);

# The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the filter reference genome fasta file.
$whitelist_infile = "" unless defined $whitelist_infile;

sub usage {
    
    die <<"USAGE";
    
    
Usage: $0 -i stacks_fasta_infile -f stacks_refgen_infile -w whitelist_infile -o output_dir

VERSION 1.21

DESCRIPTION - A Program that generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file. You can also used a whitelist of locus ids to get only the reference genome scaffolds that pass the present data filter in convert_fasta2nexus.pl. Also prints out the populations fasta sequences that passed the present data filter to an output file.
    
OPTIONS:
    
    -i stacks_fasta_infile - The Population Stacks write_single_snp fasta formatted input file.
    
    -r stacks_refgen_infile - The Stacks reference genome fasta formatted input file.
    
    -w whitelist_infile - The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the filter reference genome fasta file.
    
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
my %filtered_fasta_seqs = ();
my ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
while(<INFILE>){
    chomp $_;
    
    if($_ =~ /^>/){
        my $fasta_header = $_;
        # Parse stacks fasta output file for the locus, allele number, individual name id, reference genome sequence id, SNP position, and reference genome sequence strand orientation.
        if($fasta_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+); (\d+), (\d+), (\+|\-)\]/){
            warn $fasta_header . "\n";
            ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = ($1, $2, $3, $4, $5, $6);
            
            $fasta_refgen_ids{$locus_id} = $refgen_seq_id;
            
            push(@{$filtered_fasta_seqs{$locus_id}{$individual_id}{$allele_id}}, "$fasta_header");
            #die join("\t", $locus_id, $individual_id, $allele_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) . "\n";
        }else{
            die "Error: $fasta_header is not in the correct format!";
        }
    }elsif($_ =~ m/[ACGTN]+/){ # The fasta sequence entry.
        my $sequence = $_;
        push(@{$filtered_fasta_seqs{$locus_id}{$individual_id}{$allele_id}}, $sequence);
        ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
    }
}
close(INFILE) or die "Couldn't close file $stacks_fasta_infile";

# Adding locus sequences in a hash array in order to concatenate the sequences.
my %filtered_refgen_ids = ();
my @locus_id_list = ();
if($whitelist_infile eq ""){
    foreach my $locus_id (sort {$a <=> $b} keys %fasta_refgen_ids){
        my $refgen_seq_id = $fasta_refgen_ids{$locus_id};
        $filtered_refgen_ids{$refgen_seq_id} = $refgen_seq_id;
        
    }
}elsif(-s $whitelist_infile){
    
    open(INFILE, "<$whitelist_infile") or die "Couldn't open file $whitelist_infile for reading, $!";
    while(<INFILE>){
        chomp $_;
        
        my $locus_id = $_;
        
        if(defined($fasta_refgen_ids{$locus_id})){
            my $refgen_seq_id = $fasta_refgen_ids{$locus_id};
            $filtered_refgen_ids{$refgen_seq_id} = $refgen_seq_id;
            push(@locus_id_list, $locus_id);
        }
        
    }
    close(INFILE) or die "Couldn't close file $whitelist_infile";
    
    # Print out the fasta sequences to the filtered populations stacks outfile.
    my $filtered_pop_fasta_filename = fileparse($stacks_fasta_infile, qr/\.fa|\.fna|\.fasta|\.fna\.fasta/);
    my $filtered_pop_fasta_outfile = join('/', $output_dir, join("-", $filtered_pop_fasta_filename, "filtered.fasta"));
    open(OUTFILE, ">$filtered_pop_fasta_outfile") or die "Couldn't open file $filtered_pop_fasta_outfile for writting, $!";
    foreach my $locus_id (sort {$a <=> $b} @locus_id_list){
        foreach my $individual_id (sort {$a cmp $b} keys %{$filtered_fasta_seqs{$locus_id}}){
            foreach my $allele_id (sort {$a <=> $b} keys %{$filtered_fasta_seqs{$locus_id}{$individual_id}}){
                print OUTFILE join("\n", @{$filtered_fasta_seqs{$locus_id}{$individual_id}{$allele_id}}) . "\n";
            }
        }
    }
    close(OUTFILE) or die "Couldn't close file $filtered_pop_fasta_outfile";
    
}else{
    die "$whitelist_infile does not exist";
}


# Parse the contents of the Stacks reference genome formatted fasta file to filter sequences containing SNPs the from refgen_stacks_analysis_pipeline.pl script.
open(INFILE, "<$stacks_refgen_infile") or die "Couldn't open file $stacks_refgen_infile for reading, $!";
my $filtered_refgen_filename = fileparse($stacks_refgen_infile, qr/\.fna|\.fasta|\.fna\.fasta/);
my $filtered_refgen_outfile = join('/', $output_dir, "$filtered_refgen_filename-filtered.fasta");
# Print out the filtered sequence to the new reference genome filtered outfile.
open(OUTFILE, ">$filtered_refgen_outfile") or die "Couldn't open file $filtered_refgen_outfile for writting, $!";
my $refgen_id = "";
while(<INFILE>){
    chomp $_;
    # Parse stacks reference genome fasta input file for the reference genome.
    if($_ =~ /^>/){
        # Get the header line without the ">" character.
        my $header_length = length($_);
        $refgen_id = substr($_, 1, $header_length);

    }elsif($_ =~ m/[ACGTNRYSWKM]+/){
        my $sequence = $_;
        print OUTFILE join("\n", ">$refgen_id", $sequence) . "\n" if(defined($filtered_refgen_ids{$refgen_id}));
        $refgen_seq_id = "";
    }else{
        die "Error: $stacks_refgen_infile is not in the correct format! $_";
    }
}
close(INFILE) or die "Couldn't close file $stacks_refgen_infile";
close(OUTFILE) or die "Couldn't close file $filtered_refgen_outfile";


