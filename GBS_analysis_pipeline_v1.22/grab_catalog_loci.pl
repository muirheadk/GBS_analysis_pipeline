#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

# VERSION 1.22

#### PROGRAM NAME ####
# grab_catalog_loci.pl - Program that generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file. You can also use a whitelist of locus ids to get only the reference genome scaffolds that pass the present data filter in convert_fasta2nexus.pl. Also prints out the populations fasta sequences that passed the present data filter to an output file.

#### DESCRIPTION ####
# This program generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file. You can also use a whitelist of locus ids to get only the reference genome scaffolds that pass the present data filter in convert_fasta2nexus.pl. Also prints out the populations fasta sequences that passed the present data filter to an output file.

#### SAMPLE COMMAND ####
#
my ($stacks_catalog_infile, $whitelist_infile, $output_dir);
GetOptions(
    'i=s'    => \$stacks_catalog_infile, # The Population Stacks write_single_snp fasta formatted input file.
    'w=s'    => \$whitelist_infile, # The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the filter reference genome fasta file.
    'o=s'    => \$output_dir, # The directory to contain all the output files.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
    defined $stacks_catalog_infile
    and defined $whitelist_infile
    and defined $output_dir
);

sub usage {
    
    die <<"USAGE";
    
    
Usage: $0 -i stacks_catalog_infile -w whitelist_infile -o output_dir
    
VERSION 1.22

DESCRIPTION - A Program that generates a reference genome fasta file filtered to only contain sequences that contain SNPs generated from the refgen_stacks_analysis_pipeline.pl and populations_stacks.pl scripts. Requires a reference genome that was formatted by the refgen_stacks_analysis_pipeline.pl script so that we can grab the SNPs by the reference sequence ID found in the populations Stacks write_single_snp fasta formatted input file. You can also used a whitelist of locus ids to get only the reference genome scaffolds that pass the present data filter in convert_fasta2nexus.pl. Also prints out the populations fasta sequences that passed the present data filter to an output file.
    
OPTIONS:
    
    -i stacks_catalog_infile - The Population Stacks write_single_snp fasta formatted input file.
    
    -w whitelist_infile - The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the filter reference genome fasta file.
    
    -o output_dir - The directory to contain all the output files.
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
    mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Parse the contents of the Stacks fasta file to obtain the reference genome sequence ids.
open(INFILE, "<$stacks_catalog_infile") or die "Couldn't open file $stacks_catalog_infile for reading, $!";
my %catalog_entries = ();
my ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
while(<INFILE>){
    chomp $_;
    if($_ !~ m/\#/){
        my @split_entries = split(/\t/, $_);
    
        #>CLocus_3_Sample_58_Locus_4_Allele_0 [F11; 0037, 275, -]
        my ($sql_id, $sample_id, $locus_id, $chromosome, $basepair, $strand, $sequence_type, $stack_component, $sequence_id,
            $sequence, $deleveraged_flag, $blacklisted_flag, $lumberjackstack_flag, $log_likelihood) = @split_entries;
        $catalog_entries{$locus_id}{$chromosome}{$basepair}{$strand} = join("\n", join(" ", join("_", ">CLocus", $locus_id, $sequence_type, "sequence"), join("; ", "[" . $chromosome, $basepair, $strand . "]")), $sequence);
    
    }
}
close(INFILE) or die "Couldn't close file $stacks_catalog_infile";

# Adding locus sequences in a hash array in order to concatenate the sequences.
my %filtered_refgen_ids = ();
my @locus_id_list = ();
if(-s $whitelist_infile){
    
    open(INFILE, "<$whitelist_infile") or die "Couldn't open file $whitelist_infile for reading, $!";
    while(<INFILE>){
        chomp $_;
        
        my $locus_id = $_;
        push(@locus_id_list, $locus_id);
        
    }
    close(INFILE) or die "Couldn't close file $whitelist_infile";
    
}else{
    die "$whitelist_infile does not exist";
}

# Print out the fasta sequences to the filtered populations stacks outfile.
my $filtered_pop_fasta_filename = fileparse($whitelist_infile, qr/\.txt|\.csv|\.tsv/);
my $filtered_pop_fasta_outfile = join('/', $output_dir, join("", $filtered_pop_fasta_filename, ".fasta"));
open(OUTFILE, ">$filtered_pop_fasta_outfile") or die "Couldn't open file $filtered_pop_fasta_outfile for writting, $!";
foreach my $locus_id (sort {$a <=> $b} @locus_id_list){
    foreach my $chromosome (sort {$a <=> $b} keys %{$catalog_entries{$locus_id}}){
        foreach my $basepair (sort {$a <=> $b} keys %{$catalog_entries{$locus_id}{$chromosome}}){
            foreach my $strand (sort {$a cmp $b} keys %{$catalog_entries{$locus_id}{$chromosome}{$basepair}}){
                print OUTFILE $catalog_entries{$locus_id}{$chromosome}{$basepair}{$strand} . "\n";
            }
        }
    }
}
close(OUTFILE) or die "Couldn't close file $filtered_pop_fasta_outfile";
