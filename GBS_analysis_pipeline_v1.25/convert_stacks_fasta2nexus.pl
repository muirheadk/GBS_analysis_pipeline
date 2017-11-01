#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

# VERSION 1.25

#### PROGRAM NAME ####
# convert_stacks_fasta2nexus.pl - Program that generates a file in NEXUS file format from a Stacks fasta formatted file from population stacks.

#### DESCRIPTION ####
# This program generates a file in NEXUS file format from a Stacks fasta formatted file from population stacks.

#### SAMPLE COMMAND ####
# perl convert_stacks_fasta2nexus.pl -m ~/Desktop/sbwphylo_popmap.txt -n sbwphylo -i ~/SBWPHYLO_TRIMMED_OFFSET_5_POPULATIONS/POPULATIONS_FASTA_DIR/batch_1.fa -p 0.75 -o ~/SBWPHYLO_TRIMMED_OFFSET_5_NEXUS
my ($stacks_fasta_infile, $stacks_popmap_infile, $project_name, $whitelist_infile, $percent_present_locus_data, $percent_present_indiv_data, $output_dir);
GetOptions(
	'i=s'    => \$stacks_fasta_infile, # The Population Stacks write_single_snp fasta formatted input file.
	'm=s'    => \$stacks_popmap_infile, # The Stacks population map formatted input file.
	'n=s'    => \$project_name,  # The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output file names.
	'w=s'    => \$whitelist_infile, # The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the nexus file.
	'p=s'    => \$percent_present_locus_data, # The percent of present data at a locus across all individuals. percent_present_locus_data = (present_data_locus/(present_data_locus + missing_data_locus)) Default: 0.75
	'd=s'    => \$percent_present_indiv_data, # The percent of present data for an individual across all characters. percent_present_indiv_data = (present_data_indiv/(present_data_indiv + missing_data_indiv)) Default: 0.75
	'o=s'    => \$output_dir, # The directory to contain all the output files.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
	defined $stacks_fasta_infile
	and defined $stacks_popmap_infile
	and defined $project_name
	and defined $output_dir
);

# The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the nexus file.
$whitelist_infile = "" unless defined $whitelist_infile;

# The percent of present data at a locus across all individuals. percent_present_locus_data = (present_data_locus/(present_data_locus + missing_data_locus)) Default: 0.75
$percent_present_locus_data = 0.75 unless defined $percent_present_locus_data;

# The percent of present data for an individual across all characters. percent_present_indiv_data = (present_data_indiv/(present_data_indiv + missing_data_indiv)) Default: 0.75
$percent_present_indiv_data = 0.75 unless defined $percent_present_indiv_data;

sub usage {

die <<"USAGE";
Usage: $0 -i stacks_input_dir -f stacks_fasta_infile -n project_name -w whitelist_infile -p percent_present_locus_data -d percent_present_indiv_data -o output_dir

VERSION 1.25

DESCRIPTION - This program generates a file in NEXUS file format from a Stacks fasta formatted file from population stacks.

OPTIONS:

 -i stacks_fasta_infile - The Population Stacks write_single_snp fasta formatted input file.
 
-m stacks_popmap_infile - The Stacks population map formatted input file.
 
-n project_name - The name of the Genotyping by Sequencing (GBS) project, which is used to generate the output file names.

-w whitelist_infile - The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the nexus file.    

-p percent_present_locus_data - The percent of present data at a locus across all individuals. percent_present_locus_data = (present_data_locus/(present_data_locus + missing_data_locus)) Default: 0.75
    
-d percent_present_indiv_data - The percent of present data for an individual across all characters. percent_present_indiv_data = (present_data_indiv/(present_data_indiv + missing_data_indiv)) Default: 0.75
 
-o output_dir - The directory to contain all the output files.

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Parse the sample names from the population map to have a list of individuals so that we can figure out missing data.
my %sample_names = ();
my $sample_counts = 0;
open(INFILE, "<$stacks_popmap_infile") or die "Couldn't open file $stacks_popmap_infile for reading, $!";
while(<INFILE>){
    chomp $_;
    if($_ !~ m/^$/){
        my @split_popmap_entries = split(/\t/, $_);
    
        my ($sample_name, $population_name) = @split_popmap_entries;
        $sample_names{$sample_name} = $sample_name;
        $sample_counts++;
    }
}
close(INFILE) or die "Couldn't close file $stacks_popmap_infile";


# Parse the contents of the Stacks fasta file to generate the nexus files.
open(INFILE, "<$stacks_fasta_infile") or die "Couldn't open file $stacks_fasta_infile for reading, $!";
my (%locus_metadata, %locus_seq_lengths, %fasta_locus_seq_counts, %fasta_locus_refgen) = ();
my ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
while(<INFILE>){
    chomp $_;
    
    if($_ =~ /^>/){
        my $fasta_header = $_;
        # Parse stacks fasta output file for the locus, allele number, individual name id, reference genome sequence id, SNP position, and reference genome sequence strand orientation.
        if($fasta_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+); (\d+), (\d+), (\+|\-)\]/){
            #warn $fasta_header . "\n";
            ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = ($1, $2, $3, $4, $5, $6);

            $fasta_locus_refgen{$locus_id} = join("\t", $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand);
            $fasta_locus_seq_counts{$locus_id}++;
            #die join("\t", $locus_id, $individual_id, $allele_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) . "\n";
        }elsif($fasta_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+)\]/){
            #warn $fasta_header . "\n";
            ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = ($1, $2, $3, "N/A", "N/A", "N/A");
            
            $fasta_locus_refgen{$locus_id} = join("\t", $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand);
            $fasta_locus_seq_counts{$locus_id}++;
            #die join("\t", $locus_id, $individual_id, $allele_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) . "\n";
            
        }else{
            die "Error: $fasta_header is not in the correct format!";
        }
    }elsif($_ =~ m/[ACGTN]+/){
        my $sequence = $_;
        $locus_metadata{$locus_id}{$individual_id}{$allele_id} = $sequence;
        push(@{$locus_seq_lengths{$locus_id}}, length($sequence));
        ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
    }
}
close(INFILE) or die "Couldn't close file $stacks_fasta_infile";

# Iterate through each individual at a given locus id to obtain the sequence length consensus for that locus.
my %locus_id_seq_length = ();
foreach my $locus_id (sort {$a <=> $b} keys %locus_seq_lengths){
    my %seen = ();
    my @sequence_length = grep { !$seen{$_}++ } @{$locus_seq_lengths{$locus_id}};
    
    if(scalar(@sequence_length) eq 1){
        $locus_id_seq_length{$locus_id} = $sequence_length[0];
        #warn join("\t", $locus_id, $locus_id_seq_length{$locus_id}) . "\n";
    }else{
        my $sequence_lengths = join(", ", @sequence_length);
        die "Error: locus id: $locus_id: sequence length not unique!\n$sequence_lengths\n";
    }
} 
%locus_seq_lengths = ();

# If the nexus locus catolog is already generated then skip, otherwise generate the nexus locus catalog.
my $locus_catalog_outfile = join('/', $output_dir, join("_", $project_name, "nexus_locus_catalog.txt"));
unless(-s $locus_catalog_outfile){

    # Iterate through each locus id to obtain the metadata for that locus and print out the contents in nexus file format.
    my (%processed_loci, %snp_loci_positions) = ();
    foreach my $locus_id (sort {$a <=> $b} keys %locus_metadata){
        # warn $locus_id . "\n";
            my @nexus_sequence_list = ();
            
            # Obtain the consensus locus sequence length.
            my $gbs_sequence_length = $locus_id_seq_length{$locus_id};
            
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
                        
                        #warn join("\t", $locus_id, $individual_id, $i, $allele_sequences[$i]) . "\n";
                        
                        my @allele_i = split('', $allele_sequences[$i]);
                        die "Error: Allele_0 is not the same length as Allele_$i" if(length($allele_sequences[0]) ne length($allele_sequences[$i]));
                        die "Error: Allele_0 is not the same length as gbs_sequence_length=$gbs_sequence_length" if(length($allele_sequences[0]) ne $gbs_sequence_length);
                        
                        my %nucleotide_counts = ();
                        for(my $nuc_position = 0; $nuc_position < length($allele_sequences[$i]); $nuc_position++){
                            if($allele_0[$nuc_position] ne $allele_i[$nuc_position]){
                                
                                #warn "$allele_0[$nuc_position] ne $allele_i[$nuc_position] $nuc_position" . "\n";
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
                                #warn join("\t", $nuc_position, $concat_snps) . "\n";
                                my @split_concat_snps = split(/,/, $concat_snps);
                                push(@snp_nuc_chars, @split_concat_snps);
                            }
                            my @unique_snp_nuc_chars = do { my %seen; grep { !$seen{$_}++ } sort {$a cmp $b} @snp_nuc_chars };
                            my $unique_concat_snps = join(",", @unique_snp_nuc_chars);
                            #warn "complex SNP $unique_concat_snps at locus $locus_id and individual $individual_id" if($unique_concat_snps eq "C,G,T");
                            #warn "complex SNP $unique_concat_snps at locus $locus_id and individual $individual_id" if($unique_concat_snps eq "A,G,T");
                            #warn "complex SNP $unique_concat_snps at locus $locus_id and individual $individual_id" if($unique_concat_snps eq "A,C,T");
                            #warn "complex SNP $unique_concat_snps at locus $locus_id and individual $individual_id" if($unique_concat_snps eq "A,C,G");
                            
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

    # Print out the locus metadata contents to the nexus locus catalog.
    my $locus_catalog_outfile = join('/', $output_dir, join("_", $project_name, "nexus_locus_catalog.txt"));
    open(OUTFILE, ">$locus_catalog_outfile") or die "Couldn't open file $locus_catalog_outfile for writting, $!";
    print OUTFILE join("\t", "locus_id", "refgen_seq_id", "refgen_seq_start", "strand", "individual_id", "snp_position", "sequence") . "\n";
    foreach my $locus_id (sort {$a <=> $b} keys %locus_metadata){
        my @locus_catalog_sequence_list = ();
        my $gbs_sequence_length = $locus_id_seq_length{$locus_id};
        foreach my $individual_id (sort {$a cmp $b} keys %sample_names){
            #warn $individual_id . "\n";
            if(!defined($processed_loci{$locus_id}{$individual_id})){
                
                $processed_loci{$locus_id}{$individual_id} = '?' x $gbs_sequence_length;
                $snp_loci_positions{$locus_id}{$individual_id} = "N/A";
            }
            
            # Grab the reference genome sequence id and strand orientation and the SNP position in the sequence.
            push(@locus_catalog_sequence_list, join("\t", $locus_id, $fasta_locus_refgen{$locus_id}, $sample_names{$individual_id}, $snp_loci_positions{$locus_id}{$individual_id}, $processed_loci{$locus_id}{$individual_id}));
        }
        print OUTFILE join("\n", @locus_catalog_sequence_list) . "\n";
    }
    close(OUTFILE) or die "Couldn't close file $locus_catalog_outfile";
}

# Count the proportion for present and missing data at a given locus based on whether or not a sequence is present for that individual.
open(INFILE, "<$locus_catalog_outfile") or die "Couldn't open file $locus_catalog_outfile for reading, $!";
my $i = 0;
my (%locus_catalog_data, %locus_catalog_counts) = ();
while(<INFILE>){
    chomp $_;
    
    if(($i ne 0) and ($_ !~ /^$/)){
        
        my @split_locus_catalog_entry = split(/\t/, $_);
        my ($locus_id, $refgen_seq_id, $refgen_seq_start, $strand, $individual_id, $snp_position, $sequence) = @split_locus_catalog_entry;
		
        my $gbs_sequence_length = $locus_id_seq_length{$locus_id};
        if($sequence =~ m/^[ACGTNRYSWKM]+$/){
            
            $locus_catalog_data{$locus_id}{$individual_id} = $sequence;
            #warn join("\t", $locus_id, $individual_id, "$sequence is DNA characters");
            
            $locus_catalog_counts{$locus_id}{"PRESENT_DATA"}++;
            
        }elsif($sequence =~ m/^\?{$gbs_sequence_length}$/){
            
            $locus_catalog_data{$locus_id}{$individual_id} = $sequence;
            #warn join("\t", $locus_id, $individual_id, "$sequence is ?s");
            
            $locus_catalog_counts{$locus_id}{"MISSING_DATA"}++;
        }
        
    }
    
    $i++;
}
close(INFILE) or die "Couldn't close file $locus_catalog_outfile";

# Calculating the percentage of missing data to present data. Using the amount of present data we calculate percent present data as $perc_present_locus_data = ($present_data_locus_count/$total_data_locus_count).
my %nexus_loci_list = ();
my $present_locus_list_outfile = join('/', $output_dir, join("_", $project_name, "filtered_locus_list.txt"));
open(OUTFILE1, ">$present_locus_list_outfile") or die "Couldn't open file $present_locus_list_outfile for writting, $!";

my $missing_locus_list_outfile = join('/', $output_dir, join("_", $project_name, "removed_locus_list.txt"));
open(OUTFILE2, ">$missing_locus_list_outfile") or die "Couldn't open file $missing_locus_list_outfile for writting, $!";

my $present_locus_data_outfile = join('/', $output_dir, join("_", $project_name, "present_locus_data.txt"));
open(OUTFILE3, ">$present_locus_data_outfile") or die "Couldn't open file $present_locus_data_outfile for writting, $!";
print OUTFILE3 join("\t", "locus_id", "percent_present_data", "percent_missing_data", "present_data_count", "missing_data_count", "total_data_count") . "\n";
foreach my $locus_id (sort {$a <=> $b} keys %locus_catalog_counts){
    
    my $present_data_locus_count = 0;
    if(defined($locus_catalog_counts{$locus_id}{"PRESENT_DATA"})){
        $present_data_locus_count = $locus_catalog_counts{$locus_id}{"PRESENT_DATA"};
    }
    
    my $missing_data_locus_count = 0;
    if(defined($locus_catalog_counts{$locus_id}{"MISSING_DATA"})){
        $missing_data_locus_count = $locus_catalog_counts{$locus_id}{"MISSING_DATA"};
    }
    
    my $total_data_locus_count = ($present_data_locus_count + $missing_data_locus_count);
    die "locus id: $locus_id total data: $total_data_locus_count ne sample counts: $sample_counts" if($total_data_locus_count ne $sample_counts);
    my $perc_present_locus_data = ($present_data_locus_count/$total_data_locus_count);
    my $perc_missing_locus_data = ($missing_data_locus_count/$total_data_locus_count);
    
    #warn "$locus_id ($present_data_locus_count/$total_data_locus_count) * 100" . " " . ($present_data_locus_count/$total_data_locus_count) . " " . $present_data_locus_count . " missing data = $missing_data_locus_count";

    # Grab the locus ids that pass the present data filter test.
    if(($perc_present_locus_data >= $percent_present_locus_data) and ($perc_present_locus_data <= 1)){
        print OUTFILE1 $locus_id . "\n";
        $nexus_loci_list{$locus_id} = $locus_id;
    }else{
        print OUTFILE2 $locus_id . "\n";
    }
    
    print OUTFILE3 join("\t", $locus_id, ($perc_present_locus_data * 100), ($perc_missing_locus_data * 100), $present_data_locus_count, $missing_data_locus_count, $total_data_locus_count) . "\n";
}
close(OUTFILE1) or die "Couldn't close file $present_locus_list_outfile";
close(OUTFILE2) or die "Couldn't close file $missing_locus_list_outfile";
close(OUTFILE3) or die "Couldn't close file $present_locus_data_outfile";

# Adding locus sequences in a hash array in order to concatenate the sequences.
my %nexus_seqs = ();
# Get the number of characters based on the unique length of each locus sequence.
my $nexus_nchar_length = 0;
if($whitelist_infile eq ""){
    foreach my $locus_id (sort {$a <=> $b} keys %nexus_loci_list){
        foreach my $individual_id (sort {$a cmp $b} keys %sample_names){
            push(@{$nexus_seqs{$individual_id}}, $locus_catalog_data{$locus_id}{$individual_id});
        }
        $nexus_nchar_length += $locus_id_seq_length{$locus_id};
    }
}elsif(-s $whitelist_infile){
    
    open(INFILE, "<$whitelist_infile") or die "Couldn't open file $whitelist_infile for reading, $!";
    my @nexus_loci_list = (); 
    while(<INFILE>){
        chomp $_;
        
        my $locus_id = $_;
        if(defined($nexus_loci_list{$locus_id})){
            push(@nexus_loci_list, $locus_id);
        }
    }
    close(INFILE) or die "Couldn't close file $whitelist_infile";
    
    
    foreach my $locus_id (sort {$a <=> $b} @nexus_loci_list){
        
        foreach my $individual_id (sort {$a cmp $b} keys %sample_names){
            push(@{$nexus_seqs{$individual_id}}, $locus_catalog_data{$locus_id}{$individual_id});
        }
        $nexus_nchar_length += $locus_id_seq_length{$locus_id};
    }
}else{
    die "$whitelist_infile does not exist";
}

# Get the present, missing, and total number of characaters in the data set for each individual.
my %locus_counts_indivs = ();
foreach my $individual_id (sort {$a cmp $b} keys %sample_names){

    if(defined(@{$nexus_seqs{$individual_id}})){
        my $nexus_concat_seq = join("", @{$nexus_seqs{$individual_id}});
        my $concat_seq_length = length($nexus_concat_seq);
        my @sequence_characters = split('', $nexus_concat_seq);
        foreach my $character (@sequence_characters) {
                if($character =~ m/[ACGTNRYSWKM]/){

                        #warn join("\t", $locus_id, $individual_id, "$sequence is DNA characters");

                        $locus_counts_indivs{$individual_id}{"PRESENT_DATA"}++;

                }elsif($character eq '?'){

                        #warn join("\t", $locus_id, $individual_id, "$sequence is ?s");

                        $locus_counts_indivs{$individual_id}{"MISSING_DATA"}++;
                }
        }
    }else{
        die "Cannot process number of characters at an individual due to an insufficient number of locus sequences present. There were no locus ids that passed the percent_present_locus_data filter constraint. Please see the \$present_locus_list_outfile and rerun the script using at least the minimum percent_present_locus_data value in the file.\n$present_locus_list_outfile";
    }
}

# Calculating the percentage of missing data to present data. Using the amount of present data we calculate percent present data as $perc_present_indiv_data = ($present_data_indiv_count/$total_data_indiv_count).
my %nexus_indiv_list = ();
my $present_indiv_list_outfile = join('/', $output_dir, join("_", $project_name, "filtered_indivs_list.txt"));
open(OUTFILE1, ">$present_indiv_list_outfile") or die "Couldn't open file $present_indiv_list_outfile for writting, $!";

my $missing_indiv_list_outfile = join('/', $output_dir, join("_", $project_name, "removed_indivs_list.txt"));
open(OUTFILE2, ">$missing_indiv_list_outfile") or die "Couldn't open file $missing_indiv_list_outfile for writting, $!";

my $present_indiv_data_outfile = join('/', $output_dir, join("_", $project_name, "present_indivs_data.txt"));
open(OUTFILE3, ">$present_indiv_data_outfile") or die "Couldn't open file $present_indiv_data_outfile for writting, $!";
print OUTFILE3 join("\t", "individual_id", "percent_present_data", "percent_missing_data", "present_data_count", "missing_data_count", "total_data_count") . "\n";
my $num_indivs = 0;
foreach my $individual_id (sort {$a cmp $b} keys %locus_counts_indivs){

    my $present_data_indiv_count = 0;
    if(defined($locus_counts_indivs{$individual_id}{"PRESENT_DATA"})){
        $present_data_indiv_count = $locus_counts_indivs{$individual_id}{"PRESENT_DATA"};
    }

    my $missing_data_indiv_count = 0;
    if(defined($locus_counts_indivs{$individual_id}{"MISSING_DATA"})){
        $missing_data_indiv_count = $locus_counts_indivs{$individual_id}{"MISSING_DATA"};
    }

    my $total_data_indiv_count = ($present_data_indiv_count + $missing_data_indiv_count);
    die "individual id: $individual_id total data: $total_data_indiv_count ne number of characters: $nexus_nchar_length" if($total_data_indiv_count ne $nexus_nchar_length);
    my $perc_present_indiv_data = ($present_data_indiv_count/$total_data_indiv_count);
    my $perc_missing_indiv_data = ($missing_data_indiv_count/$total_data_indiv_count);
    #warn "$individual_id ($present_data_indiv_count/$total_data_indiv_count) * 100" . " " . ($present_data_indiv_count/$total_data_indiv_count) . " " . $present_data_indiv_count . " missing data = $missing_data_indiv_count";

    # Grab the individual ids that pass the present data filter test.
    if(($perc_present_indiv_data >= $percent_present_indiv_data) and ($perc_present_indiv_data <= 1)){
        #warn $individual_id . "\n";
        print OUTFILE1 $individual_id . "\n";
        $nexus_indiv_list{$individual_id} = $individual_id;
        $num_indivs++;
    }else{
        print OUTFILE2 $individual_id . "\n";
    }

    print OUTFILE3 join("\t", $individual_id, ($perc_present_indiv_data * 100), ($perc_missing_indiv_data * 100), $present_data_indiv_count, $missing_data_indiv_count, $total_data_indiv_count) . "\n";
}
close(OUTFILE1) or die "Couldn't close file $present_indiv_list_outfile";
close(OUTFILE2) or die "Couldn't close file $missing_indiv_list_outfile";
close(OUTFILE3) or die "Couldn't close file $present_indiv_data_outfile";

# Print out the locus metadata contents in nexus file format.
my $nexus_outfile = join('/', $output_dir, join("_", $project_name, "stacks_nexus_file.nex"));
open(OUTFILE, ">$nexus_outfile") or die "Couldn't open file $nexus_outfile for writting, $!";
print OUTFILE "#NEXUS" . "\n";
print OUTFILE "Begin data;" . "\n";
print OUTFILE "Dimensions ntax=$num_indivs nchar=$nexus_nchar_length;" . "\n";
print OUTFILE "Format datatype=dna symbols=\"ACGTNRYSWKM\" missing=? gap=-;" . "\n";
print OUTFILE "Matrix" . "\n";
foreach my $individual_id (sort {$a cmp $b} keys %nexus_indiv_list){
    my $nexus_concat_seq = join("", @{$nexus_seqs{$individual_id}});
    my $concat_seq_length = length($nexus_concat_seq);
    die "nexus nchar length: $nexus_nchar_length ne concatenated sequence length: $concat_seq_length" if($nexus_nchar_length ne $concat_seq_length);
    print OUTFILE join("  ", $individual_id, $nexus_concat_seq) . "\n";
}
print OUTFILE ";" . "\n";
print OUTFILE "End;" . "\n";
close(OUTFILE) or die "Couldn't close file $nexus_outfile";





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
