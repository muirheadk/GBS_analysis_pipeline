#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

#### PROGRAM NAME ####
# convert_stacks_fasta2nexus.pl - Program that generates a file in NEXUS file format from a Stacks fasta formatted file from population stacks.

#### DESCRIPTION ####
# This program generates a file in NEXUS file format from a Stacks fasta formatted file from population stacks.

#### SAMPLE COMMAND ####
# perl convert_stacks_fasta2nexus.pl -m ~/Desktop/sbwphylo_popmap.txt -i ~/SBWPHYLO_TRIMMED_OFFSET_5_POPULATIONS/POPULATIONS_FASTA_DIR/batch_1.fa -l 92 -p 0.75 -o ~/SBWPHYLO_TRIMMED_OFFSET_5_NEXUS
my ($stacks_fasta_infile, $stacks_popmap_infile, $whitelist_infile,$gbs_sequence_length, $percent_present_data, $output_dir);
GetOptions(
'i=s'    => \$stacks_fasta_infile, # The Population Stacks write_single_snp fasta formatted input file.
'm=s'    => \$stacks_popmap_infile, # The Stacks population map formatted input file.
'w=s'    => \$whitelist_infile, # The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the nexus file.
'l=s'    => \$gbs_sequence_length, # The GBS fasta sequence length in base pairs (bps) common to all GBS fasta sequences. Default: 92
'p=s'    => \$percent_present_data, # The percent of present data at a locus. percent_present_data = (present_data/total_data) Default: 0.75
'o=s'    => \$output_dir, # The directory to contain all the output files.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
defined $stacks_fasta_infile
and defined $stacks_popmap_infile
and defined $output_dir
);

# The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the nexus file.
$whitelist_infile = "" unless defined $whitelist_infile;

# The GBS fasta sequence length in base pairs (bps) common to all GBS fasta sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

# The percent of present data at a locus. percent_present_data = (present_data/(present_data + missing_data)) Default: 0.75
$percent_present_data = 0.75 unless defined $percent_present_data;

sub usage {

die <<"USAGE";


Usage: $0 -i stacks_input_dir -f stacks_fasta_infile -w whitelist_infile -l gbs_sequence_length -p percent_present_data -o output_dir

DESCRIPTION - This program generates a file in NEXUS file format from a Stacks fasta formatted file from population stacks.

OPTIONS:

-i stacks_fasta_infile - The Population Stacks write_single_snp fasta formatted input file.

-m stacks_popmap_infile - The Stacks population map formatted input file.

-w whitelist_infile - The whitelist input file containing the locus ids to keep after passing the present data filter. If not specified all locus ids that passed the present data filter will be used to generate the nexus file.
    
-l gbs_sequence_length - The GBS fasta sequence length in base pairs (bps) common to all GBS fasta sequences. Default: 92

-p percent_present_data - The percent of present data at a locus. percent_present_data = (present_data/total_data) Default: 0.75
    
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

# If the nexus locus catolog is already generated then skip, otherwise generate the nexus locus catalog.
my $locus_catalog_outfile = join('/', $output_dir, "nexus_locus_catalog.txt");
unless(-s $locus_catalog_outfile){
    
    # Parse the contents of the Stacks fasta file to generate the nexus files.
    open(INFILE, "<$stacks_fasta_infile") or die "Couldn't open file $stacks_fasta_infile for reading, $!";
    my (%locus_metadata, %fasta_locus_seq_counts, %fasta_locus_refgen) = ();
    my ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
    while(<INFILE>){
        chomp $_;
        
        if($_ =~ /^>/){
            my $fasta_header = $_;
            # Parse stacks fasta output file for the locus, allele number, individual name id, reference genome sequence id, SNP position, and reference genome sequence strand orientation.
            if($fasta_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+); (\d+), (\d+), (\+|\-)\]/){
                warn $fasta_header . "\n";
                ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = ($1, $2, $3, $4, $5, $6);

                $fasta_locus_refgen{$locus_id} = join("\t", $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand);
                $fasta_locus_seq_counts{$locus_id}++;
                #die join("\t", $locus_id, $individual_id, $allele_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) . "\n";
            }elsif($fasta_header =~ m/^>CLocus_(\d+)_Sample_\d+_Locus_\d+_Allele_(\d+) \[(.+)\]/){
                warn $fasta_header . "\n";
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
            ($locus_id, $allele_id, $individual_id, $refgen_seq_id, $refgen_seq_start, $refgen_seq_strand) = "";
        }
    }
    close(INFILE) or die "Couldn't close file $stacks_fasta_infile";

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
                            warn "complex SNP $unique_concat_snps at locus $locus_id and individual $individual_id" if($unique_concat_snps eq "C,G,T");
                            warn "complex SNP $unique_concat_snps at locus $locus_id and individual $individual_id" if($unique_concat_snps eq "A,G,T");
                            warn "complex SNP $unique_concat_snps at locus $locus_id and individual $individual_id" if($unique_concat_snps eq "A,C,T");
                            warn "complex SNP $unique_concat_snps at locus $locus_id and individual $individual_id" if($unique_concat_snps eq "A,C,G");
                            
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
    open(OUTFILE, ">$locus_catalog_outfile") or die "Couldn't open file $locus_catalog_outfile for writting, $!";
    print OUTFILE join("\t", "locus_id", "refgen_seq_id", "refgen_seq_start", "strand", "individual_id", "snp_position", "sequence") . "\n";
    foreach my $locus_id (sort {$a <=> $b} keys %locus_metadata){
        my @locus_catalog_sequence_list = ();
        foreach my $individual_id (sort {$a cmp $b} keys %sample_names){
            warn $individual_id . "\n";
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


open(INFILE, "<$locus_catalog_outfile") or die "Couldn't open file $locus_catalog_outfile for reading, $!";
my $i = 0;
my (%locus_catalog_data, %locus_catalog_counts) = ();
while(<INFILE>){
    chomp $_;
    
    if(($i ne 0) and ($_ !~ /^$/)){
        
        my @split_locus_catalog_entry = split(/\t/, $_);
        my ($locus_id, $refgen_seq_id, $refgen_seq_start, $strand, $individual_id, $snp_position, $sequence) = @split_locus_catalog_entry;
        if($sequence =~ m/^[ACGTNRYSWKM]+$/){
            
            $locus_catalog_data{$locus_id}{$individual_id} = $sequence;
            warn join("\t", $locus_id, $individual_id, "$sequence is DNA characters");
            $locus_catalog_counts{$locus_id}{"PRESENT_DATA"}++;
        }elsif($sequence =~ m/^\?{$gbs_sequence_length}$/){
            
            $locus_catalog_data{$locus_id}{$individual_id} = $sequence;
            warn join("\t", $locus_id, $individual_id, "$sequence is ?s");
            $locus_catalog_counts{$locus_id}{"MISSING_DATA"}++;
        }
        
    }
    
    $i++;
}
close(INFILE) or die "Couldn't close file $locus_catalog_outfile";

# Calculating the percentage of missing data to present data. Using the amount of present data we calculate percent present data as $percent_present_data = ($present_data/$total_data_count).
my %percent_missing_data = ();
my %nexus_loci_list = ();
my $present_loci_count = 0;
my $locus_list_outfile = join('/', $output_dir, "filtered_locus_list.txt");
open(OUTFILE, ">$locus_list_outfile") or die "Couldn't open file $locus_list_outfile for writting, $!";
foreach my $locus_id (sort {$a <=> $b} keys %locus_catalog_counts){
    
    my $present_data_count = 0;
    if(defined($locus_catalog_counts{$locus_id}{"PRESENT_DATA"})){
        $present_data_count = $locus_catalog_counts{$locus_id}{"PRESENT_DATA"};
    }
    
    my $missing_data_count = 0;
    if(defined($locus_catalog_counts{$locus_id}{"MISSING_DATA"})){
        $missing_data_count = $locus_catalog_counts{$locus_id}{"MISSING_DATA"};
    }
    
    my $total_data_count = ($present_data_count + $missing_data_count);
    die "locus id: $locus_id total data: $total_data_count ne sample counts: $sample_counts" if($total_data_count ne $sample_counts);
    my $perc_present_data = ($present_data_count/$total_data_count);
    $percent_missing_data{$locus_id} = $present_data_count;
    warn "$locus_id ($present_data_count/$total_data_count) * 100" . " " . ($present_data_count/$total_data_count) . " " . $percent_missing_data{$locus_id};

    # Grab the locus ids that pass the present data filter test.
    if($perc_present_data >= $percent_present_data){
        print OUTFILE $locus_id . "\n";
        $nexus_loci_list{$locus_id} = $locus_id;
        $present_loci_count++;
    }
}
close(OUTFILE) or die "Couldn't close file $locus_list_outfile";

# Adding locus sequences in a hash array in order to concatenate the sequences.
my %nexus_seqs = ();
if($whitelist_infile eq ""){
    foreach my $locus_id (sort {$a <=> $b} keys %nexus_loci_list){
        foreach my $individual_id (sort {$a cmp $b} keys %sample_names){
            push(@{$nexus_seqs{$individual_id}}, $locus_catalog_data{$locus_id}{$individual_id});
        }
    }
}elsif(-s $whitelist_infile){
    
    open(INFILE, "<$whitelist_infile") or die "Couldn't open file $whitelist_infile for reading, $!";
    my @nexus_loci_list = ();
    $present_loci_count = 0;
    while(<INFILE>){
        chomp $_;
        
        my $locus_id = $_;
        if(defined($nexus_loci_list{$locus_id})){
            push(@nexus_loci_list, $locus_id);
            $present_loci_count++;
        }
    }
    close(INFILE) or die "Couldn't close file $whitelist_infile";
    
    
    foreach my $locus_id (sort {$a <=> $b} @nexus_loci_list){
        
        foreach my $individual_id (sort {$a cmp $b} keys %sample_names){
            push(@{$nexus_seqs{$individual_id}}, $locus_catalog_data{$locus_id}{$individual_id});
        }
    }
}else{
    die "$whitelist_infile does not exist";
}

# Print out the locus metadata contents in nexus file format.
my $nexus_nchar_length = ($present_loci_count * $gbs_sequence_length);

my $nexus_outfile = join('/', $output_dir, "stacks_nexus_file.nex");
open(OUTFILE, ">$nexus_outfile") or die "Couldn't open file $nexus_outfile for writting, $!";
print OUTFILE "#NEXUS" . "\n";
print OUTFILE "Begin data;" . "\n";
print OUTFILE "Dimensions ntax=$sample_counts nchar=$nexus_nchar_length;" . "\n";
print OUTFILE "Format datatype=dna symbols=\"ACGTNRYSWKM\" missing=? gap=-;" . "\n";
print OUTFILE "Matrix" . "\n";
foreach my $individual_id (sort {$a cmp $b} keys %sample_names){
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
