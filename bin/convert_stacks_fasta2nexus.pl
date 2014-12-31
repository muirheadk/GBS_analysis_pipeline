#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

# perl ~/Desktop/GBS_analysis_pipeline-master/bin/convert_stacks_fasta2nexus.pl -i ~/Desktop/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -f ~/Desktop/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/batch_1.fa -o ~/Desktop
my ($stacks_input_dir, $stacks_fasta_infile, $gbs_sequence_length, $output_dir);
GetOptions(
    'i=s'    => \$stacks_input_dir,
    'f=s'    => \$stacks_fasta_infile,
    'l=s'    => \$gbs_sequence_length, # The GBS fasta sequence length in base pairs (bps) common to all GBS fasta sequences. Default: 92
	'o=s'    => \$output_dir,
);

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

DESCRIPTION - 

OPTIONS:

-i stacks_input_dir -

-f stacks_fasta_infile -

-l gbs_sequence_length -

-o output_dir

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $nexus_output_dir = join('/', $output_dir, "NEXUS_OUTFILES");
unless(-d $nexus_output_dir){
    mkdir($nexus_output_dir, 0777) or die "Can't make directory: $!";
}

my ($stacks_alleles_files, $stacks_alleles_file_count) = find_files($stacks_input_dir, "alleles.tsv");

my %sample_names = ();
foreach my $file_name (sort keys %{$stacks_alleles_files}){
	if($file_name !~ m/batch_\d+\.catalog\.alleles\.tsv/){
		my $stacks_alleles_infile = $stacks_alleles_files->{$file_name};
		
		# Get the basename of the alleles filename without the .alleles.tsv extension.
		my $stacks_alleles_filename = fileparse($stacks_alleles_infile, qr/\.alleles\.tsv/);
		
		warn "Processing " . $stacks_alleles_filename . ".....\n";
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

my $seqio = Bio::SeqIO->new(-file => $stacks_fasta_infile, '-format' => 'Fasta');
my (%fasta_locus_nexus, %fasta_locus_seq_counts, %fasta_locus_refgen) = ();

while(my $seq_entry = $seqio->next_seq) {
    
    my $seq_id = $seq_entry->id;
    my $sequence = $seq_entry->seq;
    my $seq_desc = $seq_entry->desc;
    
    my $fasta_header = join(" ", $seq_id, $seq_desc);
    
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

foreach my $locus_id (sort {$a <=> $b} keys %fasta_locus_nexus){
#    warn $locus_id . "\n";
    if($fasta_locus_seq_counts{$locus_id} > 1){
        my @nexus_sequence_list = ();
        foreach my $sample_id (sort {$a <=> $b} keys %{$fasta_locus_nexus{$locus_id}}){
            
            foreach my $allele_id (sort {$a <=> $b} keys %{$fasta_locus_nexus{$locus_id}{$sample_id}}){
                push(@nexus_sequence_list, join("  ", join("_", $sample_names{$sample_id}, "Allele", $allele_id), $fasta_locus_nexus{$locus_id}{$sample_id}{$allele_id}));
            }
        }
        
        my ($refgen_seq_id, $snp_position, $refgen_seq_strand) = split(/\t/, $fasta_locus_refgen{$locus_id});
        
        my $refgen_strand = "";
        $refgen_strand = "Plus" if($refgen_seq_strand eq "+");
        $refgen_strand = "Minus" if($refgen_seq_strand eq "-");
        
        my $nexus_outfile = join('/', $nexus_output_dir, join("_", "locus$locus_id", "refSeq$refgen_seq_id", "pos$snp_position", "strand$refgen_strand") . ".nex");
        open(OUTFILE, ">$nexus_outfile") or die "Couldn't open file $nexus_outfile for writting, $!";
        print OUTFILE "#NEXUS" . "\n";
        print OUTFILE "Begin data;" . "\n";
        print OUTFILE "Dimensions ntax=$fasta_locus_seq_counts{$locus_id} nchar=$gbs_sequence_length" . "\n";
        print OUTFILE "Format datatype=dna symbols=\"ACTG\" missing=N gap=-;" . "\n";
        print OUTFILE "Matrix;" . "\n";
        print OUTFILE join("\n", @nexus_sequence_list) . "\n";
        print OUTFILE ";" . "\n";
        print OUTFILE "End;" . "\n";
        close(OUTFILE) or die "Couldn't close file $nexus_outfile";
    }
}

sub find_files{
    
	my $infile_dir = shift;
	die "Error lost input file directory" unless defined $infile_dir;
	
	my $suffix = shift;
	die "Error lost file extension suffix directory" unless defined $suffix;
	
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
}
