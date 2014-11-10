#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

# perl reformat_refgen_to_bwatassel.pl -i ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/mountain_pine_beetle/DendPond_male_1.0/Primary_Assembly/unplaced_scaffolds/FASTA/unplaced.scaf.fa -p MPB_MALE_GBS -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_GBS_ANALYSIS/REFERENCE_GENOMEmy ($fasta_infile, $project_name, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'p=s'    => \$project_name,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $fasta_infile
      and defined $project_name
      and defined $output_dir
);

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -p project_name -o output_dir
    
Description - 
    
OPTIONS:

      -i fasta_infile - 
      
      -p project_name - 
      
      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $fasta_filename = fileparse($fasta_infile);

my $fasta_outfile = join('/', $output_dir, join("_", $project_name, $fasta_filename . ".fasta"));
open(OUTFILE1, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!";

my $desc_outfile = join('/', $output_dir, join("_", $project_name, $fasta_filename, "toc.txt"));
open(OUTFILE2, ">$desc_outfile") or die "Couldn't open file $desc_outfile for writting, $!"; 
print OUTFILE2 join("\t", "bwatassel_fasta_header", "original_fasta_header") . "\n";
warn "Converting $fasta_filename to bwatassel format....\n\n";
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
my $i = 1;
while(my $seq_entry = $seqio->next_seq) {

	my $seq_id = $seq_entry->id;
	my $sequence = $seq_entry->seq;
	my $seq_desc = $seq_entry->desc;
 
#  	warn join("\n", ">$seq_id", $sequence) . "\n";
	my $fasta_header = join(" ", $seq_id, $seq_desc);
	warn $fasta_header . "\n";
	print OUTFILE1 join("\n", ">$i", $sequence) . "\n";
	print OUTFILE2 join("\t", $i, $fasta_header) . "\n";
	
	$i++;
     
}
close(OUTFILE1) or die "Couldn't close file $fasta_outfile"; 
close(OUTFILE2) or die "Couldn't close file $desc_outfile";

# Clean out the sequence I/O object.
$seqio = ();