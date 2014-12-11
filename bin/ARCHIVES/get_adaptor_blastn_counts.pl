#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long;

# perl get_adaptor_blastn_counts.pl -i /home/cookeadmin/workspace/gbs_adaptor_blastn_counts_plate5.txt -o ~/workspace/GBS_data-08-10-2013/GBS_ADAPTOR_BLASTN_COUNTS

my ($infile, $output_dir);
GetOptions(
	'i=s'    => \$infile,
	'o=s'    => \$output_dir,
);

usage() unless (
	defined $infile
	and defined $output_dir
);

sub usage {
    
    die <<"USAGE";
    
    
Usage: $0 -i infile_dir -o output_dir
    
    Description -
    
OPTIONS:
    
	-i infile_dir -
	
    
	-o output_dir -
    
    
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %adaptor_blastn_hits = ();
open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
my $i = 0;
while(<INFILE>){
	chomp $_;
	if($i ne 0){
		my @split_adaptor_blastn_hit =  split(/\t/, $_);
		my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
		$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, $gylph_colour) = @split_adaptor_blastn_hit;
		
        push(@{$adaptor_blastn_hits{$target_name}}, join("\t", $query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
		$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, $gylph_colour));
     
        
		
	}
	$i++
	
}
close(INFILE) or die "Couldn't close file $infile";

my %adaptor_align_length_counts = ();
my %adaptor_target_end_counts = ();
foreach my $target_id (keys %adaptor_blastn_hits){
    my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
    $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, $gylph_colour);
    
    warn $target_id . "\n";    
    my $optimal_adaptor_blastn_hit = "";
    if(scalar(@{$adaptor_blastn_hits{$target_id}}) > 1){
        my $current_max_adaptor_length = 0;
        foreach my $adaptor_blastn_hit (@{$adaptor_blastn_hits{$target_id}}){
            my @split_adaptor_blastn_hit =  split(/\t/, $adaptor_blastn_hit);
            ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
            $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, $gylph_colour) = @split_adaptor_blastn_hit;
            if($align_length > $current_max_adaptor_length){
                $current_max_adaptor_length = $align_length;
                $optimal_adaptor_blastn_hit = $adaptor_blastn_hit;
            }
        }
    }elsif(@{$adaptor_blastn_hits{$target_id}} eq 1){
        $optimal_adaptor_blastn_hit = @{$adaptor_blastn_hits{$target_id}}[0];
        my @split_adaptor_blastn_hit =  split(/\t/, $optimal_adaptor_blastn_hit);
        ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
        $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score, $gylph_colour) = @split_adaptor_blastn_hit;
        
    }
    
    my $adaptor_concatenated_sequence;
    if($query_name =~ m/([\w_]+)_([ACGT]+)/){
        my ($query_id, $adaptor_sequence) = ($1, $2);
        my $adaptor_sequence_length = length($adaptor_sequence);
        my $adaptor_end = length($adaptor_sequence);
        my $adaptor_sub_sequence = get_subseq($adaptor_sequence, $query_start, $query_end);
        my $adaptor_start_sequence = get_subseq($adaptor_sequence, 1, ($query_start - 1));
        my $adaptor_end_sequence = get_subseq($adaptor_sequence, ($query_end + 1), $adaptor_end);
        if($align_length eq $adaptor_sequence_length){
            $adaptor_concatenated_sequence  = join("-", $query_start, $adaptor_sub_sequence, $query_end, $adaptor_end_sequence);
        }elsif($adaptor_start_sequence eq ""){
            $adaptor_concatenated_sequence  = join("-", $query_start, $adaptor_sub_sequence, $query_end);
        }elsif($adaptor_end_sequence eq ""){
            $adaptor_concatenated_sequence  = join("-", $adaptor_start_sequence, $query_start, $adaptor_sub_sequence, $query_end);
        }else{
            $adaptor_concatenated_sequence  = join("-", $adaptor_start_sequence, $query_start, $adaptor_sub_sequence, $query_end, $adaptor_end_sequence);
        }
    }
    
    $adaptor_align_length_counts{$align_length}{$adaptor_concatenated_sequence}++;
    
    if($align_length eq 64){
        $adaptor_target_end_counts{$target_end}++;
    }
}

my $adaptor_length_counts_outfile = join("/", $output_dir, "all_adaptor_length_counts.txt");
open(OUTFILE, ">$adaptor_length_counts_outfile") or die "Couldn't open file $adaptor_length_counts_outfile for writting, $!";
print OUTFILE join("\t", "adapter_sequence_id", "adapter_length", "adaptor_sequence_count") . "\n";
foreach my $adaptor_length (sort {$b <=> $a} keys %adaptor_align_length_counts){
	foreach my $adaptor_concatenated_sequence (keys %{$adaptor_align_length_counts{$adaptor_length}}){
		my $adaptor_sequence_count = $adaptor_align_length_counts{$adaptor_length}{$adaptor_concatenated_sequence};
		print OUTFILE join("\t", $adaptor_concatenated_sequence, $adaptor_length, $adaptor_sequence_count) . "\n";
		warn join("\t", $adaptor_concatenated_sequence, $adaptor_length, $adaptor_sequence_count) . "\n";
	}
	
}
close(OUTFILE) or die "Couldn't close file $adaptor_length_counts_outfile";

my $adaptor_target_end_counts_outfile = join("/", $output_dir, "all_target_end_counts_adaptor_length_64.txt");
open(OUTFILE, ">$adaptor_target_end_counts_outfile") or die "Couldn't open file $adaptor_target_end_counts_outfile for writting, $!";
print OUTFILE join("\t", "adapter_sequence_id", "adapter_length", "adaptor_sequence_count") . "\n";
foreach my $target_end (sort {$b <=> $a} keys %adaptor_target_end_counts){

	my $adaptor_target_end_count = $adaptor_target_end_counts{$target_end};
	print OUTFILE join("\t", $target_end, $adaptor_target_end_count) . "\n";
	
}
close(OUTFILE) or die "Couldn't close file $adaptor_target_end_counts_outfile";

# my $seq = get_subseq("AGCTTGCGTT", 3, 8);
# warn $seq . "\n";
sub get_subseq{

        my $sequence = shift;
        die "Error lost sequence" unless defined $sequence;

        my $seq_start = shift;
        die "Error lost start of sequence" unless defined $seq_start;

        my $seq_end = shift;
        die "Error lost end of sequence" unless defined $seq_end;

        $seq_start = $seq_start - 1;
        $seq_end = $seq_end;

        my $length = ($seq_end - $seq_start);

        my $trimmed_seq = substr($sequence, $seq_start, $length);

        return uc($trimmed_seq);
}
