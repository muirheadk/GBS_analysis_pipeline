#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use Bio::Graphics;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;
use File::Basename;

my ($adaptor_blastn_infile, $output_dir);
GetOptions(
      'i=s'    => \$adaptor_blastn_infile,
      'o=s'    => \$output_dir,
);
# Needs either fastafile parser or gbs adaptor blast to add lengths so that this file becomes more generic. Needs to have adaptor sequence trimmed 2-AGCT-6 cut correctly. Needs to be incorporated twice in pipeline. Once after blast and once after the trimming to see which sequences ended up getting trimmed.
usage() unless (
      defined $adaptor_blastn_infile
      and defined $output_dir
);

sub usage {

die <<"USAGE";

Usage: $0 -i adaptor_blastn_infile -o output_dir

Description - 

OPTIONS:


      -i adaptor_blastn_infile -

      -o output_dir -

USAGE
}

my $gbs_sequence_length;
$gbs_sequence_length = 100 unless defined $gbs_sequence_length;

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %blast_entries = ();
open(INFILE, "<$adaptor_blastn_infile") or die "Couldn't open file $adaptor_blastn_infile for reading, $!";
my $i = 0;
while(<INFILE>){
      chomp $_;
      if($i ne 0){
#  	    warn $_ . "\n";
	    
	    my @split_row_entry = split(/\t/, $_);
	    my ($query_id,$target_id,$query_coverage,$percent_identity,$align_length,$num_mismatch,
		  $num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score, $glyph_colour) = @split_row_entry;
	    my $blast_entry = join("\t",$query_id,$target_id,$query_coverage,$percent_identity,$align_length,$num_mismatch,
		  $num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score,$glyph_colour);

	    push(@{$blast_entries{$target_id}}, $blast_entry);
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $adaptor_blastn_infile";

# length needs actual sequence length from blast.
my $panel = Bio::Graphics::Panel->new(
    -length    => $gbs_sequence_length,
    -width     => 1500,
    -pad_left  => 10,
    -pad_right => 800,
);

foreach my $target_name (keys %blast_entries){

      warn "Processing blast results for " . $target_name . "....\n";


	my ($fastq_sequence_length, $full_length);

	my $i = 0;
	foreach my $blast_entry (@{$blast_entries{$target_name}}){
    my @split_blast_entry = split(/\t/, $blast_entry);
    my ($query_id,$target_id,$query_coverage,$percent_identity,$align_length,$num_mismatch,$num_gaps,$query_start,
    $query_end,$target_start,$target_end,$e_value,$bit_score,$glyph_colour) = @split_blast_entry;

    if($query_start > $query_end){
        warn "Found hit in antisense direction for query_id, so skipping because we aren't interested in this alignment....\n";
        next;
    }
    if($target_start > $target_end){
        warn "Found hit in antisense direction for target_id, so skipping because we aren't interested in this alignment....\n";
        next;
    }
    if($num_mismatch > 0){
        warn "Found hit with number of mismatches greater than 0, so skipping because we aren't interested in this alignment....\n";
        next;
    }
    if($num_gaps > 0){
        warn "Found hit with number of gaps greater than 0, so skipping because we aren't interested in this alignment....\n";
        next;
    }
	if($i eq 0){
		   if($target_name =~ m/length=(\d+)/){
		$fastq_sequence_length = $1;
    }
      $full_length = Bio::SeqFeature::Generic->new(
	    -start        => 1,
	    -end          => $fastq_sequence_length,
	    -display_name => $target_name,
      );

      $panel->add_track(
	    $full_length,
	    -glyph   => 'arrow',
	    -tick    => 2,
	    -fgcolor => 'black',
	    -double  => 1,
	    -label   => 1,
      );
	}

    my ($query_name, $adaptor_sequence, $description);  
    if($query_id =~ m/([\w_]+)_([ACGT]+)/) {
	($query_name, $adaptor_sequence) = ($1, $2);
	my $adaptor_end = length($adaptor_sequence);
	my $adaptor_sub_sequence = get_subseq($adaptor_sequence, $query_start, $query_end);
	my $adaptor_start_sequence = get_subseq($adaptor_sequence, 1, ($query_start - 1));
	my $adaptor_end_sequence = get_subseq($adaptor_sequence, ($query_end + 1), $adaptor_end);


	my $adaptor_concatenated_sequence  = join("-", $adaptor_start_sequence, $query_start, $adaptor_sub_sequence, $query_end, $adaptor_end_sequence);
    	$description = join("; ", join("=", "sequence", $adaptor_concatenated_sequence), join("=", "query_coverage", $query_coverage), join("=", "percent_idenity", $percent_identity), join("=", "sequence_start", $target_start), join("=", "sequence_end", $target_end), join("=", "sequence_length", $align_length));
    }else{
	$query_name = $query_id;
	$description = join("; ", join("=", "sequence_start", $target_start), join("=", "sequence_end", $target_end), join("=", "sequence_length", $align_length));
    }

      my $track = $panel->add_track(
	    -glyph       => 'graded_segments',
	    -label       => 1,
	    -connector   => 'dashed',
	    -bgcolor     => $glyph_colour,
	    -fgcolor     => 'black',
	    -font2color  => 'red',
	    -sort_order  => 'high_score',
	    -description => sub {
		  my $feature = shift;
		  return unless $feature->has_tag('description');
		  my $description = $feature->each_tag_value('description');
		  "$description";
	    },
      );
    my $feature = Bio::SeqFeature::Generic->new(
      -start        => $target_start,
      -end          => $target_end,
      -score        => $bit_score,
      -display_name => $query_name,
      -tag          => {
        description => $description
      },
    );

    $track->add_feature($feature);
	
	$i++;

	}
}

my $adaptor_blastn_filename = fileparse($adaptor_blastn_infile);
warn "Generating alignment image for " . $adaptor_blastn_infile . "....\n";
my $outfile = join('/', $output_dir, join("_",$adaptor_blastn_filename, "alignment" . ".png"));
open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print OUTFILE $panel->png;
close(OUTFILE) or die "Couldn't close file $outfile";

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

