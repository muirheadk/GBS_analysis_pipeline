#!/usr/bin/perl
use warnings;
use strict;

use Bio::SeqIO;
use List::Util qw(min max);
use File::Basename;

my $query_fasta_output_dir = "/home/cookeadmin/workspace/GBS_data-08-10-2013/Papilio_GBS_Data/Papilio_UNEAK_GBS-2014-12-04/UNEAK_FASTA_QUERY_FILES";
my ($fasta_query_files, $fasta_query_file_counter) = find_files($query_fasta_output_dir, "fasta");
my %fasta_query_regex = ();
my %fasta_query_range_lengths = ();
foreach my $filename (sort keys %{$fasta_query_files}){
	warn $filename . "\n";
	my $fasta_query_infile = $fasta_query_files->{$filename};
	
	my $fasta_filename = fileparse($fasta_query_infile, qr/\.fasta/);
	my @split_query_file = split(/_/, $fasta_filename);
	my $individual_name = $split_query_file[0];
	
	my @fasta_query_lengths = ();
	my $seqio = Bio::SeqIO->new(-file => $fasta_query_infile, '-format' => 'Fasta');
	while(my $seq_entry = $seqio->next_seq) {

		my $fasta_query_seq_id = $seq_entry->id;
		my $fasta_query_sequence = $seq_entry->seq;
		# FS034_TP10005_hit_64
		my ($individual_id, $snp_id, $fasta_query_type, $fasta_query_length) = split(/_/, $fasta_query_seq_id);
# 		$fasta_query_regex{$fasta_query_sequence}{$fasta_query_length} = $fasta_query_seq_id;
		$fasta_query_regex{$individual_id}{$fasta_query_sequence} = $fasta_query_seq_id;
		push(@fasta_query_lengths, $fasta_query_length);
	}
	
	my $fasta_query_min = min(@fasta_query_lengths);
	$fasta_query_range_lengths{$individual_name}{"min"} = $fasta_query_min;
	
	my $fasta_query_max = max(@fasta_query_lengths);
	$fasta_query_range_lengths{$individual_name}{"max"} = $fasta_query_max;
	
	# Clean out the sequence I/O object.
	$seqio = ();
	
}

# foreach my $individual_name (sort keys %fasta_query_range_lengths){
# 	if(defined($fasta_query_range_lengths{$individual_name}{"min"}) and defined($fasta_query_range_lengths{$individual_name}{"max"})){
# 		my ($fasta_query_min, $fasta_query_max) = ($fasta_query_range_lengths{$individual_name}{"min"}, $fasta_query_range_lengths{$individual_name}{"max"});
# 		warn join(" ", $individual_name, $fasta_query_min, $fasta_query_max) . "\n";
# 	}
# }

my $blastn_output_dir = "/home/cookeadmin/workspace/GBS_data-08-10-2013/Papilio_GBS_Data/Papilio_UNEAK_GBS-2014-12-04/UNEAK_BLASTN_FILES";
my $target_fasta_output_dir = "/home/cookeadmin/workspace/GBS_data-08-10-2013/Papilio_GBS_Data/Papilio_UNEAK_GBS-2014-12-04/UNEAK_FASTA_TARGET_FILES";
my ($fasta_target_files, $fasta_target_file_counter) = find_files($target_fasta_output_dir, "fasta");
foreach my $filename (sort keys %{$fasta_target_files}){
	warn $filename . "\n";
	my $fasta_target_infile = $fasta_target_files->{$filename};
	my $fasta_filename = fileparse($fasta_target_infile, qr/\.fasta/);
	my @split_target_file = split(/_/, $fasta_filename);
	my $individual_name = $split_target_file[0];
	
	my $regex_alignment_outfile = join('/', $blastn_output_dir, $individual_name . ".uneak_blastn.tsv");
	open(OUTFILE, ">$regex_alignment_outfile") or die "Couldn't open file $regex_alignment_outfile for writting, $!";
	print OUTFILE join("\t", "query_name", "target_name", "align_length", "query_start", "query_end", "target_start", "target_end") . "\n";
	my %regex_align_hits = ();
	
	if(defined($fasta_query_range_lengths{$individual_name}{"min"}) and defined($fasta_query_range_lengths{$individual_name}{"max"})){
		
		my ($fasta_query_min, $fasta_query_max) = ($fasta_query_range_lengths{$individual_name}{"min"}, $fasta_query_range_lengths{$individual_name}{"max"});
		my $seqio = Bio::SeqIO->new(-file => $fasta_target_infile, '-format' => 'Fasta');
		while(my $seq_entry = $seqio->next_seq) {

			my $fasta_target_seq_id = $seq_entry->id;
			my $fasta_target_sequence = $seq_entry->seq;

			my ($query_name, $align_length, $query_start, $query_end, $target_start, $target_end);
			my $fasta_target_length_count = 1;
			my $align_found = "false";
			for(my $i = $fasta_query_max; $i > $fasta_query_min; $i--){
				if(defined($fasta_query_regex{$individual_name}{$fasta_target_sequence})){
					$query_name = $fasta_query_regex{$individual_name}{$fasta_target_sequence};
					print "$i eq $fasta_query_max\n";
					$target_start = 1;
					$target_end = $i;
					$align_length = $i;
					print join("\t", $align_length, $target_start, $target_end) . "\n";
					$align_found = "true";
					$query_start = 1;
					$query_end = $i;
					warn join("\t", $query_name, $i, $fasta_target_sequence) . "\n";
					last;
				}
				
				$query_start = 1;
				$query_end = ($fasta_query_max - $fasta_target_length_count);
				
				$fasta_target_sequence = get_subseq($fasta_target_sequence, $query_start, $query_end);
				warn join("\t", ($fasta_query_max - $fasta_target_length_count), $fasta_target_sequence) . "\n";
				$fasta_target_length_count++;
			}
			if($align_found eq "true"){
				
				my $regex_alignment = join("\t", $query_name, $fasta_target_seq_id, $align_length, $query_start, $query_end, $target_start, $target_end);
				print OUTFILE $regex_alignment . "\n";
			}
		}
		# Clean out the sequence I/O object.
		$seqio = ();
	}
	close(OUTFILE) or die "Couldn't close file $regex_alignment_outfile";
	

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