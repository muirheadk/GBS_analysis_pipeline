#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;
use Bio::SeqIO;

# perl uneak_post_analysis.pl -i ~/workspace/GBS_data-08-10-2013/PROJECT_LEADER_DIR/JULIAN_DUPUIS -p Papilio_UNEAK_GBS -w ~/workspace/GBS_data-08-10-2013/Papilio_GBS_Data/Papilio_notrim-2014-11-18 -o ~/workspace/GBS_data-08-10-2013/Papilio_GBS_Data/uneak_snps_cleanup
my ($fastq_file_dir, $project_name, $uneak_project_dir, $gbs_sequence_length, $uneak_sequence_length, $blast_num_cpu, $output_dir);
GetOptions(
	'i=s'    => \$fastq_file_dir,
	'p=s'    => \$project_name,
	'w=s'    => \$uneak_project_dir,
	'l=s'    => \$gbs_sequence_length,
	's=s'    => \$uneak_sequence_length,
	'c=s'    => \$blast_num_cpu,
	'o=s'    => \$output_dir,
);


usage() unless (
	defined $fastq_file_dir
	and defined $project_name
	and defined $uneak_project_dir
	and defined $output_dir
);

$gbs_sequence_length = 100 unless defined $gbs_sequence_length;
$uneak_sequence_length = 64 unless defined $uneak_sequence_length;
$blast_num_cpu = 2 unless defined $blast_num_cpu;



sub usage {

die <<"USAGE";


Usage: $0 -i fastq_file_dir -p project_name -w uneak_project_dir -l gbs_sequence_length -s uneak_sequence_length -c blast_num_cpu -o output_dir

DESCRIPTION - 

OPTIONS:

-i fastq_file_dir - 

-p project_name -

-w uneak_project_dir -

-l gbs_sequence_length -

-s uneak_sequence_length -

-c blast_num_cpu -

-o output_dir -

USAGE
}

# The cut-site for MspI.
# 5' ---C   CGG--- 3'
# 3' ---GGC   C--- 5'
my $restriction_enzyme_cut_site = qr(CCGG);

# The first seven bases of the GBS common adapter sequencefor MspI.
my $common_adapter_prefix = qr(CCGAGAT);
					
# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create the project output directory if it doesn't already exist.
my $project_dir = join('/', $output_dir, $project_name);
unless(-d $project_dir){
	mkdir($project_dir, 0777) or die "Can't make directory: $!";
}

# Find all files in the specified directory with the extension *.fastq.
my ($fastq_files, $file_count) = find_fastq_files($fastq_file_dir);

# Create output directory if it doesn't already exist.
my $fasta_target_output_dir = join('/', $project_dir, "UNEAK_FASTA_TARGET_FILES");
unless(-d $fasta_target_output_dir){
	mkdir($fasta_target_output_dir, 0777) or die "Can't make directory: $!";
}

# Iterate through the files with the extension *.fastq.
my %fasta_filenames = ();
foreach my $fastq_filename (sort keys %{$fastq_files}){
	my $fastq_infile = $fastq_files->{$fastq_filename};
	my $fasta_filename = fileparse($fastq_infile, qr/\.fastq/);
	my ($individual_id, $barcode, $plate_num, $well_num) = split(/_/, $fasta_filename);
	my $fasta_target_outfile = join("/", $fasta_target_output_dir, join("_", $individual_id, "target") . ".fasta");
	$fasta_filenames{$individual_id} = $fasta_target_outfile;
	unless(-s $fasta_target_outfile){
		warn "Processing " . $fastq_infile . ".....\n";
		open(OUTFILE, ">$fasta_target_outfile") or die "Couldn't open file $fasta_target_outfile for writting, $!";
		open(INFILE, "<$fastq_infile") or die "Couldn't open file $fastq_infile for reading, $!";
		
		my $barcode_length = length($barcode);
		
		my ($fastq_header, $fastq_sequence, $fastq_plus, $fastq_quality_scores);
		my $i = 1;
		# Parse the fastq files for the fastq header, sequence, plus, and quality scores and reformat to *.fasta format so that we can use adaptor blastn.
		while(<INFILE>){
			chomp $_;
			#warn $_ . "\n";
			if(($_ =~ m/^\@[A-Za-z0-9-_]+:\d+:[A-Za-z0-9]+:\d+:\d+:\d+:\d+ \d:[A-Z]:\d:[ACGTRYKMSWBDHVN]*$/)
				and ($i eq 1)){ # @HWI-ST767:215:C30VBACXX:8:1101:1801:1484 1:N:0:
				$fastq_header = $_;
	#  				die $fastq_header;
			}elsif(($_ =~ m/^[ACGTRYKMSWBDHVN]+$/i) and ($i eq 2)){
				$fastq_sequence = $_;
	#  				die $fastq_sequence;
			}elsif(($_ =~ m/^\+$/) and ($i eq 3)){
				$fastq_plus = $_;
	#  				die $fastq_plus;
			}elsif(($_ =~ m/^.+$/) and (($i % 4) eq 0)){
				$fastq_quality_scores = $_;
	#   				die $fastq_quality_scores;
			}
			
			if(($i % 4) eq 0){
				
				die "Error: fastq_header is undefined" unless(defined($fastq_header));
				die "Error: fastq_sequence is undefined" unless(defined($fastq_sequence));
				die "Error: fastq_plus is undefined" unless(defined($fastq_plus));
				die "Error: fastq_quality_scores is undefined" unless(defined($fastq_quality_scores));
				
				my $fastq_sequence_length = length($fastq_sequence);
				my $fastq_quality_scores_length = length($fastq_quality_scores);
				
				die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($fastq_sequence_length ne $gbs_sequence_length);
				die "Error: $fastq_header: fastq_sequence_length=$fastq_sequence_length ne fastq_quality_scores_length=$fastq_quality_scores_length" if($fastq_sequence_length ne $fastq_quality_scores_length);

				my $tassel_trimmed_sequence = get_subseq(get_subseq($fastq_sequence, ($barcode_length + 1), $gbs_sequence_length), 1, $uneak_sequence_length);
				
				# need to see if there are Ns within the barcode + trimmed sequence.
				my $barcode_sequence = get_subseq($fastq_sequence, 1, $barcode_length);
# 				die join("\t", $barcode, $barcode_sequence);
				my $barcode_tassel_trimmed_sequence = join("", $barcode_sequence, $tassel_trimmed_sequence);
				
				if(($barcode_tassel_trimmed_sequence !~ m/N/) and ($barcode_sequence =~ m/$barcode/)){
					my $new_tassel_trimmed_sequence = $tassel_trimmed_sequence;
					# The cut-site for MspI.
					# 5' ---C   CGG--- 3'
					# 3' ---GGC   C--- 5'
					if($tassel_trimmed_sequence =~ m/$restriction_enzyme_cut_site/g){
						my $target_start = ($-[0] + 1);
						my $target_end = $+[0];
						my $tassel_trimmed_subsequence = get_subseq($tassel_trimmed_sequence, 1, $target_start);
						my $tassel_trimmed_subseq_length = length($tassel_trimmed_subsequence);
						my $padded_poly_A_length =  ($uneak_sequence_length - $tassel_trimmed_subseq_length);
 						my $padded_poly_A_seq = 'A' x $padded_poly_A_length;
						$new_tassel_trimmed_sequence = join("", $tassel_trimmed_subsequence, $padded_poly_A_seq);
						
						my $new_fastq_header = join("\001", $fastq_header, $fasta_filename, join("=", "length", $tassel_trimmed_subseq_length));
						warn join("\t", join("\n", join("\t", $new_fastq_header, "MspI cut-site"), $new_tassel_trimmed_sequence, $tassel_trimmed_sequence), $target_start, $target_end) . "\n";
						print OUTFILE join("\n", join("", ">", $new_fastq_header, $new_tassel_trimmed_sequence)) . "\n" if($tassel_trimmed_subseq_length > 32);
					}
					elsif($tassel_trimmed_sequence =~ m/$common_adapter_prefix/g){
						# CCGAGAT
						my $target_start = ($-[0] + 1);
						my $target_end = $+[0];
						my $tassel_trimmed_subsequence = get_subseq($tassel_trimmed_sequence, 1, ($target_start + 3));
						my $tassel_trimmed_subseq_length = length($tassel_trimmed_subsequence);
						my $padded_poly_A_length =  ($uneak_sequence_length - $tassel_trimmed_subseq_length);
 						my $padded_poly_A_seq = 'A' x $padded_poly_A_length;
						$new_tassel_trimmed_sequence = join("", $tassel_trimmed_subsequence, $padded_poly_A_seq);
						
						my $new_fastq_header = join("\001", $fastq_header, $fasta_filename, join("=", "length", $tassel_trimmed_subseq_length));
						warn join("\t", join("\n", join("\t", $new_fastq_header, "MspI common adapter"), $new_tassel_trimmed_sequence, $tassel_trimmed_sequence), $target_start, $target_end) . "\n";
						print OUTFILE join("\n", join("", ">", $new_fastq_header, $new_tassel_trimmed_sequence)) . "\n" if($tassel_trimmed_subseq_length > 32);
					}
				}
				$i = 1;
				
			}else{
				$i++;
			}
		}
		close(INFILE) or die "Couldn't close file $fastq_infile";
		close(OUTFILE) or die "Couldn't close file $fasta_target_outfile";
	}
}

# Parsing HapMap.fas.txt to obtain SNP position within UNEAK sequence.
warn "Parsing HapMap.fas.txt to obtain SNP position within UNEAK sequence.....\n";
my $hap_map_fasta_infile = join('/', $uneak_project_dir, "hapMap", "HapMap.fas.txt");
my $seqio = Bio::SeqIO->new(-file => $hap_map_fasta_infile, '-format' => 'Fasta');
my %snp_fasta_seqs = ();
my %hap_map_fasta_seqs = ();
while(my $seq_entry = $seqio->next_seq) {

	my $seq_id = $seq_entry->id;
	my $sequence = $seq_entry->seq;
	my $sequence_desc = $seq_entry->desc;

	my ($hap_map_reference_id, $type, $length) = split(/_/, $seq_id);
	my $sequence_type_length = join("_", $type, $length);
	
	# sequence for blastn runs.
	my $hap_map_fasta_id = join("_", $hap_map_reference_id, $type);
	$hap_map_fasta_seqs{$hap_map_fasta_id} = get_subseq($sequence, 1, $length);
	
	# sequence for snp position search.
	push(@{$snp_fasta_seqs{$hap_map_reference_id}}, $sequence);

}

# Clean out the sequence I/O object.
$seqio = ();

my %uneak_snp_positions = ();
foreach my $hap_map_reference_id (sort keys %snp_fasta_seqs){
	my @sequence1 = split('', @{$snp_fasta_seqs{$hap_map_reference_id}}[0]);
	my @sequence2 = split('', @{$snp_fasta_seqs{$hap_map_reference_id}}[1]);
	
	for(my $i = 0; $i < $uneak_sequence_length; $i++){
		if($sequence1[$i] ne $sequence2[$i]){
			my $snp_position = ($i + 1);
			my $snp_allele1 = $sequence1[$i];
			my $snp_allele2 = $sequence2[$i];
			my $snp_alleles = join("/", $snp_allele1, $snp_allele2);
			
			# Getting the SNP id from the hap map reference id so that we can sort in ascending order easier.
			my $snp_id = "";
			if($hap_map_reference_id =~ /^[A-Z]+([0-9]+)$/){
				$snp_id = $1;
			}
			$uneak_snp_positions{$snp_id} = join("\t", $hap_map_reference_id, $snp_alleles, $snp_position);
			last;
		}
	}
}

# Generating the UNEAK snp_positions.txt output file
warn "Generating the UNEAK snp_positions.txt output file.....\n";
my $snp_positions_outfile = join('/', $project_dir, join("_", $project_name, "snp_positions.txt"));
unless(-s $snp_positions_outfile){
	open(OUTFILE, ">$snp_positions_outfile") or die "Couldn't open file $snp_positions_outfile for writting, $!";
	print OUTFILE join("\t", "SNP_ID", "alleles", "position") . "\n";
	foreach my $snp_id (sort {$a <=> $b} keys %uneak_snp_positions){
		print OUTFILE $uneak_snp_positions{$snp_id} . "\n";
	}
	close(OUTFILE) or die "Couldn't close file $snp_positions_outfile";
}


# Parsing HapMap.hmc.txt to obtain sequence counts
warn "Parsing HapMap.hmc.txt to obtain sequence counts.....\n";
my $hap_map_hmc_infile = join('/', $uneak_project_dir, "hapMap", "HapMap.hmc.txt");
open(INFILE, "<$hap_map_hmc_infile") or die "Couldn't open file $hap_map_hmc_infile for reading, $!";
my $i = 0;
my %individual_ids = ();
my %hap_map_hmc_counts = ();
while(<INFILE>){
	chomp $_;
#  	warn $_ . "\n";
	if($i eq 0){
		my @split_hmc_entries = split(/\t/, $_);
		my @hap_map_hmc_header = ();
		push(@hap_map_hmc_header, $split_hmc_entries[0]);
		for(my $j = 1; $j < (scalar(@split_hmc_entries) - 5); $j++){
# 			warn $split_hmc_entries[$j] . "\n";
			my @split_individual_name = split(/_/, $split_hmc_entries[$j]);
			my $individual_id = $split_individual_name[0];
			$individual_ids{$j} = $individual_id;
			push(@hap_map_hmc_header, $individual_id);
		}
		@hap_map_hmc_header = ();
	}elsif($i ne 0){
		my @split_hmc_entries = split(/\t/, $_);
		my $hap_map_reference_id = $split_hmc_entries[0];
		my @hap_map_hmc_entries = ();
		push(@hap_map_hmc_entries, $hap_map_reference_id);
		for(my $j = 1; $j < (scalar(@split_hmc_entries) - 5); $j++){
			my $individual_id = $individual_ids{$j};
			my @split_sequence_counts = split(/\|/, $split_hmc_entries[$j]);
			
			push(@hap_map_hmc_entries, $split_hmc_entries[$j]);
			my $query_count = $split_sequence_counts[0];
			my $hit_count = $split_sequence_counts[1];
			my $query_id = join("_", $hap_map_reference_id, "query");
			my $hit_id = join("_", $hap_map_reference_id, "hit");
			$hap_map_hmc_counts{$individual_id}{$query_id} = $query_count;
			$hap_map_hmc_counts{$individual_id}{$hit_id} = $hit_count;
		}
		@hap_map_hmc_entries = ();
	}
	$i++;
}
close(INFILE) or die "Couldn't close file $hap_map_hmc_infile";

# Create output directory if it doesn't already exist.
my $blastn_output_dir = join('/', $project_dir, "UNEAK_BLASTN_FILES");
unless(-d $blastn_output_dir){
	mkdir($blastn_output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $query_fasta_output_dir = join('/', $project_dir, "UNEAK_FASTA_QUERY_FILES");
unless(-d $query_fasta_output_dir){
	mkdir($query_fasta_output_dir, 0777) or die "Can't make directory: $!";
}



my @individual_ids = ();
foreach my $individual_id (sort {$a cmp $b} keys %hap_map_hmc_counts){
	my $max_target_seqs = 0;
	my $num_query_seqs = 0;
	my %query_fasta_sequences = ();
	foreach my $sequence_id (sort keys %{$hap_map_hmc_counts{$individual_id}}){
		# number of sequences present with this allele added together to use for the max_target_seqs parameter for blastn.
		$max_target_seqs += $hap_map_hmc_counts{$individual_id}{$sequence_id};
		
		# number of sequences present with this allele.
		my $num_snp_seqs = $hap_map_hmc_counts{$individual_id}{$sequence_id};
		if($num_snp_seqs ne 0){
			
			my ($hap_map_reference_id, $hap_map_reference_type) = split(/_/, $sequence_id);
			
			# Get the length of the sequence
			my $hap_map_fasta_sequence_length = length($hap_map_fasta_seqs{$sequence_id});
			
			my $snp_id = "";
			if($hap_map_reference_id =~ /^[A-Z]+([0-9]+)$/){
				$snp_id = $1;
			}
			
			my $padded_poly_A_length =  ($uneak_sequence_length - $hap_map_fasta_sequence_length);
 			my $padded_poly_A_seq = 'A' x $padded_poly_A_length;
			my $hap_map_padded_sequence = join("", $hap_map_fasta_seqs{$sequence_id}, $padded_poly_A_seq);
			$query_fasta_sequences{$snp_id}{$hap_map_reference_type} = join("\n", join("", ">", join("_", $individual_id, $sequence_id, $hap_map_fasta_sequence_length)), $hap_map_padded_sequence);
			
		}
		$num_query_seqs++;
	}
	
	
	my $fasta_query_outfile = join("/", $query_fasta_output_dir, join("_", $individual_id, "query") . ".fasta");
	open(OUTFILE, ">$fasta_query_outfile") or die "Couldn't open file $fasta_query_outfile for writting, $!";
	foreach my $snp_id (sort {$a <=> $b} keys %query_fasta_sequences){
		foreach my $hap_map_reference_type (sort {$a cmp $b} keys %{$query_fasta_sequences{$snp_id}}){
			print OUTFILE $query_fasta_sequences{$snp_id}{$hap_map_reference_type} . "\n";
		}
	}
	close(OUTFILE) or die "Couldn't close file $fasta_query_outfile";
	
	$max_target_seqs += $num_query_seqs;
	
	push(@individual_ids, $individual_id);
}

# Parsing uneak_blastn.tsv files to get sequence counts for each SNP and extracting the GBS fastq sequence headers to grab the trimming data
warn "Parsing uneak_blastn.tsv files to get sequence counts for each SNP and extracting the GBS fastq sequence headers to grab the trimming data.....\n";
my ($blastn_files, $blastn_file_counter) = find_files($blastn_output_dir, "uneak_blastn.tsv");
my $blastn_outfile = join('/', $project_dir, join("_", $project_name, "uneak_blastn.tsv"));
open(OUTFILE, ">$blastn_outfile") or die "Couldn't open file $blastn_outfile for writting, $!";
print OUTFILE join("\t", "query_name", "target_name", "align_length", "query_start", "query_end", "target_start", "target_end") . "\n";
my %blastn_snps_counter = ();
my %blastn_snps_gbs_references = ();
foreach my $blastn_filename (sort keys %{$blastn_files}){
	my $blastn_infile = $blastn_files->{$blastn_filename};
	open(INFILE, "<$blastn_infile") or die "Couldn't open file $blastn_infile for reading, $!";
	my $i = 0;
	while(<INFILE>){
		chomp $_;
		warn $_ . "\n";
		if($i ne 0){
			my @split_blastn_hit =  split(/\t/, $_);
			my ($query_name, $target_name, $align_length, $query_start, $query_end, $target_start, $target_end) = @split_blastn_hit;
			#JRD229_TP1731_hit
			my ($individual_id, $hap_map_reference_id, $query_type) = split(/_/, $query_name);
			
			# Getting the SNP id from the hap map reference id so that we can sort in ascending order easier.
			my $snp_id = "";
			if($hap_map_reference_id =~ /^[A-Z]+([0-9]+)$/){
				$snp_id = $1;
			}
			$blastn_snps_counter{$snp_id}{$individual_id}{$query_type}++;
			push(@{$blastn_snps_gbs_references{$snp_id}{$individual_id}{$query_type}}, $target_name);
			print OUTFILE $_ . "\n";
		}
		$i++;
	}
	close(INFILE) or die "Couldn't close file $blastn_infile";
}
close(OUTFILE) or die "Couldn't close file $blastn_outfile";

# Generating the pretrim hap map counts file so that we can use it for generating the trimmed hap map counts file.
my $pretrim_hap_map_hmc_outfile = join('/', $project_dir, join("_", $project_name, "pretrim_hap_map_hmc.txt"));
open(OUTFILE, ">$pretrim_hap_map_hmc_outfile") or die "Couldn't open file $pretrim_hap_map_hmc_outfile for writting, $!";
print OUTFILE join("\t", "SNP_ID", @individual_ids) . "\n";
foreach my $snp_id (sort {$a <=> $b} keys %uneak_snp_positions){
	my @pretrim_hap_map_hmc_counts = ();
	foreach my $individual_id (sort {$a cmp $b} @individual_ids){
		my $snp_query_count;
		if(defined($blastn_snps_counter{$snp_id}{$individual_id}{"query"})){
			$snp_query_count = $blastn_snps_counter{$snp_id}{$individual_id}{"query"};
		}
		$snp_query_count = 0 unless(defined($snp_query_count));
		
		my $snp_hit_count;
		if(defined($blastn_snps_counter{$snp_id}{$individual_id}{"hit"})){
			$snp_hit_count = $blastn_snps_counter{$snp_id}{$individual_id}{"hit"};
		}
		$snp_hit_count = 0 unless(defined($snp_hit_count));
		
		my $snp_individual_counts = join("|", $snp_query_count, $snp_hit_count);
		push(@pretrim_hap_map_hmc_counts, $snp_individual_counts);
	}
	print OUTFILE join("\t", join("", "TP", $snp_id), @pretrim_hap_map_hmc_counts) . "\n";
	@pretrim_hap_map_hmc_counts = ();
}
close(OUTFILE) or die "Couldn't close file $pretrim_hap_map_hmc_outfile";


sub find_fastq_files{
    
	my $fastq_file_dir = shift;
	die "Error lost fastq file directory" unless defined $fastq_file_dir;

	my (%fastq_files, $file_count);
	$file_count = 0;
	opendir(DIR, $fastq_file_dir) || die "Error in opening dir $fastq_file_dir\n";
	while( my $file_name = readdir(DIR)){
		my $fastq_file_name = join('/', $fastq_file_dir, $file_name) if ($file_name =~ m/\.fastq$/);
		warn "$fastq_file_name\n" if ($file_name =~ m/\.fastq$/);
		$fastq_files{$file_name} = $fastq_file_name if ($file_name =~ m/\.fastq$/);
		$file_count++ if ($file_name =~ m/\.fastq$/);
	}
	closedir(DIR);
	return (\%fastq_files, $file_count);
}

sub find_blastdb_files{
    
	my $blastdb_dir = shift;
	die "Error lost blastdb file directory" unless defined $blastdb_dir;
	
	my @blastdb_files = ();
	my $blastdb_files_counter = 0;
	opendir(DIR, $blastdb_dir) || die "Error in opening dir $blastdb_dir\n";
	while( my $file_name = readdir(DIR)){
		my $blastdb_file_name = join('/', $blastdb_dir, $file_name) if ($file_name =~ m/\.fasta$/ or $file_name =~ m/\.nin$/ or $file_name =~ m/\.nsq$/ or $file_name =~ m/\.nhr$/);
		warn "$blastdb_file_name\n" if ($file_name =~ m/\.fasta$/ or $file_name =~ m/\.nin$/ or $file_name =~ m/\.nsq$/ or $file_name =~ m/\.nhr$/);
		push(@blastdb_files, $blastdb_file_name) if ($file_name =~ m/\.fasta$/ or $file_name =~ m/\.nin$/ or $file_name =~ m/\.nsq$/ or $file_name =~ m/\.nhr$/);
		$blastdb_files_counter++ if ($file_name =~ m/\.fasta$/ or $file_name =~ m/\.nin$/ or $file_name =~ m/\.nsq$/ or $file_name =~ m/\.nhr$/);
	}
	closedir(DIR);
	return (\@blastdb_files, $blastdb_files_counter);
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
