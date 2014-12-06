#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;
use Bio::SeqIO;

# perl uneak_snps_cleanup.pl -i ~/workspace/GBS_data-08-10-2013/PROJECT_LEADER_DIR/JULIAN_DUPUIS -p Papilio_UNEAK_GBS -w ~/workspace/GBS_data-08-10-2013/Papilio_GBS_Data/Papilio_notrim-2014-11-18 -o ~/workspace/GBS_data-08-10-2013/Papilio_GBS_Data/uneak_snps_cleanup
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

my ($makeblastdb, $blastn, $cdhit_est, $cdbfasta, $cdbyank);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$blastn				= '/usr/local/bin/blastn';
$cdhit_est				= '/usr/bin/cdhit-est';
$cdbfasta				= '/usr/bin/cdbfasta';
$cdbyank				= '/usr/bin/cdbyank';

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

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create the project output directory if it doesn't already exist.
my $project_dir = join('/', $output_dir, $project_name);
unless(-d $project_dir){
	mkdir($project_dir, 0777) or die "Can't make directory: $!";
}

# my $parsed_tags_by_taxa_outfile = join('/', $project_dir, join("_", $project_name, "parsed_tags_by_taxa.txt"));
# unless(-s $parsed_tags_by_taxa_outfile){
# 	my $tags_by_taxa_infile = join('/', $uneak_project_dir, "tagsByTaxa", "tbt.bin.txt");
# 	open(OUTFILE, ">$parsed_tags_by_taxa_outfile") or die "Couldn't open file $parsed_tags_by_taxa_outfile for writting, $!";
# 	open(INFILE, "<$tags_by_taxa_infile") or die "Couldn't open file $tags_by_taxa_infile for reading, $!";
# 	my $i = 0;
# 	while(<INFILE>){
# 		chomp $_;
# 		warn $_ . "\n";
# 		if($i eq 1){
# 			my @split_tbt_entries = split(/\t/, $_);
# 			my $j = 0;
# 			my @individual_ids = ();
# 			foreach my $tbt_entry (@split_tbt_entries){
# 				my @split_individual_name = split(/_/, $tbt_entry);
# 				my $individual_id = $split_individual_name[0];
# 				$individual_ids[$j] = $individual_id;
# 				$j++;
# 			}
# 			print OUTFILE join("\t", "sequence", "length", @individual_ids) . "\n";
# 		}elsif($i > 1){
# 			print OUTFILE $_ . "\n";
# 		}
# 		$i++;
# 	}
# 	close(INFILE) or die "Couldn't close file $tags_by_taxa_infile";
# 	close(OUTFILE) or die "Couldn't close file $parsed_tags_by_taxa_outfile";
# }

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
	# Create output directory if it doesn't already exist.
	my $individual_target_output_dir = join('/', $fasta_target_output_dir, $individual_id);
	unless(-d $individual_target_output_dir){
		mkdir($individual_target_output_dir, 0777) or die "Can't make directory: $!";
	}
	my $fasta_target_outfile = join("/", $individual_target_output_dir, join("_", $individual_id, "target") . ".fasta");
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
				my $barcode_tassel_trimmed_sequence = join("", $barcode_sequence, $tassel_trimmed_sequence);
				
				if(($barcode_tassel_trimmed_sequence !~ m/N/) and ($barcode_sequence =~ m/$barcode/)){
					print OUTFILE join("\n", join("", ">", join("_", $fasta_filename, $fastq_header)), $tassel_trimmed_sequence) . "\n";
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
	my $snp_id = "";
	if($hap_map_reference_id =~ /^[A-Z]+([0-9]+)$/){
		$snp_id = $1;
	}

	$hap_map_fasta_seqs{$snp_id}{$type} = join("\t", join("_", $hap_map_reference_id, $type, $length), get_subseq($sequence, 1, $length));
	
	# sequence for snp position search.
	push(@{$snp_fasta_seqs{$hap_map_reference_id}}, $sequence);

}

# Clean out the sequence I/O object.
$seqio = ();

# Get the SNP position
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
			my $snp_id = "";
			if($hap_map_reference_id =~ /^[A-Z]+([0-9]+)$/){
				$snp_id = $1;
			}
			$hap_map_hmc_counts{$individual_id}{$snp_id}{"query"} = $query_count;
			$hap_map_hmc_counts{$individual_id}{$snp_id}{"hit"} = $hit_count;
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
	my $fasta_query_outfile = join("/", $query_fasta_output_dir, join("_", $individual_id, "query") . ".fasta");
	warn join(" ", "Generating", join("_", $individual_id, "query") . ".fasta", "query output file.....\n");
	open(OUTFILE, ">$fasta_query_outfile") or die "Couldn't open file $fasta_query_outfile for writting, $!";
	my %cdhit_est_seq_lengths = ();
	foreach my $snp_id (sort {$a <=> $b} keys %{$hap_map_hmc_counts{$individual_id}}){
	
		if(defined($hap_map_hmc_counts{$individual_id}{$snp_id}{"query"})){
			# number of snp query sequences present with this allele.
			my $num_snp_query_seqs = $hap_map_hmc_counts{$individual_id}{$snp_id}{"query"};
			if($num_snp_query_seqs ne 0){
				my ($sequence_query_id, $query_sequence) = split(/\t/, $hap_map_fasta_seqs{$snp_id}{"query"});
				my ($hap_map_reference_id, $type, $length) = split(/_/, $sequence_query_id);
				$cdhit_est_seq_lengths{$length}++;
				print OUTFILE join("\n", join("", ">", join("_", $individual_id, $sequence_query_id)), $query_sequence) . "\n";
				
			}
		}
		
		if(defined($hap_map_hmc_counts{$individual_id}{$snp_id}{"hit"})){
			# number of snp hit sequences present with this allele.
			my $num_snp_hit_seqs = $hap_map_hmc_counts{$individual_id}{$snp_id}{"hit"};
			if($num_snp_hit_seqs ne 0){
				my ($sequence_hit_id, $hit_sequence) = split(/\t/, $hap_map_fasta_seqs{$snp_id}{"hit"});
				my ($hap_map_reference_id, $type, $length) = split(/_/, $sequence_hit_id);
				$cdhit_est_seq_lengths{$length}++;
				print OUTFILE join("\n", join("", ">", join("_", $individual_id, $sequence_hit_id)), $hit_sequence) . "\n";
				
			}
		}
	}
	close(OUTFILE) or die "Couldn't close file $fasta_query_outfile";
	
	my $fasta_target_infile = $fasta_filenames{$individual_id};
	my $num_iterations = 0;
	
 	my (@cdhit_est_clstr_fasta_outfile_list, @full_clustered_list);
 	my ($cdhit_est_clustered_fasta_outfile, $cdhit_outfile, $clustered_list, $recluster_list);
	foreach my $query_seq_length (sort {$b <=> $a} keys %cdhit_est_seq_lengths){
		print $query_seq_length . "\n";
		if($num_iterations >= 1){
			
			my $cdhit_est_reclustered_fasta_outfile = cdbfasta_yank($cdhit_outfile, $query_seq_length, $recluster_list);
			($cdhit_est_clustered_fasta_outfile, $cdhit_outfile, $clustered_list, $recluster_list) = cluster_cdhit_est($cdhit_est_reclustered_fasta_outfile, $query_seq_length);
			push(@cdhit_est_clstr_fasta_outfile_list, $cdhit_est_clustered_fasta_outfile);
			push(@full_clustered_list, @{$clustered_list});
		}else{
			
			($cdhit_est_clustered_fasta_outfile, $cdhit_outfile, $clustered_list, $recluster_list) = cluster_cdhit_est($fasta_target_infile, $query_seq_length);
			push(@cdhit_est_clstr_fasta_outfile_list, $cdhit_est_clustered_fasta_outfile);
			push(@full_clustered_list, @{$clustered_list});
		}
		$num_iterations++;
	}
	my ($fasta_filename, $fasta_dir) = fileparse($fasta_target_infile, qr/\.fasta/);
	my $cdhit_est_cluster_report_outfile = join("/", $fasta_dir, join("", $fasta_filename . ".fasta.cd-hit.all.report"));
	open(OUTFILE, ">$cdhit_est_cluster_report_outfile") or die "Couldn't open file $cdhit_est_cluster_report_outfile for writting, $!";
	my $num_clusters = 1;
	foreach my $cluster (@full_clustered_list){
		print OUTFILE join("\t", "Cluster$num_clusters", $cluster) . "\n"
	}
	close(OUTFILE) or die "Couldn't close file $cdhit_est_cluster_report_outfile";
# 	generate_uneak_blastn($individual_id, $fasta_query_outfile, $fasta_filenames{$individual_id}, $max_target_seqs, $uneak_sequence_length, $blast_num_cpu, $blastn_output_dir);
 	push(@individual_ids, $individual_id);
 	%cdhit_est_seq_lengths = ();
}


# # Parsing uneak_blastn.tsv files to get sequence counts for each SNP and extracting the GBS fastq sequence headers to grab the trimming data
# warn "Parsing uneak_blastn.tsv files to get sequence counts for each SNP and extracting the GBS fastq sequence headers to grab the trimming data.....\n";
# my ($blastn_files, $blastn_file_counter) = find_files($blastn_output_dir, "uneak_blastn.tsv");
# my $blastn_outfile = join('/', $project_dir, join("_", $project_name, "uneak_blastn.tsv"));
# open(OUTFILE, ">$blastn_outfile") or die "Couldn't open file $blastn_outfile for writting, $!";
# print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch",
# 	"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
# my %blastn_snps_counter = ();
# my %blastn_snps_gbs_references = ();
# foreach my $blastn_filename (sort keys %{$blastn_files}){
# 	my $blastn_infile = $blastn_files->{$blastn_filename};
# 	open(INFILE, "<$blastn_infile") or die "Couldn't open file $blastn_infile for reading, $!";
# 	my $i = 0;
# 	while(<INFILE>){
# 		chomp $_;
# 		warn $_ . "\n";
# 		if($i ne 0){
# 			my @split_blastn_hit =  split(/\t/, $_);
# 			my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
# 			$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @split_blastn_hit;
# 			#JRD229_TP1731_hit
# 			my ($individual_id, $hap_map_reference_id, $query_type) = split(/_/, $query_name);
# 			
# 			# Getting the SNP id from the hap map reference id so that we can sort in ascending order easier.
# 			my $snp_id = "";
# 			if($hap_map_reference_id =~ /^[A-Z]+([0-9]+)$/){
# 				$snp_id = $1;
# 			}
# 			$blastn_snps_counter{$snp_id}{$individual_id}{$query_type}++;
# 			push(@{$blastn_snps_gbs_references{$snp_id}{$individual_id}{$query_type}}, $target_name);
# 			print OUTFILE $_ . "\n";
# 		}
# 		$i++;
# 	}
# 	close(INFILE) or die "Couldn't close file $blastn_infile";
# }
# close(OUTFILE) or die "Couldn't close file $blastn_outfile";
# 
# # Generating the pretrim hap map counts file so that we can use it for generating the trimmed hap map counts file.
# my $pretrim_hap_map_hmc_outfile = join('/', $project_dir, join("_", $project_name, "pretrim_hap_map_hmc.txt"));
# open(OUTFILE, ">$pretrim_hap_map_hmc_outfile") or die "Couldn't open file $pretrim_hap_map_hmc_outfile for writting, $!";
# print OUTFILE join("\t", "SNP_ID", @individual_ids) . "\n";
# foreach my $snp_id (sort {$a <=> $b} keys %uneak_snp_positions){
# 	my @pretrim_hap_map_hmc_counts = ();
# 	foreach my $individual_id (sort {$a cmp $b} @individual_ids){
# 		my $snp_query_count;
# 		if(defined($blastn_snps_counter{$snp_id}{$individual_id}{"query"})){
# 			$snp_query_count = $blastn_snps_counter{$snp_id}{$individual_id}{"query"};
# 		}
# 		$snp_query_count = 0 unless(defined($snp_query_count));
# 		
# 		my $snp_hit_count;
# 		if(defined($blastn_snps_counter{$snp_id}{$individual_id}{"hit"})){
# 			$snp_hit_count = $blastn_snps_counter{$snp_id}{$individual_id}{"hit"};
# 		}
# 		$snp_hit_count = 0 unless(defined($snp_hit_count));
# 		
# 		my $snp_individual_counts = join("|", $snp_query_count, $snp_hit_count);
# 		push(@pretrim_hap_map_hmc_counts, $snp_individual_counts);
# 	}
# 	print OUTFILE join("\t", join("", "TP", $snp_id), @pretrim_hap_map_hmc_counts) . "\n";
# 	@pretrim_hap_map_hmc_counts = ();
# }
# close(OUTFILE) or die "Couldn't close file $pretrim_hap_map_hmc_outfile";


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

# makeblastdb -in ncbi_nr_db_2014-05-30_organisms.fasta -dbtype 'nucl' -out ncbi_nr_db_2014-05-30_organisms.fasta
sub makeblastdb_nuc{
	my $fastadb = shift;
	die "Error lost fastadb to makeblastdb" unless defined $fastadb;

	# format the database file into .nin .nsq .nhr files.
	my ($fastadbNIN, $fastadbNSQ, $fastadbNHR);
	$fastadbNIN = $fastadb . '.nin';
	$fastadbNSQ = $fastadb . '.nsq';
	$fastadbNHR = $fastadb . '.nhr';
	unless(-s $fastadbNIN and -s $fastadbNSQ and -s $fastadbNHR){
		warn "Calling makeblastdb for $fastadb....\n";
		my $makeblastdbCmd = "$makeblastdb -in $fastadb -dbtype nucl";
		warn $makeblastdbCmd . "\n\n";
		system($makeblastdb, 
			'-in', $fastadb, 
			'-dbtype', 'nucl'
		) == 0 or die "Error calling $makeblastdbCmd: $?";
	}

}

sub generate_uneak_blastn{
	my $individual_id = shift;
	die "Error lost name of the individual" unless defined $individual_id;
	
	my $fasta_query_infile = shift;
	die "Error lost fasta query input file" unless defined $fasta_query_infile;
	
	my $fasta_target_infile = shift;
	die "Error lost fasta database target input file" unless defined $fasta_target_infile;
	
	my $max_target_seqs = shift;
	die "Error lost maximum number of sequences" unless defined $max_target_seqs;
	
	my $uneak_sequence_length = shift;
	die "Error lost uneak sequence length" unless defined $uneak_sequence_length;
	
	my $blast_num_cpu = shift;
	die "Error lost number of cpus to allocate" unless defined $blast_num_cpu;
	
	my $blastn_output_dir = shift;
	die "Error lost uneak blastn output directory" unless defined $blastn_output_dir;
	
	my $fasta_target_name = fileparse($fasta_target_infile);
	
	makeblastdb_nuc($fasta_target_infile);

	my $uneak_blastn_tsv_outfile = join('/', $blastn_output_dir, $individual_id . ".uneak_blastn.tsv");
	unless(-s $uneak_blastn_tsv_outfile){
		warn join(" ", "Generating blast tab-delimited file", join("", $individual_id, ".uneak_blastn", ".tsv"), "using uneak sequence blastn.....") . "\n";
		my $uneakSequenceBlastnCmd  = "$blastn -query $fasta_query_infile -db $fasta_target_infile -task blastn -strand plus -word_size 11 -dust yes -evalue 1e-6 -max_target_seqs 10000 -perc_identity 100 -outfmt '6 qseqid salltitles qcovhsp pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
		warn $uneakSequenceBlastnCmd . "\n\n";
		open(OUTFILE, ">$uneak_blastn_tsv_outfile") or die "Couldn't open file $uneak_blastn_tsv_outfile for writting, $!";
                print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch",
                "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
		local (*UNEAK_BLASTN_OUT, *UNEAK_BLASTN_IN);
		my $pid = open2(\*UNEAK_BLASTN_OUT,\*UNEAK_BLASTN_IN, $uneakSequenceBlastnCmd) or die "Error calling open2: $!";
		close UNEAK_BLASTN_IN or die "Error closing STDIN to uneak sequence blastn process: $!";
		while(<UNEAK_BLASTN_OUT>){
			chomp $_;
			my @split_uneak_blastn_hit =  split(/\t/, $_);
			my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
			$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @split_uneak_blastn_hit;
			
			my ($individual, $hap_map_reference_id, $type, $sequence_length) = split(/_/, $query_name);
			if((($target_start eq 1) and ($target_end eq $sequence_length)) and ($sequence_length eq $align_length)){
				print OUTFILE join("\t", $query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch, 
					$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
			}
		}
		close UNEAK_BLASTN_OUT or die "Error closing STDOUT from uneak sequence blastn process: $!";
		wait;
		close(OUTFILE) or die "Couldn't close file $uneak_blastn_tsv_outfile";
	}
	return $uneak_blastn_tsv_outfile;
}

sub cluster_cdhit_est {
	my $fasta_infile = shift;
	die "Error lost fasta input file" unless defined $fasta_infile;

	my $cluster_sequence_length = shift;
	die "Error lost cluster sequence length" unless defined $cluster_sequence_length;
	
	my ($fasta_filename, $fasta_dir) = fileparse($fasta_infile, qr/\.fasta/);
	# Clustering via cdhit-est
	
	my $cdhit_outfile = join("/", $fasta_dir, join("_", $fasta_filename, $cluster_sequence_length . ".fasta.cd-hit"));
 	unless(-s $cdhit_outfile){
		my $cdhitESTCmd = "$cdhit_est -c 1 -G 1 -n 10 -d 0 -M 0 -T 7 -i $fasta_infile -o $cdhit_outfile";
		warn "$cdhitESTCmd\n\n";
		system($cdhit_est, 
			'-c',	1,
			'-G',	1,
			'-n',	10,
			'-d',	0,
			'-M',	0,
			'-T',	7,
			'-i',	$fasta_infile, 
			'-o',	$cdhit_outfile
		) == 0 or die "Error calling $cdhitESTCmd: $?";
 	}

	# Create a summary report from the CD-hit step
	my $cdhit_est_clstr_infile = join("", $cdhit_outfile, ".clstr");
	open(INFILE, "<$cdhit_est_clstr_infile") or die "Couldn't open file $cdhit_est_clstr_infile for reading, $!";
	my %cdhit_clusters = ();
	my %cdhit_clustered_seq = ();
	my $cluster_id;
	while(<INFILE>){
		chomp $_;
		warn $_ . "\n";
		if($_ =~ m/>Cluster\s(\d+)/){
			$cluster_id = $1;
		}elsif($_ =~ m/\d+\t(\d+)nt,\s>(.+)\.\.\. \*/){
			my ($cluster_seq_length, $sequence_id) = ($1, $2);
			$cdhit_clustered_seq{$cluster_id}{"sequence_id"} = join(" ", $sequence_id, "1:N:0:");
			$cdhit_clustered_seq{$cluster_id}{"length"} = $cluster_seq_length;
			$cdhit_clustered_seq{$cluster_id}{"counts"}++;
			push(@{$cdhit_clusters{$cluster_id}}, join(" ", $sequence_id, "1:N:0:"));
		}elsif($_ =~ m/\d+\t(\d+)nt,\s>(.+)\.\.\.\sat\s\+\/100\.00\%/){
			my ($cluster_seq_length, $sequence_id) = ($1, $2);
			push(@{$cdhit_clusters{$cluster_id}}, join(" ", $sequence_id, "1:N:0:"));
			$cdhit_clustered_seq{$cluster_id}{"counts"}++;
		}
		
	}
	close(INFILE) or die "Couldn't close file $cdhit_est_clstr_infile";
	

	my @recluster_list = ();
	my @clustered_list = ();
	my @cdbyank_clustered_list = ();
	foreach my $cluster_id (sort {$a <=> $b} keys %cdhit_clusters){
		if($cdhit_clustered_seq{$cluster_id}{"counts"} > 1){
			push(@clustered_list, join("\t", $cdhit_clustered_seq{$cluster_id}{"length"}, $cdhit_clustered_seq{$cluster_id}{"counts"}, $cdhit_clustered_seq{$cluster_id}{"sequence_id"}, join(",", @{$cdhit_clusters{$cluster_id}})));
			push(@cdbyank_clustered_list, $cdhit_clustered_seq{$cluster_id}{"sequence_id"});
		}else{
			push(@recluster_list, $cdhit_clustered_seq{$cluster_id}{"sequence_id"});
		}
	}
	
	my $cdhit_est_clustered_fasta_outfile = cdbfasta_yank($cdhit_outfile, $cdhit_clustered_seq{$cluster_id}{"length"}, \@cdbyank_clustered_list);
	return ($cdhit_est_clustered_fasta_outfile, $cdhit_outfile, \@clustered_list, \@recluster_list);
}

sub cdbfasta_yank{

	
	my $cdhit_infile = shift;
	die "Error lost cdhit input file" unless defined $cdhit_infile;

	my $cluster_sequence_length = shift;
	die "Error lost cluster sequence length" unless defined $cluster_sequence_length;
	
	my $cdbyank_clustered_list = shift;
	die "Error lost cdbyank clustered sequence id list" unless defined $cdbyank_clustered_list;
	# Make sure that the clustered results are indexed with CDBFASTA
	my $cdbfasta_infile = $cdhit_infile . '.cidx';
	unless(-s $cdbfasta_infile){
		# Create a CDBFasta index of the GBS sequences.
		my $cdbfastaCmd = "$cdbfasta $cdhit_infile";
		warn "$cdbfastaCmd\n\n";
		system($cdbfastaCmd) == 0 or die "Error calling $cdbfastaCmd: $?";
	}

	# Create a fasta file of GBS sequences which were clustered down into a non-redundant set.
	my $cdhit_est_clstr_fasta_outfile = $cdhit_infile . '.fna';
# 	unless (-s $cdhit_est_clstr_fasta_outfile){
		open(OUTFILE, ">$cdhit_est_clstr_fasta_outfile") or die "Couldn't open file $cdhit_est_clstr_fasta_outfile for writting, $!";
		foreach my $sequence_id (@{$cdbyank_clustered_list}){
			my $cdbyankCmd = "$cdbyank $cdbfasta_infile -a '" . $sequence_id . "'";
			warn "$cdbyankCmd\n";
		
			local (*CDBYANK_OUT, *CDBYANK_IN);
			my $pid = open2(\*CDBYANK_OUT,\*CDBYANK_IN, $cdbyankCmd) or die "Error calling open2: $!";
			close CDBYANK_IN or die "Error closing STDIN to cdbyank process: $!";	
			
			while (<CDBYANK_OUT>){
				chomp $_;
				if($_ =~ m/^>/){
					print OUTFILE $_ . "\n";
				}elsif($_ =~ m/^[ACGTRYKMSWBDHVN]+$/i){
					print OUTFILE get_subseq($_, 1, $cluster_sequence_length) . "\n";
				}
			}
			
			close CDBYANK_OUT or die "Error closing STDOUT from cdbyank process: $!";
			wait;
		}
		close(OUTFILE) or die "Couldn't close file $cdhit_est_clstr_fasta_outfile";
# 	}
	return $cdhit_est_clstr_fasta_outfile;
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
