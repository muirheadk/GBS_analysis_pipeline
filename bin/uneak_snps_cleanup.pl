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

my ($makeblastdb, $blastn);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$blastn				= '/usr/local/bin/blastn';

sub usage {

die <<"USAGE";


Usage: $0 -i fastq_file_dir -p project_name -w uneak_project_dir -l gbs_sequence_length -c blast_num_cpu -o output_dir

DESCRIPTION - 

OPTIONS:

-i fastq_file_dir - 

-p project_name -

-w uneak_project_dir -

-l gbs_sequence_length -

-c blast_num_cpu -

-o output_dir -

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create output directory if it doesn't already exist.
my $project_dir = join('/', $output_dir, $project_name);
unless(-d $project_dir){
	mkdir($project_dir, 0777) or die "Can't make directory: $!";
}

my $parsed_tags_by_taxa_outfile = join('/', $project_dir, join("_", $project_name, "parsed_tags_by_taxa.txt"));
unless(-s $parsed_tags_by_taxa_outfile){
	my $tags_by_taxa_infile = join('/', $uneak_project_dir, "tagsByTaxa", "tbt.bin.txt");
	open(OUTFILE, ">$parsed_tags_by_taxa_outfile") or die "Couldn't open file $parsed_tags_by_taxa_outfile for writting, $!";
	open(INFILE, "<$tags_by_taxa_infile") or die "Couldn't open file $tags_by_taxa_infile for reading, $!";
	my $i = 0;
	while(<INFILE>){
		chomp $_;
		warn $_ . "\n";
		if($i eq 1){
			my @split_tbt_entries = split(/\t/, $_);
			my $j = 0;
			my @individual_ids = ();
			foreach my $tbt_entry (@split_tbt_entries){
				my @split_individual_name = split(/_/, $tbt_entry);
				my $individual_id = $split_individual_name[0];
				$individual_ids[$j] = $individual_id;
				$j++;
			}
			print OUTFILE join("\t", "sequence", "length", @individual_ids) . "\n";
		}elsif($i > 1){
			print OUTFILE $_ . "\n";
		}
		$i++;
	}
	close(INFILE) or die "Couldn't close file $tags_by_taxa_infile";
	close(OUTFILE) or die "Couldn't close file $parsed_tags_by_taxa_outfile";
}

# Find all files in the specified directory with the extension *.fastq.
my ($fastq_files, $file_count) = find_fastq_files($fastq_file_dir);

# Create output directory if it doesn't already exist.
my $fasta_output_dir = join('/', $project_dir, "UNEAK_FASTA_FILES");
unless(-d $fasta_output_dir){
	mkdir($fasta_output_dir, 0777) or die "Can't make directory: $!";
}

# Iterate through the files with the extension *.fastq.
my %fasta_filenames = ();
foreach my $fastq_filename (sort keys %{$fastq_files}){
	my $fastq_infile = $fastq_files->{$fastq_filename};
	my $fasta_filename = fileparse($fastq_infile, qr/\.fastq/);
	my ($individual_id, $barcode, $plate_num, $well_num) = split(/_/, $fasta_filename);
	my $fasta_target_outfile = join("/", $fasta_output_dir, $individual_id . ".fasta");
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
				
				print OUTFILE join("\n", join("", ">", join("_", $fasta_filename, $fastq_header)), $tassel_trimmed_sequence) . "\n";
				
				$i = 1;
				
			}else{
				$i++;
			}
		}
		close(INFILE) or die "Couldn't close file $fastq_infile";
		close(OUTFILE) or die "Couldn't close file $fasta_target_outfile";
	}
}
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
	$hap_map_fasta_seqs{$hap_map_fasta_id} = $sequence;
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
			$uneak_snp_positions{$hap_map_reference_id} = $snp_position;
			last;
		}
	}
}

warn "Parsing HapMap.hmp.txt to store SNP positions within UNEAK parsed_hap_map_hmp.txt file.....\n";
my $parsed_hap_map_hmp_outfile = join('/', $project_dir, join("_", $project_name, "parsed_hap_map_hmp.txt"));
unless(-s $parsed_hap_map_hmp_outfile){
	my $hap_map_hmp_infile = join('/', $uneak_project_dir, "hapMap", "HapMap.hmp.txt");
	open(INFILE, "<$hap_map_hmp_infile") or die "Couldn't open file $hap_map_hmp_infile for reading, $!";
	open(OUTFILE, ">$parsed_hap_map_hmp_outfile") or die "Couldn't open file $parsed_hap_map_hmp_outfile for writting, $!";
	my $i = 0;
	while(<INFILE>){
		chomp $_;
		warn $_ . "\n";
		if($i eq 0){
			my @split_hmp_entries = split(/\t/, $_);
			for(my $j = 11; $j < scalar(@split_hmp_entries); $j++){
				my @split_individual_name = split(/_/, $split_hmp_entries[$j]);
				my $individual_id = $split_individual_name[0];
				$split_hmp_entries[$j] = $individual_id;
			}
			print OUTFILE join("\t", @split_hmp_entries) . "\n";
			
		}elsif($i ne 0){
			my @split_hmp_entries = split(/\t/, $_);
			my ($hap_map_reference_id, $hap_map_pos) = ($split_hmp_entries[0], $split_hmp_entries[3]);
			my $snp_position = $uneak_snp_positions{$hap_map_reference_id};
			$split_hmp_entries[3] = $snp_position;
			print OUTFILE join("\t", @split_hmp_entries) . "\n";
		}
		$i++;
	}
	close(INFILE) or die "Couldn't close file $hap_map_hmp_infile";
	close(OUTFILE) or die "Couldn't close file $parsed_hap_map_hmp_outfile";
}

my $hap_map_hmc_infile = join('/', $uneak_project_dir, "hapMap", "HapMap.hmc.txt");
open(INFILE, "<$hap_map_hmc_infile") or die "Couldn't open file $hap_map_hmc_infile for reading, $!";
my $i = 0;
my %individual_ids = ();
my %hap_map_hmc_counts = ();
my $total_num_seqs_count = 0;
while(<INFILE>){
	chomp $_;
# 	warn $_ . "\n";
	if($i eq 0){
		my @split_hmp_entries = split(/\t/, $_);
		for(my $j = 1; $j < (scalar(@split_hmp_entries) - 5); $j++){
# 			warn $split_hmp_entries[$j] . "\n";
			my @split_individual_name = split(/_/, $split_hmp_entries[$j]);
			my $individual_id = $split_individual_name[0];
			$individual_ids{$j} = $individual_id;
		}
	}elsif($i ne 0){
		my @split_hmp_entries = split(/\t/, $_);
		my $snp_id = $split_hmp_entries[0];
		for(my $j = 1; $j < (scalar(@split_hmp_entries) - 5); $j++){
			my $individual_id = $individual_ids{$j};
			my @split_sequence_counts = split(/\|/, $split_hmp_entries[$j]);
			my $query_count = $split_sequence_counts[0];
			my $hit_count = $split_sequence_counts[1];
			my $query_id = join("_", $snp_id, "query");
			my $hit_id = join("_", $snp_id, "hit");
			$hap_map_hmc_counts{$individual_id}{$query_id} = $query_count;
			$hap_map_hmc_counts{$individual_id}{$hit_id} = $hit_count;
			$total_num_seqs_count += ($query_count + $hit_count);
		}
	}
	$i++;
}
close(INFILE) or die "Couldn't close file $hap_map_hmc_infile";


# Create output directory if it doesn't already exist.
my $blastn_output_dir = join('/', $project_dir, "UNEAK_BLASTN_FILES");
unless(-d $blastn_output_dir){
	mkdir($blastn_output_dir, 0777) or die "Can't make directory: $!";
}

foreach my $individual_id (sort keys %hap_map_hmc_counts){
	foreach my $sequence_id (sort keys $hap_map_hmc_counts{$individual_id}){
# 		print join("\t", $individual_id, $sequence_id, $hap_map_hmc_counts{$individual_id}{$sequence_id}) . "\n";
		
		# number of sequences present with this allele.
		my $max_target_seqs = $hap_map_hmc_counts{$individual_id}{$sequence_id};
		my $blastn_filename = join("_", $individual_id, $sequence_id);
		warn join("\t", $individual_id, $sequence_id, $max_target_seqs) . "\n";
		if($max_target_seqs ne 0){
			generate_uneak_blastn($blastn_filename, $hap_map_fasta_seqs{$sequence_id}, $fasta_filenames{$individual_id}, $max_target_seqs, $uneak_sequence_length, $blast_num_cpu, $blastn_output_dir);
		}
		
	}
}

my ($blastn_files, $blastn_file_counter) = find_files($blastn_output_dir, "uneak_blastn.tsv");
my $blastn_outfile = join('/', $uneak_project_dir, $project_name . ".uneak_blastn.tsv");
open(OUTFILE, ">$blastn_outfile") or die "Couldn't open file $blastn_outfile for writting, $!";
print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch",
	"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
my $entries_counter = 0;
foreach my $blastn_filename (sort keys %{$blastn_files}){
	my $blastn_infile = $blastn_files->{$blastn_filename};
	open(INFILE, "<$blastn_infile") or die "Couldn't open file $blastn_infile for reading, $!";
	my $i = 0;
	while(<INFILE>){
		chomp $_;
		warn $_ . "\n";
		if($i ne 0){
			print OUTFILE $_ . "\n";
			$entries_counter++;
		}
		$i++;
	}
	close(INFILE) or die "Couldn't close file $blastn_infile";
	my ($individual_id, $query_id) = split(/_/, $blastn_filename, 2);
	my $hap_map_hmc_counter = $hap_map_hmc_counts{$individual_id}{$query_id};
	die "Error: hap_map_hmc_counter=$hap_map_hmc_counter ne hap_map_hmc_counter=$hap_map_hmc_counter" if($hap_map_hmc_counter ne $hap_map_hmc_counter);
	$entries_counter = 0;
}
close(OUTFILE) or die "Couldn't close file $blastn_outfile";

warn $total_num_seqs_count . "\n";

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
	my $file_counter = 0;	opendir(DIR1, $infile_dir) || die "Error in opening dir $infile_dir\n";
	while( my $dir_name = readdir(DIR1)){
		my $individual_dir = join('/', $infile_dir, $dir_name);
		opendir(DIR2, $individual_dir) || die "Error in opening dir $individual_dir\n";
		while( my $file_name = readdir(DIR2)){
			my $infile_name = join('/', $individual_dir, $file_name) if ($file_name =~ m/\.$suffix$/);
			warn "$infile_name\n" if ($file_name =~ m/\.$suffix$/);
			$files{$file_name} = $infile_name if ($file_name =~ m/\.$suffix$/);
			$file_counter++ if ($file_name =~ m/\.$suffix$/);
		}
		closedir(DIR2);
	}
	closedir(DIR1);
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

	my $blastn_filename = shift;
	die "Error lost blastn file name" unless defined $blastn_filename;
	
	my $fastq_adaptor_sequence = shift;
	die "Error lost fastq adaptor sequence" unless defined $fastq_adaptor_sequence;
	
	my $fasta_target = shift;
	die "Error lost fasta database target file" unless defined $fasta_target;
	
	my $max_target_seqs = shift;
	die "Error lost maximum number of sequences" unless defined $max_target_seqs;
	
	my $uneak_sequence_length = shift;
	die "Error lost uneak sequence length" unless defined $uneak_sequence_length;
	
	my $blast_num_cpu = shift;
	die "Error lost number of cpus to allocate" unless defined $blast_num_cpu;
	
	my $blastn_output_dir = shift;
	die "Error lost adaptor blastn output directory" unless defined $blastn_output_dir;
	
	my ($individual_id, $query_id) = split(/_/, $blastn_filename, 2);
	
	my $fasta_target_name = fileparse($fasta_target);
	
	makeblastdb_nuc($fasta_target);
	
	# Create output directory if it doesn't already exist.
	my $individual_blastn_output_dir = join('/', $blastn_output_dir, $individual_id);
	unless(-d $individual_blastn_output_dir){
		mkdir($individual_blastn_output_dir, 0777) or die "Can't make directory: $!";
	}

	my $uneak_blastn_tsv_outfile = join('/', $individual_blastn_output_dir, $blastn_filename . ".uneak_blastn.tsv");
	unless(-s $uneak_blastn_tsv_outfile){
		warn join(" ", "Generating blast tab-delimited file", join("", $blastn_filename, ".uneak_blastn", ".tsv"), "using uneak sequence blastn.....") . "\n";
		my $uneakSequenceBlastnCmd  = "echo -e \"$fastq_adaptor_sequence\" | $blastn -query - -db $fasta_target -task blastn -strand plus -word_size 11 -dust yes -evalue 1e-6 -max_target_seqs $max_target_seqs -perc_identity 100 -outfmt '6 qseqid salltitles qcovhsp pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
		warn $uneakSequenceBlastnCmd . "\n\n";
		open(OUTFILE, ">$uneak_blastn_tsv_outfile") or die "Couldn't open file $uneak_blastn_tsv_outfile for writting, $!";
                print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch",
                "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
		local (*UNEAK_BLASTN_OUT, *UNEAK_BLASTN_IN);
		my $pid = open2(\*UNEAK_BLASTN_OUT,\*UNEAK_BLASTN_IN, $uneakSequenceBlastnCmd) or die "Error calling open2: $!";
		close UNEAK_BLASTN_IN or die "Error closing STDIN to uneak sequence blastn process: $!";
		while(<UNEAK_BLASTN_OUT>){
			chomp $_;
			my @split_adaptor_blastn_hit =  split(/\t/, $_);
			my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
			$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @split_adaptor_blastn_hit;
			if($align_length eq $uneak_sequence_length){
				print OUTFILE join("\t", join("_", $blastn_filename, $fastq_adaptor_sequence), $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch, 
					$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
			}
		}
		close UNEAK_BLASTN_OUT or die "Error closing STDOUT from uneak sequence blastn process: $!";
		wait;
		close(OUTFILE) or die "Couldn't close file $uneak_blastn_tsv_outfile";
	}
	return $uneak_blastn_tsv_outfile;
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
