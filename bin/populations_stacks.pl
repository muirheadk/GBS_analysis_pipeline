#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

# perl populations_stacks.pl -i ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR/STEPHEN_TREVOY/TRIMMED_OUTPUT_FILES/STEPHEN_TREVOY_trimmed_offset_3.fastq.gz -p MPB_MALE_GBS -g ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_sequence_data/DendPond_male_1.0/Primary_Assembly/unplaced_scaffolds/FASTA/DendPond_male_1.0_unplaced.scaf.fa -c 7 -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3
# populations -P ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -M ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/mpb_poputlations.txt -b 1 -k -p 12 -r 0.75 -m 5 -a 0.1 -f p_value --p_value_cutoff 0.05 -t 7 --fasta

my ($gbs_fastq_dir, $gbs_fastq_file_type, $project_name, $refgen_infile, $gbs_sequence_length, $min_depth_coverage_pstacks, $alpha_value_pstacks, $bwa_num_cpu, $stacks_num_cpu, $output_dir);

GetOptions(
	'i=s'    => \$gbs_fastq_dir,
	't=s'    => \$gbs_fastq_file_type,
	'p=s'    => \$project_name,
	'g=s'    => \$refgen_infile,
	'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
	'd=s'    => \$min_depth_coverage_pstacks,
	'a=s'    => \$alpha_value_pstacks,
	'c=s'    => \$bwa_num_cpu,
	's=s'    => \$stacks_num_cpu,
	'o=s'    => \$output_dir,
);

usage() unless (
	defined $gbs_fastq_dir
	and defined $project_name
	and defined $refgen_infile
	and defined $output_dir
);

$gbs_fastq_file_type = 'gzfastq' unless defined $gbs_fastq_file_type;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

$min_depth_coverage_pstacks = 1 unless defined $min_depth_coverage_pstacks;

$alpha_value_pstacks = 0.05 unless defined $alpha_value_pstacks;

$bwa_num_cpu = 2 unless defined $bwa_num_cpu;

$stacks_num_cpu = 2 unless defined $stacks_num_cpu;

# Program dependencies - The stacks populations program.
my $populations			= '/usr/local/bin/populations';

sub usage {

die <<"USAGE";


Usage: $0 -i gbs_fastq_dir -p project_name -g refgen_infile -c bwa_num_cpu -o output_dir

DESCRIPTION - 

OPTIONS:

-i gbs_fastq_dir -

-p project_name - 

-g refgen_infile -  

-c bwa_num_cpu - 

-o output_dir - 

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create the reference genome output directory if it doesn't already exist.
my $refgen_output_dir = join('/', $output_dir, "REFERENCE_GENOME");
unless(-d $refgen_output_dir){
      mkdir($refgen_output_dir, 0777) or die "Can't make directory: $!";
}

# Converts the reference genome to BWA input format and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers.
my ($refgen_fasta_outfile, $refgen_toc_outfile, $num_chromosomes) = convert_refgen_bwa_input_format($refgen_infile, $project_name, $refgen_output_dir);

# Creates the BWA reference genome index file
bwa_index($refgen_fasta_outfile, $project_name, $refgen_output_dir);

# Create the BWA alignment file output directory if it doesn't already exist.
my $bwa_output_dir = join('/', $output_dir, "BWA_ALIGNMENT_FILES");
unless(-d $bwa_output_dir){
      mkdir($bwa_output_dir, 0777) or die "Can't make directory: $!";
}

# Find all files in the specified directory with the extension *.fastq or *.fastq.qz.
my ($gbs_fastq_files, $gbs_fastq_file_count);
if($gbs_fastq_file_type eq "gzfastq"){
	($gbs_fastq_files, $gbs_fastq_file_count) = find_files($gbs_fastq_dir, "fastq.gz");
}elsif($gbs_fastq_file_type eq "fastq"){
	($gbs_fastq_files, $gbs_fastq_file_count) = find_files($gbs_fastq_dir, "fastq");
}

foreach my $file_name (sort keys %{$gbs_fastq_files}){
	warn "Processing " . $file_name . ".....\n";
	my $gbs_fastq_infile = $gbs_fastq_files->{$file_name};
    
	# If the bulk fastq file is compressed, uncompress the file and set the resulting fastq filename to be the fastq infile.
	if($gbs_fastq_file_type eq "gzfastq"){
		my $uncompressed_fastq_file = gunzip_fastq_file($gbs_fastq_infile);
		$gbs_fastq_infile = $uncompressed_fastq_file;
	}

	# Creates the BWA alignment file using the GBS fastq input file to align the fastq sequence reads to the reference genome.
	my $bwa_alignment_outfile = bwa_aln($refgen_fasta_outfile, $gbs_fastq_infile, $project_name, $bwa_num_cpu, $bwa_output_dir);

	# Creates the BWA single-ended alignment file in sam format using the GBS fastq input file to align the fastq sequence reads to the reference genome and the BWA aln format file.
	my $bwa_aligned_master_outfile = bwa_samse($refgen_fasta_outfile, $gbs_fastq_infile, $project_name, $bwa_alignment_outfile, $bwa_output_dir);
}


# Create the padded sam file output directory if it doesn't already exist.
my $padded_sam_output_dir = join('/', $output_dir, "PADDED_SAM_OUTFILES");
unless(-d $padded_sam_output_dir){
	mkdir($padded_sam_output_dir, 0777) or die "Can't make directory: $!";
}
	
my ($bwa_align_files, $bwa_align_file_count) = find_files($bwa_output_dir, "sam");
foreach my $file_name (sort keys %{$bwa_align_files}){
	warn "Processing " . $file_name . ".....\n";
	my $bwa_align_infile = $bwa_align_files->{$file_name};
	bwa_pad_sam_files($bwa_align_infile, $padded_sam_output_dir);
}

# Create the stacks output directory if it doesn't already exist.
my $stacks_output_dir = join('/', $output_dir, "STACKS_OUTFILES");
unless(-d $stacks_output_dir){
	mkdir($stacks_output_dir, 0777) or die "Can't make directory: $!";
}

my ($padded_sam_files, $padded_sam_file_count) = find_files($padded_sam_output_dir, "sam");
my $sql_id = 1;
foreach my $file_name (sort keys %{$padded_sam_files}){
	warn "Processing " . $file_name . ".....\n";
	my $padded_sam_infile = $padded_sam_files->{$file_name};
	
	pstacks($padded_sam_infile, $sql_id, $min_depth_coverage_pstacks, $stacks_num_cpu, $alpha_value_pstacks, $stacks_output_dir);
	
	$sql_id++;
}

my $cstacks_file = cstacks($stacks_output_dir, $stacks_num_cpu);

my ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");
foreach my $file_name (sort keys %{$pstacks_tags_files}){
	if($file_name !~ m/batch_\d+\.catalog\.tags\.tsv/){
		my $pstacks_tags_infile = $pstacks_tags_files->{$file_name};
		
		# Get the basename of the tags filename without the .tags.tsv extension.
		my $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv/);
		
		warn "Processing " . $pstacks_filename . ".....\n";
		my $sstacks_infile = join('/', $stacks_output_dir, $pstacks_filename);
		
		sstacks($cstacks_file, $sstacks_infile, $stacks_num_cpu, $stacks_output_dir);
	}
}

sub convert_refgen_bwa_input_format{

	my $refgen_fasta_infile = shift;
	die "Error lost the reference genome fasta input file" unless defined $refgen_fasta_infile;
	
	my $project_name = shift;
	die "Error lost the project name" unless defined $project_name;
	
	my $refgen_output_dir = shift;
	die "Error lost reference genome output directory" unless defined $refgen_output_dir;

	my $refgen_fasta_filename = fileparse($refgen_fasta_infile);
	my $refgen_fasta_outfile = join('/', $refgen_output_dir, join("_", $project_name, $refgen_fasta_filename . ".fasta"));
	my $refgen_toc_outfile = join('/', $refgen_output_dir, join("_", $project_name, $refgen_fasta_filename . ".toc.txt"));
	
	unless(-s $refgen_fasta_outfile and -s $refgen_toc_outfile){
		warn "Converting $refgen_fasta_filename to BWA input format....\n\n";
		my $seqio = Bio::SeqIO->new(-file => $refgen_fasta_infile, '-format' => 'Fasta');
		my $seq_counter = 1;
		my %fasta_seqs = ();
		while(my $seq_entry = $seqio->next_seq) {

			my $seq_id = $seq_entry->id;
			my $sequence = $seq_entry->seq;
			my $seq_desc = $seq_entry->desc;
		
			my $fasta_header = join(" ", $seq_id, $seq_desc);
			$fasta_seqs{$seq_counter} = join("\t", $fasta_header, $sequence);
			$seq_counter++;
		
		}
 		my $num_digits = length($seq_counter - 1);
		open(REFGEN_FASTA_OUTFILE, ">$refgen_fasta_outfile") or die "Couldn't open file $refgen_fasta_outfile for writting, $!";
		open(REFGEN_TOC_OUTFILE, ">$refgen_toc_outfile") or die "Couldn't open file $refgen_toc_outfile for writting, $!"; 
		print REFGEN_TOC_OUTFILE join("\t", "bwa_input_fasta_header", "original_fasta_header") . "\n";
		foreach my $seq_num (sort {$a <=> $b} keys %fasta_seqs){
 			my ($fasta_header, $sequence) = split(/\t/, $fasta_seqs{$seq_num});
 			my $padded_zeros_length = ($num_digits - length($seq_num));
 			my $padded_seq_num = '0' x $padded_zeros_length . $seq_num;
			
			warn join("\t", $padded_seq_num, $fasta_header) . "\n";
			print REFGEN_FASTA_OUTFILE join("\n", ">$padded_seq_num", $sequence) . "\n";
			print REFGEN_TOC_OUTFILE join("\t", $padded_seq_num, $fasta_header) . "\n";

		}
		close(REFGEN_FASTA_OUTFILE) or die "Couldn't close file $refgen_fasta_outfile";
		close(REFGEN_TOC_OUTFILE) or die "Couldn't close file $refgen_toc_outfile";

		# Clean out the sequence I/O object.
		$seqio = ();
	}
	
	# Need to get the actual number of chromosomes or scaffolds for the -eC ($num_chromosomes option).
	open(REFGEN_TOC_INFILE, "<$refgen_toc_outfile") or die "Couldn't open file $refgen_toc_outfile for reading, $!";
	my $i = 0;
	my $num_chromosomes = 0;
	while(<REFGEN_TOC_INFILE>){
		chomp $_;
		warn $_ . "\n";
		if($i ne 0){
			$num_chromosomes++;
		}
		$i++;
	}
	close(REFGEN_TOC_INFILE) or die "Couldn't close file $refgen_toc_outfile";
	
	return ($refgen_fasta_outfile, $refgen_toc_outfile, $num_chromosomes);
}

#bwa index -a bwtsw ./refgen/renumbered_Msequence.fasta
sub bwa_index{

	my $renum_refgen_fasta_infile = shift;
	die "Error lost the renumbered reference genome input file" unless defined $renum_refgen_fasta_infile;

	my $project_name = shift;
	die "Error lost the project name" unless defined $project_name;
        
	my $refgen_output_dir = shift;
	die "Error lost reference genome output directory" unless defined $refgen_output_dir;

	# format the index file into .amb .ann .bwt .pac .sa files.
	my ($renum_refgen_fastadbAMB, $renum_refgen_fastadbANN, $renum_refgen_fastadbBWT, $renum_refgen_fastadbPAC, $renum_refgen_fastadbSA);
	$renum_refgen_fastadbAMB = $renum_refgen_fasta_infile . '.amb';
	$renum_refgen_fastadbANN = $renum_refgen_fasta_infile . '.ann';
	$renum_refgen_fastadbBWT = $renum_refgen_fasta_infile . '.bwt';
	$renum_refgen_fastadbPAC = $renum_refgen_fasta_infile . '.pac';
	$renum_refgen_fastadbSA = $renum_refgen_fasta_infile . '.sa';
	unless(-s $renum_refgen_fastadbAMB and -s $renum_refgen_fastadbANN and -s $renum_refgen_fastadbBWT 
		and -s $renum_refgen_fastadbPAC and -s $renum_refgen_fastadbSA){
		warn "Generating the bwa index file.....\n\n";
		my $bwaIndexCmd  = "$bwa index -a bwtsw $renum_refgen_fasta_infile";
		warn $bwaIndexCmd . "\n\n"; 
		system($bwaIndexCmd) == 0 or die "Error calling $bwaIndexCmd: $?";
	}
}

#bwa aln -t 4 ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/mpbGBSTags.cnt.fq > ./mergedTagCounts/AlignedGBSTags1.sai
sub bwa_aln{

	my $renum_refgen_fasta_infile = shift;
	die "Error lost the renumbered reference genome input file" unless defined $renum_refgen_fasta_infile;

	my $gbs_fastq_infile = shift;
	die "Error lost the GBS fastq input file" unless defined $gbs_fastq_infile;

	my $project_name = shift;
	die "Error lost the project name" unless defined $project_name;

	my $bwa_num_cpu = shift;
	die "Error lost the number of cores for BWA" unless defined $bwa_num_cpu;

	my $bwa_output_dir = shift;
	die "Error lost the bwa output directory" unless defined $bwa_output_dir;

	# Get the basename of the fastq filename without the .fastq extension.
	my $gbs_fastq_filename = fileparse($gbs_fastq_infile, qr/\.fastq/);
	
	# Split the GBS fastq filename so that we can get the individual id.
    	my @split_gbs_fastq_filename = split(/_/, $gbs_fastq_filename);
    	my $individual_id = $split_gbs_fastq_filename[0];
    	
	my $bwa_alignment_outfile = join('/', $bwa_output_dir, join("_", $individual_id, $project_name . ".sai"));

	unless(-s $bwa_alignment_outfile){
		warn "Generating the bwa alignment file.....\n\n";
		my $bwaAlnCmd  = "$bwa aln -t $bwa_num_cpu $renum_refgen_fasta_infile $gbs_fastq_infile > $bwa_alignment_outfile";
		warn $bwaAlnCmd . "\n\n";
		system($bwaAlnCmd) == 0 or die "Error calling $bwaAlnCmd: $?";
	}
	return $bwa_alignment_outfile;
}

#bwa samse ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/AlignedGBSTags1.sai ./mergedTagCounts/mpbGBSTags.cnt.fq > mergedTagCounts/AlignedMasterTagsMPB.sam
sub bwa_samse{

	my $renum_refgen_fasta_infile = shift;
	die "Error lost the renumbered reference genome input file" unless defined $renum_refgen_fasta_infile;

	my $gbs_fastq_infile = shift;
	die "Error lost the GBS fastq input file" unless defined $gbs_fastq_infile;

	my $project_name = shift;
	die "Error lost the project name" unless defined $project_name;

	my $bwa_alignment_infile = shift;
	die "Error lost the BWA SAI formatted alignment (.sai) file" unless defined $bwa_alignment_infile;

	my $bwa_output_dir = shift;
	die "Error lost the bwa output directory" unless defined $bwa_output_dir;

	# Get the basename of the fastq filename without the .fastq extension.
	my $gbs_fastq_filename = fileparse($gbs_fastq_infile, qr/\.fastq/);
	
	# Split the GBS fastq filename so that we can get the individual id.
    	my @split_gbs_fastq_filename = split(/_/, $gbs_fastq_filename);
    	my $individual_id = $split_gbs_fastq_filename[0];

	my $bwa_aligned_master_outfile = join('/', $bwa_output_dir, join("_", $individual_id, $project_name . ".sam"));
	unless(-s $bwa_aligned_master_outfile){
		warn "Generating the bwa sam file.....\n\n";
		my $bwaSamseCmd  = "$bwa samse $renum_refgen_fasta_infile $bwa_alignment_infile $gbs_fastq_infile > $bwa_aligned_master_outfile";
		warn $bwaSamseCmd . "\n\n";
		system($bwaSamseCmd) == 0 or die "Error calling $bwaSamseCmd: $?";
	}
	return $bwa_aligned_master_outfile;
}

sub bwa_pad_sam_files{

	my $sam_infile = shift;
	die "Error lost the sam input file" unless defined $sam_infile;

	my $padded_sam_output_dir = shift;
	die "Error lost the padded sam output file directory" unless defined $padded_sam_output_dir;
	
	# Get the basename of the sam filename without the .sam extension.
	my $sam_filename = fileparse($sam_infile, qr/\.sam/);
	
	my $padded_sam_outfile = join('/', $padded_sam_output_dir, $sam_filename . ".sam");
	unless(-s $padded_sam_outfile){
		open(PADDED_SAM_OUTFILE, ">$padded_sam_outfile") or die "Couldn't open file $padded_sam_outfile for writting, $!";
		open(SAM_INFILE, "<$sam_infile") or die "Couldn't open file $sam_infile for reading, $!";
		while(<SAM_INFILE>){
			chomp $_;
	#  		warn $_ . "\n";
			if($_ =~ m/^\@HD|^\@SQ|^\@RG|^\@PG|^\@CO/){ # Parse header lines starting with @HD, @SQ, @RG, @PG, or @CO.
				print PADDED_SAM_OUTFILE $_ . "\n";
			}else{ #8_1308_8038_19954_1LL-06_TTCTG_5_A3length=66    16      7376    784845  37      66M     *       0       0       AATAATGCGCCAGCCAAGAGCTGTTTGGTAGCATTCTGCTGACGCTCGCACTCCCGTACACCTGCA      DDCB@>>>AFEFFFHC@GHGGGEGEGGIGCDEHBGGCIGIGGIIHGHEFCE:F<GIGIHFHHDFDD       XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:66

				my @split_sam_entry = split(/\t/, $_);
				my ($fastq_header, $bam_bitwise_flag, $fastq_sequence, $fastq_quality_scores) = ($split_sam_entry[0], $split_sam_entry[1], $split_sam_entry[9], $split_sam_entry[10]);
				
				my $fastq_sequence_length = length($fastq_sequence);
				my $fastq_quality_scores_length = length($fastq_quality_scores);
# 				die join("\t", $fastq_sequence_length, $fastq_header, $bam_bitwise_flag, $fastq_sequence, $fastq_quality_scores);
				if(($fastq_sequence_length ne $gbs_sequence_length) and ($fastq_quality_scores_length ne $gbs_sequence_length)){
					my $padded_N_length =  ($gbs_sequence_length - $fastq_sequence_length);
					my $padded_N_seq = 'N' x $padded_N_length;
					
					my $padded_fastq_sequence = "";
					# Pad sequence if aligned to sense strand or unmatched sequence.
					$padded_fastq_sequence = join("", $fastq_sequence, $padded_N_seq) if(($bam_bitwise_flag eq 0) or ($bam_bitwise_flag eq 4));
					
					# Pad sequence if aligned to antisense strand.
					$padded_fastq_sequence = join("", $padded_N_seq, $fastq_sequence) if(($bam_bitwise_flag eq 16) or ($bam_bitwise_flag eq 20));
					
					my $padded_score_length =  ($gbs_sequence_length - $fastq_quality_scores_length);
					my $padded_score_seq = '#' x $padded_score_length;
					
					my $padded_fastq_quality_scores = "";
					# Pad quality scores if aligned to sense strand or unmatched sequence in sense orientation.
					$padded_fastq_quality_scores = join("", $fastq_quality_scores, $padded_score_seq) if(($bam_bitwise_flag eq 0) or ($bam_bitwise_flag eq 4));
					
					# Pad quality scores if aligned to antisense strand or unmatched sequence in antisense orientation.
					$padded_fastq_quality_scores = join("", $padded_score_seq, $fastq_quality_scores) if(($bam_bitwise_flag eq 16) or ($bam_bitwise_flag eq 20));
					
					my $padded_fastq_sequence_length = length($padded_fastq_sequence);
					my $padded_fastq_quality_scores_length = length($padded_fastq_quality_scores);
					
					die "Error: $fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($padded_fastq_sequence_length ne $gbs_sequence_length);
					die "Error: $fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length ne padded_fastq_quality_scores_length=$padded_fastq_quality_scores_length" if($padded_fastq_sequence_length ne $padded_fastq_quality_scores_length);

					$split_sam_entry[9] = $padded_fastq_sequence;
					$split_sam_entry[10] = $padded_fastq_quality_scores;
					
					print PADDED_SAM_OUTFILE join("\t", @split_sam_entry) . "\n";
				}elsif(($fastq_sequence_length eq $gbs_sequence_length) and ($fastq_quality_scores_length eq $gbs_sequence_length)){
				
					print PADDED_SAM_OUTFILE $_ . "\n";
				}
	# 			die $_;
			}
		}
		close(SAM_INFILE) or die "Couldn't close file $sam_infile";
		close(PADDED_SAM_OUTFILE) or die "Couldn't close file $padded_sam_outfile";
	}

}

# pstacks -t sam -f ../PADDED_SAM_OUTFILES/LL-06_MPB-MALE-GBS.sam -o .. -i 1 -m 1 -p 7 --model_type snp --alpha 0.05
sub pstacks{

	my $sam_infile = shift;
	die "Error lost the sam input file" unless defined $sam_infile;
	
	my $sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $sql_id;
	
	my $min_depth_coverage = shift;
	die "Error lost the minimum depth of coverage to report a stack" unless defined $min_depth_coverage;
	
	my $stacks_num_cpu = shift;
	die "Error lost the number of cores for stacks" unless defined $stacks_num_cpu;
	
	my $alpha_value = shift;
	die "Error lost the chi square significance level required to call a heterozygote or homozygote" unless defined $alpha_value; 
	
	my $stacks_output_dir = shift;
	die "Error lost the stacks output file directory" unless defined $stacks_output_dir;
	
	# Get the basename of the sam filename without the .sam extension.
	my $sam_filename = fileparse($sam_infile, qr/\.sam/);
	
	my ($pstacks_alleles_file, $pstacks_snps_file, $pstacks_tags_file);
	$pstacks_alleles_file = join('/', $stacks_output_dir, $sam_filename . '.alleles.tsv');
	$pstacks_snps_file = join('/', $stacks_output_dir, $sam_filename . '.snps.tsv');
	$pstacks_tags_file = join('/', $stacks_output_dir, $sam_filename . '.tags.tsv');
	
	unless(-s $pstacks_alleles_file and -s $pstacks_snps_file and -s $pstacks_tags_file){
		warn "Executing pstacks.....\n\n";
		my $pstacksCmd  = "$pstacks -t sam -f $sam_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -p $stacks_num_cpu --model_type snp --alpha $alpha_value";
		warn $pstacksCmd . "\n\n";
		system($pstacksCmd) == 0 or die "Error calling $pstacksCmd: $?";
	}
}

# cstacks -b 1 -o /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -g -p 7 \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/LL-06_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M004-13-01-1D-G14-DA01_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M004-13-01-2D-G05-DA01_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M004-13-01-2F-G24-DA01_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M004-13-01-5ST-G36-DA02_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M024-07-01-02-DA07_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M033-06-01-05-DA126_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M033-07-01-09-DA36_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/M033-11-02-05-DA41_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/NG-38_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/QU-48_MPB-MALE-GBS \
# -s /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES/RR-45_MPB-MALE-GBS
sub cstacks{

	my $stacks_input_dir = shift;
	die "Error lost the padded sam output file directory" unless defined $stacks_input_dir;
	
	my $stacks_num_cpu = shift;
	die "Error lost the number of cores for stacks" unless defined $stacks_num_cpu;
	
	my $cstacks_file = join('/', $stacks_output_dir, "batch_1");
	
	my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file);
	$cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
	$cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
	$cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
	
	unless(-s $cstacks_alleles_file and -s $cstacks_snps_file and -s $cstacks_tags_file){
	
		my ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_input_dir, "tags.tsv");
		my @cstacks_soptions = ();
		foreach my $file_name (sort keys %{$pstacks_tags_files}){
			my $pstacks_tags_infile = $pstacks_tags_files->{$file_name};
			
			# Get the basename of the tags filename without the .tags.tsv extension.
			my $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv/);
			
			warn "Processing " . $pstacks_filename . ".....\n";
			my $cstacks_infile = join('/', $stacks_input_dir, $pstacks_filename);
			push(@cstacks_soptions, "-s $cstacks_infile ");
		}
		
		warn "Executing cstacks.....\n\n";
		my $cstacks_joined_soptions = join("\\\n", @cstacks_soptions);
		my $cstacksCmd  = "$cstacks -b 1 -o $stacks_output_dir -g -p $stacks_num_cpu \\\n $cstacks_joined_soptions";
		warn $cstacksCmd . "\n\n";
		system($cstacksCmd) == 0 or die "Error calling $cstacksCmd: $?";
	}
	
	return $cstacks_file;
}

sub sstacks{

	my $stacks_catalog_infile = shift;
	die "Error lost the stacks catalog input file directory" unless defined $stacks_catalog_infile;
	
	my $stacks_sample_infile = shift;
	die "Error lost the stacks sample input file directory" unless defined $stacks_sample_infile;
	
	my $stacks_num_cpu = shift;
	die "Error lost the number of cores for stacks" unless defined $stacks_num_cpu;
	
	my $stacks_output_dir = shift;
	die "Error lost the stacks output file directory" unless defined $stacks_output_dir;
	
 	# Get the basename of the sam filename without the .sam extension.
 	my $stacks_sample_filename = fileparse($stacks_sample_infile, qr//);

 	my $sstacks_matches_file = join('/', $stacks_output_dir, $stacks_sample_filename . '.matches.tsv');
	
	unless(-s $sstacks_matches_file){
		warn "Executing sstacks.....\n\n";
		my $sstacksCmd  = "$sstacks -b 1 -c $stacks_catalog_infile -s $stacks_sample_infile -o $stacks_output_dir -p $stacks_num_cpu";
		warn $sstacksCmd . "\n\n";
		system($sstacksCmd) == 0 or die "Error calling $sstacksCmd: $?";
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

# execute the gunzip program to uncompress the compressed fastq file.
sub gunzip_fastq_file{
	
	my $fastq_file = shift;
	die "Error lost the fastq file to compress using gunzip" unless defined $fastq_file;
	
	my ($fastq_filename, $fastq_dir) = fileparse($fastq_file, ".gz");
	
	my $uncompressed_fastq_file = join('/', $fastq_dir, $fastq_filename);
	unless(-s $uncompressed_fastq_file){
		warn "Calling gunzip for $fastq_file....\n";
		my $gunzipCmd  = "$gunzip -c $fastq_file > $uncompressed_fastq_file";
		warn $gunzipCmd . "\n\n";
		system($gunzipCmd) == 0 or die "Error calling $gunzipCmd: $?";
	}
	return $uncompressed_fastq_file;
}
