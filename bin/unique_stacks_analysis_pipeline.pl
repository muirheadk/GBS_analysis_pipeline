#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

#### PROGRAM NAME ####
# unique_stacks_analysis_pipeline.pl - Program that aligns quality filtered, demultiplexed, and adapter trimmed GBS data sequences (padded) into exactly-matching stacks. Comparing the stacks it will form a set of loci and detect SNPs at each locus using a maximum likelihood framework. Runs the Stacks unique pipeline using the ustacks, cstacks, and sstacks programs of the Stacks Software Suite.

#### DESCRIPTION ####
# This program takes the quality filtered, demultiplexed, and adapter trimmed GBS *.fastq files (padded) as input. It executes the ustacks program, which extracts sequence stacks using a denovo assembly approach to form exact matching stacks. Comparing the stacks it will form a set of loci and detect SNPs at each locus using a maximum likelihood framework. These sequence stacks are then processed using cstacks and sstacks to obtain the filtered SNP stacks output files.

#### SAMPLE COMMAND ####
# perl unique_stacks_analysis_pipeline.pl -i ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS/GBS_TRIMMED_ADAPTER_DIR/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR/CHRISTIANNE_MCDONALD/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -c 7 -o ~/workspace/GBS_data-08-10-2013/CHRISTIANNE_MCDONALD_POLYGONIA
my ($gbs_fastq_dir, $gbs_fastq_file_type, $stacks_sql_id, $min_depth_coverage_ustacks, $max_nuc_distance_ustacks, $max_align_distance_ustacks, $alpha_value_ustacks, $max_locus_stacks, $num_mismatches_tag, $num_threads, $output_dir);

GetOptions(
	'i=s'    => \$gbs_fastq_dir, # The absolute path to the quality filtered, demultiplexed, and adapter trimmed *.fastq input file (padded) directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
	't=s'    => \$gbs_fastq_file_type, # The fastq input file type. Default: gzfastq
	'b=s'    => \$stacks_sql_id, # The SQL ID to insert into the output to identify this sample. Default: 1
	'd=s'    => \$min_depth_coverage_ustacks, # The minimum depth of coverage to report a stack. Default: 2
	'm=s'    => \$max_nuc_distance_ustacks, # The maximum distance (in nucleotides) allowed between stacks. Default: 2
	'n=s'    => \$max_align_distance_ustacks, # The maximum distance allowed to align secondary reads to primary stacks. Default: ($max_nuc_distance_ustacks + 2)
	'a=s'    => \$alpha_value_ustacks, # The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
	'l=s'    => \$max_locus_stacks, # The maximum number of stacks at a single de novo locus. Default: 3
	's=s'    => \$num_mismatches_tag, # The number of mismatches allowed between sample tags when generating the catalog. Default: 1
	'c=s'    => \$num_threads, # The number of cpu cores to use for the stacks programs. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2
	'o=s'    => \$output_dir, # The absolute path to the output directory to contain the Stacks output files and directories.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
	defined $gbs_fastq_dir
	and defined $output_dir
);

# The fastq input file type. Default: gzfastq
$gbs_fastq_file_type = 'gzfastq' unless defined $gbs_fastq_file_type;

# The SQL ID to insert into the output to identify this sample.
$stacks_sql_id = 1 unless defined $stacks_sql_id;

# The minimum depth of coverage to report a stack. Default: 2
$min_depth_coverage_ustacks = 2 unless defined $min_depth_coverage_ustacks;

# The maximum distance (in nucleotides) allowed between stacks. Default: 2
$max_nuc_distance_ustacks = 2 unless defined $max_nuc_distance_ustacks;

# The maximum distance allowed to align secondary reads to primary stacks. Default: ($max_nuc_distance_ustacks + 2)
$max_align_distance_ustacks = ($max_nuc_distance_ustacks + 2) unless defined $max_align_distance_ustacks;

# The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
$alpha_value_ustacks = 0.05 unless defined $alpha_value_ustacks;

# The maximum number of stacks at a single de novo locus. Default: 3
$max_locus_stacks = 3 unless defined $max_locus_stacks;

# The number of mismatches allowed between sample tags when generating the catalog. Default: 1
$num_mismatches_tag = 1 unless defined $num_mismatches_tag;

# The number of cpu cores to use for the stacks programs. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2
$num_threads = 2 unless defined $num_threads;

# Program dependencies - The absolute paths to gunzip to uncompress bulk compressed fastq.gz input file if present and stacks related programs.
my ($gunzip, $ustacks, $cstacks, $sstacks);
$gunzip				= '/bin/gunzip';
$ustacks			= '/usr/local/bin/ustacks';
$cstacks			= '/usr/local/bin/cstacks';
$sstacks			= '/usr/local/bin/sstacks';

sub usage {

die <<"USAGE";

Usage: $0 -i gbs_fastq_dir -t gbs_fastq_file_type -p project_name -b stacks_sql_id -d min_depth_coverage_pstacks -m max_nuc_distance_ustacks -n max_align_distance_ustacks -a alpha_value_ustacks -l max_locus_stacks -s num_mismatches_tag -c num_threads -o output_dir

DESCRIPTION - This program takes the quality filtered, demultiplexed, and adapter trimmed GBS *.fastq files (padded) as input. It executes the ustacks program, which extracts sequence stacks using a denovo assembly approach to form exact matching stacks. Comparing the stacks it will form a set of loci and detect SNPs at each locus using a maximum likelihood framework. These sequence stacks are then processed using cstacks and sstacks to obtain the filtered SNP stacks output files.

OPTIONS:

-i gbs_fastq_dir - The absolute path to the quality filtered, demultiplexed, and adapter trimmed *.fastq input file (padded) directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.

-t gbs_fastq_file_type - The fastq input file type. Default: gzfastq

-b stacks_sql_id - The SQL ID to insert into the output to identify this sample.

-d min_depth_coverage_pstacks - The minimum depth of coverage to report a stack. Default: 2

-m max_nuc_distance_ustacks - Maximum distance (in nucleotides) allowed between stacks. Default: 2
    
-n max_align_distance_ustacks - The maximum distance allowed to align secondary reads to primary stacks. Default: (max_nuc_distance_ustacks + 2)
    
-a alpha_value_ustacks - The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05

-l max_locus_stacks - The maximum number of stacks at a single de novo locus. Default: 3

-s num_mismatches_tag - The number of mismatches allowed between sample tags when generating the catalog. Default: 1

-c num_threads - The number of cpu cores to use for the stacks programs. You should choose a number so that this parameter is at most the total number of cpu cores on your system minus 1. Default: 2

-o output_dir - The absolute path to the output directory to contain the Stacks output files and directories.
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create the stacks output directory if it doesn't already exist.
my $stacks_output_dir = join('/', $output_dir, "STACKS_OUTFILES");
unless(-d $stacks_output_dir){
	mkdir($stacks_output_dir, 0777) or die "Can't make directory: $!";
}

# Find all files in the specified directory with the extension *.fastq or *.fastq.qz.
my ($gbs_fastq_files, $gbs_fastq_file_count);
if($gbs_fastq_file_type eq "gzfastq"){
	($gbs_fastq_files, $gbs_fastq_file_count) = find_files($gbs_fastq_dir, "fastq.gz");
}elsif($gbs_fastq_file_type eq "fastq"){
	($gbs_fastq_files, $gbs_fastq_file_count) = find_files($gbs_fastq_dir, "fastq");
}else{
	die "Error: $gbs_fastq_file_type is not a recognized input file type. Use the -t option and specify gzfastq for *.fastq.gz files or fastq for *.fastq files.";
}

die "Error: The GBS fastq file count was $gbs_fastq_file_count. The input file type was specified as $gbs_fastq_file_type. Either you specified the wrong directory or the files in $gbs_fastq_dir are a different format than $gbs_fastq_file_type. Use the -t option and specify gzfastq for *.fastq.gz files or fastq for *fastq files." unless($gbs_fastq_file_count > 0);

# Iterate through each GBS fastq file and execute the ustacks program, which extracts exact-matching stacks and detects SNPs at each locus using a maximum likelihood framework.
my $sql_id = 1;
foreach my $file_name (sort keys %{$gbs_fastq_files}){
	warn "Processing " . $file_name . ".....\n";
	my $gbs_fastq_infile = $gbs_fastq_files->{$file_name};
    
	# If the bulk fastq file is compressed, uncompress the file and set the resulting fastq filename to be the fastq infile.
	if($gbs_fastq_file_type eq "gzfastq"){
		my $uncompressed_fastq_file = gunzip_fastq_file($gbs_fastq_infile);
		$gbs_fastq_infile = $uncompressed_fastq_file;
	}
	
	# Execute the ustacks program, which extracts exact-matching stacks and detects SNPs at each locus using a maximum likelihood framework.
	ustacks($gbs_fastq_infile, $sql_id, $min_depth_coverage_ustacks, $max_nuc_distance_ustacks, $max_align_distance_ustacks, 
		$num_threads, $alpha_value_ustacks, $max_locus_stacks, $stacks_output_dir);
	
	# unlink uncompressed file if $gbs_fastq_file_type is gzfastq.
	unlink($gbs_fastq_infile) or die "Could not unlink $gbs_fastq_infile: $!" if($gbs_fastq_file_type eq "gzfastq");
	
	$sql_id++;
}

# Execute the cstacks program to build a catalog from a set of samples processed by the ustacks program. The cstacks program creates a set of consensus loci, merging alleles together.
my $cstacks_file = cstacks($stacks_output_dir, $stacks_sql_id, $num_mismatches_tag, $num_threads);

# Find all ustacks tags output files from the stacks output directory with the extension *.tags.tsv.
my ($ustacks_tags_files, $ustacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");

# Iterate through each tags output file with extension *.tags.tsv and execute the sstacks program.
foreach my $file_name (sort keys %{$ustacks_tags_files}){
	if($file_name !~ m/batch_\d+\.catalog\.tags\.tsv/){ # If *.tags.tsv file does not match batch_*.catalog.tags.tsv.
		my $ustacks_tags_infile = $ustacks_tags_files->{$file_name};
		
		# Get the basename of the tags filename without the .tags.tsv extension.
		my $ustacks_filename = fileparse($ustacks_tags_infile, qr/\.tags.tsv/);
		warn "Processing " . $ustacks_filename . ".....\n";
		my $sstacks_infile = join('/', $stacks_output_dir, $ustacks_filename);
		
		# Execute the sstacks program. Sets of stacks constructed by the pstacks program is searched against the catalog produced by cstacks.
		sstacks($cstacks_file, $sstacks_infile, $stacks_sql_id, $num_threads, $stacks_output_dir);
	}
}

# ustacks($fastq_infile, $sql_id, $min_depth_coverage, $max_nuc_distance_ustacks, $max_align_distance_ustacks, $num_threads, $alpha_value, $max_locus_stacks, $stacks_output_dir) - Executes the ustacks program in the Stacks Software Suite. Extracts exact-matching stacks and detects SNPs at each locus using a maximum likelihood framework.
#
# Input paramater(s):
#
# $fastq_infile - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
#
# $sql_id - The SQL ID to insert into the output to identify this sample.
#
# $min_depth_coverage - The minimum depth of coverage to report a stack.
#
# $max_nuc_distance_ustacks - The maximum distance (in nucleotides) allowed between stacks.
#
# $max_align_distance_ustacks - The maximum distance (in nucleotides) allowed between stacks.
#
# $num_threads - The number of threads to use for ustacks.
#
# $alpha_value - The chi square significance level required to call a heterozygote or homozygote.
#
# $max_locus_stacks - The maximum number of stacks at a single de novo locus.
#
# $stacks_output_dir - The stacks output directory that contains the results from the pstacks program.

# ustacks -t fastq -f ./samples/f0_male.fq -o ./stacks -i 1 -d -r -m 3 -p 15
sub ustacks{

	# The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
	my $fastq_infile = shift;
	die "Error lost the GBS fastq input file" unless defined $fastq_infile;
	
	# The SQL ID to insert into the output to identify this sample.
	my $sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $sql_id;
	
	# The minimum depth of coverage to report a stack.
	my $min_depth_coverage = shift;
	die "Error lost the minimum depth of coverage to report a stack" unless defined $min_depth_coverage;
	
	# The maximum distance (in nucleotides) allowed between stacks.
	my $max_nuc_distance_ustacks = shift;
	die "Error lost the maximum distance (in nucleotides) allowed between stacks" unless defined $max_nuc_distance_ustacks;

	# The maximum distance allowed to align secondary reads to primary stacks.
	my $max_align_distance_ustacks = shift;
	die "Error lost the maximum distance allowed to align secondary reads to primary stacks" unless defined $max_align_distance_ustacks;

	# The number of threads to use for ustacks.
	my $num_threads = shift;
	die "Error lost the number of threads to use for ustacks" unless defined $num_threads;
	
	# The chi square significance level required to call a heterozygote or homozygote.
	my $alpha_value = shift;
	die "Error lost the chi square significance level required to call a heterozygote or homozygote" unless defined $alpha_value;
	
	# The maximum number of stacks at a single de novo locus.
	my $max_locus_stacks = shift;
	die "Error lost the maximum number of stacks at a single de novo locus" unless defined $max_locus_stacks; 
	
	# The stacks output directory that contains the results from the pstacks program.
	my $stacks_output_dir = shift;
	die "Error lost the stacks output file directory" unless defined $stacks_output_dir;
	
	# Get the basename of the fastq filename without the .fastq extension.
	my $fastq_filename = fileparse($fastq_infile, qr/\.fastq/);
	
	# Format the ustacks individual alleles, snps, and tags output files.
	my ($ustacks_alleles_file, $ustacks_snps_file, $ustacks_tags_file);
	$ustacks_alleles_file = join('/', $stacks_output_dir, $fastq_filename . '.alleles.tsv');
	$ustacks_snps_file = join('/', $stacks_output_dir, $fastq_filename . '.snps.tsv');
	$ustacks_tags_file = join('/', $stacks_output_dir, $fastq_filename . '.tags.tsv');
	
	# Create the stacks log output directory if it doesn't already exist.
	my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
	unless(-d $stacks_log_output_dir){
		mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
	}

	# The standard out log file for the ustacks program.
	my $ustacks_log_outfile = join('/', $stacks_log_output_dir, join("_", $fastq_filename, "ustacks.log"));
        
	# Execute the ustacks program if the pstacks alleles, snps, and tags output files are not already generated.
	unless(-s $ustacks_alleles_file and -s $ustacks_snps_file and -s $ustacks_tags_file){
		warn "Executing ustacks.....\n\n";
		my $ustacksCmd  = "$ustacks -t fastq -f $fastq_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -M $max_nuc_distance_ustacks -N $max_align_distance_ustacks -p $num_threads -d -r -R --model_type snp --alpha $alpha_value --max_locus_stacks $max_locus_stacks 2> $ustacks_log_outfile";
		warn $ustacksCmd . "\n\n";
		system($ustacksCmd) == 0 or die "Error calling $ustacksCmd: $?";
	}
}

# $cstacks_file = cstacks($stacks_output_dir, $stacks_sql_id, $num_mismatches_tag, $num_threads) - Executes the cstacks program in the Stacks Software Suite. Build a catalog from a set of samples processed by the ustacks program. Creates a set of consensus loci, merging alleles together.
# 
# Input paramater(s):
# 
# $stacks_output_dir - The stacks directory that contains the results from the ustacks program.
#
# $stacks_sql_id - The SQL ID to insert into the output to identify this sample.
#
# $num_mismatches_tag - The number of mismatches allowed between sample tags when generating the catalog.
#
# $num_threads - The number of threads to use for sstacks.
# 
# Output paramater(s):
# 
# $cstacks_file - The cstacks output file prefix for the tab-delmited catalog alleles, snps, and tags files.

# cstacks -b 1 -o /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -n 1 -p 7 \
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

	# The stacks directory that contains the results from the pstacks program.
	my $stacks_output_dir = shift;
	die "Error lost the stacks output directory" unless defined $stacks_output_dir;
	
	# The SQL ID to insert into the output to identify this sample.
	my $stacks_sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $stacks_sql_id;
	
	# The number of mismatches allowed between sample tags when generating the catalog.
	my $num_mismatches_tag = shift;
	die "Error lost the number of mismatches allowed between sample tags when generating the catalog" unless defined $num_mismatches_tag;
	
	# The number of threads to use for cstacks.
	my $num_threads = shift;
	die "Error lost the number of threads to use for cstacks" unless defined $num_threads;
	
	# The cstacks batch file prefix file path.
	my $cstacks_file = join('/', $stacks_output_dir, join("_", "batch", $stacks_sql_id));
	
	# The cstacks catalog alleles, snps, and tags file paths.
	my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file);
	$cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
	$cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
	$cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
	
	# Create the stacks log output directory if it doesn't already exist.
	my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
	unless(-d $stacks_log_output_dir){
		mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
	}

	# The standard out log file for the cstacks program.
	my $cstacks_log_outfile = join('/', $stacks_log_output_dir, "cstacks.log");
	
	# Execute the cstacks program if the following files are not already generated.
	unless(-s $cstacks_alleles_file and -s $cstacks_snps_file and -s $cstacks_tags_file){
	
		# Find all ustacks tags output files from the stacks output directory with the extension *.tags.tsv.
		my ($ustacks_tags_files, $ustacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");
		my @cstacks_soptions = ();
		# Iterate through each ustacks tags output file with extension *.tags.tsv and execute the cstacks program.
		foreach my $file_name (sort keys %{$ustacks_tags_files}){
			my $ustacks_tags_infile = $ustacks_tags_files->{$file_name};
			
			# Get the basename of the tags filename without the .tags.tsv extension.
			my $ustacks_filename = fileparse($ustacks_tags_infile, qr/\.tags.tsv/);
			
			# Obtain a list of all the file prefix paths specifed by the -s option for cstacks.
			warn "Processing " . $ustacks_filename . ".....\n";
			my $cstacks_infile = join('/', $stacks_output_dir, $ustacks_filename);
			push(@cstacks_soptions, "-s $cstacks_infile ");
		}
		
		#### USED PARAMETERS ####
		# b — MySQL ID of this batch.
		# o — output path to write results.
		# s — TSV file from which to load radtags.
		# p — enable parallel execution with num_threads threads.
		# n — number of mismatches allowed between sample tags when generating the catalog.

		#### UNUSED PARAMETERS ####
		# g — base catalog matching on genomic location, not sequence identity.
		# m — include tags in the catalog that match to more than one entry.
		# h — display this help messsage.

		# Catalog editing:
		# --catalog [path] — provide the path to an existing catalog. cstacks will add data to this existing catalog.

		# Advanced options:
		# --report_mmatches — report query loci that match more than one catalog locus.
		
		warn "Executing cstacks.....\n\n";
		my $cstacks_joined_soptions = join("\\\n", @cstacks_soptions);
		my $cstacksCmd  = "$cstacks -b $stacks_sql_id -o $stacks_output_dir -n $num_mismatches_tag -p $num_threads \\\n $cstacks_joined_soptions 2> $cstacks_log_outfile";
		warn $cstacksCmd . "\n\n";
		system($cstacksCmd) == 0 or die "Error calling $cstacksCmd: $?";
	}
	
	# Returns the cstacks batch file prefix file path.
	return $cstacks_file;
}

# sstacks($stacks_catalog_infile, $stacks_sample_infile, $stacks_sql_id, $num_threads, $stacks_output_dir) - Executes the sstacks program in the Stacks Software Suite. Sets of stacks constructed by the pstacks program is searched against the catalog produced by cstacks.
# 
# Input paramater(s):
# 
# $stacks_catalog_infile - The stacks catalog input batch_*.catalog.*.tsv file from which to load the catalog GBS-Tags.
# 
# $stacks_sample_infile - The stacks sample input *.tags.tsv file from which to load sample GBS-Tags.
#
# $stacks_sql_id - The SQL ID to insert into the output to identify this sample.
# 
# $num_threads - The number of threads to use for sstacks.
# 
# $stacks_output_dir - The stacks output file directory that contains all the results files generated by sstacks.
sub sstacks{
    
	# The stacks catalog input batch_*.catalog.*.tsv file from which to load the catalog GBS-Tags.
	my $stacks_catalog_infile = shift;
	die "Error lost the stacks catalog input file directory" unless defined $stacks_catalog_infile;
    
	# The stacks sample input *.tags.tsv file from which to load sample GBS-Tags.
	my $stacks_sample_infile = shift;
	die "Error lost the stacks sample input file directory" unless defined $stacks_sample_infile;
	
	# The SQL ID to insert into the output to identify this sample.
	my $stacks_sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $stacks_sql_id;
	
	# The number of threads to use for sstacks.
	my $num_threads = shift;
	die "Error lost the number of cores for stacks" unless defined $num_threads;
	
	# The stacks output file directory that contains all the results files generated by sstacks.my $stacks_output_dir = shift;
	die "Error lost the stacks output file directory" unless defined $stacks_output_dir;
	
	# Get the basename of the sam filename without the .sam extension.
	my $stacks_sample_filename = fileparse($stacks_sample_infile, qr//);
    
	# The sstacks output matches file.
	my $sstacks_matches_file = join('/', $stacks_output_dir, $stacks_sample_filename . '.matches.tsv');
	
	# Create the stacks log output directory if it doesn't already exist.
	my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
	unless(-d $stacks_log_output_dir){
		mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
	}

	# The standard out log file for the cstacks program.
	my $sstacks_log_outfile = join('/', $stacks_log_output_dir, join("_", $stacks_sample_filename, "sstacks.log"));
    
	#### USED PARAMETERS ####
	# b — MySQL ID of this batch.
	# c — TSV file from which to load the catalog RAD-Tags.
	# s — TSV file from which to load sample RAD-Tags.
	# o — output path to write results.
	# p — enable parallel execution with num_threads threads.

	#### UNUSED PARAMETERS ####
	# g — base matching on genomic location, not sequence identity.
	# r — Load the TSV file of a single sample instead of a catalog.
	# x — don’t verify haplotype of matching locus.
	# v — print program version.
	# h — display this help messsage.

	# Execute the sstacks program if the matches sstacks results files are not already generated.
	unless(-s $sstacks_matches_file){
		warn "Executing sstacks.....\n\n";
		my $sstacksCmd  = "$sstacks -b $stacks_sql_id -c $stacks_catalog_infile -s $stacks_sample_infile -o $stacks_output_dir -p $num_threads 2> $sstacks_log_outfile";
		warn $sstacksCmd . "\n\n";
		system($sstacksCmd) == 0 or die "Error calling $sstacksCmd: $?";
	}
}

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
sub find_files{
    
	# The input file directory.
	my $infile_dir = shift;
	die "Error lost input file directory" unless defined $infile_dir;
	
	# The file extension suffix.
	my $suffix = shift;
	die "Error lost file extension suffix directory" unless defined $suffix;
	
	if(-d $infile_dir){ # Check if $infile_dir is a directory.
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
	}else{
		die "Error $infile_dir does not exist!\n";
	}
}

# $output_dir = gunzip_fastq_file($fastq_file) - Execute the gunzip program to uncompress the compressed fastq file.
#
# Input paramater(s):
#
# $fastq_file - The fastq file to uncompress using gunzip.
#
# Output paramater(s):
#
# $uncompressed_fastq_file - The uncompressed fastq file path.
sub gunzip_fastq_file{
	
	# The fastq file to uncompress using gunzip.
	my $fastq_file = shift;
	die "Error lost the fastq file to uncompress using gunzip" unless defined $fastq_file;
	
	my ($fastq_filename, $fastq_dir) = fileparse($fastq_file, qr/\.fastq.gz/);

	# Split the fastq filename so that we can get the individual id.
	my @split_fastq_filename = split(/_/, $fastq_filename);
	my $individual_id = $split_fastq_filename[0];
	
	my $uncompressed_fastq_file = join('/', $fastq_dir, $individual_id . ".fastq");
	
	unless(-s $uncompressed_fastq_file){
		warn "Calling gunzip for $fastq_file....\n";
		my $gunzipCmd  = "$gunzip -c $fastq_file > $uncompressed_fastq_file";
		warn $gunzipCmd . "\n\n";
		system($gunzipCmd) == 0 or die "Error calling $gunzipCmd: $?";
	}
	return $uncompressed_fastq_file;
}
