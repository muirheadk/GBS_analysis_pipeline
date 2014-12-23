#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
use File::Copy;
use Switch;

# perl unique_stacks_analysis_pipeline.pl -i ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR_UNPADDED/CHRISTIANNE_MCDONALD/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -p CHRISTIANNE_MCDONALD -o ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS/STACKS_DENOVO_ASSEMBLY
my ($gbs_fastq_dir, $gbs_fastq_file_type, $project_name, $min_depth_coverage_ustacks, $max_nuc_distance_ustacks, $max_align_distance_ustacks, $alpha_value_ustacks, $max_locus_stacks, $num_cpu_cores, $output_dir);

GetOptions(
	'i=s'    => \$gbs_fastq_dir,
	't=s'    => \$gbs_fastq_file_type,
	'p=s'    => \$project_name,
	'd=s'    => \$min_depth_coverage_ustacks,
	'm=s'    => \$max_nuc_distance_ustacks,
	'n=s'    => \$max_align_distance_ustacks,
	'a=s'    => \$alpha_value_ustacks,
	'l=s'    => \$max_locus_stacks,
	'c=s'    => \$num_cpu_cores,
	'o=s'    => \$output_dir,
);

usage() unless (
	defined $gbs_fastq_dir
	and defined $project_name
	and defined $output_dir
);

$gbs_fastq_file_type = 'gzfastq' unless defined $gbs_fastq_file_type;

$min_depth_coverage_ustacks = 2 unless defined $min_depth_coverage_ustacks;

$max_nuc_distance_ustacks = 2 unless defined $max_nuc_distance_ustacks;

$max_align_distance_ustacks = ($max_nuc_distance_ustacks + 2) unless defined $max_align_distance_ustacks;

$alpha_value_ustacks = 0.05 unless defined $alpha_value_ustacks;

$max_locus_stacks = 3 unless defined $max_locus_stacks;

$num_cpu_cores = 2 unless defined $num_cpu_cores;

# Program dependencies - The absolute paths to gunzip to uncompress bulk compressed fastq.gz input file if present and stacks related programs.
my ($gunzip, $ustacks, $cstacks, $sstacks);
$gunzip				= '/bin/gunzip';
$ustacks			= '/usr/local/bin/ustacks';
$cstacks			= '/usr/local/bin/cstacks';
$sstacks			= '/usr/local/bin/sstacks';

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
}

my $sql_id = 1;
foreach my $file_name (sort keys %{$gbs_fastq_files}){
	warn "Processing " . $file_name . ".....\n";
	my $gbs_fastq_infile = $gbs_fastq_files->{$file_name};
    
	# If the bulk fastq file is compressed, uncompress the file and set the resulting fastq filename to be the fastq infile.
	if($gbs_fastq_file_type eq "gzfastq"){
		my $uncompressed_fastq_file = gunzip_fastq_file($gbs_fastq_infile);
		$gbs_fastq_infile = $uncompressed_fastq_file;
	}
	
	ustacks($gbs_fastq_infile, $sql_id, $min_depth_coverage_ustacks, $max_nuc_distance_ustacks, $max_align_distance_ustacks, 
		$num_cpu_cores, $alpha_value_ustacks, $max_locus_stacks, $stacks_output_dir);
	
	$sql_id++;
}

my $cstacks_file = cstacks($stacks_output_dir, $num_cpu_cores);

my ($ustacks_tags_files, $ustacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");
foreach my $file_name (sort keys %{$ustacks_tags_files}){
	if($file_name !~ m/batch_\d+\.catalog\.tags\.tsv/){
		my $ustacks_tags_infile = $ustacks_tags_files->{$file_name};
		
		# Get the basename of the tags filename without the .tags.tsv extension.
		my $ustacks_filename = fileparse($ustacks_tags_infile, qr/\.tags.tsv/);
		
		warn "Processing " . $ustacks_filename . ".....\n";
		my $sstacks_infile = join('/', $stacks_output_dir, $ustacks_filename);
		
		sstacks($cstacks_file, $sstacks_infile, $num_cpu_cores, $stacks_output_dir);
	}
}

# ustacks -t fastq -f ./samples/f0_male.fq -o ./stacks -i 1 -d -r -m 3 -p 15
sub ustacks{

	my $fastq_infile = shift;
	die "Error lost the input file" unless defined $fastq_infile;
	
	my $sql_id = shift;
	die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $sql_id;
	
	my $min_depth_coverage = shift;
	die "Error lost the minimum depth of coverage to report a stack" unless defined $min_depth_coverage;
	
	my $max_nuc_distance_ustacks = shift;
	die "Error lost the maximum distance (in nucleotides) allowed between stacks" unless defined $max_nuc_distance_ustacks;

	my $max_align_distance_ustacks = shift;
	die "Error lost the maximum distance allowed to align secondary reads to primary stacks" unless defined $max_align_distance_ustacks;

	my $num_cpu_cores = shift;
	die "Error lost the number of cores for stacks" unless defined $num_cpu_cores;
	
	my $alpha_value = shift;
	die "Error lost the chi square significance level required to call a heterozygote or homozygote" unless defined $alpha_value;
	
	my $max_locus_stacks = shift;
	die "Error lost the maximum number of stacks at a single de novo locus" unless defined $max_locus_stacks; 
	
	my $stacks_output_dir = shift;
	die "Error lost the stacks output file directory" unless defined $stacks_output_dir;
	
	# Get the basename of the fastq filename without the .fastq extension.
	my $fastq_filename = fileparse($fastq_infile, qr/\.fastq/);
	
	my ($ustacks_alleles_file, $ustacks_snps_file, $ustacks_tags_file);
	$ustacks_alleles_file = join('/', $stacks_output_dir, $fastq_filename . '.alleles.tsv');
	$ustacks_snps_file = join('/', $stacks_output_dir, $fastq_filename . '.snps.tsv');
	$ustacks_tags_file = join('/', $stacks_output_dir, $fastq_filename . '.tags.tsv');
	
	unless(-s $ustacks_alleles_file and -s $ustacks_snps_file and -s $ustacks_tags_file){
		warn "Executing ustacks.....\n\n";
		my $ustacksCmd  = "$ustacks -t fastq -f $fastq_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -M $max_nuc_distance_ustacks -N $max_align_distance_ustacks -p $num_cpu_cores -d -r -R --model_type snp --alpha $alpha_value --max_locus_stacks $max_locus_stacks";
		warn $ustacksCmd . "\n\n";
		system($ustacksCmd) == 0 or die "Error calling $ustacksCmd: $?";
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
	
	my $num_cpu_cores = shift;
	die "Error lost the number of cores for stacks" unless defined $num_cpu_cores;
	
	my $cstacks_file = join('/', $stacks_output_dir, "batch_1");
	
	my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file);
	$cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
	$cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
	$cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
	
	unless(-s $cstacks_alleles_file and -s $cstacks_snps_file and -s $cstacks_tags_file){
	
		my ($ustacks_tags_files, $ustacks_tags_file_count) = find_files($stacks_input_dir, "tags.tsv");
		my @cstacks_soptions = ();
		foreach my $file_name (sort keys %{$ustacks_tags_files}){
			my $ustacks_tags_infile = $ustacks_tags_files->{$file_name};
			
			# Get the basename of the tags filename without the .tags.tsv extension.
			my $ustacks_filename = fileparse($ustacks_tags_infile, qr/\.tags.tsv/);
			
			warn "Processing " . $ustacks_filename . ".....\n";
			my $cstacks_infile = join('/', $stacks_input_dir, $ustacks_filename);
			push(@cstacks_soptions, "-s $cstacks_infile ");
		}
		
		warn "Executing cstacks.....\n\n";
		my $cstacks_joined_soptions = join("\\\n", @cstacks_soptions);
		my $cstacksCmd  = "$cstacks -b 1 -o $stacks_output_dir -p $num_cpu_cores \\\n $cstacks_joined_soptions";
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
	
	my $num_cpu_cores = shift;
	die "Error lost the number of cores for stacks" unless defined $num_cpu_cores;
	
	my $stacks_output_dir = shift;
	die "Error lost the stacks output file directory" unless defined $stacks_output_dir;
	
 	# Get the basename of the sam filename without the .sam extension.
 	my $stacks_sample_filename = fileparse($stacks_sample_infile, qr//);

 	my $sstacks_matches_file = join('/', $stacks_output_dir, $stacks_sample_filename . '.matches.tsv');
	
	unless(-s $sstacks_matches_file){
		warn "Executing sstacks.....\n\n";
		my $sstacksCmd  = "$sstacks -b 1 -c $stacks_catalog_infile -s $stacks_sample_infile -o $stacks_output_dir -p $num_cpu_cores";
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
	
	my $project_name = shift;
	die "Error lost the project name" unless defined $project_name;
	
	my ($fastq_filename, $fastq_dir) = fileparse($fastq_file, qr/\.fastq.gz/);
	
	#9315_CGCCTTAT_5_E11_trimmed_offset_3.fastq.gz
	$fastq_filename =~ s/\_[ACGT]+\_\d+\_[A-Z]+\d+\_trimmed\_offset\_\d+//g;
	my $uncompressed_fastq_file = join('/', $fastq_dir, join("_", $fastq_filename, $project_name) . ".fastq");
	
	unless(-s $uncompressed_fastq_file){
		warn "Calling gunzip for $fastq_file....\n";
		my $gunzipCmd  = "$gunzip -c $fastq_file > $uncompressed_fastq_file";
		warn $gunzipCmd . "\n\n";
		system($gunzipCmd) == 0 or die "Error calling $gunzipCmd: $?";
	}
	return $uncompressed_fastq_file;
}
