#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use IPC::Open2;
use File::Basename;
use File::Copy;

# VERSION 1.25

#### PROGRAM NAME ####
# refgen_stacks_analysis.pl - Program that runs the Stacks reference genome pipeline using the pstacks, cstacks, and sstacks programs of the Stacks Software Suite.

#### DESCRIPTION ####
# This program takes sorted bam or sam alignment file directory as input and executes the pstacks program, which extracts sequence stacks that were aligned to a reference genome using the BWA alignment program and identifies SNPs. These sequence stacks are then processed using cstacks and sstacks to obtain the filtered SNP stacks output files.

#### SAMPLE COMMANDS ####

## Standard example
# perl refgen_stacks_analysis_pipeline.pl -i ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR_UNPADDED/STEPHEN_TREVOY/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_sequence_data/DendPond_male_1.0/Primary_Assembly/unplaced_scaffolds/FASTA/DendPond_male_1.0_unplaced.scaf.fa -c 7 -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3

## Using the cstacks catalog prefix path.
# perl ~/refgen_stacks_analysis_pipeline.pl -i ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_TRIMMED_GBS/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_MALE_GBS_STACKS/REFERENCE_GENOME/DendPond_male_1.0_unplaced.scaf.fa.fasta -p ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_MALE_GBS_STACKS/STACKS_OUTFILES/batch_1 -t fastq -f bam -l 92 -b 1 -d 5 -a 0.05 -s 1 -m 20 -c 24 -o ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/test_output_dir/
my ($bwa_input_dir, $pstacks_infile_type, $pstacks_model_type, $cstacks_catalog_prefix, $use_existing_catalog, $gbs_sequence_length, $stacks_sql_id, $min_depth_coverage_pstacks, $alpha_value_pstacks, $bound_low_value_pstacks, $bound_high_value_pstacks, $num_mismatches_tag, $num_threads, $output_dir);
GetOptions(
'i=s'    => \$bwa_input_dir, # The absolute path to the padded alignment file directory that contains files with the extension .sam or .bam for each individual within the Genotyping by Sequencing (GBS) project.
'f=s'    => \$pstacks_infile_type, # The type of alignment file used for input into the pstacks program. Can either be "bam" or "sam". Default: bam.
'm=s'    => \$pstacks_model_type, # The type of snp model used for input into the pstacks program. Can either be "snp" or "bounded" or "fixed". Default: snp
'p=s'    => \$cstacks_catalog_prefix, # The cstacks catalog prefix file path consisting of a set of consensus loci built from a set of samples processed by the pstacks program. i.e. /path/to/catalog/batch_1
'e=s'    => \$use_existing_catalog, # The use existing catalog bitwise flag. If true add data to the existing catalog, Otherwise run the default cstacks program. Default: false
'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
'b=s'    => \$stacks_sql_id, # The SQL ID to insert into the output to identify this sample. Default: 1
'd=s'    => \$min_depth_coverage_pstacks, # The minimum depth of coverage to report a stack. Default: 5
'a=s'    => \$alpha_value_pstacks, # The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
'w=s'    => \$bound_low_value_pstacks, # Lower bound for epsilon, the error rate, between 0 and 1.0. Only used if pstacks_model_type is set to "bounded". Default: 0
'h=s'    => \$bound_high_value_pstacks, # Upper bound for epsilon, the error rate, between 0 and 1.0. Only used if pstacks_model_type is set to "bounded". Default: 1
's=s'    => \$num_mismatches_tag, # The number of mismatches allowed between sample tags when generating the catalog. Default: 1
'c=s'    => \$num_threads, # The number of cpu threads to use for the stacks programs. Default: 2
'o=s'    => \$output_dir, # The absolute path to the output directory to contain the renumerated reference genome, BWA sam, padded sam, and Stacks output files and directories.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
    defined $bwa_input_dir
    and defined $output_dir
);

# The type of alignment file used for input into the pstacks program. Can either be "bam" or "sam". If "sam" is specified the keep_sam_output_files cannot be set to false. Default: bam
$pstacks_infile_type = "bam" unless defined $pstacks_infile_type;

# The type of snp model used for input into the pstacks program. Can either be "snp" or "bounded". Default: snp
$pstacks_model_type = "snp" unless defined $pstacks_model_type;

# The cstacks catalog prefix file path consisting of a set of consensus loci built from a set of samples processed by the pstacks program. i.e. /path/to/catalog/batch_1
$cstacks_catalog_prefix = "" unless defined $cstacks_catalog_prefix;

# The use existing catalog bitwise flag. If true add data to the existing catalog, Otherwise run the default cstacks program. Default: false
$use_existing_catalog = "false" unless defined $use_existing_catalog;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

# The SQL ID to insert into the output to identify this sample. Default: 1
$stacks_sql_id = 1 unless defined $stacks_sql_id;

# The minimum depth of coverage to report a stack. Default: 5
$min_depth_coverage_pstacks = 5 unless defined $min_depth_coverage_pstacks;

# The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
$alpha_value_pstacks = 0.05 unless defined $alpha_value_pstacks;

# Lower bound for epsilon, the error rate, between 0 and 1.0. Only used if pstacks_model_type is set to "bounded". Default: 0
$bound_low_value_pstacks = 0 unless defined $bound_low_value_pstacks;

# Upper bound for epsilon, the error rate, between 0 and 1.0. Only used if pstacks_model_type is set to "bounded". Default: 1
$bound_high_value_pstacks = 1 unless defined $bound_high_value_pstacks;

# The number of mismatches allowed between sample tags when generating the catalog. Default: 1
$num_mismatches_tag = 1 unless defined $num_mismatches_tag;

# The number of cpu threads to use for the stacks programs. Default: 2
$num_threads = 2 unless defined $num_threads;

# Program dependencies - The absolute paths to the Stacks related programs.
my ($pstacks, $cstacks, $sstacks);
$pstacks			= $ENV{"HOME"} . '/software/stacks-1.47/pstacks';
$cstacks			= $ENV{"HOME"} . '/software/stacks-1.47/cstacks';
$sstacks			= $ENV{"HOME"} . '/software/stacks-1.47/sstacks';

sub usage {
    
    die <<"USAGE";
    
Usage: $0 -i bwa_input_dir -f pstacks_input_file_type -m pstacks_model_type -p cstacks_catalog_prefix -e use_existing_catalog -l gbs_sequence_length -b stacks_sql_id -d min_depth_coverage_pstacks -a alpha_value_pstacks -s num_mismatches_tag -c num_threads -o output_dir
    
    
    VERSION 1.25
    
    DESCRIPTION - This program takes sorted bam or sam alignment file directory as input and executes the pstacks program, which extracts sequence stacks that were aligned to a reference genome using the BWA alignment program and identifies SNPs. These sequence stacks are then processed using cstacks and sstacks to obtain the filtered SNP stacks output files.
    
    OPTIONS:

    -i bwa_input_dir - The absolute path to the padded alignment file directory that contains files with the extension .sam or .bam for each individual within the Genotyping by Sequencing (GBS) project.

    -f pstacks_input_file_type - The type of alignment file used for input into the pstacks program. Can either be "bam" or "sam". Default: bam
    
    -m pstacks_model_type - The type of snp model used for input into the pstacks program. Can either be "snp" or "bounded" or "fixed". Default: snp
        
    -p cstacks_catalog_prefix - The cstacks catalog prefix file path consisting of a set of consensus loci built from a set of samples processed by the pstacks program. i.e. /path/to/catalog/batch_1

    -e use_existing_catalog - The use existing catalog bitwise flag. If true add data to the existing catalog, Otherwise run the default cstacks program. Default: false

    -l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92

    -b stacks_sql_id - The SQL ID to insert into the output to identify this sample. Default: 1

    -d min_depth_coverage_pstacks - The minimum depth of coverage to report a stack. Default: 5

    -a alpha_value_pstacks - The chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05, 0.01, or 0.001. Default: 0.05
    
    -w bound_low_value_pstacks - Lower bound for epsilon, the error rate, between 0 and 1.0. Only used if pstacks_model_type is set to "bounded". Default: 0
    
    -h bound_high_value_pstacks - Upper bound for epsilon, the error rate, between 0 and 1.0. Only used if pstacks_model_type is set to "bounded". Default: 1
        
    -s num_mismatches_tag - The number of mismatches allowed between sample tags when generating the catalog. Default: 1

    -c num_threads - The number of cpu cores to use for the stacks programs. Default: 2

    -o output_dir - The absolute path to the output directory to contain the renumerated reference genome, BWA, and Stacks output files and directories.

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

# Pstacks parameters
my %params = ();
if($pstacks_model_type eq "snp"){
    $params{"MODEL_TYPE"} = $pstacks_model_type;
    $params{"ALPHA_VALUE"} = $alpha_value_pstacks;
}elsif($pstacks_model_type eq "bounded"){
    $params{"MODEL_TYPE"} = $pstacks_model_type;
    $params{"ALPHA_VALUE"} = $alpha_value_pstacks;
    $params{"LOW_VALUE"} = $bound_low_value_pstacks;
    $params{"HIGH VALUE"} = $bound_high_value_pstacks;
}

if($pstacks_infile_type eq "bam"){
    # Find all padded bam alignment output files from the padded bam output output directory with the extension *.bam.
    my ($padded_bam_files, $padded_bam_file_count) = find_files($bwa_input_dir, "bam");
    
    # Iterate through each padded bam alignment output file with extension *.sam and execute the pstacks program.
    my $pstacks_sql_id = 1;
    foreach my $file_name (sort keys %{$padded_bam_files}){
        warn "Processing " . $file_name . ".....\n";
        my $padded_bam_infile = $padded_bam_files->{$file_name};
        
        # Execute the pstacks program, which extract stacks that have been aligned to a reference genome and identify SNPs. These stacks can then be processed with cstacks and sstacks.
        pstacks($padded_bam_infile, "bam", $pstacks_sql_id, $min_depth_coverage_pstacks, $num_threads, \%params, $stacks_output_dir);
        
        $pstacks_sql_id++;
    }
}elsif($pstacks_infile_type eq "sam"){
    # Find all padded sam alignment output files from the padded sam output output directory with the extension *.sam.
    my ($padded_sam_files, $padded_sam_file_count) = find_files($bwa_input_dir, "sam");
    
    # Iterate through each padded sam alignment output file with extension *.sam and execute the pstacks program.
    my $pstacks_sql_id = 1;
    foreach my $file_name (sort keys %{$padded_sam_files}){
        warn "Processing " . $file_name . ".....\n";
        my $padded_sam_infile = $padded_sam_files->{$file_name};
        
        # Execute the pstacks program, which extract stacks that have been aligned to a reference genome and identify SNPs. These stacks can then be processed with cstacks and sstacks.
        pstacks($padded_sam_infile, "sam", $pstacks_sql_id, $min_depth_coverage_pstacks, $num_threads, \%params, $stacks_output_dir);
        
        $pstacks_sql_id++;
    }
}

# Execute the cstacks program to build a catalog from a set of samples processed by the pstacks program if a catalog tags file was not specified. The cstacks program creates a set of consensus loci, merging alleles together. If a catalog tags file was specified and the use_existing_catalog bitwise flag is set to true, then execute the cstacks program to add data to an existing catalog. If a catalog tags file was specified and the use_existing_catalog bitwise flag is set to false, then skip the execution of the cstacks program and use the catalog when executing sstacks.
my $cstacks_file = "";
if($cstacks_catalog_prefix eq ""){
    $cstacks_file = cstacks($stacks_output_dir, $pstacks_infile_type, $cstacks_catalog_prefix, $stacks_sql_id, $num_mismatches_tag, $num_threads);
    
}elsif($cstacks_catalog_prefix =~ m/batch_\d+$/){
    
    if($use_existing_catalog eq "true"){
        $cstacks_file = cstacks($stacks_output_dir, $pstacks_infile_type, $cstacks_catalog_prefix, $stacks_sql_id, $num_mismatches_tag, $num_threads);
    }elsif($use_existing_catalog eq "false"){
        
        # The existing cstacks catalog alleles, snps, and tags file paths.
        my ($existing_cstacks_alleles_file, $existing_cstacks_snps_file, $existing_cstacks_tags_file) = "";
        if($pstacks_infile_type eq "bam"){
            $existing_cstacks_alleles_file = $cstacks_catalog_prefix . '.catalog.alleles.tsv.gz';
            $existing_cstacks_snps_file = $cstacks_catalog_prefix . '.catalog.snps.tsv.gz';
            $existing_cstacks_tags_file = $cstacks_catalog_prefix . '.catalog.tags.tsv.gz';
        }elsif($pstacks_infile_type eq "sam"){
            $existing_cstacks_alleles_file = $cstacks_catalog_prefix . '.catalog.alleles.tsv';
            $existing_cstacks_snps_file = $cstacks_catalog_prefix . '.catalog.snps.tsv';
            $existing_cstacks_tags_file = $cstacks_catalog_prefix . '.catalog.tags.tsv';
        }
        
        # The cstacks batch file prefix file path.
        $cstacks_file = join('/', $stacks_output_dir, join("_", "batch", $stacks_sql_id));
        
        # The cstacks catalog alleles, snps, and tags file paths.
        my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file) = "";
        if($pstacks_infile_type eq "bam"){
            $cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv.gz';
            $cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv.gz';
            $cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv.gz';
        }elsif($pstacks_infile_type eq "sam"){
            $cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
            $cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
            $cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
        }
        warn "Copying $existing_cstacks_alleles_file to $cstacks_alleles_file.....";
        copy($existing_cstacks_alleles_file, $cstacks_alleles_file) or die "Copy failed: $!";
        
        warn "Copying $existing_cstacks_snps_file to $cstacks_snps_file.....";
        copy($existing_cstacks_snps_file, $cstacks_snps_file) or die "Copy failed: $!";
        
        warn "Copying $existing_cstacks_tags_file to $cstacks_tags_file.....";
        copy($existing_cstacks_tags_file, $cstacks_tags_file) or die "Copy failed: $!";
    }
}

my ($pstacks_tags_files, $pstacks_tags_file_count) = "";
if($pstacks_infile_type eq "bam"){
    
    # Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.gz.
    ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv.gz");
    
}elsif($pstacks_infile_type eq "sam"){
    
    # Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.
    ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");
    
}

# Iterate through each tags output file with extension *.tags.tsv and execute the sstacks program.
foreach my $file_name (sort keys %{$pstacks_tags_files}){
    if(($file_name !~ m/batch_\d+\.catalog\.tags\.tsv/)
        or $file_name !~ m/batch_\d+\.catalog\.tags\.tsv\.gz/){ # If *.tags.tsv file does not match batch_*.catalog.tags.tsv.
            
            my $pstacks_tags_infile = $pstacks_tags_files->{$file_name};
            
            my $pstacks_filename = "";
            if($pstacks_infile_type eq "bam"){
                
                # Get the basename of the tags filename without the .tags.tsv.gz extension.
                $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv.gz/);
                
            }elsif($pstacks_infile_type eq "sam"){
                
                # Get the basename of the tags filename without the .tags.tsv extension.
                $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv/);
                
            }
            
            warn "Processing " . $pstacks_filename . ".....\n";
            
            my $sstacks_infile = join('/', $stacks_output_dir, $pstacks_filename);
            
            # Execute the sstacks program. Sets of stacks constructed by the pstacks program is searched against the catalog produced by cstacks.
            sstacks($cstacks_file, $sstacks_infile, $pstacks_infile_type, $stacks_sql_id, $num_threads, $stacks_output_dir);
        }
}

# pstacks($alignment_infile, $alignment_file_type, $sql_id, $min_depth_coverage, $num_threads, $alpha_value) - Executes the pstacks program in the Stacks Software Suite. Extracts stacks that have been aligned to a reference genome and identify SNPs. These stacks can then be processed with cstacks and sstacks.
#
# Input paramater(s):
#
# $alignment_infile - The padded alignment input file.
#
# $alignment_file_type - The padded alignment input file type.
#
# $sql_id - The SQL ID to insert into the output to identify this sample.
#
# $min_depth_coverage - The minimum depth of coverage to report a stack.
#
# $num_threads - The number of threads to use for pstacks.
#
# $alpha_value - The chi square significance level required to call a heterozygote or homozygote.

# pstacks -t bam -f ../PADDED_BAM_OUTFILES/LL-06_MPB-MALE-GBS.bam -o .. -i 1 -m 1 -p 7 --model_type snp --alpha 0.05
# pstacks -t sam -f ../PADDED_SAM_OUTFILES/LL-06_MPB-MALE-GBS.sam -o .. -i 1 -m 1 -p 7 --model_type snp --alpha 0.05
sub pstacks{
    
    # The padded alignment input file.
    my $alignment_infile = shift;
    die "Error lost the alignment input file" unless defined $alignment_infile;
    
    # The padded alignment input file type.
    my $alignment_file_type = shift;
    die "Error lost the alignment input file type" unless defined $alignment_file_type;
    
    # The SQL ID to insert into the output to identify this sample.
    my $sql_id = shift;
    die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $sql_id;
    
    # The minimum depth of coverage to report a stack.
    my $min_depth_coverage = shift;
    die "Error lost the minimum depth of coverage to report a stack" unless defined $min_depth_coverage;
    
    # The number of threads to use for pstacks.
    my $num_threads = shift;
    die "Error lost the number of cores for stacks" unless defined $num_threads;
    
#    # The chi square significance level required to call a heterozygote or homozygote.
#    my $alpha_value = shift;
#    die "Error lost the chi square significance level required to call a heterozygote or homozygote" unless defined $alpha_value;
#
    my $params = shift;
    die "Error lost the parameters for pstacks" unless defined $params;
    
    # The stacks output directory that contains the results from the pstacks program.
    my $stacks_output_dir = shift;
    die "Error lost the stacks output file directory" unless defined $stacks_output_dir;
    
    
    # Format the pstacks individual alleles, snps, and tags output files.
    my ($alignment_filename, $pstacks_alleles_file, $pstacks_snps_file, $pstacks_tags_file) = "";
    if($alignment_file_type eq "bam"){
        
        # Get the basename of the filename without the extension.
        $alignment_filename = fileparse($alignment_infile, qr/\.bam/);
        $pstacks_alleles_file = join('/', $stacks_output_dir, $alignment_filename . '.alleles.tsv.gz');
        $pstacks_snps_file = join('/', $stacks_output_dir, $alignment_filename . '.snps.tsv.gz');
        $pstacks_tags_file = join('/', $stacks_output_dir, $alignment_filename . '.tags.tsv.gz');
        
    }elsif($alignment_file_type eq "sam"){
        
        # Get the basename of the filename without the extension.
        $alignment_filename = fileparse($alignment_infile, qr/\.sam/);
        $pstacks_alleles_file = join('/', $stacks_output_dir, $alignment_filename . '.alleles.tsv');
        $pstacks_snps_file = join('/', $stacks_output_dir, $alignment_filename . '.snps.tsv');
        $pstacks_tags_file = join('/', $stacks_output_dir, $alignment_filename . '.tags.tsv');
        
    }
    
    #### USED PARAMETERS ####
    # t — input file Type. Supported types: bowtie, sam, or bam.
    # f — input file path.
    # o — output path to write results.
    # i — SQL ID to insert into the output to identify this sample.
    # m — minimum depth of coverage to report a stack (default 1).
    # p — enable parallel execution with num_threads threads.
    
    # Model options:
    # --model_type [type] — either 'snp' (default), 'bounded', or 'fixed'
    
    # For the SNP or Bounded SNP model:
    # --alpha [num] — chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.
    
    #### UNUSED PARAMETERS ####
    # h — display this help messsage.
    
    # For the Bounded SNP model:
    # --bound_low [num] — lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).
    # --bound_high [num] — upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).
    
    # For the Fixed model:
    # --bc_err_freq [num] — specify the barcode error frequency, between 0 and 1.0.
    
    # Create the stacks log output directory if it doesn't already exist.
    my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
    unless(-d $stacks_log_output_dir){
        mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
    }
    
    # The standard output file for the pstacks program.
    my $pstacks_log_outfile = join('/', $stacks_log_output_dir, join("_", $alignment_filename, "pstacks.log"));
#
    # Initialize Pstacks parameters from params hash reference
    my ($pstacks_model_type, $alpha_value_pstacks, $bound_low_value_pstacks, $bound_high_value_pstacks) = "";
    
    $pstacks_model_type = $params->{"MODEL_TYPE"};
    if($pstacks_model_type eq "snp"){
        $alpha_value_pstacks = $params->{"ALPHA_VALUE"};
    }elsif($pstacks_model_type eq "bounded"){
        $alpha_value_pstacks = $params->{"ALPHA_VALUE"};
        $bound_low_value_pstacks = $params->{"LOW_VALUE"};
        $bound_high_value_pstacks = $params->{"HIGH VALUE"};
    }
    
    # Execute the pstacks program if the pstacks alleles, snps, and tags output files are not already generated.
    unless(-s $pstacks_alleles_file and -s $pstacks_snps_file and -s $pstacks_tags_file){
        warn "Executing pstacks.....\n\n";
        my $pstacksCmd = "";
        if($alignment_file_type eq "bam"){
            if($pstacks_model_type eq "snp"){
                $pstacksCmd  = "$pstacks -t bam -f $alignment_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -p $num_threads --model_type snp --alpha $alpha_value_pstacks 2> $pstacks_log_outfile";
            }elsif($pstacks_model_type eq "bounded"){
                $pstacksCmd  = "$pstacks -t bam -f $alignment_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -p $num_threads --model_type bounded --alpha $alpha_value_pstacks --bound_low $bound_low_value_pstacks --bound_high $bound_high_value_pstacks 2> $pstacks_log_outfile";
            }
        }elsif($alignment_file_type eq "sam"){
            
            if($pstacks_model_type eq "snp"){
                $pstacksCmd  = "$pstacks -t sam -f $alignment_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -p $num_threads --model_type snp --alpha $alpha_value_pstacks 2> $pstacks_log_outfile";
            }elsif($pstacks_model_type eq "bounded"){
                $pstacksCmd  = "$pstacks -t sam -f $alignment_infile -o $stacks_output_dir -i $sql_id -m $min_depth_coverage -p $num_threads --model_type bounded --alpha $alpha_value_pstacks --bound_low $bound_low_value_pstacks --bound_high $bound_high_value_pstacks 2> $pstacks_log_outfile";
            }
        }
        
        warn $pstacksCmd . "\n\n";
        system($pstacksCmd) == 0 or die "Error calling $pstacksCmd: $?";
    }
}

# $cstacks_file = cstacks($stacks_output_dir, $cstacks_catalog_prefix, $stacks_sql_id, $num_mismatches_tag, $num_threads) - Executes the cstacks program in the Stacks Software Suite. Build a catalog from a set of samples processed by the ustacks program. Creates a set of consensus loci, merging alleles together.
#
# Input paramater(s):
#
# $stacks_output_dir - The stacks directory that contains the results from the pstacks program.
#
# $cstacks_catalog_prefix - The absolute path to an existing catalog. If path is not defined execute the default cstacks program, otherwise add data to an existing catalog.
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

# cstacks -b 1 -o /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -g -n 1 -p 7 \
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
    
    # The padded alignment input file type.
    my $alignment_file_type = shift;
    die "Error lost the padded alignment input file type" unless defined $alignment_file_type;
    
    # The absolute path to an existing catalog. If path is not defined execute the default cstacks program, otherwise add data to an existing catalog.
    my $cstacks_catalog_prefix = shift;
    die "Error lost he absolute path to an existing catalog. If path is not defined execute the default cstacks program, otherwise add data to an existing catalog" unless defined $cstacks_catalog_prefix;
    
    # The SQL ID to insert into the output to identify this sample.
    my $stacks_sql_id = shift;
    die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $stacks_sql_id;
    
    # The number of mismatches allowed between sample tags when generating the catalog.
    my $num_mismatches_tag = shift;
    die "Error lost the number of mismatches allowed between sample tags when generating the catalog" unless defined $num_mismatches_tag;
    
    # The number of threads to use for cstacks.
    my $num_threads = shift;
    die "Error lost the number of threads to use for cstacks" unless defined $num_threads;
    
    # Create the stacks log output directory if it doesn't already exist.
    my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
    unless(-d $stacks_log_output_dir){
        mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
    }
    
    # The standard out log file for the cstacks program.
    my $cstacks_log_outfile = join('/', $stacks_log_output_dir, "cstacks.log");
    
    # The cstacks input file
    my $cstacks_file = "";
    if($cstacks_catalog_prefix eq ""){
        
        # The cstacks batch file prefix file path.
        $cstacks_file = join('/', $stacks_output_dir, join("_", "batch", $stacks_sql_id));
        
        # The cstacks catalog alleles, snps, and tags file paths.
        my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file) = "";
        if($alignment_file_type eq "bam"){
            $cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv.gz';
            $cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv.gz';
            $cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv.gz';
        }elsif($alignment_file_type eq "sam"){
            $cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
            $cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
            $cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
        }
        
        
        # Execute the cstacks program if the following files are not already generated.
        unless(-s $cstacks_alleles_file and -s $cstacks_snps_file and -s $cstacks_tags_file){
            
            my ($pstacks_tags_files, $pstacks_tags_file_count) = "";
            if($alignment_file_type eq "bam"){
                
                # Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.gz.
                ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv.gz");
                
            }elsif($alignment_file_type eq "sam"){
                
                # Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.
                ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");
                
            }
            
            my @cstacks_soptions = ();
            # Iterate through each pstacks tags output file with extension *.tags.tsv and execute the cstacks program.
            foreach my $file_name (sort keys %{$pstacks_tags_files}){
                my $pstacks_tags_infile = $pstacks_tags_files->{$file_name};
                
                my $pstacks_filename = "";
                if($alignment_file_type eq "bam"){
                    
                    # Get the basename of the tags filename without the .tags.tsv.gz extension.
                    $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv.gz/);
                    
                }elsif($alignment_file_type eq "sam"){
                    
                    # Get the basename of the tags filename without the .tags.tsv extension.
                    $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv/);
                    
                }
                
                
                # Obtain a list of all the file prefix paths specifed by the -s option for cstacks.
                warn "Processing " . $pstacks_filename . ".....\n";
                my $cstacks_infile = join('/', $stacks_output_dir, $pstacks_filename);
                push(@cstacks_soptions, "-s $cstacks_infile ");
            }
            
            #### USED PARAMETERS ####
            # b — MySQL ID of this batch.
            # o — output path to write results.
            # s — TSV file from which to load radtags.
            # g — base catalog matching on genomic location, not sequence identity.
            # n — number of mismatches allowed between sample tags when generating the catalog.
            # p — enable parallel execution with num_threads threads.
            
            #### UNUSED PARAMETERS ####
            # m — include tags in the catalog that match to more than one entry.
            # h — display this help messsage.
            
            # Catalog editing:
            # --catalog [path] — provide the path to an existing catalog. cstacks will add data to this existing catalog.
            
            # Advanced options:
            # --report_mmatches — report query loci that match more than one catalog locus.
            
            warn "Executing cstacks.....\n\n";
            my $cstacks_joined_soptions = join("\\\n", @cstacks_soptions);
            my $cstacksCmd  = "$cstacks -b $stacks_sql_id -o $stacks_output_dir -g -n $num_mismatches_tag -p $num_threads \\\n $cstacks_joined_soptions 2> $cstacks_log_outfile";
            warn $cstacksCmd . "\n\n";
            system($cstacksCmd) == 0 or die "Error calling $cstacksCmd: $?";
        }
    }elsif($cstacks_catalog_prefix =~ m/batch_\d+$/){
        
        
        # The existing cstacks catalog alleles, snps, and tags file paths.
        my ($existing_cstacks_alleles_file, $existing_cstacks_snps_file, $existing_cstacks_tags_file);
        if($alignment_file_type eq "bam"){
            $existing_cstacks_alleles_file = $cstacks_catalog_prefix . '.catalog.alleles.tsv.gz';
            $existing_cstacks_snps_file = $cstacks_catalog_prefix . '.catalog.snps.tsv.gz';
            $existing_cstacks_tags_file = $cstacks_catalog_prefix . '.catalog.tags.tsv.gz';
        }elsif($alignment_file_type eq "sam"){
            $existing_cstacks_alleles_file = $cstacks_catalog_prefix . '.catalog.alleles.tsv';
            $existing_cstacks_snps_file = $cstacks_catalog_prefix . '.catalog.snps.tsv';
            $existing_cstacks_tags_file = $cstacks_catalog_prefix . '.catalog.tags.tsv';
        }
        
        # The cstacks batch file prefix file path.
        $cstacks_file = join('/', $stacks_output_dir, join("_", "batch", $stacks_sql_id));
        
        # The cstacks catalog alleles, snps, and tags file paths.
        my ($cstacks_alleles_file, $cstacks_snps_file, $cstacks_tags_file);
        if($alignment_file_type eq "bam"){
            $cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv.gz';
            $cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv.gz';
            $cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv.gz';
        }elsif($alignment_file_type eq "sam"){
            $cstacks_alleles_file = $cstacks_file . '.catalog.alleles.tsv';
            $cstacks_snps_file = $cstacks_file . '.catalog.snps.tsv';
            $cstacks_tags_file = $cstacks_file . '.catalog.tags.tsv';
        }
        
        warn "Copying $existing_cstacks_alleles_file to $cstacks_alleles_file.....";
        copy($existing_cstacks_alleles_file, $cstacks_alleles_file) or die "Copy failed: $!";
        
        warn "Copying $existing_cstacks_snps_file to $cstacks_snps_file.....";
        copy($existing_cstacks_snps_file, $cstacks_snps_file) or die "Copy failed: $!";
        
        warn "Copying $existing_cstacks_tags_file to $cstacks_tags_file.....";
        copy($existing_cstacks_tags_file, $cstacks_tags_file) or die "Copy failed: $!";
        
        my ($pstacks_tags_files, $pstacks_tags_file_count) = "";
        if($alignment_file_type eq "bam"){
            # Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.gz.
            ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv.gz");
        }elsif($alignment_file_type eq "sam"){
            # Find all pstacks tags output files from the stacks output directory with the extension *.tags.tsv.
            ($pstacks_tags_files, $pstacks_tags_file_count) = find_files($stacks_output_dir, "tags.tsv");
        }
        
        my @cstacks_soptions = ();
        # Iterate through each pstacks tags output file with extension *.tags.tsv and execute the cstacks program.
        foreach my $file_name (sort keys %{$pstacks_tags_files}){
            my $pstacks_tags_infile = $pstacks_tags_files->{$file_name};
            
            my $pstacks_filename = "";
            if($alignment_file_type eq "bam"){
                
                # Get the basename of the tags filename without the .tags.tsv.gz extension.
                $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv.gz/);
                
            }elsif($alignment_file_type eq "sam"){
                
                # Get the basename of the tags filename without the .tags.tsv extension.
                $pstacks_filename = fileparse($pstacks_tags_infile, qr/\.tags.tsv/);
                
            }
            
            # Obtain a list of all the file prefix paths specifed by the -s option for cstacks.
            warn "Processing " . $pstacks_filename . ".....\n";
            my $cstacks_infile = join('/', $stacks_output_dir, $pstacks_filename);
            push(@cstacks_soptions, "-s $cstacks_infile ");
        }
        
        #### USED PARAMETERS ####
        # b — MySQL ID of this batch.
        # o — output path to write results.
        # s — TSV file from which to load radtags.
        # g — base catalog matching on genomic location, not sequence identity.
        # n — number of mismatches allowed between sample tags when generating the catalog.
        # p — enable parallel execution with num_threads threads.
        
        #### UNUSED PARAMETERS ####
        # m — include tags in the catalog that match to more than one entry.
        # h — display this help messsage.
        
        # Catalog editing:
        # --catalog [path] — provide the path to an existing catalog. cstacks will add data to this existing catalog.
        
        # Advanced options:
        # --report_mmatches — report query loci that match more than one catalog locus.
        
        warn "Executing cstacks.....\n\n";
        my $cstacks_joined_soptions = join("\\\n", @cstacks_soptions);
        my $cstacksCmd  = "$cstacks --catalog $cstacks_file -b $stacks_sql_id -o $stacks_output_dir -g -n $num_mismatches_tag -p $num_threads \\\n $cstacks_joined_soptions 2> $cstacks_log_outfile";
        warn $cstacksCmd . "\n\n";
        system($cstacksCmd) == 0 or die "Error calling $cstacksCmd: $?";
    }
    
    # Returns the cstacks batch file prefix file path.
    return $cstacks_file;
}

# sstacks($stacks_catalog_infile, $stacks_sample_infile, $alignment_file_type, $stacks_sql_id, $num_threads, $stacks_output_dir) - Executes the sstacks program in the Stacks Software Suite. Sets of stacks constructed by the pstacks program is searched against the catalog produced by cstacks.
#
# Input paramater(s):
#
# $stacks_catalog_infile - The stacks catalog input batch_*.catalog.*.tsv file from which to load the catalog GBS-Tags.
#
# $stacks_sample_infile - The stacks sample input *.tags.tsv file from which to load sample GBS-Tags.
#
# $alignment_file_type - The padded alignment input file type.
#
# $stacks_sql_id - The SQL ID to insert into the output to identify this sample.
#
# $num_threads - The number of threads to use for sstacks.
#
# $stacks_output_dir - The stacks output file directory that contains all the results files generated by sstacks.
sub sstacks{
    
    # The stacks catalog input batch_*.catalog.*.tsv file from which to load the catalog GBS-Tags.
    my $stacks_catalog_infile = shift;
    die "Error lost the stacks catalog input file" unless defined $stacks_catalog_infile;
    
    # The stacks sample input *.tags.tsv file from which to load sample GBS-Tags.
    my $stacks_sample_infile = shift;
    die "Error lost the stacks sample input file" unless defined $stacks_sample_infile;
    
    # The padded alignment input file type.
    my $alignment_file_type = shift;
    die "Error lost the padded alignment input file type" unless defined $alignment_file_type;
    
    # The SQL ID to insert into the output to identify this sample.
    my $stacks_sql_id = shift;
    die "Error lost the SQL ID to insert into the output to identify this sample" unless defined $stacks_sql_id;
    
    # The number of threads to use for sstacks.
    my $num_threads = shift;
    die "Error lost the number of threads to use for sstacks" unless defined $num_threads;
    
    # The stacks output file directory that contains all the results files generated by sstacks.
    my $stacks_output_dir = shift;
    die "Error lost the stacks output file directory that contains all files generated by sstacks" unless defined $stacks_output_dir;
    
    # Get the basename of the sam filename without the .sam extension.
    my $stacks_sample_filename = fileparse($stacks_sample_infile, qr//);
    
    # The sstacks output matches file.
    my $sstacks_matches_file = "";
    if($alignment_file_type eq "bam"){
        $sstacks_matches_file = join('/', $stacks_output_dir, $stacks_sample_filename . '.matches.tsv.gz');
    }elsif($alignment_file_type eq "sam"){
        $sstacks_matches_file = join('/', $stacks_output_dir, $stacks_sample_filename . '.matches.tsv');
    }
    
    # Create the stacks log output directory if it doesn't already exist.
    my $stacks_log_output_dir = join('/', $stacks_output_dir, "STACKS_LOG_FILES");
    unless(-d $stacks_log_output_dir){
        mkdir($stacks_log_output_dir, 0777) or die "Can't make directory: $!";
    }
    
    # The standard out log file for the sstacks program.
    my $sstacks_log_outfile = join('/', $stacks_log_output_dir, join("_", $stacks_sample_filename, "sstacks.log"));
    
    #### USED PARAMETERS ####
    # b — MySQL ID of this batch.
    # c — TSV file from which to load the catalog RAD-Tags.
    # s — TSV file from which to load sample RAD-Tags.
    # o — output path to write results.
    # g — base matching on genomic location, not sequence identity.
    # p — enable parallel execution with num_threads threads.
    
    #### UNUSED PARAMETERS ####
    # r — Load the TSV file of a single sample instead of a catalog.
    # x — don’t verify haplotype of matching locus.
    # v — print program version.
    # h — display this help messsage.
    
    # Execute the sstacks program if the matches sstacks results files are not already generated.
    unless(-s $sstacks_matches_file){
        warn "Executing sstacks.....\n\n";
        my $sstacksCmd  = "$sstacks -b $stacks_sql_id -c $stacks_catalog_infile -s $stacks_sample_infile -o $stacks_output_dir -g -p $num_threads 2> $sstacks_log_outfile";
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
