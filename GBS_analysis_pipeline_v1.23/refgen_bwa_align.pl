#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use IPC::Open2;
use File::Basename;
use File::Copy;

# VERSION 1.23

#### PROGRAM NAME ####
# refgen_bwa_align.pl - Program that aligns quality filtered, demultiplexed, and adapter trimmed GBS data sequences (unpadded) against a reference genome using the BWA alignment program.

#### DESCRIPTION ####
# This program takes the quality filtered, demultiplexed, and adapter trimmed *.fastq input files (unpadded) and a reference genome fasta input file as input. Converts a reference genome fasta file to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers. It then performs a BWA alignment to align the GBS fastq sequences to a reference genome.

#### SAMPLE COMMANDS ####

## Standard example
# perl refgen_bwa_align.pl -i ~/workspace/GBS_data-08-10-2013/PROCESSED_RADTAGS/TRIMMED_OFFSET_3_ADAPTOR_REGEX_PARALLEL_FASTQ_DIR_UNPADDED/STEPHEN_TREVOY/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_sequence_data/DendPond_male_1.0/Primary_Assembly/unplaced_scaffolds/FASTA/DendPond_male_1.0_unplaced.scaf.fa -c 7 -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3

## Using the cstacks catalog prefix path.
# perl refgen_bwa_align.pl -i ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_TRIMMED_GBS/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_MALE_GBS_STACKS/REFERENCE_GENOME/DendPond_male_1.0_unplaced.scaf.fa.fasta -p ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/MPB_MALE_GBS_STACKS/STACKS_OUTFILES/batch_1 -t fastq -f bam -w samse -k false -l 92 -m 20 -c 24 -o ~/scratch/JASMINE_JANES_GBS_Data-2015-09-15/MPB_GBS_Data/test_output_dir/
my ($gbs_fastq_dir, $gbs_fastq_file_type, $refgen_infile, $cstacks_catalog_prefix, $use_existing_catalog, $bwa_algorithm, $align_outfile_type, $keep_sam_output_files, $gbs_sequence_length, $sam_mapq_threshold, $num_threads, $output_dir);
GetOptions(
    'i=s'    => \$gbs_fastq_dir, # The absolute path to the quality filtered, demultiplexed, and adapter trimmed *.fastq input file (unpadded) directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.
    't=s'    => \$gbs_fastq_file_type, # The fastq input file type. Can be either fastq or gzfastq. Default: gzfastq
    'g=s'    => \$refgen_infile, # The absolute path to the reference genome input fasta file to align GBS fastq sequences.
    'w=s'    => \$bwa_algorithm, # The alignment algorithm that will be used in the BWA program. Can either be samse or mem. Default: samse
    'f=s'    => \$align_outfile_type, # The type of alignment output file generated. Can either be "bam" or "sam". If "sam" is specified the keep_sam_output_files cannot be set to false. Default: bam
    'k=s'    => \$keep_sam_output_files, # The bitwise flag to keep the sam outfiles from the BWA program using the samse or mem alignment algorithm. If true keep the sam output files, Otherwise delete the sam output files to save space. Cannot be used if the align_outfile_type option is set to "sam". Default: true
    'l=s'    => \$gbs_sequence_length, # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
    'm=s'    => \$sam_mapq_threshold, # The mapping quality score threshold filter for sam alignments. Any mapping quality score below this value is filtered out of the sam alignment files. Default: 20
    'c=s'    => \$num_threads, # The number of cpu threads to use for the stacks programs. Default: 2
    'o=s'    => \$output_dir, # The absolute path to the output directory to contain the renumerated reference genome, BWA sam, padded sam, and Stacks output files and directories.
);

# Print usage message if the following input parameters are not specified.
usage() unless (
    defined $gbs_fastq_dir
    and defined $refgen_infile
    and defined $output_dir
);

# The fastq input file type. Can be either fastq or gzfastq. Default: gzfastq
$gbs_fastq_file_type = 'gzfastq' unless defined $gbs_fastq_file_type;

# The alignment algorithm that will be used in the BWA program. Can either be samse or mem. Default: samse
$bwa_algorithm = "samse" unless defined $bwa_algorithm;

# The type of alignment output file generated. Can either be "bam" or "sam". If "sam" is specified the keep_sam_output_files cannot be set to false. Default: bam
$align_outfile_type = "bam" unless defined $align_outfile_type;

# The bitwise flag to keep the sam outfiles from the BWA program using the samse or mem alignment algorithm. If true keep the sam output files, Otherwise delete the sam output files to save space. Cannot be used if the align_outfile_type option is set to "sam". Default: false
$keep_sam_output_files = "true" unless defined $keep_sam_output_files;

# The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92
$gbs_sequence_length = 92 unless defined $gbs_sequence_length;

# The mapping quality score threshold filter for sam alignments. Any mapq score below this value is filtered out of the sam alignment files. Default: 20
$sam_mapq_threshold = 20 unless defined $sam_mapq_threshold;

# The number of cpu threads to use for the stacks programs. Default: 2
$num_threads = 2 unless defined $num_threads;

# Program dependencies - The absolute paths to gunzip to uncompress fastq.gz input files, BWA alignment program, and Stacks related programs.
my ($gunzip, $bwa, $samtools);
$gunzip				= '/bin/gunzip';
$bwa				= $ENV{"HOME"} . '/software/bwa-0.7.15/bwa';
$samtools			= $ENV{"HOME"} . '/software/samtools-1.3.1/samtools';

sub usage {
    
    die <<"USAGE";
    
    Usage: $0 -i gbs_fastq_dir -t gbs_fastq_file_type -g refgen_infile -w bwa_algorithm -f align_outfile_type -k keep_original_sam_output_files -l gbs_sequence_length -m sam_mapq_threshold -c num_threads -o output_dir
    
    
    VERSION 1.23
    
    DESCRIPTION - This program takes the quality filtered, demultiplexed, and adapter trimmed *.fastq input files (unpadded) and reference genome fasta input files as input. Converts the reference genome fasta file to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers. It then performs a BWA alignment to align the GBS fastq sequences to the reference genome.
    
    OPTIONS:
    
    -i gbs_fastq_dir - The absolute path to the quality filtered, demultiplexed, and adapter trimmed *.fastq input file (unpadded) directory that contains files with the extension .fastq for each individual within the Genotyping by Sequencing (GBS) project.

    -t gbs_fastq_file_type - The fastq input file type. Can be either fastq or gzfastq. Default: gzfastq

    -g refgen_infile - The absolute path to the reference genome input fasta file to align GBS fastq sequences using the BWA alignment program.

    -w bwa_algorithm - The alignment algorithm that will be used in the BWA program. Can either be samse or mem. Default: samse

    -f align_outfile_type - The type of alignment file used for input into the pstacks program. Can either be "bam" or "sam". If "sam" is specified the keep_sam_output_files cannot be set to false. Default: bam

    -k keep_sam_output_files - The bitwise flag to keep the sam outfiles from the BWA program using the samse or mem alignment algorithm. If true keep the sam output files, Otherwise delete the sam output files to save space. Cannot be used if the align_outfile_type option is set to "sam". Default: true

    -l gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences. Default: 92

    -m sam_mapq_threshold - The mapping quality score threshold filter for sam alignments. Any mapq score below this value is filtered out of the sam alignment files. Default: 20

    -c num_threads - The number of cpu cores to use for the stacks programs. Default: 2

    -o output_dir - The absolute path to the output directory to contain the renumerated reference genome, BWA, and Stacks output files and directories.

USAGE
}

die "Error: keep_sam_output_files option cannot be set to \"false\" when align_outfile_type is specified" if(($align_outfile_type eq "sam") and ($keep_sam_output_files eq "false"));

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
    mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Create the reference genome output directory if it doesn't already exist.
my $refgen_output_dir = join('/', $output_dir, "REFERENCE_GENOME");
unless(-d $refgen_output_dir){
    mkdir($refgen_output_dir, 0777) or die "Can't make directory: $!";
}

# Converts the reference genome to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers.
my ($refgen_fasta_outfile, $refgen_toc_outfile, $num_chromosomes) = convert_refgen_bwa_input_format($refgen_infile, $refgen_output_dir);

# Creates the BWA reference genome index file
bwa_index($refgen_fasta_outfile, $refgen_output_dir);

# Create the BWA alignment file output directory if it doesn't already exist.
my $bwa_output_dir = join('/', $output_dir, "BWA_ALIGNMENT_FILES");
unless(-d $bwa_output_dir){
    mkdir($bwa_output_dir, 0777) or die "Can't make directory: $!";
}

# Find all files in the processed GBS fastq file directory with the extension *.fastq or *.fastq.qz.
my ($gbs_fastq_files, $gbs_fastq_file_count);
if($gbs_fastq_file_type eq "gzfastq"){
    ($gbs_fastq_files, $gbs_fastq_file_count) = find_files($gbs_fastq_dir, "fastq.gz");
}elsif($gbs_fastq_file_type eq "fastq"){
    ($gbs_fastq_files, $gbs_fastq_file_count) = find_files($gbs_fastq_dir, "fastq");
}else{
    die "Error: $gbs_fastq_file_type is not a recognized input file type. Use the -t option and specify gzfastq for *.fastq.gz files or fastq for *.fastq files.";
}

die "Error: The GBS fastq file count was $gbs_fastq_file_count. The input file type was specified as $gbs_fastq_file_type. Either you specified the wrong directory or the files in $gbs_fastq_dir are a different format than $gbs_fastq_file_type. Use the -t option and specify gzfastq for *.fastq.gz files or fastq for *fastq files." unless($gbs_fastq_file_count > 0);

# Iterate through each GBS fastq file and perform an alignment to the reference genome using the BWA alignment program.
foreach my $file_name (sort keys %{$gbs_fastq_files}){
    warn "Processing " . $file_name . ".....\n";
    my $gbs_fastq_infile = $gbs_fastq_files->{$file_name};
    
    # Get the fastq file directory and file name without the extension.
    my ($fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/_trimmed_offset_\d+\.fastq/);
    
    # Process fastq files for BWA alignments or skip if sam files were previously generated.
    my $sam_outfile = join('/', $bwa_output_dir, $fastq_filename . ".sam");
    unless(-s $sam_outfile){
        # If the fastq file is compressed, uncompress the file and set the resulting fastq filename to be the fastq infile.
        if($gbs_fastq_file_type eq "gzfastq"){
            my $uncompressed_fastq_file = gunzip_fastq_file($gbs_fastq_infile);
            $gbs_fastq_infile = $uncompressed_fastq_file;
        }elsif(($gbs_fastq_file_type eq "fastq") and ($file_name =~ m/trimmed_offset_\d+/)){ # If the fastq file is not compressed set the resulting fastq filename to be the fastq infile.
            my $fastq_infile = join('/', $fastq_dir, $fastq_filename . ".fastq");
            warn "Copying $gbs_fastq_infile to $fastq_infile.....";
            copy($gbs_fastq_infile, $fastq_infile) or die "Copy failed: $!";
            $gbs_fastq_infile = $fastq_infile;
        }
        
        if($bwa_algorithm eq "samse"){
            # Creates the BWA alignment file using the GBS fastq input file to align the fastq sequence reads to the reference genome.
            my $bwa_alignment_outfile = bwa_aln($refgen_fasta_outfile, $gbs_fastq_infile, $num_threads, $bwa_output_dir);
            
            # Creates the BWA single-ended alignment file in sam format using the GBS fastq input file to align the fastq sequence reads to the reference genome and the BWA aln format file.
            my $bwa_aligned_master_outfile = bwa_samse($refgen_fasta_outfile, $gbs_fastq_infile, $bwa_alignment_outfile, $sam_mapq_threshold, $bwa_output_dir);
            
        }elsif($bwa_algorithm eq "mem"){
            
            # Creates the BWA single-ended alignment file in sam format using the GBS fastq input file to align the fastq sequence reads to the reference genome.
            my $bwa_alignment_master_outfile = bwa_mem($refgen_fasta_outfile, $gbs_fastq_infile, $num_threads, $bwa_output_dir);
            
        }
        
        # Remove the uncompressed fastq file to save space as we do not need file after this point.
        unlink($gbs_fastq_infile) or die "Could not unlink $gbs_fastq_infile: $!" if(($gbs_fastq_file_type eq "gzfastq") and ($gbs_fastq_infile =~ m/\.fastq$/));
        unlink($gbs_fastq_infile) or die "Could not unlink $gbs_fastq_infile: $!" if(($gbs_fastq_file_type eq "fastq") and ($file_name =~ m/trimmed_offset_\d+/));
    }
}

# Create the padded sam alignment file output directory if it doesn't already exist.
my $padded_sam_output_dir = join('/', $output_dir, "PADDED_SAM_OUTFILES");
unless(-d $padded_sam_output_dir){
    mkdir($padded_sam_output_dir, 0777) or die "Can't make directory: $!";
}


# Create the padded bam alignment file output directory if it doesn't already exist.
my $padded_bam_output_dir = join('/', $output_dir, "PADDED_BAM_OUTFILES");
if($align_outfile_type eq "bam"){
    unless(-d $padded_bam_output_dir){
        mkdir($padded_bam_output_dir, 0777) or die "Can't make directory: $!";
    }
}

# Find all the BWA sam alignment output files from the BWA alignment directory with the extension *.sam.
my ($bwa_align_files, $bwa_align_file_count) = find_files($bwa_output_dir, "sam");

# Iterate through each BWA sam alignment output file with extension *.sam and padd sequences with poly-Ns. If the sequence is aligned to the reference genome in the sense strand orientation padd the sequence on the 3' end. Otherwise the sequence is aligned to the reference genome in the antisense strand orientation padd the sequence on the 5' end.
foreach my $file_name (sort keys %{$bwa_align_files}){
    warn "Processing " . $file_name . ".....\n";
    my $bwa_align_infile = $bwa_align_files->{$file_name};
    my $padded_sam_outfile = bwa_pad_sam_files($bwa_align_infile, $gbs_sequence_length, $sam_mapq_threshold, $padded_sam_output_dir);
    
    my $filename_prefix = fileparse($padded_sam_outfile, qr/\.sam/);
    my $bwa_sai_file = join("/", $bwa_output_dir, "$filename_prefix.sai");
    if((-s $padded_sam_outfile) and ($keep_sam_output_files eq "false")){
        
        # Remove the BWA sam file to save space as we do not need the file after this point.
        if(-s $bwa_align_infile){
            unlink($bwa_align_infile) or die "Could not unlink $bwa_align_infile: $!";
        }
        # Remove the BWA sai file to save space as we do not need the file after this point.
        if(-s $bwa_sai_file){
            unlink($bwa_sai_file) or die "Could not unlink $bwa_sai_file: $!";
        }
    }
    
    if($align_outfile_type eq "bam"){
        
        my $padded_bam_outfile = join("/", $padded_bam_output_dir, "$filename_prefix.bam");
        my $padded_sam_tempfile = join("/", $padded_bam_output_dir, "$filename_prefix-temp.sam");
        unless(-s $padded_bam_outfile){
            my $samtoolsCmd1  = "$samtools sort -@ $num_threads $padded_sam_outfile > $padded_sam_tempfile";
            warn $samtoolsCmd1 . "\n\n";
            system($samtoolsCmd1) == 0 or die "Error calling $samtoolsCmd1: $?";
            
            my $samtoolsCmd2  = "$samtools view -h -b -@ $num_threads $padded_sam_tempfile > $padded_bam_outfile";
            warn $samtoolsCmd2 . "\n\n";
            system($samtoolsCmd2) == 0 or die "Error calling $samtoolsCmd2: $?";
        }
        
        if(-s $padded_bam_outfile){
            # Remove the padded sam file to save space as we do not need the file after this point.
            if(-s $padded_sam_outfile){
                unlink($padded_sam_outfile) or die "Could not unlink $padded_sam_outfile: $!" if($keep_sam_output_files eq "false");
            }
        }
        
        if(-s $padded_sam_tempfile){
            # Remove the padded sam temp file to save space as we do not need the file after this point.
            unlink($padded_sam_tempfile) or die "Could not unlink $padded_sam_tempfile: $!"
        }
    }
    
}
# Remove the BWA sam file directory as it is empty and we do not need the file after this point.
rmdir($bwa_output_dir) or die "Could not unlink $bwa_output_dir: $!" if($keep_sam_output_files eq "false");

# Remove the padded sam file directory as it is empty and we do not need the file after this point.
rmdir($padded_sam_output_dir) or die "Could not unlink $padded_sam_output_dir: $!" if($keep_sam_output_files eq "false");

# convert_refgen_bwa_input_format($renum_refgen_fasta_infile, $refgen_output_dir) - Converts the reference genome to BWA input format by renumerating the fasta headers and generates a table of contents file referencing the sequence headers to the new BWA input format sequence headers.
#
# Input paramater(s):
#
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
#
# $refgen_output_dir - The reference genome output directory.
sub convert_refgen_bwa_input_format{
    
    # The reference genome input fasta file.
    my $refgen_fasta_infile = shift;
    die "Error lost the reference genome input fasta file." unless defined $refgen_fasta_infile;
    
    # The reference genome output directory.
    my $refgen_output_dir = shift;
    die "Error lost reference genome output directory" unless defined $refgen_output_dir;
    
    # Format the renumerated fasta and table of contents .toc files.
    my $refgen_fasta_filename = fileparse($refgen_fasta_infile);
    my $refgen_fasta_outfile = join('/', $refgen_output_dir, $refgen_fasta_filename . ".fasta");
    my $refgen_toc_outfile = join('/', $refgen_output_dir, $refgen_fasta_filename . ".toc.txt");
    
    # Generate the renumerated fasta and table of contents .toc files if the files are not already generated.
    unless(-s $refgen_fasta_outfile and -s $refgen_toc_outfile){
        warn "Converting $refgen_fasta_filename to BWA input format....\n\n";
        my $seq_counter = 0;
        my %fasta_header = ();
        my %fasta_seqs = ();
        
        # Parse the contents of the reference genome fasta file to filter.
        open(INFILE, "<$refgen_fasta_infile") or die "Couldn't open file $refgen_fasta_infile for reading, $!";
        while(<INFILE>){
            chomp $_;
            
            if($_ !~ /^$/){
                if($_ =~ /^>/){
                    # Parse stacks reference genome fasta input file for the reference genome.
                    $seq_counter++;
                    my $refgen_seq_id = substr($_, 1, length $_);
                    $fasta_header{$seq_counter}{"HEADER"} = $refgen_seq_id;
                    print $fasta_header{$seq_counter}{"HEADER"} . "\n";
                }elsif($_ =~ m/^[ACGTNRYSWKM]+/i){
                    my $sequence = $_;
                    push(@{$fasta_seqs{$seq_counter}{"SEQUENCE"}}, $sequence);
                }else{
                    die "Error: $refgen_fasta_infile is not in the correct format!\n$_";
                }
            }
        }
        close(INFILE) or die "Couldn't close file $refgen_fasta_infile";
        
        my $num_digits = length($seq_counter - 1);
        open(REFGEN_FASTA_OUTFILE, ">$refgen_fasta_outfile") or die "Couldn't open file $refgen_fasta_outfile for writting, $!";
        open(REFGEN_TOC_OUTFILE, ">$refgen_toc_outfile") or die "Couldn't open file $refgen_toc_outfile for writting, $!";
        print REFGEN_TOC_OUTFILE join("\t", "bwa_input_fasta_header", "original_fasta_header") . "\n";
        foreach my $seq_num (sort {$a <=> $b} keys %fasta_seqs){
            #warn "Printing $seq_num";
            my $fasta_header = $fasta_header{$seq_num}{"HEADER"};
            
            my $sequence = join("", @{$fasta_seqs{$seq_num}{"SEQUENCE"}});
            my $padded_zeros_length = ($num_digits - length($seq_num));
            my $padded_seq_num = '0' x $padded_zeros_length . $seq_num;
            
            # 			warn join("\t", $padded_seq_num, $fasta_header) . "\n";
            print REFGEN_FASTA_OUTFILE join("\n", ">$padded_seq_num", $sequence) . "\n";
            print REFGEN_TOC_OUTFILE join("\t", $padded_seq_num, $fasta_header) . "\n";
            
        }
        close(REFGEN_FASTA_OUTFILE) or die "Couldn't close file $refgen_fasta_outfile";
        close(REFGEN_TOC_OUTFILE) or die "Couldn't close file $refgen_toc_outfile";
        
    }
    
    # Need to get the actual number of chromosomes or scaffolds for the -eC ($num_chromosomes option).
    open(REFGEN_TOC_INFILE, "<$refgen_toc_outfile") or die "Couldn't open file $refgen_toc_outfile for reading, $!";
    my $i = 0;
    my $num_chromosomes = 0;
    while(<REFGEN_TOC_INFILE>){
        chomp $_;
        # 		warn $_ . "\n";
        if($i ne 0){
            $num_chromosomes++;
        }
        $i++;
    }
    close(REFGEN_TOC_INFILE) or die "Couldn't close file $refgen_toc_outfile";
    
    return ($refgen_fasta_outfile, $refgen_toc_outfile, $num_chromosomes);
}


# bwa_index($renum_refgen_fasta_infile, $refgen_output_dir) - Executes the BWA alignment program using the index option to generate BWA index .amb .ann .bwt .pac .sa files from a renumerated reference genome file.
#
# Input paramater(s):
#
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
#
# $refgen_output_dir - The reference genome output directory that contains the .amb .ann .bwt .pac .sa index files.

#bwa index -a bwtsw ./refgen/renumbered_Msequence.fasta
sub bwa_index{
    
    # The renumerated reference genome input fasta file.
    my $renum_refgen_fasta_infile = shift;
    die "Error lost the renumbered reference genome input fasta file" unless defined $renum_refgen_fasta_infile;
    
    # The reference genome output directory that contains the .amb .ann .bwt .pac .sa index files.
    my $refgen_output_dir = shift;
    die "Error lost reference genome output directory" unless defined $refgen_output_dir;
    
    # Format the index file into .amb .ann .bwt .pac .sa files.
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

# $bwa_alignment_outfile = bwa_samse($renum_refgen_fasta_infile, $gbs_fastq_infile, $num_threads, $bwa_output_dir) - Executes the BWA alignment program using the aln option to generate BWA sai files from a renumerated reference genome and quality filtered and trimmed adapter GBS fastq input files.
#
# Input paramater(s):
#
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
#
# $gbs_fastq_infile - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
#
# $num_threads - The number of threads to use for BWA.
#
# $bwa_output_dir - The BWA alignment output directory that contains the sai alignment files.
#
# Output paramater(s):
#
# $bwa_alignment_outfile - The BWA sai alignment output file.

#bwa aln -t 4 ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/mpbGBSTags.cnt.fq > ./mergedTagCounts/AlignedGBSTags1.sai
sub bwa_aln{
    
    # The renumerated reference genome input fasta file.
    my $renum_refgen_fasta_infile = shift;
    die "Error lost the renumbered reference genome input fasta file" unless defined $renum_refgen_fasta_infile;
    
    # The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
    my $gbs_fastq_infile = shift;
    die "Error lost the GBS fastq input file" unless defined $gbs_fastq_infile;
    
    # The number of threads to use for BWA.
    my $num_threads = shift;
    die "Error lost the number of cores for BWA" unless defined $num_threads;
    
    # The BWA alignment output directory that contains the sai alignment files.
    my $bwa_output_dir = shift;
    die "Error lost the bwa output directory" unless defined $bwa_output_dir;
    
    my $individual_id = "";
    if($gbs_fastq_infile =~ m/trimmed_offset_\d+\.fastq/){
        # Get the basename of the fastq filename without the _trimmed_offset_\d+\.fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/_trimmed_offset_\d+\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }else{
        # Get the basename of the fastq filename without the .fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }
    
    # Format the sai file.
    my $bwa_alignment_outfile = join('/', $bwa_output_dir, $individual_id . ".sai");
    
    # Execute the BWA alignment program if the SAI file is not already generated.
    unless(-s $bwa_alignment_outfile){
        warn "Generating the bwa alignment file.....\n\n";
        my $bwaAlnCmd  = "$bwa aln -t $num_threads $renum_refgen_fasta_infile $gbs_fastq_infile > $bwa_alignment_outfile";
        warn $bwaAlnCmd . "\n\n";
        system($bwaAlnCmd) == 0 or die "Error calling $bwaAlnCmd: $?";
    }
    
    # Return the BWA sai alignment output file.
    return $bwa_alignment_outfile;
}

# $bwa_aligned_master_outfile = bwa_samse($renum_refgen_fasta_infile, $gbs_fastq_infile, $bwa_alignment_infile, $sam_mapq_threshold, $bwa_output_dir) - Executes the BWA alignment program using the samse option to generate sam alignment files from a renumerated reference genome, BWA sai, and quality filtered and trimmed adapter GBS fastq input files.
#
# Input paramater(s):
#
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
#
# $gbs_fastq_infile - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
#
# $bwa_alignment_infile - The BWA alignment index input file.
#
# $sam_mapq_threshold - The sam alignment file mapping quality score threshold.
#
# $bwa_output_dir - The BWA alignment output directory that contains the sam alignment files.
#
# Output paramater(s):
#
# $bwa_aligned_master_outfile - The BWA sam alignment output file.

#bwa samse ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/AlignedGBSTags1.sai ./mergedTagCounts/mpbGBSTags.cnt.fq > mergedTagCounts/AlignedMasterTagsMPB.sam
sub bwa_samse{
    
    # The renumerated reference genome input fasta file.
    my $renum_refgen_fasta_infile = shift;
    die "Error lost the renumbered reference genome input file" unless defined $renum_refgen_fasta_infile;
    
    # The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
    my $gbs_fastq_infile = shift;
    die "Error lost the GBS fastq input file" unless defined $gbs_fastq_infile;
    
    # The BWA alignment index input file.
    my $bwa_alignment_infile = shift;
    die "Error lost the BWA SAI formatted alignment (.sai) file" unless defined $bwa_alignment_infile;
    
    # The sam alignment file mapping quality score threshold.
    my $sam_mapq_threshold = shift;
    die "Error lost the sam alignment file mapping quality score threshold" unless defined $sam_mapq_threshold;
    
    # The BWA alignment output directory that contains the sam alignment files.
    my $bwa_output_dir = shift;
    die "Error lost the bwa output directory" unless defined $bwa_output_dir;
    
    my $individual_id = "";
    if($gbs_fastq_infile =~ m/trimmed_offset_\d+\.fastq/){
        # Get the basename of the fastq filename without the _trimmed_offset_\d+\.fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/_trimmed_offset_\d+\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }else{
        # Get the basename of the fastq filename without the .fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }
    
    # Execute the BWA alignment program if the sam alignment file is not already generated.
    my $bwa_aligned_master_outfile = join('/', $bwa_output_dir, $individual_id . ".sam");
    unless(-s $bwa_aligned_master_outfile){
        warn "Generating the bwa sam file.....\n\n";
        my $bwaSamseCmd  = "$bwa samse $renum_refgen_fasta_infile $bwa_alignment_infile $gbs_fastq_infile";
        warn $bwaSamseCmd . "\n\n";
        
        open(SAM_OUTFILE, ">$bwa_aligned_master_outfile") or die "Couldn't open file $bwa_aligned_master_outfile for writting, $!";
        local (*BWA_OUT, *BWA_IN);
        my $pid = open2(\*BWA_OUT,\*BWA_IN, $bwaSamseCmd) or die "Error calling open2: $!";
        close BWA_IN or die "Error closing STDIN to bwa process: $!";
        while(<BWA_OUT>){
            chomp $_;
            if($_ =~ m/^\@HD|^\@SQ|^\@RG|^\@PG|^\@CO/){ # Parse header lines starting with @HD, @SQ, @RG, @PG, or @CO.
                print SAM_OUTFILE $_ . "\n";
            }else{ #8_1308_8038_19954_1LL-06_TTCTG_5_A3length=66    16      7376    784845  37      66M     *       0       0       AATAATGCGCCAGCCAAGAGCTGTTTGGTAGCATTCTGCTGACGCTCGCACTCCCGTACACCTGCA      DDCB@>>>AFEFFFHC@GHGGGEGEGGIGCDEHBGGCIGIGGIIHGHEFCE:F<GIGIHFHHDFDD       XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:66
                
                my @split_sam_entry = split(/\t/, $_);
                my ($fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext,
                $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, @optional_fields) = @split_sam_entry;
                
                my $optional_fields = join("\t", @optional_fields);
                
                # die join("\t", $fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext, $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, $optional_fields);
                
                # Filter alignment entry if the read is unmapped in the sense and antisense orientation.
                next if(($bam_bitwise_flag eq 4) or ($bam_bitwise_flag eq 20));
                
                # Filter alignment entry if the mapping quality of the read is below this threshold.
                #next if($mapq < $sam_mapq_threshold);
                
                # Filter alignment entry if the cigar string contains an insertion/deletion in the alignment.
                # next if($cigar !~ m/^\d+M$/);
                
                # Filter alignment entry if the read is not uniquely mapped to the reference.
                #next if($optional_fields !~ m/XT:A:U/);
                
                print SAM_OUTFILE $_ . "\n";
            }
            
        }
        close BWA_OUT or die "Error closing STDOUT from bwa process: $!";
        wait;
        close(SAM_OUTFILE) or die "Couldn't close file $bwa_aligned_master_outfile";
    }
    
    # The BWA sam alignment output file.
    return $bwa_aligned_master_outfile;
}

# $bwa_alignment_master_outfile = bwa_mem($renum_refgen_fasta_infile, $gbs_fastq_infile, $num_threads, $bwa_output_dir) - Executes the BWA alignment program using the aln option to generate BWA sai files from$
#
# Input paramater(s):
#
# $renum_refgen_fasta_infile - The renumerated reference genome input fasta file.
#
# $gbs_fastq_infile - The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
#
# $num_threads - The number of threads to use for BWA.
#
# $bwa_output_dir - The BWA alignment output directory that contains the sai alignment files.
#
# Output paramater(s):
#
# $bwa_alignment_outfile - The BWA sai alignment output file.

#bwa mem -t 4 ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/mpbGBSTags.cnt.fq > ./mergedTagCounts/AlignedGBSTags1.sam
sub bwa_mem{
    
    # The renumerated reference genome input fasta file.
    my $renum_refgen_fasta_infile = shift;
    die "Error lost the renumbered reference genome input fasta file" unless defined $renum_refgen_fasta_infile;
    
    # The quality filtered and adapter trimmed GBS fastq input file for an individual within the Genotyping by Sequencing (GBS) project.
    my $gbs_fastq_infile = shift;
    die "Error lost the GBS fastq input file" unless defined $gbs_fastq_infile;
    
    # The number of threads to use for BWA.
    my $num_threads = shift;
    die "Error lost the number of cores for BWA" unless defined $num_threads;
    
    # The BWA alignment output directory that contains the sai alignment files.
    my $bwa_output_dir = shift;
    die "Error lost the bwa output directory" unless defined $bwa_output_dir;
    
    my $individual_id = "";
    if($gbs_fastq_infile =~ m/\.fastq/){
        # Get the basename of the fastq filename without the .fastq extension.
        my ($gbs_fastq_filename, $fastq_dir) = fileparse($gbs_fastq_infile, qr/\.fastq/);
        $individual_id = $gbs_fastq_filename;
    }
    
   	# Execute the BWA alignment program if the sam alignment file is not already generated.
    my $bwa_aligned_master_outfile = join('/', $bwa_output_dir, $individual_id . ".sam");
    unless(-s $bwa_aligned_master_outfile){
        warn "Generating the bwa sam file.....\n\n";
        my $bwaMemSamseCmd  = "$bwa mem -M -t $num_threads $renum_refgen_fasta_infile $gbs_fastq_infile";
        warn $bwaMemSamseCmd . "\n\n";
        
        open(SAM_OUTFILE, ">$bwa_aligned_master_outfile") or die "Couldn't open file $bwa_aligned_master_outfile for writting, $!";
        local (*BWA_OUT, *BWA_IN);
        my $pid = open2(\*BWA_OUT,\*BWA_IN, $bwaMemSamseCmd) or die "Error calling open2: $!";
        close BWA_IN or die "Error closing STDIN to bwa process: $!";
        while(<BWA_OUT>){
            chomp $_;
            if($_ =~ m/^\@HD|^\@SQ|^\@RG|^\@PG|^\@CO/){ # Parse header lines starting with @HD, @SQ, @RG, @PG, or @CO.
                print SAM_OUTFILE $_ . "\n";
            }else{ #8_1308_8038_19954_1LL-06_TTCTG_5_A3length=66    16	7376    784845  37	66M     *	0	0	AATAATGCGCCAGCCAAGAGCTGTTTGGTAGCATTCTGCTGACGCTCGCACTCCCGTAC$
                
                my @split_sam_entry = split(/\t/, $_);
                my ($fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext,
                $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, @optional_fields) = @split_sam_entry;
                
                my $optional_fields = join("\t", @optional_fields);
                
                # die join("\t", $fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext, $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, $optional_fields);
                
                # Filter alignment entry if the read is unmapped in the sense and antisense orientation.
                next if(($bam_bitwise_flag eq 4) or ($bam_bitwise_flag eq 20));
                
                # Filter alignment entry if the mapping quality of the read is below this threshold.
                #next if($mapq < $sam_mapq_threshold);
                
                # Filter alignment entry if the cigar string contains an insertion/deletion in the alignment.
                # next if($cigar !~ m/^\d+M$/);
                
                # Filter alignment entry if the read is not uniquely mapped to the reference.
                #next if($optional_fields =~ m/XA:Z:\d+/);
                #next if($optional_fields =~ m/SA:Z:\d+/);
                print SAM_OUTFILE $_ . "\n";
            }
            
        }
        close BWA_OUT or die "Error closing STDOUT from bwa process: $!";
        wait;
        close(SAM_OUTFILE) or die "Couldn't close file $bwa_aligned_master_outfile";
    }
    
    # The BWA sam alignment output file.
    return $bwa_aligned_master_outfile;
}


# $padded_sam_file = bwa_pad_sam_files($sam_infile, $gbs_sequence_length, $sam_mapq_threshold, $padded_sam_output_dir) - Pads the BWA sam alignment files with poly-Ns so that all sequences are of uniform length.
#
# Input paramater(s):
#
# $sam_infile - The unpadded BWA sam alignment input file.
#
# $gbs_sequence_length - The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences.
#
# $sam_mapq_threshold - The sam alignment file mapping quality score threshold.
#
# $padded_sam_output_dir - The output directory that contains all the BWA sam alignment files.
#
# Output parameter(s):
#
# $padded_sam_file - The padded BWA sam alignment output file.
sub bwa_pad_sam_files{
    
    # The unpadded BWA sam alignment input file.
    my $sam_infile = shift;
    die "Error lost the sam input file" unless defined $sam_infile;
    
    # The GBS fastq sequence length in base pairs (bps) common to all GBS fastq sequences.
    my $gbs_sequence_length = shift;
    die "Error lost the GBS fastq sequence length in base pairs (bps)" unless defined $gbs_sequence_length;
    
    # The sam alignment file mapping quality score threshold.
    my $sam_mapq_threshold = shift;
    die "Error lost the sam alignment file mapping quality score threshold" unless defined $sam_mapq_threshold;
    
    # The output directory that contains all the BWA sam alignment files.
    my $padded_sam_output_dir = shift;
    die "Error lost the padded sam output file directory" unless defined $padded_sam_output_dir;
    
    # Get the basename of the sam filename without the .sam extension.
    my $sam_filename = fileparse($sam_infile, qr/\.sam/);
    
    # Generate the padded sam output file if the file is not already generated.
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
                my ($fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext,
                $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, @optional_fields) = @split_sam_entry;
                
                my $optional_fields = join("\t", @optional_fields);
                
                my $fastq_sequence_length = length($fastq_sequence);
                my $fastq_quality_scores_length = length($fastq_quality_scores);
                # 				die join("\t", $fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $cigar, $rnext, $pnext, $tlen, $fastq_sequence, $fastq_quality_scores, $optional_fields);
                
                # Filter alignment entry if the read is unmapped in the sense and antisense orientation.
                next if(($bam_bitwise_flag eq 4) or ($bam_bitwise_flag eq 20));
                
                # Filter alignment entry if the mapping quality of the read is below this threshold.
                next if($mapq < $sam_mapq_threshold);
                
                # Filter alignment entry if the cigar string contains an insertion/deletion in the alignment.
                next if($cigar !~ m/^\d+M$/);
                
                # Filter alignment entry if the read is not uniquely mapped to the reference.
                if($bwa_algorithm eq "samse"){
                    next if($optional_fields !~ m/XT:A:U/);
                }
                
                if($bwa_algorithm eq "mem"){
                    next if($optional_fields =~ m/XA:Z:\d+/);
                    next if($optional_fields =~ m/SA:Z:\d+/);
                }
                
                if(($fastq_sequence_length ne $gbs_sequence_length) and ($fastq_quality_scores_length ne $gbs_sequence_length)){
                    
                    # Get the number of poly-Ns to pad.
                    my $padded_nucleotide_length = ($gbs_sequence_length - $fastq_sequence_length);
                    my $padded_nucleotide_seq = 'N' x $padded_nucleotide_length;
                    
                    # Pad sequence based on alignment type.
                    my $padded_fastq_sequence = "";
                    
                    # Pad sequence if aligned to sense strand.
                    $padded_fastq_sequence = join("", $fastq_sequence, $padded_nucleotide_seq) if($bam_bitwise_flag eq 0);
                    
                    # Pad sequence if aligned to antisense strand.
                    $padded_fastq_sequence = join("", $padded_nucleotide_seq, $fastq_sequence) if($bam_bitwise_flag eq 16);
                    
                    # Get the number of quality scores to pad.
                    my $padded_score_length =  ($gbs_sequence_length - $fastq_quality_scores_length);
                    my $padded_score_seq = '#' x $padded_score_length;
                    
                    # Pad quality scores based on alignment type.
                    my $padded_fastq_quality_scores = "";
                    
                    # Pad quality scores if aligned to sense strand sequence in sense orientation.
                    $padded_fastq_quality_scores = join("", $fastq_quality_scores, $padded_score_seq) if($bam_bitwise_flag eq 0);
                    
                    # Pad quality scores if aligned to antisense strand in antisense orientation.
                    $padded_fastq_quality_scores = join("", $padded_score_seq, $fastq_quality_scores) if($bam_bitwise_flag eq 16);
                    
                    # Get the number of poly-Ns to pad.
                    #                    my $extended_cigar_string = $padded_nucleotide_length . "P";
                    
                    # Pad quality scores based on alignment type.
                    #                    my $padded_extended_cigar_string = "";
                    my $padded_extended_cigar_string = $gbs_sequence_length . "M";
                    # Pad quality scores if aligned to sense strand sequence in sense orientation.
                    #                    $padded_extended_cigar_string = join("", $cigar, $extended_cigar_string) if($bam_bitwise_flag eq 0);
                    
                    # Pad quality scores if aligned to antisense strand in antisense orientation.
                    #                    $padded_extended_cigar_string = join("", $extended_cigar_string, $cigar) if($bam_bitwise_flag eq 16);
                    
                    my $padded_fastq_sequence_length = length($padded_fastq_sequence);
                    my $padded_fastq_quality_scores_length = length($padded_fastq_quality_scores);
                    
                    die "Error: $fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length bp ne gbs_sequence_length=$gbs_sequence_length bp" if($padded_fastq_sequence_length ne $gbs_sequence_length);
                    die "Error: $fastq_header: padded_fastq_sequence_length=$padded_fastq_sequence_length ne padded_fastq_quality_scores_length=$padded_fastq_quality_scores_length" if($padded_fastq_sequence_length ne $padded_fastq_quality_scores_length);
                    
                    #print PADDED_SAM_OUTFILE join("\t", $fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, join("", $gbs_sequence_length, "M"), $rnext, $pnext, $tlen, $padded_fastq_sequence, $padded_fastq_quality_scores, $optional_fields) . "\n";
                    print PADDED_SAM_OUTFILE join("\t", $fastq_header, $bam_bitwise_flag, $rname, $r_pos, $mapq, $padded_extended_cigar_string, $rnext, $pnext, $tlen, $padded_fastq_sequence, $padded_fastq_quality_scores, $optional_fields) . "\n";
                    
                }elsif(($fastq_sequence_length eq $gbs_sequence_length) and ($fastq_quality_scores_length eq $gbs_sequence_length)){
                    
                    print PADDED_SAM_OUTFILE $_ . "\n";
                }
                # 			die $_;
            }
        }
        close(SAM_INFILE) or die "Couldn't close file $sam_infile";
        close(PADDED_SAM_OUTFILE) or die "Couldn't close file $padded_sam_outfile";
    }
    
    # The padded BWA sam alignment output file.
    return $padded_sam_outfile;
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
    
    my ($fastq_filename, $fastq_dir, $individual_id) = "";
    if($fastq_file =~ m/trimmed_offset_\d+\.fastq\.gz/){
        # get the basename of the fastq filename without the _trimmed_offset_\d+\.fastq extension.
        ($fastq_filename, $fastq_dir) = fileparse($fastq_file, qr/_trimmed_offset_\d+\.fastq\.gz/);
        $individual_id = $fastq_filename;
    }else{
        # get the basename of the fastq filename without the \.fastq extension.
        ($fastq_filename, $fastq_dir) = fileparse($fastq_file, qr/\.fastq/);
        $individual_id = $fastq_filename;
    }         	
    my $uncompressed_fastq_file = join('/', $fastq_dir, $individual_id . ".fastq");
    
    unless(-s $uncompressed_fastq_file){
        warn "calling gunzip for $fastq_file....\n";
        my $gunzipCmd  = "$gunzip -c $fastq_file > $uncompressed_fastq_file";
        warn $gunzipCmd . "\n\n";
        system($gunzipCmd) == 0 or die "Error calling $gunzipCmd: $?";
    }
    return $uncompressed_fastq_file;
}
