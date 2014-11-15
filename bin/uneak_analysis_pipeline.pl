#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Copy;
use File::Basename;

# perl uneak_analysis_pipeline.pl -i ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/GQ03122013_5_fastq.txt.gz -n MPB_MALE_UNEAK_GBS -k ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/mpb_barcodes.txt -o ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/UNEAK_MPB_GBS_ANALYSIS

my ($fastq_infile, $project_name, $fastq_barcodes_infile, $restriction_enzymes, $max_num_barcode_reads, $min_num_present_ufastq_to_tag, $min_num_present_merge_multiple_tag, 
$error_rate, $min_minor_allele_freq_hap_map, $max_minor_allele_freq_hap_map, $min_call_rate, $max_call_rate, $tassel_num_ram, $uneak_project_dir);

GetOptions(
	'i=s'    => \$fastq_infile,
	'n=s'    => \$project_name,
	'k=s'    => \$fastq_barcodes_infile,
	'e=s'    => \$restriction_enzymes,
	's=s'    => \$max_num_barcode_reads,
	'c=s'    => \$min_num_present_ufastq_to_tag,
	'u=s'    => \$min_num_present_merge_multiple_tag,
	'r=s'    => \$error_rate,
	'f=s'    => \$min_minor_allele_freq_hap_map,
	'a=s'    => \$max_minor_allele_freq_hap_map,
	'y=s'    => \$min_call_rate,
	'x=s'    => \$max_call_rate,
	'm=s'    => \$tassel_num_ram,
	'o=s'    => \$uneak_project_dir,
);

usage() unless (
	defined $fastq_infile
	and defined $project_name
 	and defined $fastq_barcodes_infile
	and defined $uneak_project_dir
);

$restriction_enzymes = 'PstI-MspI' unless defined $restriction_enzymes;
$max_num_barcode_reads = 250000000 unless defined $max_num_barcode_reads;
$min_num_present_ufastq_to_tag = 1 unless defined $min_num_present_ufastq_to_tag;
$min_num_present_merge_multiple_tag = 5 unless defined $min_num_present_merge_multiple_tag;
$error_rate = 0.03 unless defined $error_rate;
$min_minor_allele_freq_hap_map = 0.00 unless defined $min_minor_allele_freq_hap_map;
$max_minor_allele_freq_hap_map = 1.00 unless defined $max_minor_allele_freq_hap_map;
$min_call_rate = 0 unless defined $min_call_rate;
$max_call_rate = 1 unless defined $max_call_rate;
$tassel_num_ram = 12 unless defined $tassel_num_ram;

my ($run_pipeline, $bwa);
$run_pipeline		= '/home/cookeadmin/tassel-tassel3-standalone/run_pipeline.pl';
$bwa 			= '/usr/local/bin/bwa';

sub usage {

die <<"USAGE";


Full Usage: $0 -i fastq_infile -n project_name -k fastq_barcodes_infile -e restriction_enzymes -s max_num_barcode_reads -c min_num_present_ufastq_to_tag -u min_num_present_merge_multiple_tag -r error_rate -f min_minor_allele_freq_hap_map -a max_minor_allele_freq_hap_map -y min_call_rate -x max_call_rate -m tassel_num_ram -o uneak_project_dir

Minimum Usage: $0 -i fastq_infile -n project_name -k fastq_barcodes_infile -o uneak_project_dir

DESCRIPTION - 

OPTIONS:

-i fastq_infile

-n project_name

-k fastq_barcodes_infile

-e restriction_enzymes

-s max_num_barcode_reads

-c min_num_present_ufastq_to_tag

-u min_num_present_merge_multiple_tag

-r error_rate

-f min_minor_allele_freq_hap_map

-a max_minor_allele_freq_hap_map

-y min_call_rate

-x max_call_rate

-m tassel_num_ram

-o uneak_project_dir

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $uneak_project_dir){
      mkdir($uneak_project_dir, 0777) or die "Can't make directory: $!";
}


# ./Illumina/ (original raw data files, one file per flowcell lane)
# Create the Illumina output directory if it doesn't already exist.
my $illumina_output_dir = join('/', $uneak_project_dir, "Illumina");
unless(-d $illumina_output_dir){
      mkdir($illumina_output_dir, 0777) or die "Can't make directory: $!";
}

# Copying fastq_filename to the Illumina directory 
my $fastq_filename = basename($fastq_infile);
warn "Copying $fastq_filename to the Illumina directory $illumina_output_dir.....\n\n";
my $uneak_fastq_infile = join('/',$illumina_output_dir, $fastq_filename);
unless(-s $uneak_fastq_infile){
	copy($fastq_infile, $uneak_fastq_infile) or die "Copying $fastq_infile to $uneak_fastq_infile failed: $!";
}

# ./key/ (Barcode key file of original raw data files)
# Create the fastq barcodes key output directory if it doesn't already exist.
my $barcodes_key_output_dir = join('/', $uneak_project_dir, "key");
unless(-d $barcodes_key_output_dir){
      mkdir($barcodes_key_output_dir, 0777) or die "Can't make directory: $!";
}

# Copying fastq_barcodes_infile to the barcodes key directory
my $fastq_barcodes_filename = basename($fastq_barcodes_infile);
warn "Copying $fastq_barcodes_filename to the barcodes key directory $barcodes_key_output_dir.....\n\n";
my $uneak_fastq_barcodes_infile = join('/', $barcodes_key_output_dir, $fastq_barcodes_filename);
unless(-s $uneak_fastq_barcodes_infile){
	copy($fastq_barcodes_infile, $uneak_fastq_barcodes_infile) or die "Copying $fastq_barcodes_infile to $uneak_fastq_barcodes_infile failed: $!";
}

# ./tagCounts/ (for output from UQseqToTagCountPlugin OR UFastqToTagCountPlugin OR UMergeTaxaTagCountPlugin)
# Create the tag counts output directory if it doesn't already exist.
my $tag_counts_output_dir = join('/', $uneak_project_dir, "tagCounts");
unless(-d $tag_counts_output_dir){
      mkdir($tag_counts_output_dir, 0777) or die "Can't make directory: $!";
}

ufastq_to_tag_counts($uneak_project_dir, $project_name, $fastq_barcodes_infile, $restriction_enzymes, $max_num_barcode_reads, $min_num_present_ufastq_to_tag, 
 	$tassel_num_ram, $tag_counts_output_dir);
 	
#run binarytotext plugin to make the ufastqtotagcount result from above human readable
#change GQ25032013_5 in the following line to the actual name of your files (should be GQ25032013_yourlanenumber 
# run_pipeline.pl -Xmx4g -fork1 -BinaryToTextPlugin -i ./tagCounts/GQ25032013_7.cnt -o ./tagCounts/GQ25032013_7_cnt.txt -t TagCounts -endPlugin -runfork1 > ./bintotex.log
binary_to_text($uneak_project_dir, $project_name, "TagCounts", $tassel_num_ram, $tag_counts_output_dir);

# ./mergedTagCounts/ (for output fromUMergeTaxaTagCountPlugin)
# Create the merged tag counts output directory if it doesn't already exist.
my $merged_tag_counts_output_dir = join('/', $uneak_project_dir, "mergedTagCounts");
unless(-d $merged_tag_counts_output_dir){
      mkdir($merged_tag_counts_output_dir, 0777) or die "Can't make directory: $!";
}

my $merged_tag_count_outfile = umerge_taxa_tag_counts($uneak_project_dir, $project_name, $min_num_present_merge_multiple_tag, $tassel_num_ram, $merged_tag_counts_output_dir);
#run binarytotext plugin to make the umergetaxatagcount result from above human readable
# run_pipeline.pl -Xmx4g -fork1 -BinaryToTextPlugin -i ./mergedTagCounts/mergeAll.cnt -o ./mergedTagCounts/mergeAll_cnt.txt -t TagCounts -endPlugin -runfork1 >> ./bintotex.log
binary_to_text($uneak_project_dir, $project_name, "TagCounts", $tassel_num_ram, $merged_tag_counts_output_dir);

# ./tagPair/ (for output fromUTagCountToTagPairPlugin)
# Create the tag pair output directory if it doesn't already exist.
my $tag_pair_output_dir = join('/', $uneak_project_dir, "tagPair");
unless(-d $tag_pair_output_dir){
      mkdir($tag_pair_output_dir, 0777) or die "Can't make directory: $!";
}

my $tag_pair_outfile = utag_count_to_tag_pair($uneak_project_dir, $project_name, $error_rate, $tassel_num_ram, $tag_pair_output_dir);

# ./tagsByTaxa/ (for output fromUTBTToMapInfoPlugin)
# Create the tags by taxa output directory if it doesn't already exist.
my $tags_by_taxa_output_dir = join('/', $uneak_project_dir, "tagsByTaxa");
unless(-d $tags_by_taxa_output_dir){
      mkdir($tags_by_taxa_output_dir, 0777) or die "Can't make directory: $!";
}

my $tbt_outfile = utag_pair_to_tbt_counts($uneak_project_dir, $project_name, $tassel_num_ram, $tags_by_taxa_output_dir);
        
# ./mapInfo/ (for output from UTBTToMapInfoPlugin)
# Create the map info output directory if it doesn't already exist.
my $map_info_output_dir = join('/', $uneak_project_dir, "mapInfo");
unless(-d $map_info_output_dir){
      mkdir($map_info_output_dir, 0777) or die "Can't make directory: $!";
} 

# my $tbt_outfile = 
utag_by_taxa_to_map_info($uneak_project_dir, $project_name, $tassel_num_ram, $map_info_output_dir);

# ./hapMap/ (for output from TagsToSNPByAlignmentPlugin)
# Create the hap map output directory if it doesn't already exist.
my $hap_map_output_dir = join('/', $uneak_project_dir, "hapMap");
unless(-d $hap_map_output_dir){
      mkdir($hap_map_output_dir, 0777) or die "Can't make directory: $!";
}

umap_info_to_hap_map($uneak_project_dir, $project_name, $min_minor_allele_freq_hap_map, $max_minor_allele_freq_hap_map, $min_call_rate, 
	$max_call_rate, $tassel_num_ram, $tassel_num_ram, $hap_map_output_dir);
        
#run ufastqtotagcount plugin from the current directory (-w) using PstI-MspI combo (-e), retaining all reads with atleast 1 copy (-c). Expected number of good barcoded reads to be less than 250 million (-s). 
#don't change the -c option setting here since you'll be able to filter by read depth at the umergetaxatagcount step later on.
# run_pipeline.pl -Xmx4g -fork1 -UFastqToTagCountPlugin -w ./ -e PstI-MspI -s 250000000 -c 1 -endPlugin -runfork1 > ./UtagCount.log
sub ufastq_to_tag_counts{

	my $uneak_project_dir = shift;
        die "Error lost the UNEAK project directory" unless defined $uneak_project_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $fastq_barcodes_infile = shift;
        die "Error lost the key file listing barcodes for each sample and plate layout" unless defined $fastq_barcodes_infile;

        my $restriction_enzymes = shift;
        die "Error lost the enzyme(s) used to create the GBS library" unless defined $restriction_enzymes;
        
        my $max_num_barcode_reads = shift;
        die "Error lost the maximum number of good barcoded reads per lane" unless defined $max_num_barcode_reads;
        
        my $min_num_present_tag = shift;
        die "Error lost the minimum number of times a tag must be present to be output" unless defined $min_num_present_tag;
        
        my $tassel_num_ram = shift;
        die "Error lost the amount of tassel memory" unless defined $tassel_num_ram;
        
        my $tag_count_output_dir = shift;
        die "Error lost the output directory to contain output .cnt (tag count) files, one per input FASTQ file" unless defined $tag_count_output_dir;
        
        my $uneak_fastq_to_tag_counts_stdout_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "FastqToTagCounts.stdout.log"));
        my $uneak_fastq_to_tag_counts_stderr_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "FastqToTagCounts.stderr.log"));
        
        my ($tag_count_files, $tag_count_file_counter) = find_files($tag_counts_output_dir, "cnt");
	my $non_zero_tag_count_files = 0;
	foreach my $file_name (sort keys %{$tag_count_files}){
#  		warn $file_name . "\n";
		if(-s $tag_count_files->{$file_name}){
			$non_zero_tag_count_files++;
		}
	}
	
	unless(($non_zero_tag_count_files eq $tag_count_file_counter) and ($tag_count_file_counter ne 0)){
		warn "Generating the UNEAK tag counts list file using the UFastqToTagCountPlugin.....\n";
		my $fastqToTagCountsCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -UFastqToTagCountPlugin -w $uneak_project_dir -e $restriction_enzymes -s $max_num_barcode_reads -c $min_num_present_tag -endPlugin -runfork1 1> $uneak_fastq_to_tag_counts_stdout_log_outfile 2> $uneak_fastq_to_tag_counts_stderr_log_outfile");
		warn $fastqToTagCountsCmd . "\n\n";
		system($fastqToTagCountsCmd) == 0 or die "Error calling $fastqToTagCountsCmd: $?";
	}
}

#run binarytotext plugin to make the ufastqtotagcount result from above human readable
#change GQ25032013_5 in the following line to the actual name of your files (should be GQ25032013_yourlanenumber 
# run_pipeline.pl -Xmx4g -fork1 -BinaryToTextPlugin -i ./tagCounts/GQ25032013_7.cnt -o ./tagCounts/GQ25032013_7_cnt.txt -t TagCounts -endPlugin -runfork1 > ./bintotex.log
sub binary_to_text{

	my $uneak_project_dir = shift;
        die "Error lost the UNEAK project directory" unless defined $uneak_project_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;

        my $file_type = shift;
        die "Error lost the project name" unless defined $file_type;
        
        my $tassel_num_ram = shift;
        die "Error lost the amount of tassel memory" unless defined $tassel_num_ram;
        
        my $output_dir = shift;
        die "Error lost the output directory to contain output .txt (text) files" unless defined $output_dir;

        my ($tag_count_txt_files, $tag_count_txt_file_counter) = find_files($output_dir, "cnt.txt");
        
	my $non_zero_tag_count_txt_files = 0;
	foreach my $file_name (sort keys %{$tag_count_txt_files}){
# 		warn $file_name . "\n";
		if(-s $tag_count_txt_files->{$file_name}){
			$non_zero_tag_count_txt_files++;
		}
	}
	
	unless(($non_zero_tag_count_txt_files eq $tag_count_txt_file_counter) and ($tag_count_txt_file_counter ne 0)){
		# Create the fastq barcodes key output directory if it doesn't already exist.
		my $binary_to_text_log_output_dir = join('/', $uneak_project_dir, "binaryToTextLogFiles");
		unless(-d $binary_to_text_log_output_dir){
			mkdir($binary_to_text_log_output_dir, 0777) or die "Can't make directory: $!";
		}
		
		my ($tag_count_files, $tag_count_file_counter) = find_files($output_dir, "cnt");
		foreach my $file_name (sort keys %{$tag_count_files}){
# 			warn $file_name . "\n";
			my $binary_infile = $tag_count_files->{$file_name};
			
			my $binary_filename = basename($binary_infile);
			my $text_outfile = join('/', $output_dir, $binary_filename . ".txt");

			my $binary_to_text_stdout_log_outfile = join('/', $binary_to_text_log_output_dir, join("_", $project_name, $binary_filename, "BinaryToText.stdout.log"));
			my $binary_to_text_stderr_log_outfile = join('/', $binary_to_text_log_output_dir, join("_", $project_name, $binary_filename, "BinaryToText.stderr.log"));
			
			warn "Generating the UNEAK $file_type text file using the BinaryToTextPlugin on $binary_filename.....\n";
			my $binaryToTextCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -BinaryToTextPlugin -i $binary_infile -o $text_outfile -t $file_type -endPlugin -runfork1 1> $binary_to_text_stdout_log_outfile 2> $binary_to_text_stderr_log_outfile");
			warn $binaryToTextCmd . "\n\n";
			system($binaryToTextCmd) == 0 or die "Error calling $binaryToTextCmd: $?";
		}
	}
}

#run umergetaxatagcounts plugin from the current directory (-w), using a read depth cutoff of 5 (-c). You can change the cutoff to whatever you like.
# run_pipeline.pl -Xmx4g -fork1 -UMergeTaxaTagCountPlugin -w ./ -c 5 -endPlugin -runfork1 > UmergeTagCount.log
sub umerge_taxa_tag_counts{

	my $uneak_project_dir = shift;
        die "Error lost the UNEAK project directory" unless defined $uneak_project_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $min_num_present_merge_multiple_tag = shift;
        die "minimum number of times a tag must be present to be output" unless defined $min_num_present_merge_multiple_tag;
        
        my $tassel_num_ram = shift;
        die "Error lost the amount of tassel memory" unless defined $tassel_num_ram;
        
        my $merged_tag_counts_output_dir = shift;
        die "Error lost the output directory to contain output .cnt (merged tag count) files" unless defined $merged_tag_counts_output_dir;
        
        my $umerge_taxa_tag_counts_stdout_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UmergeTagCount.stdout.log"));
        my $umerge_taxa_tag_counts_stderr_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UmergeTagCount.stderr.log"));
        
        # The merged tag counts output file name.
	my $merged_tag_count_outfile = join('/', $merged_tag_counts_output_dir, "mergedAll.cnt");
	unless(-s $merged_tag_count_outfile){
		warn "Generating the UNEAK merged taxa tag counts list file using the UMergeTaxaTagCountPlugin.....\n";
		my $umergeTaxaTagCountCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -UMergeTaxaTagCountPlugin -w $uneak_project_dir -c $min_num_present_merge_multiple_tag -endPlugin -runfork1 1> $umerge_taxa_tag_counts_stdout_log_outfile 2> $umerge_taxa_tag_counts_stderr_log_outfile");
		warn $umergeTaxaTagCountCmd . "\n\n";
		system($umergeTaxaTagCountCmd) == 0 or die "Error calling $umergeTaxaTagCountCmd: $?";
	}
	return $merged_tag_count_outfile;
}

#run utagcounttotagpair plugin from current directory (-w), with a sequencing error threshold of 3% (this should not be higher than 5% based on known Illumina 
#sequencing error rates.
# run_pipeline.pl -Xmx4g -fork1 -UTagCountToTagPairPlugin -w ./ -e 0.03 -endPlugin -runfork1 > utagtotp.log
sub utag_count_to_tag_pair{

	my $uneak_project_dir = shift;
        die "Error lost the UNEAK project directory" unless defined $uneak_project_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $error_rate = shift;
        die "Error lost error tolerance rate in the network filter. " unless defined $error_rate;
        
        my $tassel_num_ram = shift;
        die "Error lost the amount of tassel memory" unless defined $tassel_num_ram;
        
        my $tag_pair_output_dir = shift;
        die "Error lost the output directory to contain output .tps (tag pair) files" unless defined $tag_pair_output_dir;
        
        my $utag_count_to_tag_pair_stdout_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UTagCountToTagPair.stdout.log"));
        my $utag_count_to_tag_pair_stderr_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UTagCountToTagPair.stderr.log"));
        
        # The tag pair output file name.
	my $tag_pair_outfile = join('/', $tag_pair_output_dir, "tagPair.tps");
	unless(-s $tag_pair_outfile){
		warn "Generating the UNEAK tag count to tag pairs file using the UTagCountToTagPairPlugin.....\n";
		my $utagCountToTagPairCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -UTagCountToTagPairPlugin -w $uneak_project_dir -e $error_rate -endPlugin -runfork1 1> $utag_count_to_tag_pair_stdout_log_outfile 2> $utag_count_to_tag_pair_stderr_log_outfile");
		warn $utagCountToTagPairCmd . "\n\n";
		system($utagCountToTagPairCmd) == 0 or die "Error calling $utagCountToTagPairCmd: $?";
	}
	return $tag_pair_outfile;
}

#run utagpairtoTBTplugin from current directory (-w).
# run_pipeline.pl -Xmx4g -fork1 -UTagPairToTBTPlugin -w ./ -endPlugin -runfork1 > utagtotbt.log
sub utag_pair_to_tbt_counts{

	my $uneak_project_dir = shift;
        die "Error lost the UNEAK project directory" unless defined $uneak_project_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $tassel_num_ram = shift;
        die "Error lost the amount of tassel memory" unless defined $tassel_num_ram;
        
        my $tags_by_taxa_output_dir = shift;
        die "Error lost the output directory to contain output .tbt.bin (tag by taxa) files" unless defined $tags_by_taxa_output_dir;
        
        my $utag_pair_to_tbt_file_stdout_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UTagPairToTBT.stdout.log"));
        my $utag_pair_to_tbt_file_stderr_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UTagPairToTBT.stderr.log"));
        
        # The tag by taxa output file name.
	my $tbt_outfile = join('/', $tags_by_taxa_output_dir, "tbt.bin");
	unless(-s $tbt_outfile){
		warn "Generating the UNEAK tags by taxa file using the UTagPairToTBTPlugin.....\n";
		my $utagPairToTBTCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -UTagPairToTBTPlugin -w $uneak_project_dir -endPlugin -runfork1 1> $utag_pair_to_tbt_file_stdout_log_outfile 2> $utag_pair_to_tbt_file_stderr_log_outfile");
		warn $utagPairToTBTCmd . "\n\n";
		system($utagPairToTBTCmd) == 0 or die "Error calling $utagPairToTBTCmd: $?";
	}
	return $tbt_outfile;
}

#run uTBTtomapinfoplugin from current directory (-w).
# run_pipeline.pl -Xmx4g -fork1 -UTBTToMapInfoPlugin -w ./ -endPlugin -runfork1 > utbttomapinfo.log
sub utag_by_taxa_to_map_info{

	my $uneak_project_dir = shift;
        die "Error lost the UNEAK project directory" unless defined $uneak_project_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
        
        my $tassel_num_ram = shift;
        die "Error lost the amount of tassel memory" unless defined $tassel_num_ram;
        
        my $map_info_output_dir = shift;
        die "Error lost the output directory to contain output .tbt.bin (map info) files" unless defined $map_info_output_dir;
        
        my $utag_by_taxa_to_map_info_stdout_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UTBTToMapInfo.stdout.log"));
        my $utag_by_taxa_to_map_info_stderr_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UTBTToMapInfo.stderr.log"));
        
        # The tag by taxa output file name.
	my $map_info_outfile = join('/', $map_info_output_dir, "mapInfo.bin");
	unless(-s $map_info_outfile){
		warn "Generating the UNEAK map info file using the UTBTToMapInfoPlugin.....\n";
		my $utbtToMapInfoCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -UTBTToMapInfoPlugin -w $uneak_project_dir -endPlugin -runfork1 1> $utag_by_taxa_to_map_info_stdout_log_outfile 2> $utag_by_taxa_to_map_info_stderr_log_outfile");
		warn $utbtToMapInfoCmd . "\n\n";
		system($utbtToMapInfoCmd) == 0 or die "Error calling $utbtToMapInfoCmd: $?";
	}
	return $map_info_outfile;
}

#run umapinfotohapmapplugin from current directory (-w), with a minor allele frequency of 1%, a major allele frequency of 99%, minimum call rate of 0, maximum 
#callrate of 1.
# run_pipeline.pl -Xmx4g -fork1 -UMapInfoToHapMapPlugin -w ./ -mnMAF 0.00 -mxMAF 1.00 -mnC 0 -mxC 1 -endPlugin -runfork1 > umapinfotohapmap.log
sub umap_info_to_hap_map{

	my $uneak_project_dir = shift;
        die "Error lost the UNEAK project directory" unless defined $uneak_project_dir;
        
        my $project_name = shift;
        die "Error lost the project name" unless defined $project_name;
	
	my $min_minor_allele_freq_hap_map = shift;
	die "Error lost the minimum minor allele frequency" unless defined $min_minor_allele_freq_hap_map;
	
	my $max_minor_allele_freq_hap_map = shift;
	die "Error lost the minimum site coverage" unless defined $max_minor_allele_freq_hap_map;
	
	my $min_call_rate = shift;
	die "Error lost minimum call rate" unless defined $min_call_rate;
	
	my $max_call_rate = shift;
	die "Error lost maximum call rate" unless defined $max_call_rate;

        my $tassel_num_ram = shift;
        die "Error lost amount of tassel memory" unless defined $tassel_num_ram;
        
        my $hap_map_output_dir = shift;
        die "Error lost the output directory to contain output HapMap SNPs files." unless defined $hap_map_output_dir;

        my $umap_info_to_hap_map_stdout_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UMapInfoToHapMap.stdout.log"));
        my $umap_info_to_hap_map_stderr_log_outfile = join('/', $uneak_project_dir, join("_", $project_name, "UMapInfoToHapMap.stderr.log"));
        
        # The HapMap files
        my ($hap_map_fas_txt_outfile, $hap_map_hmc_txt_outfile, $hap_map_hmp_txt_outfile);
	$hap_map_fas_txt_outfile = join('/', $hap_map_output_dir, "HapMap.fas.txt");
	$hap_map_hmc_txt_outfile = join('/', $hap_map_output_dir, "HapMap.hmc.txt");
	$hap_map_hmp_txt_outfile = join('/', $hap_map_output_dir, "HapMap.hmp.txt");
	
	unless(-s $hap_map_fas_txt_outfile and -s $hap_map_hmc_txt_outfile and -s $hap_map_hmp_txt_outfile){
		warn "Generating hap map SNPs file using the UMapInfoToHapMapPlugin.....\n\n";
		my $umapInfoToHapMapCmd  = join("", "$run_pipeline -Xmx", $tassel_num_ram, "g -fork1 -UMapInfoToHapMapPlugin -w $uneak_project_dir -mnMAF $min_minor_allele_freq_hap_map -mxMAF $max_minor_allele_freq_hap_map -mnC $min_call_rate -mxC $max_call_rate -endPlugin -runfork1 1> $umap_info_to_hap_map_stdout_log_outfile 2> $umap_info_to_hap_map_stderr_log_outfile");
		warn $umapInfoToHapMapCmd . "\n\n";
		system($umapInfoToHapMapCmd) == 0 or die "Error calling $umapInfoToHapMapCmd: $?";
	}
	return ($hap_map_fas_txt_outfile, $hap_map_hmc_txt_outfile, $hap_map_hmp_txt_outfile);
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