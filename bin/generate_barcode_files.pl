#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

#### PROGRAM NAME ####
# generate_barcode_files.pl - Program to generate barcode input files for input into the Stacks fastq_quality_barcode_splitter.pl script and TASSEL GBS variant pipelines. Takes a GBS shipment and barcode legend file as input and generates a barcode file for each GBS plate listed in the GBS shipment file.

#### DESCRIPTION ####
# This program generates barcode input files for input into the Stacks fastq_quality_barcode_splitter.pl script and TASSEL GBS variant pipelines. Takes a GBS shipment and barcode legend file as input and generates a barcode file for each GBS plate listed in the GBS shipment file.

#### SAMPLE COMMAND ####
# perl generate_barcode_files.pl -i ~/workspace/GBS_data-08-10-2013/GBS_shipment_9Oct2013.csv -l ~/workspace/GBS_data-08-10-2013/Pstl_MspI_barcodes_legend.csv -o ~/workspace/GBS_data-08-10-2013/GBS_BARCODES
my ($GBS_shipment_infile, $barcodes_legend_infile, $flowcell_name, $barcodes_output_dir);
GetOptions(
	'i=s'    => \$GBS_shipment_infile, # The GBS shipment file in tab-delimited format (i.e. GBS_shipment_9Oct2013.csv)
	'l=s'    => \$barcodes_legend_infile, # The barcodes legend file in tab-delimited format (i.e. Pstl_MspI_barcodes_legend.csv)
	'f=s'    => \$flowcell_name, # The flow cell name used for the TASSEL GBS variant pipeline barcodes file format. Default: GQ03122013
	'o=s'    => \$barcodes_output_dir,  # The absolute path to the output directory to contain the barcodes output files
);

# Display a usage message if the following parameters are not specified.
usage() unless (
	defined $GBS_shipment_infile
	and defined $barcodes_legend_infile
	and defined $barcodes_output_dir
);

# The flow cell name used for the TASSEL GBS variant pipeline barcodes file format.
$flowcell_name = "GQ03122013" unless defined $flowcell_name;

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i GBS_shipment_infile -l barcodes_legend_infile -o barcodes_output_dir
    
DESCRIPTION - This program generates barcode input files for input into the Stacks fastq_quality_barcode_splitter.pl script and TASSEL GBS variant pipelines. Takes a GBS shipment and barcode legend file as input and generates a barcode file for each GBS plate listed in the GBS shipment file.

OPTIONS:

-i GBS_shipment_infile - The GBS shipment file in tab-delimited format (i.e. GBS_shipment_9Oct2013.csv)

	e.g. /path/to/GBS_shipment_9Oct2013.csv

-l barcodes_legend_infile - The barcodes legend file in tab-delimited format (i.e. Pstl_MspI_barcodes_legend.csv)

	e.g. /path/to/Pstl_MspI_barcodes_legend.csv

-f flowcell_name - The flow cell name used for the TASSEL GBS variant pipeline barcodes file format. Default: GQ03122013

-o barcodes_output_dir - The absolute path to the output directory to contain the barcodes output files

	e.g. /path/to/barcodes_output_dir

USAGE
}

# Create barcodes output directory if it doesn't already exist.
unless(-d $barcodes_output_dir){
      mkdir($barcodes_output_dir, 0777) or die "Can't make directory: $!";
}

# Parse the GBS shipment file in tab-delimited format (i.e. GBS_shipment_9Oct2013.csv) and store the contents into a hash array variable that indexes
# well number, sample/individual id, and the name of the project leader based on the plate number.
# GBS shipment file format from GBS_shipment_9Oct2013.csv. Taken from plate 1 of the GBS_shipment_9Oct2013.xlsx excel file.
# Plate #	Well #	ID #	Nanodrop ng/ul	260/280	260/230	Qubit ng/ul	Project Leader
# 1	A1	SCP-9	60.6	1.74	1.67	53.6	Greg Breed
# 1	B1	BUL-4	63.64	1.75	1.54	31	Greg Breed
# 1	C1	VT4-4	53.07	1.78	1.67	44.8	Greg Breed
# 1	D1	VT1-7	51.11	1.65	1.44	22.6	Greg Breed
# 1	E1	PIS-7	25.12	1.8	1.92	15.7	Greg Breed
# 1	F1	PIS-12	36.28	1.74	1.61	33	Greg Breed
# 1	G1	HARV-17L	22.81	1.72	1.71	17.3	Greg Breed
my %GBS_samples = ();
open(INFILE, "<$GBS_shipment_infile") or die "Couldn't open file $GBS_shipment_infile for reading, $!";
my $i = 0;
while(<INFILE>){
	chomp $_;
	#    warn "$_\n";
	if($i ne 0){ # Skip the GBS shipment header because we know the positions of the contents within the file.
	
		# Split the contents of the GBS shipment file into the plate number, well number, sample/individual id, and the name of the project leader.
		# We are only interested in these values of the file so that we can concatenate based on the well alphanumeric code. 
		my @split_entries = split(/\t/, $_);
		my ($plate_num, $well_num, $sample_id, $project_leader_name) = ($split_entries[0], $split_entries[1], $split_entries[2], $split_entries[7]);
		
		# Store each well number, sample/individual id, and project leader entry into a hash array indexed by the plate number. 
		# We want to generate barcode output files based on plate number. Once we have the files we can sort by project leader or grep based 
		# on the sample/individual ids of interest.
		push(@{$GBS_samples{$plate_num}}, join("\t", $sample_id, $well_num, $project_leader_name));
	}
	$i++;
    
}
close(INFILE) or die "Couldn't close file $GBS_shipment_infile";

# Parse the barcodes legend file in tab-delimited format (i.e. Pstl_MspI_barcodes_legend.csv) to generate a hash variable that indexes 
# barcodes based on the well alphanumeric code. 
# i.e. {A}{1} is well number A1 and stores barcode sequence "CTCG".
#
# Barcodes legend format from Pstl_MspI_barcodes_legend.csv. Taken from the PstI/MspI portion of the Identifiants RAD.xlsx excel file.
# 	1	2	3	4	5	6	7	8	9	10	11	12
# A	CTCG	AGCG	TTCTG	ATTGA	TCGTT	GAGATA	CTATTA	CTTGCTT	AATATGG	GCGGAAT	TGCAAGGA	CCATGGGT
# B	TGCA	GATG	AGCCG	CATCT	GGTTGT	ATGCCT	GCCAGT	ATGAAAG	ACGTGTT	TAGCGGA	TGGTACGT	CGCGGAGA
# C	ACTA	TCAG	GTATT	CCTAG	CCACGT	AGTGGA	GGAAGA	AAAAGTT	ATTAATT	TCGAAGA	TCTCAGTG	CGTGTGGT
# D	CAGA	TGCGA	CTGTA	GAGGA	TTCAGA	ACCTAA	GTACTT	GAATTCA	ATTGGAT	TCTGTGA	CGCGATAT	GCTGTGGA
# E	AACT	CGCTT	ACCGT	GGAAG	TAGGAA	ATATGT	GTTGAA	GAACTTG	CATAAGT	TGCTGGA	CGCCTTAT	GGATTGGT
# F	GCGT	TCACG	GCTTA	GTCAA	GCTCTA	ATCGTA	TAACGA	GGACCTA	CGCTGAT	ACGACTAG	AACCGAGA	GTGAGGGT
# G	CGAT	CTAGG	GGTGT	TAATA	CCACAA	CATCGT	TGGCTA	GTCGATT	CGGTAGA	TAGCATGG	ACAGGGA	TATCGGGA
# H	GTAA	ACAAA	AGGAT	TACAT	CTTCCA	CGCGGT	TATTTTT	AACGCCT	CTACGGA	TAGGCCAT	ACGTGGTA	TTCCTGGA
my %GBS_barcodes_legend = ();
open(INFILE, "<$barcodes_legend_infile") or die "Couldn't open file $barcodes_legend_infile for reading, $!";
$i = 0;
while(<INFILE>){
	chomp $_;
	#    warn "$_\n";
	if($i ne 0){ # Skip legend header because we know the positions of the alphabetical code and where the barcodes start.
	
		# Split entries into well alphabetical code and a list of barcode sequences.
		my @split_entries = split(/\t/, $_);
		my ($well_alpha_code, @barcodes) = @split_entries;

		# Well alphabetical code starts at 0 and barcode sequences start at array index 1-11.
		my $well_num_counter = 1;
		
		# Iterate through each barcode sequence and index by well alphanumeric code.
		foreach my $barcode (@barcodes){
			$GBS_barcodes_legend{$well_alpha_code}{$well_num_counter} = $barcode;
			$well_num_counter++;
		}
	}
	$i++;
    
}
close(INFILE) or die "Couldn't close file $barcodes_legend_infile";

# Create the TASSEL barcode output directory if it doesn't already exist.
my $tassel_barcodes_output_dir = join('/', $barcodes_output_dir, "TASSEL_BARCODE_OUTFILES");
unless(-d $tassel_barcodes_output_dir){
    mkdir($tassel_barcodes_output_dir, 0777) or die "Can't make directory: $!";
}

# Generate TASSEL GBS variant pipeline barcodes file format.
# Concatenates the contents of the GBS shipment file and barcodes legend file based on well alphanumeric code and separates based on plate number.
# An example of the TASSEL GBS variant pipeline barcodes file format.
# Flowcell        Lane    Barcode Sample  PlateName       Row     Column
# GQ03122013      5       AGCG    M033-11-02-05-DA41      plate5  A       2
# Iterate through each GBS plate number and print the barcode entry based on the GBS shipment and barcode legend files.
foreach my $plate_num (sort keys %GBS_samples){

	# Get the barcode output file path and open the file for writting the barcode entries for each GBS plate.
	my $plate_name = join("_", "plate", $plate_num);
	my $tassel_barcodes_outfile = join('/', $tassel_barcodes_output_dir, join("-", "tassel-barcodes", $plate_name . ".txt"));
	open(OUTFILE, ">$tassel_barcodes_outfile") or die "Couldn't open file $tassel_barcodes_outfile for writting, $!";
	print OUTFILE join("\t", "Flowcell", "Lane", "Barcode", "Sample", "PlateName", "Row", "Column") . "\n";
	
	# Iterate through each GBS sample entry and print the barcode entry.
	foreach my $GBS_sample (@{$GBS_samples{$plate_num}}){
	
		# Split the GBS sample entry and parse for the sample id, alpha-numeric well number, and project leader.
		my @split_GBS_samples = split(/\t/, $GBS_sample);
		my ($sample_id, $well_num, $project_leader) = @split_GBS_samples;
		
		# Split the well alpha numeric code and obtain the alpha and number codes separately.
		my ($well_alpha_code, $well_num_code) = split('', $well_num);

		# Print the barcode entry.
		print OUTFILE join("\t", $flowcell_name, $plate_num, $GBS_barcodes_legend{$well_alpha_code}{$well_num_code}, $sample_id, $plate_name, $well_alpha_code, $well_num_code) . "\n";
	}
	close(OUTFILE) or die "Couldn't close file $tassel_barcodes_outfile";
    
}

# Create the Stacks GBS variant pipeline barcodes output directory if it doesn't already exist.
my $stacks_barcodes_output_dir = join('/', $barcodes_output_dir, "STACKS_BARCODE_OUTFILES");
unless(-d $stacks_barcodes_output_dir){
    mkdir($stacks_barcodes_output_dir, 0777) or die "Can't make directory: $!";
}

# Generate the Stacks barcodes file format for input into the fastq_quality_barcode_splitter.pl script.
# Concatenates the contents of the GBS shipment file and barcodes legend file based on well alphanumeric code and separates based on plate number.
# An example of the barcodes file format.
# plate_number    well_row        well_column     run_id  project_leader_name  barcode_sequence
# 5       A       2       M033_11_02_05_DA41      JASMINE_JANES   AGCG
# Iterate through each GBS plate number and print the barcode entry based on the GBS shipment and barcode legend files.
foreach my $plate_num (sort keys %GBS_samples){

	# Get the barcode output file path and open the file for writting the barcode entries for each GBS plate.
	my $plate_name = join("_", "plate", $plate_num);
	my $barcodes_outfile = join('/', $stacks_barcodes_output_dir, join("-", "stacks-barcodes", $plate_name . ".txt"));
	open(OUTFILE, ">$barcodes_outfile") or die "Couldn't open file $barcodes_outfile for writting, $!";
	print OUTFILE join("\t", "plate_number", "well_row", "well_column", "run_id", "project_leader_name", "barcode_sequence") . "\n";
	
	# Iterate through each GBS sample entry and print the barcode entry.
	foreach my $GBS_sample (@{$GBS_samples{$plate_num}}){
		# Split the GBS sample entry and parse for the sample id, alpha-numeric well number, and project leader.
		my @split_GBS_samples = split(/\t/, $GBS_sample);
		my ($sample_id, $well_num, $project_leader_name) = @split_GBS_samples;

		# Grab the well alpha numeric code and obtain the alpha and number codes separately.
		my ($well_alpha_code, $well_num_code);
		if($well_num =~ m/([A-Z]{1})([0-9]+)/){
			($well_alpha_code, $well_num_code) = ($1, $2);
		}

		# Replace any spaces in the sample id with an underscore.
		$sample_id =~ s/ /-/g;

		# Replace any spaces in the project leader name with an underscore and convert to uppercase if not already done.
		$project_leader_name =~ s/ /_/g;
		$project_leader_name = uc($project_leader_name);

		# Print the barcode entry.
		print OUTFILE join("\t", $plate_num, $well_alpha_code, $well_num_code, $sample_id, $project_leader_name, $GBS_barcodes_legend{$well_alpha_code}{$well_num_code}) . "\n";
	}
	close(OUTFILE) or die "Couldn't close file $barcodes_outfile";
}
