#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

# perl generate_barcode_files.pl -i ~/workspace/GBS_data-08-10-2013/GBS_BARCODES/GBS_shipment_9Oct2013.csv -l ~/workspace/GBS_data-08-10-2013/GBS_BARCODES/Pstl_MspI_barcodes_legend.csv -o ~/workspace/GBS_data-08-10-2013/GBS_BARCODES/BARCODE_OUTFILES
my ($gbs_shipment_infile, $barcodes_legend_infile, $flowcell_name, $output_dir);
GetOptions(
    'i=s'    => \$gbs_shipment_infile, # The GBS shipment file in tab-delimited format (i.e. GBS_shipment_9Oct2013.csv)
    'l=s'    => \$barcodes_legend_infile, # The barcodes legend file in tab-delimited format (i.e. Pstl_MspI_barcodes_legend.csv)
    'f=s'    => \$flowcell_name, # The flow cell name used for the TASSEL GBS variant pipeline barcodes file format. Default: GQ03122013
    'o=s'    => \$output_dir,  # The absolute path to the output directory to contain the barcodes output files
);

# Display usage message if the following parameters are not specified.
usage() unless (
      defined $gbs_shipment_infile
      and defined $barcodes_legend_infile
      and defined $output_dir
);

# The flow cell name used for the TASSEL GBS variant pipeline barcodes file format.
$flowcell_name = "GQ03122013" unless defined $flowcell_name;

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i gbs_shipment_infile -l barcodes_legend_infile -f flowcell_name -o output_dir
    
DESCRIPTION - Generate TASSEL GBS variant pipeline and FASTX toolkit barcodes file formats based on a GBS shipment and Barcode legend file.
    
OPTIONS:

-i gbs_shipment_infile - The GBS shipment file in tab-delimited format (i.e. GBS_shipment_9Oct2013.csv)

	e.g. /path/to/GBS_shipment_9Oct2013.csv

-l barcodes_legend_infile - The barcodes legend file in tab-delimited format (i.e. Pstl_MspI_barcodes_legend.csv)

	e.g. /path/to/Pstl_MspI_barcodes_legend.csv
	
-f flowcell_name - The flow cell name used for the TASSEL GBS variant pipeline barcodes file format. Default: GQ03122013
      
-o output_dir - The absolute path to the output directory to contain the barcodes output files

	e.g. /path/to/output_dir

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
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
my %gbs_samples = ();
open(INFILE, "<$gbs_shipment_infile") or die "Couldn't open file $gbs_shipment_infile for reading, $!";
my $i = 0;
while(<INFILE>){
	chomp $_;
	#    warn "$_\n";
	if($i ne 0){ # Skip the GBS shipment header because we know the positions of the contents within the file.
		# Split the contents of the GBS shipment file into the plate number, well number, sample/individual id, and the name of the project leader.
		# We are only interested in these values of the file so that we can concatenate based on the well alphanumeric code. 
		my @split_entries = split(/\t/, $_);
		my ($plate_num, $well_num, $sample_id, $project_leader) = ($split_entries[0], $split_entries[1], $split_entries[2], $split_entries[7]);
		# Store each well number, sample/individual id, and project leader entry into a hash array indexed by the plate number. 
		# We want to generate barcode output files based on plate number. Once we have the files we can sort by project leader or grep based 
		# on the sample/individual ids of interest.
		push(@{$gbs_samples{$plate_num}}, join("\t", $sample_id, $well_num, $project_leader));
	}
	$i++;
    
}
close(INFILE) or die "Couldn't close file $gbs_shipment_infile";

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
my %gbs_barcodes_legend = ();
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
            $gbs_barcodes_legend{$well_alpha_code}{$well_num_counter} = $barcode;
            $well_num_counter++;
        }
    }
    $i++;
    
}
close(INFILE) or die "Couldn't close file $barcodes_legend_infile";

# Create the TASSEL barcode output directory if it doesn't already exist.
my $tassel_barcodes_output_dir = join('/', $output_dir, "TASSEL_BARCODE_OUTFILES");
unless(-d $tassel_barcodes_output_dir){
    mkdir($tassel_barcodes_output_dir, 0777) or die "Can't make directory: $!";
}

# Generate TASSEL GBS variant pipeline barcodes file format.
# Concatenates the contents of the GBS shipment file and barcodes legend file based on well alphanumeric code and separates based on plate number.
# An example of the TASSEL GBS variant pipeline barcodes file format.
# Flowcell        Lane    Barcode Sample  PlateName       Row     Column
# GQ03122013      5       AGCG    M033-11-02-05-DA41      plate5  A       2
foreach my $plate_num (sort keys %gbs_samples){
    
    my $plate_name = join("", "plate", $plate_num);
    my $tassel_barcodes_outfile = join('/', $tassel_barcodes_output_dir, join("_", "tassel_barcodes", $plate_name . ".txt"));
    open(OUTFILE, ">$tassel_barcodes_outfile") or die "Couldn't open file $tassel_barcodes_outfile for writting, $!";
    print OUTFILE join("\t", "Flowcell", "Lane", "Barcode", "Sample", "PlateName", "Row", "Column") . "\n";
    foreach my $gbs_sample (@{$gbs_samples{$plate_num}}){
        
        my @split_gbs_samples = split(/\t/, $gbs_sample);
        my ($sample_id, $well_num, $project_leader) = @split_gbs_samples;
        my ($well_alpha_code, $well_num_code) = split('', $well_num);
        print OUTFILE join("\t", $flowcell_name, $plate_num, $gbs_barcodes_legend{$well_alpha_code}{$well_num_code}, $sample_id, $plate_name, $well_alpha_code, $well_num_code) . "\n";
    }
    close(OUTFILE) or die "Couldn't close file $tassel_barcodes_outfile";
    
}

# Create the fastx barcode splitter cascade barcode file output directory if it doesn't already exist.
my $split_barcodes_output_dir = join('/', $output_dir, "FASTX_BARCODE_SPLITTER_OUTFILES");
unless(-d $split_barcodes_output_dir){
    mkdir($split_barcodes_output_dir, 0777) or die "Can't make directory: $!";
}

# Generate FASTX toolkit barcodes file format for input into the fastx_barcode_splitter_cascade.pl script.
# Concatenates the contents of the GBS shipment file and barcodes legend file based on well alphanumeric code and separates based on plate number.
# An example of the FASTX toolkit barcodes file format.
# plate_number    well_row        well_column     run_id  project_leader  barcode_sequence
# 5       A       2       M033_11_02_05_DA41      JASMINE_JANES   AGCG
foreach my $plate_num (sort keys %gbs_samples){
    
    my $plate_name = join("", "plate", $plate_num);
    my $split_barcodes_outfile = join('/', $split_barcodes_output_dir, join("_", "fastx_split_barcodes", $plate_name . ".txt"));
    open(OUTFILE, ">$split_barcodes_outfile") or die "Couldn't open file $split_barcodes_outfile for writting, $!";
    print OUTFILE join("\t", "plate_number", "well_row", "well_column", "run_id", "project_leader", "barcode_sequence") . "\n";

    foreach my $gbs_sample (@{$gbs_samples{$plate_num}}){
        
        my @split_gbs_samples = split(/\t/, $gbs_sample);
        my ($sample_id, $well_num, $project_leader) = @split_gbs_samples;
        my ($well_alpha_code, $well_num_code);
        if($well_num =~ m/([A-Z]{1})([0-9]+)/){
        	($well_alpha_code, $well_num_code) = ($1, $2);
        }
        $sample_id =~ s/ /_/g;
        $sample_id =~ s/-/_/g;
        $project_leader =~ s/ /_/g;
        $project_leader = uc($project_leader);
        print OUTFILE join("\t", $plate_num, $well_alpha_code, $well_num_code, $sample_id, $project_leader, $gbs_barcodes_legend{$well_alpha_code}{$well_num_code}) . "\n";
    }
    close(OUTFILE) or die "Couldn't close file $split_barcodes_outfile";
    
}

