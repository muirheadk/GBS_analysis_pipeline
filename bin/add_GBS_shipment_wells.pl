#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

my ($GBS_shipment_infile, $output_dir);
GetOptions(
	'i=s'    => \$GBS_shipment_infile, # The GBS shipment file in tab-delimited format (i.e. GBS_shipment_9Oct2013.csv)
	'o=s'    => \$output_dir,  # The absolute path to the output directory to contain the barcodes output files
);

# Display a usage message if the following parameters are not specified.
usage() unless (
	defined $GBS_shipment_infile
	and defined $output_dir
);

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i GBS_shipment_infile -o output_dir
    
Description - 
    
OPTIONS:
      -i GBS_shipment_infile - 

      -o output_dir -
    
USAGE
}

# Create barcodes output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %sbw_entries = ();
my $header = "";
open(INFILE, "<$GBS_shipment_infile") or die "Couldn't open file $GBS_shipment_infile for reading, $!";
my $i = 0;
while(<INFILE>){
	chomp $_;
	
	if($i eq 0){
		my @split_header = split(/\t/, $_);
		
		$header = join("\t", @split_header[1..$#split_header]);
	}elsif($i ne 0){
		my @split_sbw_entry = split(/\t/, $_);
		my ($flowcell_file, $plate_num, $well_num, $run_id, $nandrop, $f260_280, $f260_230, $qubit, $project_leader) = @split_sbw_entry;
		my $flowcell_filename = fileparse($flowcell_file, qr/.fastq.gz/);
		
		my $new_project_leader = join("_", $project_leader, $flowcell_filename);
  		$sbw_entries{$flowcell_filename}{$well_num} = join("\t", $flowcell_filename, $well_num, $run_id, $nandrop, $f260_280, $f260_230, $qubit, $new_project_leader);
	}
	$i++;
}
close(INFILE) or die "Couldn't close file $GBS_shipment_infile";

my @alpha_well_codes = ("A","B","C","D","E","F","G","H");

my @alpha_numeric_well_codes = ();
for(my $i = 1; $i <= 12; $i++){
	foreach my $alpha_well_code (@alpha_well_codes){
# 		warn join("", $alpha_well_code, $i);
		push(@alpha_numeric_well_codes, join("", $alpha_well_code, $i));
	}
}

my $sbw_filename = fileparse($GBS_shipment_infile, qr/\.csv/);
my $sbw_outfile = join('/', $output_dir, join("_", $sbw_filename, "all_wells") . ".csv");
open(OUTFILE, ">$sbw_outfile") or die "Couldn't open file $sbw_outfile for writting, $!";
print OUTFILE $header . "\n";
foreach my $flowcell_filename (sort keys %sbw_entries){
	warn $flowcell_filename . "\n";
	my $i = 1;
	foreach my $well_num (@alpha_numeric_well_codes){
		if(defined($sbw_entries{$flowcell_filename}{$well_num})){
 			print OUTFILE $sbw_entries{$flowcell_filename}{$well_num} . "\n";
		}else{
			print OUTFILE join("\t", $flowcell_filename, $well_num, join("_", "miscellaneous", $well_num) , 0, 0, 0, 0, join("_", "miscellaneous", $flowcell_filename)) . "\n";
		}
	}
}
close(OUTFILE) or die "Couldn't close file $sbw_outfile";
