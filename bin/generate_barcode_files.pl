#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

# perl generate_barcode_files.pl -i ~/workspace/GBS_data-08-10-2013/GBS_BARCODES/GBS_shipment_9Oct2013.csv -l ~/workspace/GBS_data-08-10-2013/GBS_BARCODES/Pstl_MspI_barcodes_legend.csv -o ~/workspace/GBS_data-08-10-2013/GBS_BARCODES/BARCODE_OUTFILES

my ($gbs_shipment_infile, $barcodes_legend_infile, $flowcell_name, $output_dir);
GetOptions(
    'i=s'    => \$gbs_shipment_infile,
    'l=s'    => \$barcodes_legend_infile,
    'f=s'    => \$flowcell_name,
    'o=s'    => \$output_dir,
);
usage() unless (
      defined $gbs_shipment_infile
      and defined $barcodes_legend_infile
      and defined $output_dir
);

$flowcell_name = "GQ03122013" unless defined $flowcell_name;

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i gbs_shipment_infile -l barcodes_legend_infile -f flowcell_name -o output_dir
    
Description - 
    
OPTIONS:
      -i gbs_shipment_infile -
    
      -l barcodes_legend_infile -
    
      -f flowcell_name -
      
      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %gbs_samples = ();
open(INFILE, "<$gbs_shipment_infile") or die "Couldn't open file $gbs_shipment_infile for reading, $!";
my $i = 0;
while(<INFILE>){
    chomp $_;
#    warn "$_\n";
    if($i ne 0){
        my @split_entries = split(/\t/, $_);
        my ($plate_num, $well_num, $sample_id, $project_leader) = ($split_entries[0], $split_entries[1], $split_entries[2], $split_entries[7]);
        push(@{$gbs_samples{$plate_num}}, join("\t", $sample_id, $well_num, $project_leader));
        
    }
    $i++;
    
}
close(INFILE) or die "Couldn't close file $gbs_shipment_infile";

#foreach my $sample_id (sort keys %gbs_samples){
#    warn join("\t", $sample_id, $gbs_samples{$sample_id}) . "\n";
#    
#}
my %gbs_barcodes_legend = ();
open(INFILE, "<$barcodes_legend_infile") or die "Couldn't open file $barcodes_legend_infile for reading, $!";
$i = 0;
while(<INFILE>){
    chomp $_;
#    warn "$_\n";
    if($i ne 0){
        my @split_entries = split(/\t/, $_);
        my ($well_alpha_code, @barcodes) = @split_entries;
        
#        die $well_alpha_code . "\n";
        
        my $well_num_counter = 1;
        foreach my $barcode (@barcodes){
            $gbs_barcodes_legend{$well_alpha_code}{$well_num_counter} = $barcode;
            $well_num_counter++;
        }
    }
    $i++;
    
}
close(INFILE) or die "Couldn't close file $barcodes_legend_infile";

#foreach my $well_alpha_code (sort {$a cmp $b} keys %gbs_barcodes_legend){
#    foreach my $well_num_counter (sort {$a <=> $b} keys $gbs_barcodes_legend{$well_alpha_code}){
#        warn join("\t", $well_alpha_code, $well_num_counter,$gbs_barcodes_legend{$well_alpha_code}{$well_num_counter}) . "\n";
#    }
#
#}

# Create output directory if it doesn't already exist.
my $tassel_barcodes_output_dir = join('/', $output_dir, "TASSEL_BARCODE_OUTFILES");
unless(-d $tassel_barcodes_output_dir){
    mkdir($tassel_barcodes_output_dir, 0777) or die "Can't make directory: $!";
}

#Flowcell        Lane    Barcode Sample  PlateName       Row     Column
#GQ03122013      5       AGCG    M033-11-02-05-DA41      plate5  A       2

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

# Create output directory if it doesn't already exist.
my $split_barcodes_output_dir = join('/', $output_dir, "FASTX_BARCODE_SPLITTER_OUTFILES");
unless(-d $split_barcodes_output_dir){
    mkdir($split_barcodes_output_dir, 0777) or die "Can't make directory: $!";
}

#plate_number    well_row        well_column     run_id  project_leader  barcode_sequence
foreach my $plate_num (sort keys %gbs_samples){
    
    my $plate_name = join("", "plate", $plate_num);
    my $split_barcodes_outfile = join('/', $split_barcodes_output_dir, join("_", "fastx_split_barcodes", $plate_name . ".txt"));
    open(OUTFILE, ">$split_barcodes_outfile") or die "Couldn't open file $split_barcodes_outfile for writting, $!";
    print OUTFILE join("\t", "plate_number", "well_row", "well_column", "run_id", "project_leader", "barcode_sequence") . "\n";
#    5       A       2       M033_11_02_05_DA41      JASMINE_JANES   AGCG
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

