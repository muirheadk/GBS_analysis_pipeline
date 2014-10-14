#!/usr/bin/perl
use warnings;
use strict;

my $i = 0;
while(<>){
	chomp $_;
	
	if($i eq 0){
		print $_ . "\n";
	}else{
		my @split_row_entry = split(/\t/, $_);
		my ($fastq_plate_num, $fastq_well_num, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) = ($split_row_entry[0], $split_row_entry[1], $split_row_entry[2], $split_row_entry[3], $split_row_entry[4]);
		$fastq_run_id =~ s/\s|-/_/g;
		$fastq_project_leader =~ s/ /_/g;
		$fastq_project_leader = uc($fastq_project_leader);
		print join("\t", $fastq_plate_num, $fastq_well_num, $fastq_run_id, $fastq_project_leader, $fastq_barcode_seq) . "\n";
	}
	$i++;
}