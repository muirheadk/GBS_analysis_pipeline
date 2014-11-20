#!/usr/bin/perl
use warnings;
use strict;

my $uneak_project_dir = "/home/cookeadmin/workspace/GBS_data-08-10-2013/Papilio_GBS_Data/Papilio_notrim-2014-11-18";
my $uneak_sequence_length = 64;


my $hap_map_hmc_infile = join('/', $uneak_project_dir, "hapMap", "HapMap.hmc.txt");
open(INFILE, "<$hap_map_hmc_infile") or die "Couldn't open file $hap_map_hmc_infile for reading, $!";
my $i = 0;
my %individual_ids = ();
my %hap_map_hmc_counts = ();
while(<INFILE>){
	chomp $_;
# 	warn $_ . "\n";
	if($i eq 0){
		my @split_hmp_entries = split(/\t/, $_);
		for(my $j = 1; $j < (scalar(@split_hmp_entries) - 5); $j++){
# 			warn $split_hmp_entries[$j] . "\n";
			my @split_individual_name = split(/_/, $split_hmp_entries[$j]);
			my $individual_id = $split_individual_name[0];
			$individual_ids{$j} = $individual_id;
		}
	}elsif($i ne 0){
		my @split_hmp_entries = split(/\t/, $_);
		my $snp_id = $split_hmp_entries[0];
		for(my $j = 1; $j < (scalar(@split_hmp_entries) - 5); $j++){
			my $individual_id = $individual_ids{$j};
			my @split_sequence_counts = split(/\|/, $split_hmp_entries[$j]);
			my $query_count = $split_sequence_counts[0];
			my $hit_count = $split_sequence_counts[1];
			my $query_id = join("_", $snp_id, "query", $uneak_sequence_length);
			my $hit_id = join("_", $snp_id, "hit", $uneak_sequence_length);
			$hap_map_hmc_counts{$individual_id}{$query_id} = $query_count;
			$hap_map_hmc_counts{$individual_id}{$hit_id} = $hit_count;
			
			
		}
	}
	$i++;
}
close(INFILE) or die "Couldn't close file $hap_map_hmc_infile";
my $individual_counts = 0;
my $snp_id_counts = 0;
foreach my $individual_id (sort keys %hap_map_hmc_counts){
	foreach my $sequence_id (sort keys $hap_map_hmc_counts{$individual_id}){
		print join("\t", $individual_id, $sequence_id, $hap_map_hmc_counts{$individual_id}{$sequence_id}) . "\n";
		$snp_id_counts++;
	}
	$individual_counts++;
}
print $individual_counts . "\n";

print $snp_id_counts . "\n";

my $num_snps = ($snp_id_counts/$individual_counts);
print $num_snps . "\n";