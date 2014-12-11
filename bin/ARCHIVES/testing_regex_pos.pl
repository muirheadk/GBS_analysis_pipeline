#!/usr/bin/perl
use warnings;
use strict;

# @HWI-ST767:215:C30VBACXX:8:1306:10922:63979 1:N:0:
# TGCAAGGATGCAGCCACCAACCTTGGGGCCAAGCGCCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGA

my $sequence = 'CGCGGAGATGCAGGCGGCGGGCGCCGAGGTGATCCGCTCCACCGACAACCTGGTGCGCGCCGCGCGGGACGCCATCACGTGCGACGACGAGCGCAGCCTC';
my $gbs_sequence_length = 100;
my $fastq_adaptor_sequence = qr(C);
print "Sequence: $sequence\n";
print "fastq_adaptor_sequence: $fastq_adaptor_sequence\n";



    while ($sequence =~ m/$fastq_adaptor_sequence/g) {
	if($+[0] eq $gbs_sequence_length){
	
		my $target_start = ($-[0] + 1);
		my $target_end = $+[0];
		print join("\t", $target_start, $target_end) . "\n";
# 	print ($+[0] - $-[0]) . "\n";
	}
    }
