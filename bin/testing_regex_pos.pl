#!/usr/bin/perl
use warnings;
use strict;

# @HWI-ST767:215:C30VBACXX:8:1306:10922:63979 1:N:0:
# TGCAAGGATGCAGCCACCAACCTTGGGGCCAAGCGCCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGA

my $sequence = 'TGCAAGGATGCAGCCACCAACCTTGGGGCCAAGCGCCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGA';
my $fastq_adaptor_sequence = qr(CCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG);
print "Sequence: $sequence\n";
print "fastq_adaptor_sequence: $fastq_adaptor_sequence\n";

if($sequence =~ /$fastq_adaptor_sequence/g){
	print join("\t", $-[0], $+[0]) . "\n";
	print ($+[0] - $-[0]);
}

