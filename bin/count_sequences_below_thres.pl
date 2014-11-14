#!/usr/bin/perl
use warnings;
use strict;

my $i = 0;
my $total_counter = 0;
my $below_counter = 0;
my $above_counter = 0;
while(<>){
	chomp $_;
	warn $_ . "\n";
	if($i ne 0){
		my @split_entries = split(/\t/, $_);
		my $adaptor_align_length = $split_entries[4];
		
		my $trimmed_length = (100 - $adaptor_align_length);
		
		$above_counter++ if($trimmed_length >= 72);
		$below_counter++ if($trimmed_length < 72);
		$total_counter++;
	}
	$i++;
}

print "total=$total_counter\n";
print "below=$below_counter\n";
print "above=$above_counter\n";