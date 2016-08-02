#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);

use File::Basename;
use File::Copy;

# VERSION 1.22

# populations -P ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -M ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/mpb_poputlations.txt -b 1 -k -p 10 -r 0.75 -m 5 -a 0.1 -f p_value --p_value_cutoff 0.05 -t 7 --fasta

my ($stacks_dir, $pop_map_infile, $batch_num, $min_num_pops_locus, $min_percent_indvs_pop, $min_stack_depth, $min_allele_freq, $outfile_type, $num_cpu_cores, $output_dir);
GetOptions(
	'P=s'    => \$stacks_dir, # The absolute path to the Stacks output files.
	'M=s'    => \$pop_map_infile, # The  absolute path to the population map, a tab separated file describing which individuals belong in which population.
	'b=s'    => \$batch_num, # The batch ID to examine when exporting from the catalog.
	'p=s'    => \$min_num_pops_locus, # The minimum number of populations a locus must be present in to process a locus. Default: 2
	'r=s'    => \$min_percent_indvs_pop, # The minimum percentage of individuals in a population required to process a locus for that population.  Default: 0
	'm=s'    => \$min_stack_depth, # The minimum stack depth required for individuals at a locus. Default: 10
	'a=s'    => \$min_allele_freq, # The minimum minor allele frequency required before calculating Fst at a locus (0 < a < 0.5). Default: 0.05
	'f=s'    => \$outfile_type, # The output file format type. Outputs the full sequence for each allele, from each sample locus in FASTA format, Population structure for each allele, from each sample locus in Structure format, etc. Can be either fasta, structure, fasta_single_snp, structure_single_snp, write_singleDefault: 'fasta'
	't=s'    => \$num_cpu_cores, # The number of threads to run in parallel sections of code. Default: 2
	'o=s'    => \$output_dir, # The absolute path to the output directory to contain the Stacks output files and directories.
);

usage() unless (
	defined $stacks_dir
	and defined $pop_map_infile
	and defined $output_dir
);

# The batch ID to examine when exporting from the catalog.
$batch_num = 1 unless defined $batch_num;

# The minimum number of populations a locus must be present in to process a locus. Default: 2
$min_num_pops_locus = 2 unless defined $min_num_pops_locus;

# The minimum percentage of individuals in a population required to process a locus for that population. Default: 0
$min_percent_indvs_pop = 0 unless defined $min_percent_indvs_pop;

# The minimum stack depth required for individuals at a locus. Default: 10
$min_stack_depth = 10 unless defined $min_stack_depth;

# The minimum minor allele frequency required before calculating Fst at a locus (0 < a < 0.5). Default: 0.05
$min_allele_freq = 0.05 unless defined $min_allele_freq;

# The number of threads to run in parallel sections of code. Default: 2
$num_cpu_cores = 2 unless defined $num_cpu_cores;

# The output file format type. Outputs the full sequence for each allele, from each sample locus in FASTA format, Population structure for each allele, from each sample locus in Structure format, etc. Can be either fasta, structure, fasta_single_snp, structure_single_snp, write_singleDefault: 'fasta'
$outfile_type = 'fasta' unless defined $outfile_type;

# Program dependencies - The stacks populations program.
my $populations			= '/usr/local/bin/populations';

sub usage {

die <<"USAGE";

Usage: $0 -P stacks_dir -M pop_map_infile -b batch_num -p min_num_pops_locus -r min_percent_indvs_pop -m min_stack_depth -a min_allele_freq -f outfile_type -t num_cpu_cores -o output_dir

VERSION 1.22

DESCRIPTION – Program that executes the Stacks populations program that analyzes a population of individual samples given a population map file that designates what individual belongs to what population. The populations program computes a number of population genetics statistics as well as exporting a variety of standard output formats such as fasta and structure formats.

OPTIONS:

-P stacks_dir - The absolute path to the Stacks output files.

-M pop_map_infile - The  absolute path to the population map, a tab separated file describing which individuals belong in which population.

-b batch_num - The batch ID to examine when exporting from the catalog.

-p min_num_pops_locus - The minimum number of populations a locus must be present in to process a locus. Default: 2

-r min_percent_indvs_pop - The minimum percentage of individuals in a population required to process a locus for that population.  Default: 0

-m min_stack_depth - The minimum stack depth required for individuals at a locus. Default: 10

-a min_allele_freq - The minimum minor allele frequency required before calculating Fst at a locus (0 < a < 0.5). Default: 0.05

-f outfile_type – The output file format type. Outputs the full sequence for each allele, from each sample locus in FASTA format, Population structure for each allele, from each sample locus in Structure format, etc. Can be either fasta, structure, fasta_single_snp, structure_single_snp, write_singleDefault: 'fasta'

-t num_cpu_cores - The number of threads to run in parallel sections of code. Default: 2

-o output_dir - The absolute path to the output directory to contain the populations stacks output files and directories.

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my ($batch_fasta_infile, $batch_structure_infile, $batch_vcf_infile, $batch_hapstats_infile, $batch_haplotypes_infile, $batch_populations_log_infile, $batch_sumstats_infile, $batch_sumstats_summary_infile);
$batch_fasta_infile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".fa"));
$batch_structure_infile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".structure.tsv"));
$batch_vcf_infile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".vcf"));
$batch_hapstats_infile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".hapstats.tsv"));
$batch_haplotypes_infile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".haplotypes.tsv"));
$batch_populations_log_infile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".populations.log"));
$batch_sumstats_infile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".sumstats.tsv"));
$batch_sumstats_summary_infile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".sumstats_summary.tsv"));

my %batch_files = ();
if($outfile_type eq 'fasta'){
	# Create output directory if it doesn't already exist.
	my $populations_fasta_output_dir = join('/', $output_dir, "POPULATIONS_FASTA_DIR");
	
	my $batch_fasta_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".fa"));
	unless(-s $batch_fasta_outfile){
		warn "Executing populations fasta output.....\n\n";
		my $populationsCmd  = "$populations -P $stacks_dir -M $pop_map_infile -b $batch_num -p $min_num_pops_locus -r $min_percent_indvs_pop -m $min_stack_depth -a $min_allele_freq -t $num_cpu_cores --fasta";
		warn $populationsCmd . "\n\n";
		system($populationsCmd) == 0 or die "Error calling $populationsCmd: $?";
		
		# Create output directory if it doesn't already exist.
		unless(-d $populations_fasta_output_dir){
			mkdir($populations_fasta_output_dir, 0777) or die "Can't make directory: $!";
		}
			
		my ($batch_hapstats_outfile, $batch_haplotypes_outfile, $batch_populations_log_outfile, $batch_sumstats_outfile, $batch_sumstats_summary_outfile);
		$batch_hapstats_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".hapstats.tsv"));
		$batch_haplotypes_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".haplotypes.tsv"));
		$batch_populations_log_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".populations.log"));
		$batch_sumstats_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".sumstats.tsv"));
		$batch_sumstats_summary_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".sumstats_summary.tsv"));
		
		$batch_files{$batch_fasta_infile} = $batch_fasta_outfile;
		$batch_files{$batch_hapstats_infile} = $batch_hapstats_outfile;
		$batch_files{$batch_haplotypes_infile} = $batch_haplotypes_outfile;
		$batch_files{$batch_populations_log_infile} = $batch_populations_log_outfile;
		$batch_files{$batch_sumstats_infile} = $batch_sumstats_outfile;
		$batch_files{$batch_sumstats_summary_infile} = $batch_sumstats_summary_outfile;
		
		foreach my $batch_infile (sort keys %batch_files){
			my $batch_outfile = $batch_files{$batch_infile};
			warn "Copying $batch_infile to $batch_outfile.....";
			copy($batch_infile, $batch_outfile) or die "Copy failed: $!";
			warn "Unlinking $batch_infile.....";
	 		unlink($batch_infile) or die "Could not unlink $batch_infile: $!";
		}
	}
	
}if($outfile_type eq 'fasta_single_snp'){
	# Create output directory if it doesn't already exist.
	my $populations_fasta_output_dir = join('/', $output_dir, "POPULATIONS_FASTA_DIR");
	
	my $batch_fasta_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".fa"));
	unless(-s $batch_fasta_outfile){
		warn "Executing populations fasta output.....\n\n";
		my $populationsCmd  = "$populations -P $stacks_dir -M $pop_map_infile -b $batch_num -p $min_num_pops_locus -r $min_percent_indvs_pop -m $min_stack_depth -a $min_allele_freq -t $num_cpu_cores --fasta --write_single_snp";
		warn $populationsCmd . "\n\n";
		system($populationsCmd) == 0 or die "Error calling $populationsCmd: $?";
		
		# Create output directory if it doesn't already exist.
		unless(-d $populations_fasta_output_dir){
			mkdir($populations_fasta_output_dir, 0777) or die "Can't make directory: $!";
		}
			
		my ($batch_hapstats_outfile, $batch_haplotypes_outfile, $batch_populations_log_outfile, $batch_sumstats_outfile, $batch_sumstats_summary_outfile);
		$batch_hapstats_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".hapstats.tsv"));
		$batch_haplotypes_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".haplotypes.tsv"));
		$batch_populations_log_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".populations.log"));
		$batch_sumstats_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".sumstats.tsv"));
		$batch_sumstats_summary_outfile = join('/', $populations_fasta_output_dir, join("_", "batch", $batch_num . ".sumstats_summary.tsv"));
		
		$batch_files{$batch_fasta_infile} = $batch_fasta_outfile;
		$batch_files{$batch_hapstats_infile} = $batch_hapstats_outfile;
		$batch_files{$batch_haplotypes_infile} = $batch_haplotypes_outfile;
		$batch_files{$batch_populations_log_infile} = $batch_populations_log_outfile;
		$batch_files{$batch_sumstats_infile} = $batch_sumstats_outfile;
		$batch_files{$batch_sumstats_summary_infile} = $batch_sumstats_summary_outfile;
		
		foreach my $batch_infile (sort keys %batch_files){
			my $batch_outfile = $batch_files{$batch_infile};
			warn "Copying $batch_infile to $batch_outfile.....";
			copy($batch_infile, $batch_outfile) or die "Copy failed: $!";
			warn "Unlinking $batch_infile.....";
	 		unlink($batch_infile) or die "Could not unlink $batch_infile: $!";
		}
	}
	
}elsif($outfile_type eq 'struct_single_snp'){

	# Create output directory if it doesn't already exist.
	my $populations_structure_output_dir = join('/', $output_dir, "POPULATIONS_STRUCTURE_DIR");
	
	my $batch_structure_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".structure.tsv"));
	unless(-s $batch_structure_outfile){
		warn "Executing populations structure output.....\n\n";
		my $populationsCmd  = "$populations -P $stacks_dir -M $pop_map_infile -b $batch_num -p $min_num_pops_locus -r $min_percent_indvs_pop -m $min_stack_depth -a $min_allele_freq -t $num_cpu_cores --structure --write_single_snp";
		warn $populationsCmd . "\n\n";
		system($populationsCmd) == 0 or die "Error calling $populationsCmd: $?";
		
		# Create output directory if it doesn't already exist.
		unless(-d $populations_structure_output_dir){
			mkdir($populations_structure_output_dir, 0777) or die "Can't make directory: $!";
		}
		
		my ($batch_hapstats_outfile, $batch_haplotypes_outfile, $batch_populations_log_outfile, $batch_sumstats_outfile, $batch_sumstats_summary_outfile);
		$batch_hapstats_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".hapstats.tsv"));
		$batch_haplotypes_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".haplotypes.tsv"));
		$batch_populations_log_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".populations.log"));
		$batch_sumstats_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".sumstats.tsv"));
		$batch_sumstats_summary_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".sumstats_summary.tsv"));
		
		$batch_files{$batch_structure_infile} = $batch_structure_outfile;
		$batch_files{$batch_hapstats_infile} = $batch_hapstats_outfile;
		$batch_files{$batch_haplotypes_infile} = $batch_haplotypes_outfile;
		$batch_files{$batch_populations_log_infile} = $batch_populations_log_outfile;
		$batch_files{$batch_sumstats_infile} = $batch_sumstats_outfile;
		$batch_files{$batch_sumstats_summary_infile} = $batch_sumstats_summary_outfile;
		
		foreach my $batch_infile (sort keys %batch_files){
			my $batch_outfile = $batch_files{$batch_infile};
			warn "Copying $batch_infile to $batch_outfile.....";
			copy($batch_infile, $batch_outfile) or die "Copy failed: $!";
			warn "Unlinking $batch_infile.....";
	 		unlink($batch_infile) or die "Could not unlink $batch_infile: $!";
		}
	}
	
}elsif($outfile_type eq 'structure'){

	# Create output directory if it doesn't already exist.
	my $populations_structure_output_dir = join('/', $output_dir, "POPULATIONS_STRUCTURE_DIR");
	
	my $batch_structure_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".structure.tsv"));
	unless(-s $batch_structure_outfile){
		warn "Executing populations structure output.....\n\n";
		my $populationsCmd  = "$populations -P $stacks_dir -M $pop_map_infile -b $batch_num -p $min_num_pops_locus -r $min_percent_indvs_pop -m $min_stack_depth -a $min_allele_freq -t $num_cpu_cores --structure";
		warn $populationsCmd . "\n\n";
		system($populationsCmd) == 0 or die "Error calling $populationsCmd: $?";
		
		# Create output directory if it doesn't already exist.
		unless(-d $populations_structure_output_dir){
			mkdir($populations_structure_output_dir, 0777) or die "Can't make directory: $!";
		}
		
		my ($batch_hapstats_outfile, $batch_haplotypes_outfile, $batch_populations_log_outfile, $batch_sumstats_outfile, $batch_sumstats_summary_outfile);
		$batch_hapstats_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".hapstats.tsv"));
		$batch_haplotypes_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".haplotypes.tsv"));
		$batch_populations_log_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".populations.log"));
		$batch_sumstats_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".sumstats.tsv"));
		$batch_sumstats_summary_outfile = join('/', $populations_structure_output_dir, join("_", "batch", $batch_num . ".sumstats_summary.tsv"));
		
		$batch_files{$batch_structure_infile} = $batch_structure_outfile;
		$batch_files{$batch_hapstats_infile} = $batch_hapstats_outfile;
		$batch_files{$batch_haplotypes_infile} = $batch_haplotypes_outfile;
		$batch_files{$batch_populations_log_infile} = $batch_populations_log_outfile;
		$batch_files{$batch_sumstats_infile} = $batch_sumstats_outfile;
		$batch_files{$batch_sumstats_summary_infile} = $batch_sumstats_summary_outfile;
		
		foreach my $batch_infile (sort keys %batch_files){
			my $batch_outfile = $batch_files{$batch_infile};
			warn "Copying $batch_infile to $batch_outfile.....";
			copy($batch_infile, $batch_outfile) or die "Copy failed: $!";
			warn "Unlinking $batch_infile.....";
	 		unlink($batch_infile) or die "Could not unlink $batch_infile: $!";
		}
	}
	
}elsif($outfile_type eq 'vcf'){

	# Create output directory if it doesn't already exist.
	my $populations_vcf_output_dir = join('/', $output_dir, "POPULATIONS_VCF_DIR");
	
	my $batch_vcf_outfile = join('/', $populations_vcf_output_dir, join("_", "batch", $batch_num . ".vcf"));
	unless(-s $batch_vcf_outfile){
		warn "Executing populations vcf output.....\n\n";
		my $populationsCmd  = "$populations -P $stacks_dir -M $pop_map_infile -b $batch_num -p $min_num_pops_locus -r $min_percent_indvs_pop -m $min_stack_depth -a $min_allele_freq -t $num_cpu_cores --vcf";
		warn $populationsCmd . "\n\n";
		system($populationsCmd) == 0 or die "Error calling $populationsCmd: $?";
		
		# Create output directory if it doesn't already exist.
		unless(-d $populations_vcf_output_dir){
			mkdir($populations_vcf_output_dir, 0777) or die "Can't make directory: $!";
		}
		
		my ($batch_hapstats_outfile, $batch_haplotypes_outfile, $batch_populations_log_outfile, $batch_sumstats_outfile, $batch_sumstats_summary_outfile);
		$batch_hapstats_outfile = join('/', $populations_vcf_output_dir, join("_", "batch", $batch_num . ".hapstats.tsv"));
		$batch_haplotypes_outfile = join('/', $populations_vcf_output_dir, join("_", "batch", $batch_num . ".haplotypes.tsv"));
		$batch_populations_log_outfile = join('/', $populations_vcf_output_dir, join("_", "batch", $batch_num . ".populations.log"));
		$batch_sumstats_outfile = join('/', $populations_vcf_output_dir, join("_", "batch", $batch_num . ".sumstats.tsv"));
		$batch_sumstats_summary_outfile = join('/', $populations_vcf_output_dir, join("_", "batch", $batch_num . ".sumstats_summary.tsv"));
		
		$batch_files{$batch_vcf_infile} = $batch_vcf_outfile;
		$batch_files{$batch_hapstats_infile} = $batch_hapstats_outfile;
		$batch_files{$batch_haplotypes_infile} = $batch_haplotypes_outfile;
		$batch_files{$batch_populations_log_infile} = $batch_populations_log_outfile;
		$batch_files{$batch_sumstats_infile} = $batch_sumstats_outfile;
		$batch_files{$batch_sumstats_summary_infile} = $batch_sumstats_summary_outfile;
		
		foreach my $batch_infile (sort keys %batch_files){
			my $batch_outfile = $batch_files{$batch_infile};
			warn "Copying $batch_infile to $batch_outfile.....";
			copy($batch_infile, $batch_outfile) or die "Copy failed: $!";
			warn "Unlinking $batch_infile.....";
	 		unlink($batch_infile) or die "Could not unlink $batch_infile: $!";
		}
	}
	
}
