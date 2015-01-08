#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);

use Bio::SeqIO;
use File::Basename;
use feature qw(switch);

# populations -P ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_MALE_GBS_ANALYSIS_TRIMMED_OFFSET_3/STACKS_OUTFILES -M ~/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/mpb_poputlations.txt -b 1 -k -p 10 -r 0.75 -m 5 -a 0.1 -f p_value --p_value_cutoff 0.05 -t 7 --fasta

my ($stacks_dir, $pop_map_infile, $batch_num, $min_num_pops_locus, $min_percent_indvs_pop,
$min_stack_depth, $min_allele_freq, $corr_value_fst, $p_value_cutoff, $filt_loci_lnl_lim, $enable_fstat,
$enable_kernel_smoothed_calcs, $kernel_smoothed_window_size, $output_option, $num_cpu_cores);
GetOptions(
	'P=s'    => \$stacks_dir,
	'M=s'    => \$pop_map_infile,
	'b=s'    => \$batch_num,
	'p=s'    => \$min_num_pops_locus,
	'r=s'    => \$min_percent_indvs_pop,
	'm=s'    => \$min_stack_depth,
	'a=s'    => \$min_allele_freq,
	'f=s'    => \$corr_value_fst,
	'p_value_cutoff=s'    => \$p_value_cutoff,
	'lnl_lim=s'    => \$filt_loci_lnl_lim,
	'fstats=s'    => \$enable_fstat,
	'k=s'    => \$enable_kernel_smoothed_calcs,
	'window_size=s'    => \$kernel_smoothed_window_size,
	'n=s'    => \$output_option,
	't=s'    => \$num_cpu_cores,
);

usage() unless (
    defined $stacks_dir
    and defined $pop_map_infile
);

$batch_num = 1 unless defined $batch_num;

$min_num_pops_locus = 10 unless defined $min_num_pops_locus;

$min_percent_indvs_pop = 0.75 unless defined $min_percent_indvs_pop;

$min_stack_depth = 5 unless defined $min_stack_depth;

$min_allele_freq = 0.1 unless defined $min_allele_freq;

$corr_value_fst = 'p_value' unless defined $corr_value_fst;

$p_value_cutoff = 0.05 unless defined $p_value_cutoff;

$filt_loci_lnl_lim = 0.0 unless defined $filt_loci_lnl_lim;

$enable_fstat = 'true' unless defined $enable_fstat;

$enable_kernel_smoothed_calcs = 'true' unless defined $enable_kernel_smoothed_calcs;

$output_option = 'fasta' unless defined $output_option;

$num_cpu_cores = 2 unless defined $num_cpu_cores;

# Program dependencies - The stacks populations program.
my $populations			= '/usr/local/bin/populations';

sub usage {

die <<"USAGE";


Usage: $0 -P stacks_dir -M pop_map_infile -b batch_num -p min_num_pops_locus -r min_percent_indvs_pop -m min_stack_depth -a min_allele_freq -f corr_value_fst --p_value_cutoff p_value_cutoff --lnl_lim filt_loci_lnl_lim --fstats enable_fstat -k enable_kernel_smoothed_calcs --window_size kernel_smoothed_window_size -n output_option -t num_cpu_cores

DESCRIPTION - 

OPTIONS:
    
	-P stacks_dir -
	-M pop_map_infile -
	-b batch_num -
	-p min_num_pops_locus -
	-r min_percent_indvs_pop -
	-m min_stack_depth -
	-a min_allele_freq -
	-f corr_value_fst -
	--p_value_cutoff p_value_cutoff -
	--lnl_lim filt_loci_lnl_lim -
	--fstats enable_fstat -
	-k enable_kernel_smoothed_calcs -
	--window_size kernel_smoothed_window_size -
	-n output_option -
	-t num_cpu_cores -

USAGE
}

given($output_option) {
    when("fasta") {
#        warn "Executing populations\n";
        my $batch_fasta_outfile = join('/', $stacks_dir, join("_", "batch", $batch_num . ".fa"));
        unless(-s $batch_fasta_outfile){
            warn "Executing populations fasta output.....\n\n";
            my $populationsCmd  = "$populations -P $stacks_dir -M $pop_map_infile -b $batch_num -k -p $min_num_pops_locus -r $min_percent_indvs_pop -m $min_stack_depth -a $min_allele_freq -f p_value --p_value_cutoff $p_value_cutoff -t $num_cpu_cores --fasta";
            warn $populationsCmd . "\n\n";
            system($populationsCmd) == 0 or die "Error calling $populationsCmd: $?";
        }
    }
    default {
        die "Error: $output_option is invalid.";
    }
}
