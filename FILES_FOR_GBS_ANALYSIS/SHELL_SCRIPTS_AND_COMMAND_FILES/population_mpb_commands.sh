#!/bin/sh
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_no_trim_refgen_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_no_trim_refgen_populations
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_trim_offset_3_refgen_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_trim_offset_3_refgen_populations
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_trim_offset_5_refgen_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_MALE_GBS_trim_offset_5_refgen_populations
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_no_trim_refgen_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_no_trim_refgen_populations
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_trim_offset_3_refgen_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_trim_offset_3_refgen_populations
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_trim_offset_5_refgen_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_FEMALE_GBS_trim_offset_5_refgen_populations
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_notrim_unique_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_no_trim_unique_populations
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_trim_offset_3_unique_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_trim_offset_3_unique_populations
perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_trim_offset_5_unique_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_trim_offset_5_unique_populations

