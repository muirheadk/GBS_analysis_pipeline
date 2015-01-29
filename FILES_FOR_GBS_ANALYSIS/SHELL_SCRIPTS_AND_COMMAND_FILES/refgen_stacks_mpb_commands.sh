#!/bin/sh

mkdir -p ~/GBS_analysis-2015-01-13/Stephen_MPB/REFGEN_STACKS_OUFILES

for i in {1..10}
do 

echo perl ~/GBS_analysis-2015-01-13/refgen_stacks_analysis_pipeline.pl -i ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_TRIMMED_ADAPTER_UNPADDED_DIR/TRIMMED_OFFSET_$((i))_ADAPTER_DIR/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_REFERENCE_GENOMES/DendPond_male_1.0_unplaced.scaf.fa -t gzfastq -c 10 -o ~/GBS_analysis-2015-01-13/Stephen_MPB/REFGEN_STACKS_OUFILES/MPB_MALE_GBS_trim_offset_$((i))_refgen_stacks

perl ~/GBS_analysis-2015-01-13/refgen_stacks_analysis_pipeline.pl -i ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_TRIMMED_ADAPTER_UNPADDED_DIR/TRIMMED_OFFSET_$((i))_ADAPTER_DIR/TRIMMED_OUTPUT_FILES/TRIMMED_FASTQ_FILES -g ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_REFERENCE_GENOMES/DendPond_male_1.0_unplaced.scaf.fa -t gzfastq -c 10 -o ~/GBS_analysis-2015-01-13/Stephen_MPB/REFGEN_STACKS_OUFILES/MPB_MALE_GBS_trim_offset_$((i))_refgen_stacks

echo perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/REFGEN_STACKS_OUFILES/MPB_MALE_GBS_trim_offset_$((i))_refgen_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/REFGEN_STACKS_OUFILES/MPB_MALE_GBS_trim_offset_$((i))_refgen_populations

perl ~/GBS_analysis-2015-01-13/populations_stacks.pl -p 1 -r 0 -m 2 -a 0 -P ~/GBS_analysis-2015-01-13/Stephen_MPB/REFGEN_STACKS_OUFILES/MPB_MALE_GBS_trim_offset_$((i))_refgen_stacks/STACKS_OUTFILES -M ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_population_map.txt -t 10 -f fasta -o ~/GBS_analysis-2015-01-13/Stephen_MPB/REFGEN_STACKS_OUFILES/MPB_MALE_GBS_trim_offset_$((i))_refgen_populations


done

