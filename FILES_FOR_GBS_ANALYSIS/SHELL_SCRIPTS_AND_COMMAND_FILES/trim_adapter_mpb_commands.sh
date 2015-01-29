#!/bin/sh

mkdir -p ~/GBS_analysis-2015-01-13/Stephen_MPB

for i in {1..10}
do 

echo perl ~/GBS_analysis-2015-01-13/trim_adapter_fastq_parallel_regex.pl -i ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_PROCESSED_RADTAGS -p TRIMMED_OFFSET_$((i))_ADAPTER_DIR -t $i -q 32 -m 16 -l 92 -r PstI/MspI -c 10 -n true -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_TRIMMED_ADAPTER_PADDED_DIR
perl ~/GBS_analysis-2015-01-13/trim_adapter_fastq_parallel_regex.pl -i ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_PROCESSED_RADTAGS -p TRIMMED_OFFSET_$((i))_ADAPTER_DIR -t $i -q 32 -m 16 -l 92 -r PstI/MspI -c 10 -n true -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_TRIMMED_ADAPTER_PADDED_DIR

echo perl ~/GBS_analysis-2015-01-13/trim_adapter_fastq_parallel_regex.pl -i ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_PROCESSED_RADTAGS -p TRIMMED_OFFSET_$((i))_ADAPTER_DIR -t $i -q 32 -m 16 -l 92 -r PstI/MspI -c 10 -n false -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_TRIMMED_ADAPTER_UNPADDED_DIR
perl ~/GBS_analysis-2015-01-13/trim_adapter_fastq_parallel_regex.pl -i ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_PROCESSED_RADTAGS -p TRIMMED_OFFSET_$((i))_ADAPTER_DIR -t $i -q 32 -m 16 -l 92 -r PstI/MspI -c 10 -n false -o ~/GBS_analysis-2015-01-13/Stephen_MPB/MPB_GBS_TRIMMED_ADAPTER_UNPADDED_DIR

done
