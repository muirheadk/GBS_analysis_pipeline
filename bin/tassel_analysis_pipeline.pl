#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -M strevoy@ualberta.ca
#PBS -l mem=12gb
#PBS -m bea
#PBS -l procs=4
cd $PBS_O_WORKDIR

#run_pipeline.pl -Xmx12g -fork1 -FastqToTagCountPlugin -i ./fastq -k barcodes.txt \
#-e PstI-MspI -s 300000000 -c 1 -o ./tagCounts -endPlugin -runfork1 > ./tagCounts.log

#run_pipeline.pl -Xmx12g -fork1 -MergeMultipleTagCountPlugin -i ./tagCounts -o ./mergedTagCounts/mpbGBSTags.cnt \
#-c 5 -t -endPlugin -runfork1 > ./MTC.log

#run_pipeline.pl -Xmx12g -fork1 -TagCountToFastqPlugin -i ./mergedTagCounts/mpbGBSTags.cnt -o \
#./mergedTagCounts/mpbGBSTags.fq -c 5 -endPlugin -runfork1 > ./MTCtoFQ.log

#bwa index -a bwtsw ./refgen/renumbered_Msequence.fasta

#bwa aln -t 4 ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/mpbGBSTags.cnt.fq > \
#./mergedTagCounts/AlignedGBSTags1.sai

#bwa samse ./refgen/renumbered_Msequence.fasta ./mergedTagCounts/AlignedGBSTags1.sai \
#./mergedTagCounts/mpbGBSTags.cnt.fq > mergedTagCounts/AlignedMasterTagsMPB.sam

#run_pipeline.pl -Xmx12g -fork1 -SAMConverterPlugin -i ./mergedTagCounts/AlignedMasterTagsMPB.sam \
#-o ./topm/MasterTagsMPB.topm -endPlugin -runfork1 > ./SAMconvert.log

#run_pipeline.pl -Xmx12g -fork1 -FastqToTBTPlugin -i ./fastq -k barcodes.txt -e PstI-MspI -o ./tbt \
#-y -m ./topm/MasterTagsMPB.topm -endPlugin -runfork1 > ./TBT.log

#run_pipeline.pl -Xmx12g -fork1 -MergeTagsByTaxaFilesPlugin -i ./tbt -o ./mergedTBT/mpbGBSstudy.tbt.byte \
#-s 200000000 -endPlugin -runfork1 > ./mergeTBT.log

#run_pipeline.pl -Xmx12g -fork1 -TagsToSNPByAlignmentPlugin -i ./mergedTBT/mpbGBSstudy.tbt.byte \
#-y -m ./topm/MasterTagsMPB.topm -mUpd ./topm/MasterTagsMPBwVariants.topm -o ./hapmap/raw/mpbGBSGenos_chr+.hmp.txt \
#-mxSites 500000 -mnMAF 0.02 -mnMAC 100000 -ref ./refgen/renumbered_Msequence.fasta -sC 1 -eC 8188 \
#-endPlugin -runfork1 > ./TagstoSNPAlign.log

#run_pipeline.pl -Xmx12g -fork1 -MergeDuplicateSNPsPlugin -hmp ./hapmap/raw/mpbGBSGenos_chr+.hmp.txt \
#-o ./hapmap/mergedSNPs/mpbGBSGenos_mergedSNPs_chr+.hmp.txt \
#-misMat 0.1 -callHets -sC 1 -eC 8188 -endPlugin -runfork1 > MergeDupSNP.log

run_pipeline.pl -Xmx12g -fork1 -GBSHapMapFiltersPlugin -hmp ./hapmap/mergedSNPs/mpbGBSGenos_mergedSNPs_chr+.hmp.txt \
-o ./hapmap/filt/mpbGBSGenos_mergedSNPsFilt_chr+.hmp.txt -mnSCov 0.2 MnMAF 0.01 \
-sC 1 -eC 8188 -endPlugin -runfork1 > ./hmpFilt.log

