#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l mem=4000mb
#PBS -M jrdupuis@ualberta.ca
#PBS -m bea
cd $PBS_O_WORKDIR

#The above lines tell westgrid to request 72 hrs of executiong time, 4gb of RAM (you'll notice that on the pipeline steps below the -Xmx4g option specifically 
#tells tassel that there are 4gb of RAM to work with), to email you when the run has started and finished or if there are problems with the run, and to begin 
#executing the following lines from the current directory.
#you can run this script by typing qsub run_uneak_forjrd.pbs on the cmd line in your bugaboo home or scratch directory, but first read the comments below.

#Below are the steps to run uneak, but first truncate your GBS datafile name from the HI.0803.005.GQ25032013-7_R1.fastq.gz formate to the 
#GQ25032013_7_fastq.txt.gz formate and place it within the Illumina directory created by the ucreatworkdir plugin.
#Next prepare your barcodes_forjrd.txt file with the following format: GQ25032013      1       CTCG    7001    plate1  A       1, where
#each column corresponds to the run name, lane, individual-specific barcode, individual id, plate # (you can keep this the same for all samples since they 
#were all run on a single plate), row position on plate, column position on plate, and each line in the barcode file corresponds to a different individual in 
#your dataset and place the barcode file in the key directory created by the ucreatworkdir plugin. All you should have to do is change the lane number to the 
#correct one and than change the individual ids to the correct ones.
#All steps below have a "> ./*.log" segment indicating that any info printed to standard output will be saved to a log file that you can open for your own 
#viewing pleasure. 

#run ufastqtotagcount plugin from the current directory (-w) using PstI-MspI combo (-e), retaining all reads with atleast 1 copy (-c). Expected number of good barcoded reads to be less than 250 million (-s). 
#don't change the -c option setting here since you'll be able to filter by read depth at the umergetaxatagcount step later on.
run_pipeline.pl -Xmx4g -fork1 -UFastqToTagCountPlugin -w ./ -e PstI-MspI -s 250000000 -c 1 -endPlugin -runfork1 > ./UtagCount.log

#run binarytotext plugin to make the ufastqtotagcount result from above human readable
#change GQ25032013_5 in the following line to the actual name of your files (should be GQ25032013_yourlanenumber 
run_pipeline.pl -Xmx4g -fork1 -BinaryToTextPlugin -i ./tagCounts/GQ25032013_7.cnt -o ./tagCounts/GQ25032013_7_cnt.txt -t TagCounts -endPlugin -runfork1 > ./bintotex.log

#run umergetaxatagcounts plugin from the current directory (-w), using a read depth cutoff of 5 (-c). You can change the cutoff to whatever you like.
run_pipeline.pl -Xmx4g -fork1 -UMergeTaxaTagCountPlugin -w ./ -c 5 -endPlugin -runfork1 > UmergeTagCount.log

#run binarytotext plugin to make the umergetaxatagcount result from above human readable
run_pipeline.pl -Xmx4g -fork1 -BinaryToTextPlugin -i ./mergedTagCounts/mergeAll.cnt -o ./mergedTagCounts/mergeAll_cnt.txt -t TagCounts -endPlugin -runfork1 >> ./bintotex.log

#run utagcounttotagpair plugin from current directory (-w), with a sequencing error threshold of 3% (this should not be higher than 5% based on known Illumina 
#sequencing error rates.
run_pipeline.pl -Xmx4g -fork1 -UTagCountToTagPairPlugin -w ./ -e 0.03 -endPlugin -runfork1 > utagtotp.log

#run utagpairtoTBTplugin from current directory (-w).
run_pipeline.pl -Xmx4g -fork1 -UTagPairToTBTPlugin -w ./ -endPlugin -runfork1 > utagtotbt.log

#run uTBTtomapinfoplugin from current directory (-w).
run_pipeline.pl -Xmx4g -fork1 -UTBTToMapInfoPlugin -w ./ -endPlugin -runfork1 > utbttomapinfo.log

#run umapinfotohapmapplugin from current directory (-w), with a minor allele frequency of 1%, a major allele frequency of 99%, minimum call rate of 0, maximum 
#callrate of 1.
run_pipeline.pl -Xmx4g -fork1 -UMapInfoToHapMapPlugin -w ./ -mnMAF 0.01 -mxMAF 0.99 -mnC 0 -mxC 1 -endPlugin -runfork1 > umapinfotohapmap.log

