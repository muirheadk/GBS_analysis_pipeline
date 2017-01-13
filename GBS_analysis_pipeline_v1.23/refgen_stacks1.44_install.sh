#!/bin/bash

software_dir="$HOME/software"
mkdir -p $software_dir

cd $software_dir

bwa_dir="$software_dir/bwa-0.7.15"
samtools_dir="$software_dir/samtools-1.3.1"
stacks_dir="$software_dir/stacks-1.44"

wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
tar xvjf bwa-0.7.15.tar.bz2
cd $bwa_dir
make

wget https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2
tar xvjf samtools-1.3.1.tar.bz2
cd $samtools_dir
./configure --prefix="$samtools_dir"
make
make install

wget http://catchenlab.life.illinois.edu/stacks/source/stacks-1.44.tar.gz
tar xvzf stacks-1.44.tar.gz
cd $stacks_dir
./configure --prefix="$stacks_dir"
make
make install
