#!/bin/bash

software_dir="$HOME/software"
mkdir -p $software_dir

# Download the github repository for the GBS_analysis_pipeline
cd $software_dir
git clone https://github.com/muirheadk/GBS_analysis_pipeline.git

pipeline_dir="$software_dir/GBS_analysis_pipeline"

echo "export PATH=$PATH:$pipeline_dir/GBS_analysis_pipeline_v1.25" >> $HOME/.bashrc

# Install BWA version 0.7.15
bwa_dir="$software_dir/bwa-0.7.15"
cd $software_dir
wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
tar xvjf bwa-0.7.15.tar.bz2
cd $bwa_dir
make

echo "export PATH=$PATH:$bwa_dir" >> $HOME/.bashrc

# Install Samtools version 1.3.1
samtools_dir="$software_dir/samtools-1.3.1"
cd $software_dir
wget https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2
tar xvjf samtools-1.3.1.tar.bz2
cd $samtools_dir
./configure --prefix="$samtools_dir"
make
make install

echo "export PATH=$PATH:$samtools_dir/bin" >> $HOME/.bashrc

# Install Stacks Version 1.44
stacks_dir="$software_dir/stacks-1.44"
cd $software_dir
wget http://catchenlab.life.illinois.edu/stacks/source/stacks-1.44.tar.gz
tar xvzf stacks-1.44.tar.gz
cd $stacks_dir
./configure --prefix="$stacks_dir"
make
make install

echo "export PATH=$PATH:$stacks_dir/bin" >> $HOME/.bashrc

# Install Perl Modules from CPAN

perl_module_dir="$HOME/perl5/lib/perl5"

echo "export PERL5LIB=$PERLLIB:\"$perl_module_dir\"" >> $HOME/.bashrc

source $HOME/.bashrc

cpan install Parallel::Loops


