#!/usr/bin/perl
use warnings;
use strict;

my $adaptor_blastn_dir = "/home/cookeadmin/workspace/TRIMMED_OFFSET_5_ADAPTOR_FASTQ_DIR/CHRISTIANNE_MCDONALD/ADAPTOR_BLASTN_FILES";
my ($adaptor_blastn_files, $adaptor_blastn_file_counter) = find_files($adaptor_blastn_dir, "gbs_adaptor_blastn.tsv");
foreach my $adaptor_blastn_filename (sort keys %{$adaptor_blastn_files}){
	my $adaptor_blastn_infile = $fastq_files->{$fastq_filename};
	my $fasta_filename = fileparse($adaptor_blastn_infile, qr/\.fastq/);
	my ($individual_id, $barcode, $plate_num, $well_num) = split(/_/, $fasta_filename);
	# Create output directory if it doesn't already exist.
	my $individual_target_output_dir = join('/', $fasta_target_output_dir, $individual_id);
	unless(-d $individual_target_output_dir){
		mkdir($individual_target_output_dir, 0777) or die "Can't make directory: $!";
	}
	my $fasta_target_outfile = join("/", $individual_target_output_dir, join("_", $individual_id, "target") . ".fasta");
	warn "Processing " . $adaptor_blastn_infile . ".....\n";
	open(OUTFILE, ">$fasta_target_outfile") or die "Couldn't open file $fasta_target_outfile for writting, $!";
	open(INFILE, "<$adaptor_blastn_infile") or die "Couldn't open file $adaptor_blastn_infile for reading, $!";
	while(<INFILE>){
		chomp $_;
		#warn $_ . "\n";
	}
	close(INFILE) or die "Couldn't close file $adaptor_blastn_infile";
}
close(OUTFILE) or die "Couldn't close file $fasta_target_outfile";

# my $adaptor_regex_dir = "/home/cookeadmin/workspace/GBS_data-08-10-2013/TRIMMED_ADAPTOR_FASTQ_DIR/CHRISTIANNE_MCDONALD/ADAPTOR_BLASTN_FILES";
# my ($adaptor_regex_files, $adaptor_regex_file_counter) = find_files($adaptor_regex_dir, "gbs_adaptor_blastn.tsv");

sub find_files{
    
	my $infile_dir = shift;
	die "Error lost input file directory" unless defined $infile_dir;
	
	my $suffix = shift;
	die "Error lost file extension suffix directory" unless defined $suffix;
	
	my %files = ();
	my $file_counter = 0;	
	opendir(DIR, $infile_dir) || die "Error in opening dir $infile_dir\n";
	while( my $file_name = readdir(DIR)){
		my $infile_name = join('/', $infile_dir, $file_name) if ($file_name =~ m/\.$suffix$/);
		warn "$infile_name\n" if ($file_name =~ m/\.$suffix$/);
		$files{$file_name} = $infile_name if ($file_name =~ m/\.$suffix$/);
		$file_counter++ if ($file_name =~ m/\.$suffix$/);
	}
	closedir(DIR);
	return (\%files, $file_counter);
}
