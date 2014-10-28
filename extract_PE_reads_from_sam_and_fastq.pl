#!/usr/bin/perl

## Extract Illumina reads that map to target in SAM file and extract dangling ends from Fastq file for reassembly
## Contact: kbseah@mpi-bremen.de
## Version 2 -- KBS 2014-10-23 -- Added output reporting, option to specify length of read header suffixes, changed indents to spaces
## Version 1 -- KBS 2014-09-23 -- Fixed Fastq header recognition
## Version 0 -- KBS 2014-03-18

## HEAD #############################################################################################################

use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

my $targetlistfile;		# List of scaffolds for which to extract the mapped reads (thereof)
my $samfile;			# SAM file of the mapping from the readfile to the scaffolds contained in the targets list
my $readfile;			# Fastq-formatted reads
my $is_gzip = 0;		# Flag for gzipped read input
my $suffix_length=0;		# Length of suffix (in characters) appended by user to the Illumina read headers, e.g. _1 and _2 to distinguish forward and reverse reads
my $targetscaffolds_counter=0;
my $targetheaders_counter=0;
my $scannedreads_counter=0;
my $foundreads_counter=0;
my %targets_hash;
my %reads_hash;

## MAIN #############################################################################################################

if (@ARGV == 0) { usage(); }

GetOptions ("target|t=s" => \$targetlistfile,
    "sam|s=s" => \$samfile,
    "reads|r=s" => \$readfile,
    "gzipped|g" => \$is_gzip,
    "suffix|x=i" => \$suffix_length)
or die("Error in command line arguments\n");

print STDERR "Reading list of target scaffolds\n";
read_target_file();
print STDERR "Found $targetscaffolds_counter target scaffolds in list\n";
print STDERR "Getting list of reads that map to targets from SAM file\n";
scan_sam_file();
print STDERR "Found $targetheaders_counter read headers to be extracted\n";
if ($is_gzip == 1) {
    print STDERR "The input reads are gzipped (option -g|--gzipped was specified)\n";
    filter_fq_gz_file();
} 
else { filter_fq_file(); }

## SUBROUTINES ######################################################################################################

sub usage {
    print STDERR "\n";
    print STDERR "Read extractor for reassembly\n";
    print STDERR "Version 2 - Brandon Seah - 2014-10-23 \n";
    print STDERR " \t Find reads that map to specified targets in SAM file, and then extract dangling ends from Fastq file for reassembly \n";
    print STDERR " \t Assumes that Fastq headers are produced by Illumina machine and start with \@HWI\n";  
    print STDERR "\n";
    print STDERR "Usage: \n";
    print STDERR " \$ perl extract_PE_reads_from_sam_and_fastq.pl [options] -t <target_list> -s <sam_mapping> -r <reads.fastq>\n";
    print STDERR "\n";
    print STDERR "Options: \n";
    print STDERR " \t --target|-t FILE      Plain text list of target scaffolds (one scaffold name per line, no \"\>\") for which to extract the mapped reads \n";
    print STDERR " \t --sam|-s    FILE      SAM file with mapping of reads to scaffolds\n";
    print STDERR " \t --reads|-r  FILE      Fastq or gzipped Fastq formatted reads from which to extract those which map to the target scaffolds\n";
    print STDERR " \t --gzipped|-g          Input reads are gzipped (default: no)\n";
    print STDERR " \t --suffix|-x INT       Input read header names have a suffix added of length INT (e.g. \_1 or \_2) (default = 0)\n";
    print STDERR " \n";
    print STDERR "Use hyphen if streaming to STDIN, e.g.: \n";
    print STDERR " \$ samtools view BAMFILE | perl extract_PE_reads_from_sam_and_fastq.pl --target TARGET --sam - --reads READS | gzip > output.fq.gz \n";
    print STDERR "\n";
    exit;
}

sub read_target_file {  # Read file containing list of scaffolds to extract the reads that map thereto
open (TARGETSINPUT, "< $targetlistfile") or die ("Cannot open target list file: $!\n");
while (<TARGETSINPUT>) {
    chomp;
    if (!exists $targets_hash{$_}) {
        $targets_hash{$_} = 1;
        $targetscaffolds_counter++;
    }
}
close (TARGETSINPUT);
}

sub scan_sam_file {     # Scan SAM file of mapping to find header names of reads that should be extracted
    open (SAMINPUT, "< $samfile") or die ("Cannot open SAM file: $!\n");
    while (<SAMINPUT>) {
        chomp;
        if ($_ =~ /^@/) { next; }
        my @theline = split('\t',$_);
        if ( exists $targets_hash{$theline[2]} ) {
            my $header_nosuffix = substr $theline[0], 0, 0-$suffix_length;      # Chops off the last $suffix_length characters from the read header name
            $reads_hash{$header_nosuffix} = 1;
            $targetheaders_counter++;
        }
    }
    close (SAMINPUT);
}

sub filter_fq_file {    # Filter input Fastq file to output reads that map to the target scaffolds
my $writeflag = 0;
open (FASTQ, "< $readfile") or die ("Cannot open Fastq file: $!\n");
while (<FASTQ>) {
    chomp;
    if ($_ =~ /^\@HWI/) {    # Important! @ character alone will also match Phred quality score symbol and will cause quality lines to be misread as header lines. @HWI assumes Illumina output
        $scannedreads_counter++;
        if ($scannedreads_counter > 0 && $scannedreads_counter%500000==0) {     # Output a progress message every 500000 reads
            print STDERR "Scanned $scannedreads_counter reads from input Fastq file\n";
        }
        my $readname = substr $_, 1, 0-$suffix_length;  # Remove any suffix characters that have been added to the read headers, and the @ character at beginning of Fastq header line
        if (exists $reads_hash{$readname}) {
            $writeflag = 1;
            print STDOUT $_, "\n";
            $foundreads_counter++;
        }
        else { $writeflag = 0; }
    }
    elsif ( $writeflag == 1) { print STDOUT $_, "\n"; }
    else { next; }
}
close(FASTQ);
}

sub filter_fq_gz_file { # Filter input gzipped Fastq file
my $writeflag = 0;
my $theinput = IO::Uncompress::Gunzip->new( $readfile ) or die "Cannot open Fastq.gzip file: $GunzipError\n";
while (<$theinput>) {
    chomp;
    if ($_ =~ /^\@HWI/) {
        $scannedreads_counter++;
        if ($scannedreads_counter>0 && $scannedreads_counter%500000==0) {
            print STDERR "Scanned $scannedreads_counter reads from input Fastq.gz file\n";
        }
        my $readname = substr $_, 1, 0-$suffix_length;  # Remove any suffix characters that have been added to the read headers, and the @ character at beginning of Fastq header line
        if (exists $reads_hash{$readname}) {
            $writeflag = 1;
            print STDOUT $_, "\n";
            $foundreads_counter++;
        }
        else { $writeflag = 0; }
    }
    elsif ( $writeflag == 1) { print STDOUT $_, "\n"; }
    else { next; }
}
close($theinput);
}

#sub gunzip_example {
#    my $inputfile = $ARGV[0];    # specify file to be unzipped
#    my $counter = 0;
#    my $theinput = IO::Uncompress::Gunzip->new( $inputfile ) or die "Gunzip failed: $GunzipError\n";
#    while (<$theinput>) {
#        chomp;
#        my @theline = split('\t',$_);
#        if ($theline[1] == 161) { $counter = $counter + 1; }
#    }
#    print STDOUT $counter, "\n";
#}