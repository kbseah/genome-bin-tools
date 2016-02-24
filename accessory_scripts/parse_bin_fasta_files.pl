#!/usr/bin/env perl

# Given a list of fasta files representing genome bins from a metagenome
# Parse the Fasta headers to make a table of which contig belongs to which bin (allowing overlaps)
# And this table can be imported to gbtools to automatically create bins

use strict;
use warnings;
use Getopt::Long;

my $input_table = "";
my $output_file = "";
my %fasta_hash;
my %contig_HoA;

GetOptions (
    "input|i=s" =>\$input_table,
    "output|o=s" =>\$output_file,
);

## MAIN #######################################################################

if ($input_table eq "") {
    usage();
}
else {
    read_table();
    read_fasta_write_table();
}

## SUBROUTINES ################################################################

sub usage {
    print STDERR "\nParse list of Fasta files and bin names \n";
    print STDERR "\nGiven a set of Fasta files, each representing a separate \n";
    print STDERR "genome bin, parse this into a table listing each contig in \n";
    print STDERR "the metagenome, and the corresponding bin it belongs to. \n";
    print STDERR "This table can then be imported by gbtools to automatically\n";
    print STDERR "create gbtbin objects\n\n";
    print STDERR "Usage:\n";
    print STDERR "\t \$ perl parse_bin_fasta_files.pl -i input_table\n\n";
    print STDERR "Where input_table is a tab-separated text file with Fasta file\n";
    print STDERR "paths in the first column, and the bin name in second column.\n\n";
    exit;
}

sub read_table {
    open(INPUT, "<", $input_table) or die ("$!\n");
    while (<INPUT>) {
        chomp;
        my @linesplit = split "\t", $_;
        if (!exists $fasta_hash{$linesplit[1]} ) {
            $fasta_hash{$linesplit[1]} = $linesplit[0];
        }
        else {
            print STDERR "Warning: Bin with name $linesplit[1] appears twice! Skipping... \n";
        }
    }
    close (INPUT);
}

sub read_fasta_write_table {
    open(THEOUT, "> $output_file") or *THEOUT = *STDOUT;
        # If output file not specified, write to STDOUT
    foreach my $current_bin (keys %fasta_hash) {
        print STDERR "Reading bin $current_bin\n";
        open(FASTA, "<", $fasta_hash{$current_bin}) or die ("$!\n");
        while (<FASTA>) {
            chomp;
            if ($_ =~ m/^>/) {
                my $contig_header = substr $_, 1;
                print THEOUT $current_bin."\t".$contig_header."\n";
            }
        }
        close(FASTA);
    }
    close (THEOUT);
}
