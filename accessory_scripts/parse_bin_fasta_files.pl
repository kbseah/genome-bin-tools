#!/usr/bin/env perl

=head1 NAME

parse_bin_fasta_files.pl - Parse list of Fasta files and bin names for gbtools

=head1 SYNOPSIS

perl parse_bin_fasta_files.pl -i <input> -o <output>

perl parse_bin_fasta_files.pl --help

=head1 DESCRIPTION

Parse headers from Fasta files corresponding to genome bins (e.g. output of
MetaBat or MetaWatt pipelines) to a table for import to gbtools.

For more information, refer to gbtools documentation.

Part of the gbtools package by Brandon Seah:
https://github.com/kbseah/genome-bin-tools/

=head1 ARGUMENTS

=over 8

=item --input|-i <file>

Tab-separated file. Column 1 is path to Fasta files that correspond to genome
bins, and column 2 is shortname for the genome bin.

=item --output|-o <file>

Name for output file. 

=back

=head1 OUTPUT

Tab-separated file, each line with name of contig and name of corresponding
genome bin.

=head1 COPYRIGHT AND LICENSE

gbtools - Interactive tools for metagenome visualization and binning in R
Copyright (C) 2015,2016  Brandon Seah (kbseah@mpi-bremen.de)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $input_table = "";
my $output_file = "";
my %fasta_hash;
my %contig_HoA;

if (@ARGV == 0 ) {
    pod2usage(-message => "Insufficient options were supplied", -existatus => 2);
}

GetOptions (
    "input|i=s" =>\$input_table,
    "output|o=s" =>\$output_file,
    'help|h' => sub { pod2usage( -exitstatus => 2, -verbose => 2); },
    'man|m'=> sub { pod2usage ( -exitstatus => 0, -verbose => 2) }
);

## MAIN #######################################################################

if ($input_table eq "") {
    pod2usage(-message => "Input table not supplied", -exitstatus => 2, -verbose => 2);
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
