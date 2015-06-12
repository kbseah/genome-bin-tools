#!/usr/bin/env perl

## Version 1 - 2014-10-22

## Contact: kbseah@mpi-bremen.de

## Script to parse coverage info and AMPHORA2 marker phylotyping info and generate pretty plots for evaluating assembly
## after the style of Albertsen et al.

use strict;
use warnings;
use Getopt::Long;
use FindBin;

my $phylotyping_result;             # File with results of AMPHORA2 or Phyla-AMPHORA Phylotyping.pl results
my %marker_name_hash;               # Hash for gene name, marker ID as key
my %marker_taxon_hash;              # Hash for taxon (default Class level) from phylotyping result, marker ID as key
my %marker_scaffold_hash;           # Hash for scaffold containing a given marker gene, marker ID as key
my $coverage_result;                # Output of pileup.sh from BBMap suite, giving coverage, GC, and length statistics for each contig
my $assem_name="Assembly";          # Name of the assembly (for plot title)
my $dir = $FindBin::Bin;            # Folder where this script is run from: assume that assembly_parse_plot.r is in the same folder
my $taxon_level=4;                  # Which taxonomic level should we parse the output? 1=Domain, 2=Phylum, 3=Class, 4=Order, 5=Family, 6=Genus, 7=Species

## Options

if (! @ARGV) { usage(); }
GetOptions (
    "phylotypes|p=s" => \$phylotyping_result,  
    "coverage|c=s" => \$coverage_result,
    "assembly|a=s" => \$assem_name,
    "taxlevel|t=i" => \$taxon_level
);

## Main

if ($taxon_level>7) { die "Taxon level cannot be lower than species!\n"; }  # Catch smart-asses
parse_phylotyping_result();
generate_plots();

## Subroutines

sub usage {
    print STDOUT "Usage: perl assembly_parse_plot.pl -c <pileup_result> -p <phylotype_result> -a <assembly_name> -t <taxon_level>\n";
    print STDOUT "  -c     Output file from pileup.sh, must have been generated wtih fasta reference\n";
    print STDOUT "  -p     Output from AMPHORA2 or Phyla-AMPHORA Phylotyping.pl script\n";
    print STDOUT "  -a     Name of your assembly, for plot titles\n";
    print STDOUT "  -t     Taxonomic level for coloring the marker genes, default = 4 (Order)\n";
    exit;
}

sub parse_phylotyping_result {
    my $cut_level = $taxon_level+1;
    open(PHYLOTYPING, "< $phylotyping_result") or die ("Cannot open phylotyping results file: $! \n");
    my $discardheader = <PHYLOTYPING>;
    while (<PHYLOTYPING>) {
        my @currentline = split "\t", $_;
        my $currentmarker = $currentline[0];                            # Get current marker ID
        $marker_name_hash{$currentmarker} = $currentline[1];            # Save name of marker gene
        my @temparray = split "_", $currentmarker;                      # Splitting and popping to get scaffold name from marker ID, by removing ID number tacked on by getorf
        my $discard = pop @temparray;
        $marker_scaffold_hash{$currentmarker} = join "_", @temparray;   # Save scaffold containing marker
        if ( scalar(@currentline) >= $cut_level+2 ) {                               # If this marker has been assigned at least to level of Class, extract the taxonomic class as the taxon name
            my @temparray2 = split /\(/, $currentline[$cut_level];
            $marker_taxon_hash{$currentmarker} = $temparray2[0];        # Save class-level taxon assignment of marker
        }
        elsif ( scalar(@currentline) < $cut_level+2 ) {                             # Otherwise use the next-highest taxonomic level as the taxon name
            my $pos = scalar(@currentline);
            my @temparray3 = split /\(/, $currentline[$pos-1];
            $marker_taxon_hash{$currentmarker} = $temparray3[0];        # Save lowest-level taxon-assignment of marker
        }
    }
    close(PHYLOTYPING);
    open(PHYLOTYPINGOUT, "> $phylotyping_result\.parsed") or die ("Cannot open $phylotyping_result\.parsed for writing: $!\n"); # Write file containing parsed marker details
    print PHYLOTYPINGOUT "markerid", "\t", "scaffold", "\t", "gene", "\t", "taxon", "\n";                                       # Header line
    foreach my $currentmarker (keys %marker_name_hash) {
        print PHYLOTYPINGOUT $currentmarker, "\t", $marker_scaffold_hash{$currentmarker}, "\t", $marker_name_hash{$currentmarker}, "\t", $marker_taxon_hash{$currentmarker}, "\n";
    }
    close (PHYLOTYPINGOUT);
}

sub generate_plots {
    my @rscript_command = ("Rscript", "$dir/assembly_parse_plot.r", "--args", $coverage_result, "$phylotyping_result\.parsed", $assem_name);
    system(@rscript_command);
}