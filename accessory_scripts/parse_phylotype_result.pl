#!/usr/bin/env perl

=head1 NAME

parse_phylotype_result.pl - Parse results from Amphora2 or Phyla-Amphora pipelines for gbtools

=head1 SYNOPSIS

perl parse_phylotype_result.pl -p <amphora2_phylotype> > phylotype.parsed

perl parse_phylotype_result.pl --help

=head1 DESCRIPTION

Parse phylotype result table from Amphora2 or Phyla-Amphora pipelines to format
for import to gbtools in R. 

For more information, refer to gbtools documentation.

Part of the gbtools package by Brandon Seah:
https://github.com/kbseah/genome-bin-tools/

=head1 ARGUMENTS

=over 8

=item --phylotypes|-p <file>

Result of Phylotyping.pl script from Amphora2 or Phyla-Amphora pipelines.

=back

=head1 OUTPUT

Output is written to STDOUT.

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


## Script to parse AMPHORA2 or Phyla-AMPHORA marker phylotyping info to a format that can be imported to R
## after the style of Albertsen et al.

## 2016-03-09 - Documentation to POD
## 2015-06-13 - Tidied up formatting and comment lines
## Version 2 - 2014-10-27 - Full taxon string parsed for export; no longer extracting only a single taxon level
## Version 1 - 2014-10-22
## Contact: kbseah@mpi-bremen.de

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $phylotyping_result;   # File with results of AMPHORA2 or Phyla-AMPHORA Phylotyping.pl results
my %marker_name_hash;     # Hash for gene name, marker ID as key
my %marker_taxon_hash;    # Hash for taxon (default Class level) from phylotyping result, marker ID as key
my %marker_scaffold_hash; # Hash for scaffold containing a given marker gene, marker ID as key
my $taxon_level=4;        # Which taxonomic level should we parse the output? 1=Domain, 2=Phylum, 3=Class, 4=Order, 5=Family, 6=Genus, 7=Species

## MAIN #######################################################################

if (@ARGV == 0) {
    pod2usage(-message => "Insufficient options were supplied", -existatus => 2);
}

GetOptions (
    "phylotypes|p=s" => \$phylotyping_result,
    "level|l=i" => \$taxon_level,
    'help|h' => sub { pod2usage( -exitstatus => 2, -verbose => 2); },
    'man|m'=> sub { pod2usage ( -exitstatus => 0, -verbose => 2) }
);

#if ($taxon_level > 7) {    # Catch smart-asses
#   die "Taxon level cannot be lower than species!\n";
#}
#parse_phylotyping_result();
parse_phylotyping_result_new();

## SUBROUTINES ################################################################

sub usage {
    print STDERR "\n";
    print STDERR "Parse results from Phyla-Amphora or Amphora2 Phylotyping.pl\n";
    print STDERR "script for import to R\n";
    print STDERR "\n";
    print STDERR "Usage: \n";
    print STDERR " \$ perl parse_phylotype_result.pl -p <results_file>\n";
    print STDERR "\n";
    print STDERR "Options:\n";
    print STDERR " \t -p FILE   Results from the script Phylotyping.pl\n";
    print STDERR "\n";
    print STDERR "Output:\n";
    print STDERR " \t <results_file>.parsed     Suitable for import to gbtools\n";
    print STDERR "\n";
    exit;
}

sub parse_phylotyping_result_new {
    open(PHYLOTYPING, "< $phylotyping_result")
        or die ("Cannot open file $phylotyping_result: $!\n");
    my $discardfirstline = <PHYLOTYPING>;   # Throw away header line
    print STDOUT join("\t",
                      "scaffold",
                      "markerid",
                      "gene",
                      "Superkingdom",
                      "Phylum",
                      "Class",
                      "Order",
                      "Family",
                      "Genus",
                      "Species"
                      ),
                 "\n";
    while (<PHYLOTYPING>) {
        chomp;
        my @currentline= split "\t",$_;
        ## Splitting and popping to get scaffold name from marker ID, by ######
        ## removing ID number tacked on by getorf #############################
        my @temparray = split "_", $currentline[0];
        my $discard = pop @temparray;
        ## Save scaffold containing marker
        $marker_scaffold_hash{$currentline[0]} = join "_", @temparray;
        if (scalar @currentline < 9) {
            my $num_to_add = 9 - (scalar @currentline);
            my $last_string;
            ($last_string) = ($currentline[$#currentline] =~ /(.*)\([\d\.]+\)/);
            $last_string = "\(".$last_string."\)";
            ## Fill in blank taxon levels with the lowest assigned taxon name #
            while ($num_to_add > 0) {
                push @currentline, $last_string;
                $num_to_add--;
            }
        }
        ## Strip confidence levels from taxon names ###########################
        for (my $i=2; $i<(scalar @currentline); $i++) {
                if ($currentline[$i] =~ /(.*)\([\d\.]+\)/) {
                    $currentline[$i] = $1;
                }
            }
        print STDOUT join("\t",
                          $marker_scaffold_hash{$currentline[0]},
                          @currentline),
                     "\n";
    }    
    close (PHYLOTYPING);
}

sub parse_phylotyping_result {
    my $cut_level = $taxon_level + 1;
    open(PHYLOTYPING, "< $phylotyping_result")
        or die ("Cannot open phylotyping results file: $! \n");
    my $discardheader = <PHYLOTYPING>;
    while (<PHYLOTYPING>) {
        my @currentline = split "\t", $_;
        my $currentmarker = $currentline[0];                            # Get current marker ID
        $marker_name_hash{$currentmarker} = $currentline[1];            # Save name of marker gene
        my @temparray = split "_", $currentmarker;                      # Splitting and popping to get scaffold name from marker ID, by removing ID number tacked on by getorf
        my $discard = pop @temparray;
        $marker_scaffold_hash{$currentmarker} = join "_", @temparray;   # Save scaffold containing marker
        if ( scalar(@currentline) >= $cut_level+2 ) {                   # If this marker has been assigned at least to level of $taxon_level, extract the taxon name at this level
            my @temparray2 = split /\(/, $currentline[$cut_level];
            $marker_taxon_hash{$currentmarker} = $temparray2[0];        # Save taxon assignment of marker
        }
        elsif ( scalar(@currentline) < $cut_level+2 ) {                  # Otherwise use the next-highest taxonomic level as the taxon name
            my $pos = scalar(@currentline);
            my @temparray3 = split /\(/, $currentline[$pos-1];
            $marker_taxon_hash{$currentmarker} = $temparray3[0];        # Save lowest-level taxon-assignment of marker
        }
        
    }
    close(PHYLOTYPING);
    open(PHYLOTYPINGOUT, "> $phylotyping_result\.parsed")
        or die ("Cannot open $phylotyping_result\.parsed for writing: $!\n"); # Write file containing parsed marker details
    print PHYLOTYPINGOUT "markerid",  # Header line
                         "\t",
                         "scaffold",
                         "\t",
                         "gene",
                         "\t",
                         "taxon",
                         "\n";
    foreach my $currentmarker (keys %marker_name_hash) {
        print PHYLOTYPINGOUT $currentmarker,
                             "\t",
                             $marker_scaffold_hash{$currentmarker},
                             "\t",
                             $marker_name_hash{$currentmarker},
                             "\t",
                             $marker_taxon_hash{$currentmarker},
                             "\n";
    }
    close (PHYLOTYPINGOUT);
}

