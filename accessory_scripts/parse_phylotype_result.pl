#!/usr/bin/env perl

=head1 NAME

parse_phylotype_result.pl - Parse results from Amphora2 or Phyla-Amphora pipelines for gbtools

=head1 SYNOPSIS

perl parse_phylotype_result.pl -p amphora2_phylotype > phylotype.parsed

perl parse_phylotype_result.pl --help

=head1 DESCRIPTION

Parse phylotype result table from Amphora2 or Phyla-Amphora pipelines to format
for import to gbtools in R. 

For more information, refer to gbtools documentation.

Part of the gbtools package by Brandon Seah:
https://github.com/kbseah/genome-bin-tools/

=head1 ARGUMENTS

=over 8

=item --phylotypes|-p I<FILE>

Result of Phylotyping.pl script from Amphora2 or Phyla-Amphora pipelines.

=back

=head1 OUTPUT

Output is written to STDOUT.

=head1 COPYRIGHT AND LICENSE

gbtools - Interactive tools for metagenome visualization and binning in R
Copyright (C) 2015-2018  Brandon Seah (kbseah@mpi-bremen.de)

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
    'help|h' => sub { pod2usage( -exitstatus => 2, -verbose => 1); },
    'man|m'=> sub { pod2usage ( -exitstatus => 0, -verbose => 2) }
) or pod2usage(-verbose=>0);

parse_phylotyping_result();

## SUBROUTINES ################################################################

sub parse_phylotyping_result {
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
