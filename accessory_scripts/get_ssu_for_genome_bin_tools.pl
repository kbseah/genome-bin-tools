#!/usr/bin/env perl

=head1 NAME

get_ssu_for_genome_bin_tools.pl - Extract SSU data and taxonomy and parse for gbtools

=head1 SYNOPSIS

perl get_ssu_for_genome_bin_tools.pl -d path/to/db -c num_CPUs \
                                     -a scaffolds.fasta -o output_prefix

perl get_ssu_for_genome_bin_tools.pl --help

=head1 DESCRIPTION

Extracts SSU rRNA sequence from metagenome assembly in Fasta file, classifies
taxonomically by comparison to Vsearch-indexed and curated SILVA database
(produced by the PhyloFlash v1.5d+ by Harald Gruber-Vodicka) and parses into
a table that can be imported by gbtools in R.

Dependencies: Vsearch, barrnap 0.5+, bedtools 2.21.0+

For more information, refer to gbtools documentation.

Part of the gbtools package by Brandon Seah:
https://github.com/kbseah/genome-bin-tools/

=head1 ARGUMENTS

=over 8

=item --dbpath|-d I<PATH>

Path to Vsearch-indexed, curated SILVA database of SSU rRNA sequences.

=item --cpus|-c I<INTEGER>

Number of processors for barrnap and Vsearch (default: 1)

=item --assembly|-a I<FILE>

Fasta file of genome/metagenome sequences to search.

=item --output|-o I<STRING>

Prefix for output file names.

=item --help|-h

Show this help message.

=back

=head1 OUTPUT

=over 8

=item <output_prefix>.ssu.tab

Table of predicted SSU genes and their taxonomic affiliations

=item tmp.<output_prefix>.scaffolds.gff

Concatenated output from Barrnap searches

=item <output_prefix>.barrnap.ssu.gff

Predicted SSU genes, with duplicates removed, GFF feature table

=item <output_prefix>.barrnap.ssu.fasta

Fasta formatted sequences of SSU genes predicted by Barrnap

=item <output_prefix>.barrnap.ssu.usearch.out

Usearch results of predicted SSU sequences from assembly vs. Silva database

=back

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
use File::Basename;
use Pod::Usage;

my $path_to_ssu_db = "/data/db/phyloFlash_dev/old/phyloFlash_1.5/SSURef_NR99_119_for_phyloFlash.udb"; # Deprecated, keep for testing
my $path_to_ssu_db_vsearch = "/data/db/phyloFlash_dev/phyloFlash/119/SILVA_SSU.noLSU.masked.trimmed.fasta";
    # Path to Vsearch database, identical to that used by phyloFlash, change to your own local settings
my $cpus = 1;       # How many CPUs to use for parallelization-capable tools
my $assem_file;     # Assembly (Fasta formatted) to scan for SSU reads
my $output_prefix;  # Prefix for output files and intermediate files
#my $taxon_level=2; # Pick a taxonomic level to parse from the taxon string (domain=0, phylum=1, class=2...)
my %uniq_pred_hash; # Hash to make sure that each scaffold only contains at most one SSU prediction
my %SSU_assembly;   # Hash to store output from parsing Usearch results

if (! @ARGV) {
    pod2usage(-message => "Insufficient options were supplied", -existatus => 2);
}

GetOptions (
    'dbpath|d=s' => \$path_to_ssu_db_vsearch,
    'cpus|c=i' => \$cpus,
    'assembly|a=s' => \$assem_file,
    'output|o=s' => \$output_prefix,
    'help|h' => sub { pod2usage( -exitstatus => 2, -verbose => 2); },
    'man|m'=> sub { pod2usage ( -exitstatus => 0, -verbose => 2) }
    #'taxlevel|t=i' => \$taxon_level
) or pod2usage(-message => "Error in input arguments", -existatus => 2);

my ($output_prefix_file, $output_prefix_path) = fileparse ($output_prefix);
    # Parse output prefix to filename and path, in case not in current folder
do_barrnap();
filter_barrnap_results();
#do_usearch(); # Deprecated
do_vsearch();
parse_usearch_output();

## SUBROUTINES #################################################################################################

sub do_barrnap {
    ## for microbial SSU ######################################################
    foreach ('bac','arc') {
        my $barrnap_cmd = "barrnap --kingdom $_ --threads $cpus --evalue 1e-15 --reject 0.8 $assem_file | grep '16S_rRNA' \>\> $output_prefix_path\/tmp.$output_prefix_file.scaffolds.gff";
        system ($barrnap_cmd) == 0
            or warn("Could not run [$barrnap_cmd]: $!\n");
    }
    ## for the euk SSU ########################################################
        my $barrnap_euk_cmd = "barrnap --kingdom euk --threads $cpus --evalue 1e-15 --reject 0.8 $assem_file | grep '18S_rRNA' \>\> $output_prefix_path\/tmp.$output_prefix_file.scaffolds.gff";
        system ($barrnap_euk_cmd) == 0
            or warn("Could not run [$barrnap_euk_cmd]: $!\n");
            # warn instead of die on nonzero exit status because sample may not
            # contain any eukaryotic 18S
    ## for the mitochondrial SSU which is 12S not 16S #########################
        my $barrnap_mito_cmd = "barrnap --kingdom mito --threads $cpus --evalue 1e-15 --reject 0.8 $assem_file | grep '12S_rRNA' \>\> $output_prefix_path\/tmp.$output_prefix_file.scaffolds.gff";
        system ($barrnap_mito_cmd) == 0
            or warn("Could not run [$barrnap_mito_cmd]: $!\n");
            # warn instead of die on nonzero exit status because sample may not
            # contain any mitochondrial 12S
    ##
}

sub filter_barrnap_results {
    open(SSUGFF, "< $output_prefix_path\/tmp.$output_prefix_file.scaffolds.gff")
        or die ("Cannot open GFF file $output_prefix_path\/tmp.$output_prefix_file.scaffolds.gff : $!\n");
    while (<SSUGFF>) {
        chomp;
        my @theline = split /\s+/, $_;
        ## If the scaffold does not already have an SSU prediction, add it to #
        ## our final list #####################################################
        if (!exists $uniq_pred_hash{$theline[0]}) {
            $uniq_pred_hash{$theline[0]} = $_; 
        } else {next;}
    }
    close(SSUGFF);
    open(SSUGFFWRITE, "> $output_prefix.barrnap.ssu.gff")
        or die ("Cannot open file to write final barrnap SSU predictions: $!\n");
    foreach my $thescaffold (keys %uniq_pred_hash) {
        print SSUGFFWRITE $uniq_pred_hash{$thescaffold}, "\n";
    }
    close(SSUGFFWRITE);
}

sub do_usearch { # Deprecated, still here for testing purposes
    my $get_fasta_cmd = "fastaFromBed -fi $assem_file -bed $output_prefix.barrnap.ssu.gff -fo $output_prefix.barrnap.ssu.fasta -s";
    system($get_fasta_cmd)==0 or die "Could not run [$get_fasta_cmd] : $! \n";
    my $usearch_cmd = "usearch -usearch_global $output_prefix.barrnap.ssu.fasta -db $path_to_ssu_db -id 0.7 -userout $output_prefix.barrnap.ssu.usearch.out -userfields query+target+id+alnlen+evalue -threads $cpus --strand plus --notrunclabels -notmatched $output_prefix.barrnap.ssu.fasta.usearch.notmatched.fasta -dbmatched $output_prefix.barrnap.ssu.fasta.usearch.dbhits.fasta";
    system($usearch_cmd)==0 or die "Could not run [$usearch_cmd] : $!\n";
}

sub do_vsearch {
    my $get_fasta_cmd = "fastaFromBed -fi $assem_file -bed $output_prefix.barrnap.ssu.gff -fo $output_prefix.barrnap.ssu.fasta -s";
    system($get_fasta_cmd)==0 or die "Could not run [$get_fasta_cmd] : $! \n";
    my $vsearch_cmd = "vsearch -usearch_global $output_prefix.barrnap.ssu.fasta -db $path_to_ssu_db_vsearch -id 0.7 -userout $output_prefix.barrnap.ssu.usearch.out -userfields query+target+id+alnlen+evalue -threads $cpus --strand plus --notrunclabels -notmatched $output_prefix.barrnap.ssu.fasta.usearch.notmatched.fasta -dbmatched $output_prefix.barrnap.ssu.fasta.usearch.dbhits.fasta";
    system($vsearch_cmd)==0 or die "Could not run [$vsearch_cmd] : $!\n";
}


sub parse_usearch_output {
    ## Parse the output from Usearch and store in the hash %SSU_assembly ######
    open (SSUASSEM, "< $output_prefix.barrnap.ssu.usearch.out")
        or warn ("Cannot open usearch results: $!\n");
    while (<SSUASSEM>) {
      chomp;
      ## split line to recover the taxon name of the spades hits ##############
      my @ssuassem_line = split("\t",$_);
      my $ssuassem_genbankid;
      my $ssuassem_taxstring;
      ## Split taxon name into Genbank ID and taxon string proper #############
      ($ssuassem_genbankid, $ssuassem_taxstring) = ($ssuassem_line[1]
                                                    =~ /(\w+\.\d+\.\d+)\s([\w\s;\-\.]+)/);
      my $ssuassem_scaffold;
      ($ssuassem_scaffold) = ( $ssuassem_line[0] =~ /^(.*):/ );
      my @ssuassem_taxarray = split ";", $ssuassem_taxstring;
      # kingdom phylum class order family genus species
      ## If the taxonomy string does not extend to species ####################
      if (scalar @ssuassem_taxarray < 7) {  
        my $num_levels_to_add = 7 - (scalar @ssuassem_taxarray);
        my $lowest_taxon = $ssuassem_taxarray[$#ssuassem_taxarray];
        ## Add the lowest taxonomic level assigned, in brackets, to fill in
        ## the missing spaces
        while ($num_levels_to_add > 0) {
            push @ssuassem_taxarray, "\(".$lowest_taxon."\)";
            $num_levels_to_add--;
        }
      }
      #my $current_taxon_level = $taxon_level;
      #if (scalar @ssuassem_taxarray < $taxon_level+1) {
      #  $current_taxon_level = (scalar @ssuassem_taxarray) - 1;
      #}
      ## Recombine the line with info properly split ##########################
      my $ssuassem_newline = join ("\t",
                                   $ssuassem_scaffold,
                                   $ssuassem_line[0],
                                   $ssuassem_genbankid,
                                   join ("\t",
                                         @ssuassem_taxarray[0 .. 6]),
                                   $ssuassem_line[2],
                                   $ssuassem_line[3],
                                   $ssuassem_line[4]);
      $SSU_assembly{$ssuassem_line[0]} = $ssuassem_newline;
    }
    close (SSUASSEM);
    ## Account for reads which had no Usearch hits to database ################
    open(UNASSIGNED, "< $output_prefix.barrnap.ssu.fasta.usearch.notmatched.fasta")
        or die ("Cannot read file $output_prefix.barrnap.ssu.fasta.usearch.notmatched.fasta: $!\n");
    while (<UNASSIGNED>) {
        chomp;
        if ($_ =~ /^>(.*)/) {
            my $ssuassem_id = $1;
            my $ssuassem_scaffold;
            ($ssuassem_scaffold) = ($ssuassem_id =~ /^(.*):/);
            my @filler = ("unassigned") x 8;
            my $fillerstring = join "\t",@filler;
            my $ssuassem_newline = join ("\t",
                                         $ssuassem_scaffold,
                                         $ssuassem_id,
                                         $fillerstring,
                                         "",
                                         "",
                                         "");
            $SSU_assembly{$ssuassem_id} = $ssuassem_newline;
        }
    }
    close(UNASSIGNED);
    ## Print output
    open(FINALOUT, "> $output_prefix.ssu.tab")
        or die ("Cannot write to output file $output_prefix.ssu.tab: $!\n");
    my @finaloutheader = ("scaffold",
                          "SSUid",
                          "SSUgenbankid",
                          "Superkingdom",
                          "Phylum",
                          "Class",
                          "Order",
                          "Family",
                          "Genus",
                          "Species",
                          "pid",
                          "alnlen",
                          "eval");
    print FINALOUT join ("\t", @finaloutheader), "\n";   # header line for output
    foreach my $thekey (keys %SSU_assembly) {
        print FINALOUT $SSU_assembly{$thekey}, "\n";
    }
    close (FINALOUT);
}
