#!/usr/bin/env perl

## Tool to extract SSU data and phylogenetic affiliation from new assemblies and parse output to a format to import into R
## Plays with genome_bin_tools.r functions
## Incorporates code from phyloFlash version 1.5d by Harald Gruber-Vodicka

## Prerequisites:
##      Usearch
##      Barrnap 0.5
##      Bedtools 2.21.0

## Contact: kbseah@mpi-bremen.de
## Version 1 - 2014-10-27
## Version 0 - 2014-10-23

use strict;
use warnings;
use Getopt::Long;

my $path_to_ssu_db = "/data/db/phyloFlash/SSURef_NR99_119_for_phyloFlash.udb";      # Path to Usearch database, identical to that used by phyloFlash
my $cpus = 1;       # How many CPUs to use for parallelization-capable tools
my $assem_file;     # Assembly (Fasta formatted) to scan for SSU reads
my $output_prefix;  # Prefix for output files and intermediate files
#my $taxon_level=2;  # Pick a taxonomic level to parse from the taxon string (domain=0, phylum=1, class=2...)
my %uniq_pred_hash;         # Hash to make sure that each scaffold only contains at most one SSU prediction
my %SSU_assembly;           # Hash to store output from parsing Usearch results

if (@ARGV == 0 ) { usage(); }

GetOptions (
    'dbpath|d=s' => \$path_to_ssu_db,
    'cpus|c=i' => \$cpus,
    'assembly|a=s' => \$assem_file,
    'output|o=s' => \$output_prefix,
    #'taxlevel|t=i' => \$taxon_level
);

do_barrnap();
filter_barrnap_results();
do_usearch();
parse_usearch_output();


## SUBROUTINES #################################################################################################

sub usage {
    print STDERR "\n";
    print STDERR "Extract SSU from assembly, search against Silva SSU database, and parse taxonomy output\n";
    print STDERR "\n";
    print STDERR "Usage: \n";
    print STDERR " \$ perl get_ssu_for_genome_bin_tools.pl -d <path_to_ssu_db> -c <CPUs> -a <assembly.fasta> -o <output_prefix> \n";
    print STDERR "\n";
    print STDERR "Options: \n";
    print STDERR " \t -d FILE     Location of SILVA SSU database with taxonomy string, indexed by Usearch\n";
    print STDERR " \t -c INT      Number of processors for barrnap and Usearch \(default: 1\)\n";
    print STDERR " \t -a FILE     Genome assembly in Fasta format\n";
    print STDERR " \t -o STRING   Prefix for output files\n";
    print STDERR "\n";
    print STDERR "Output: \n";
    print STDERR " \t <output_prefix>.ssu.tab                 Table of predicted SSU genes and their taxonomic affiliations\n";
    print STDERR " \t tmp.<output_prefix>.scaffolds.gff       Concatenated output from Barrnap searches\n";
    print STDERR " \t <output_prefix>.barrnap.ssu.gff         Predicted SSU genes, with duplicates removed, GFF feature table\n";
    print STDERR " \t <output_prefix>.barrnap.ssu.fasta       Fasta formatted sequences of SSU genes predicted by Barrnap\n";
    print STDERR " \t <output_prefix>.barrnap.ssu.usearch.out Usearch results of predicted SSU sequences from assembly vs. Silva database\n";
    print STDERR "\n";
    exit;
}

sub do_barrnap {
    foreach ('bac','arc') {
        my $barrnap_cmd = "barrnap --kingdom $_ --threads $cpus --evalue 1e-15 --reject 0.8 $assem_file | grep '16S_rRNA' \>\> tmp.$output_prefix.scaffolds.gff";
        system ($barrnap_cmd) == 0 or die("Could not run [$barrnap_cmd]: $!\n");
    }
    # for the euk SSU
        my $barrnap_euk_cmd = "barrnap --kingdom euk --threads $cpus --evalue 1e-15 --reject 0.8 $assem_file | grep '18S_rRNA' \>\> tmp.$output_prefix.scaffolds.gff";
        system ($barrnap_euk_cmd) == 0 or warn("Could not run [$barrnap_euk_cmd]: $!\n");       # warn instead of die on nonzero exit status because sample may not contain any eukaryotic 18S
    ##
    # for the mitochondrial SSU which is 12S not 16S
        my $barrnap_mito_cmd = "barrnap --kingdom mito --threads $cpus --evalue 1e-15 --reject 0.8 $assem_file | grep '12S_rRNA' \>\> tmp.$output_prefix.scaffolds.gff";
        system ($barrnap_mito_cmd) == 0 or warn("Could not run [$barrnap_mito_cmd]: $!\n");     # warn instead of die on nonzero exit status because sample may not contain any mitochondrial 12S
    ##
}

sub filter_barrnap_results {
    open(SSUGFF, "< tmp.$output_prefix.scaffolds.gff") or die ("Cannot open GFF file tmp.$output_prefix.scaffolds.gff : $! \n");
    while (<SSUGFF>) {
        chomp;
        my @theline = split /\s+/, $_;
        if (!exists $uniq_pred_hash{$theline[0]}) {
            $uniq_pred_hash{$theline[0]} = $_;         # If the scaffold does not already have an SSU prediction, add it to our final list
        } else {next;}
    }
    close(SSUGFF);
    open(SSUGFFWRITE, "> $output_prefix.barrnap.ssu.gff") or die ("Cannot open file to write final barrnap SSU predictions: $! \n");
    foreach my $thescaffold (keys %uniq_pred_hash) {
        print SSUGFFWRITE $uniq_pred_hash{$thescaffold}, "\n";
    }
    close(SSUGFFWRITE);
}

sub do_usearch {
    my $get_fasta_cmd = "fastaFromBed -fi $assem_file -bed $output_prefix.barrnap.ssu.gff -fo $output_prefix.barrnap.ssu.fasta -s";
    system($get_fasta_cmd)==0 or die "Could not run [$get_fasta_cmd] : $! \n";
    my $usearch_cmd = "usearch -usearch_global $output_prefix.barrnap.ssu.fasta -db $path_to_ssu_db -id 0.7 -userout $output_prefix.barrnap.ssu.usearch.out -userfields query+target+id+alnlen+evalue -threads $cpus --strand plus --notrunclabels -notmatched $output_prefix.barrnap.ssu.fasta.usearch.notmatched.fasta -dbmatched $output_prefix.barrnap.ssu.fasta.usearch.dbhits.fasta";
    system($usearch_cmd)==0 or die "Could not run [$usearch_cmd] : $!\n";
}

sub parse_usearch_output {
    open (SSUASSEM, "< $output_prefix.barrnap.ssu.usearch.out") or warn ("Cannot open usearch results: $!\n");		# Parse the output from Usearch and store in the hash %SSU_assembly
    while (<SSUASSEM>) {
      chomp;
      my @ssuassem_line = split("\t",$_);	# split line to recover the taxon name of the spades hits;
      my $ssuassem_genbankid;
      my $ssuassem_taxstring;
      ($ssuassem_genbankid, $ssuassem_taxstring) = ($ssuassem_line[1] =~ /(\w+\.\d+\.\d+)\s([\w\s;\-\.]+)/);	# Split taxon name into Genbank ID and taxon string proper
      my $ssuassem_scaffold;
      ($ssuassem_scaffold) = ( $ssuassem_line[0] =~ /^(.*):/ );
      my @ssuassem_taxarray = split ";", $ssuassem_taxstring;
      # kingdom phylum class order family genus species
      if (scalar @ssuassem_taxarray < 7) {  # If the taxonomy string does not extend to species
        my $num_levels_to_add = 7 - (scalar @ssuassem_taxarray);
        my $lowest_taxon = $ssuassem_taxarray[$#ssuassem_taxarray];
        while ($num_levels_to_add > 0) {
            push @ssuassem_taxarray, "\(".$lowest_taxon."\)";     # Add the lowest taxonomic level assigned, in brackets, to fill in the missing spaces
            $num_levels_to_add--;
        }
      }
      #my $current_taxon_level = $taxon_level;
      #if (scalar @ssuassem_taxarray < $taxon_level+1) {
      #  $current_taxon_level = (scalar @ssuassem_taxarray) - 1;
      #}
      my $ssuassem_newline = join ("\t",$ssuassem_scaffold, $ssuassem_line[0], $ssuassem_genbankid,
                                   join ("\t", @ssuassem_taxarray[0 .. 6]), $ssuassem_line[2], $ssuassem_line[3], $ssuassem_line[4]);	# Recombine the line with info properly split
      $SSU_assembly{$ssuassem_line[0]} = $ssuassem_newline;
    }
    close (SSUASSEM);
    #Account for reads which had no Usearch hits to database
    open(UNASSIGNED, "< $output_prefix.barrnap.ssu.fasta.usearch.notmatched.fasta") or die ("Cannot read file $output_prefix.barrnap.ssu.fasta.usearch.notmatched.fasta: $!\n");
    while (<UNASSIGNED>) {
        chomp;
        if ($_ =~ /^>(.*)/) {
            my $ssuassem_id = $1;
            my $ssuassem_scaffold;
            ($ssuassem_scaffold) = ($ssuassem_id =~ /^(.*):/);
            my @filler = ("unassigned") x 8;
            my $fillerstring = join "\t",@filler;
            my $ssuassem_newline = join ("\t",$ssuassem_scaffold, $ssuassem_id, $fillerstring, "","","");
            $SSU_assembly{$ssuassem_id} = $ssuassem_newline;
        }
    }
    close(UNASSIGNED);
    # Print output
    open(FINALOUT, "> $output_prefix.ssu.tab") or die ("Cannot write to output file $output_prefix.ssu.tab: $!\n");
    print FINALOUT "scaffold\t","SSUid\t","SSUgenbankid\t","Superkingdom\t","Phylum\t","Class\t","Order\t","Family\t","Genus\t","Species\t","pid\t","alnlen\t","eval\n";   # header line for output
    foreach my $thekey (keys %SSU_assembly) {
        print FINALOUT $SSU_assembly{$thekey}, "\n";
    }
    close (FINALOUT);
}
