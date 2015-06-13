#!/usr/bin/env perl

=head1 NAME
blob_annotate_mod.pl
=head1 SYNOPSIS
blob_annotate_mod.pl --blasttaxid CONTIGTAXIDFILE --assembly ASSEMBLYFASTAFILE [--taxdump TAXDUMPDIR] --out OUTPUTFILE
=head1 AUTHORS
Original author: sujai.kumar@zoo.ox.ac.uk 2013.09.15
Modified by kbseah@mpi-bremen.de to produce output for genome-bin-tools 2015-03-17
=head1 LICENSE
The MIT License (MIT)
Copyright (c) 2013 blaxterlab
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($blasttaxid_file, $taxdump_dir, $assembly_file, $output_file, $blastdb_path) = ("",".","","","");
my @tax_list;
my @bam_files;
my @cov_files;
my $num_threads = 8;

GetOptions (
    "blasttaxid=s" => \$blasttaxid_file,
    "assembly=s"   => \$assembly_file,
    "blastdb=s"    => \$blastdb_path,
    "out:s"        => \$output_file,
    "num_threads:i" => \$num_threads,
    "taxdump:s"    => \$taxdump_dir,
    "taxlist:s{,}" => \@tax_list,
);

sub usage {
    print STDERR "\nAnnotate contigs with taxonomic assignment by Blast against database\n";
    print STDERR " Script originally written by Sujai Kumar for Blobology \n";
    print STDERR " Copyright (c) 2013 blaxterlab https://github.com/blaxterlab/blobology\n";
    print STDERR " Adapted by kbseah for genome-bin-tools 2015-03-17\n\n";
    print STDERR "Arguments: \n";
    print STDERR "   --blasttaxid   (Optional) Pre-computed Blast output, tab-delim, two columns: scaffold name and NCBI taxonid\n";
    print STDERR "   --assembly     Assembly file in Fasta format, if Blast output not yet computed\n";
    print STDERR "   --blastdb      Path to Blast database, e.g. NCBI nt, for computing taxonomic assignments\n";
    print STDERR "   --out          Name of output file\n";
    print STDERR "   --num_threads  Number of threads for Blast+ (default: 8)\n";
    print STDERR "   --taxdump      Folder containing NCBI taxdump files nodes.mp and names.dmp (default: current folder)\n\n";
    print STDERR "If blasttaxid is specified, assembly and blastdb will be ignored.\n\n";
    exit;
}

usage() unless -r "$taxdump_dir/nodes.dmp" and -r "$taxdump_dir/names.dmp";

# Check that blasttaxid has been run; if not, and if assembly file has been specified, run the blast job

if ($blasttaxid_file eq "" && $assembly_file ne "" && $blastdb_path ne "") {
    print STDERR scalar localtime() . " - Running blastn of assembly file $assembly_file against Blast database at $blastdb_path. \n";
    system ("blastn -task megablast -query $assembly_file -db $blastdb_path -evalue 1e-5 -num_threads $num_threads -max_target_seqs 1 -outfmt \'6 qseqid staxids\' -out tmp.$assembly_file.blastout");
    $blasttaxid_file = "tmp.$assembly_file.blastout";
}

if (not @tax_list) {@tax_list = ("superkingdom","phylum","class","order","family","genus","species")};
my %tax_levels;
foreach (@tax_list) {$tax_levels{$_}=1}

#-----------------------------------------------------------------
# Get taxon annotation info from contig-taxid file
#-----------------------------------------------------------------

my (%taxid_has_parent, %taxid_has_taxlevel, %taxid_has_name);
print STDERR scalar localtime() . " - Loading $taxdump_dir/nodes.dmp and $taxdump_dir/names.dmp into memory ...\n";
&load_nodes_names ("$taxdump_dir/nodes.dmp","$taxdump_dir/names.dmp");
print STDERR scalar localtime() . " - Loading $taxdump_dir/nodes.dmp and $taxdump_dir/names.dmp into memory ... DONE\n";

my $blasttaxid_fh = &read_fh($blasttaxid_file);
my %contig_taxinfo;
while (<$blasttaxid_fh>) {
    die "Contig-taxid file $blasttaxid_file does not seem to have two cols with the seqid in the first col and taxid in the second col" unless 
        /^(\S+)\t(\d+)/;
    $contig_taxinfo{$1} = &taxonomy_report($2);
}
close $blasttaxid_fh;

open(OUTPUTFILE, ">", $output_file) or die ("Cannot open output file $output_file for writing: $!\n");  
# print header line
print OUTPUTFILE join( "\t", ("scaffold", "markerid", "gene", "Superkingdom","Phylum", "Class", "Order","Family","Genus","Species")). "\n";
# print marker file
foreach my $thekey (sort {$a cmp $b} (keys %contig_taxinfo)) {
    print OUTPUTFILE $thekey."\t"."NA"."\t"."NA";
    foreach my $thetax (@tax_list) {
        print OUTPUTFILE "\t";
        print OUTPUTFILE (exists ($contig_taxinfo{$thekey}{$thetax}) ? $contig_taxinfo{$thekey}{$thetax} : "(unassigned)");
    }
    print OUTPUTFILE "\n";
}
close(OUTPUTFILE);

############################################################

sub taxonomy_report {
    my $hit_taxid = shift @_;
    my @parents = &get_parents($hit_taxid);
    # convert @parents to tax names:
    my %taxinfo;
    # my $taxonomy_report_string = "";
    for my $parent (@parents) {
        if (exists $taxid_has_taxlevel{$parent} and exists $tax_levels{$taxid_has_taxlevel{$parent}}) {
            $taxinfo{$taxid_has_taxlevel{$parent}} = $taxid_has_name{$parent};
        }
    }
    return \%taxinfo;
    # for my $tax_level (keys %hit_counts) {
        # for my $tax_name (keys %{$hit_counts{$tax_level}}) {
            # $taxonomy_report_string .= "$tax_level\t$tax_name\t";
        # }
    # }
    # return $taxonomy_report_string . "\n";
}

############################################################

sub get_parents {
    my @all = @_;
    my $current_id = $all[0];
    if (exists $taxid_has_parent{$current_id} and $current_id ne $taxid_has_parent{$current_id}) {
        unshift @all, $taxid_has_parent{$current_id};
        @all = &get_parents(@all);
    }
    return @all;
}

############################################################

sub load_nodes_names {
    my $fh;
    my $nodesfile = shift @_;
    my $namesfile = shift @_;
    $fh = &read_fh($nodesfile);
    while (my $line = <$fh>) {
        # line in nodes.dmp should match the regexp below.
        # Change the regexp if NCBI changes their file format
        next if $line !~ /^(\d+)\s*\|\s*(\d+)\s*\|\s*(.+?)\s*\|/;
        $taxid_has_parent{$1} = $2;
        $taxid_has_taxlevel{$1} = $3;
    }
    close $fh;
    
    $fh = &read_fh($namesfile);
    while (my $line = <$fh>) {
        next unless $line =~ /^(\d+)\s*\|\s*(.+?)\s*\|.+scientific name/;
        $taxid_has_name{$1} = $2;
    }
}

############################################################

sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}

############################################################
