#! /usr/bin/perl

## FASTG parsing tool -- by Brandon Seah (kbseah@mpi-bremen.de)

## Version 4 -- 2014-10-15
## -- Added write_cytoscape_simple() subroutine to provide "simple" network files without 5' or 3' ends indicated
## -- Modified read_bait() subroutine to fix cases where contig names include ">" character, or are Fastg headers

## Parses Fastg output from SPAdes assembler and generates the following:
## TSV file for visualizing network, with each contig represented as one node
## TSV file with contig attributes, for import to Cytoscape together with the TSV network file
## SIF file for visualizing network with each contig represented as two nodes connected by an edge of class "intra"; the two nodes represent 3' and 5' ends of the contig
## TSV file for visualizing network with each contig represented as two nodes; contains node attribute information too (use Cytoscape -- import network from table to visualize)
## Fasta and Fastg files containing all contigs connected to a list of specified "bait" contigs
## Log file giving summary of fishing expedition

## PACKAGES USED

use strict;
use warnings;
use Getopt::Long;


## DEFINE GLOBAL VARIABLES 

my $input_fastg;		# input fastg filename
my $input_bait_file;		# list of contigs to use as bait for connectivity fishing; names of contigs should be identical to those in headers, and without > or ; characters
my $output_prefix;		# Prefix for all the output files for each run
my $intra_nodes_counter = 0;	# Counter for number of contigs, for logfile reporting
my $inter_nodes_counter = 0;	# Counter for number of connections between contigs, for logfile reporting
my $baits_counter = 0;		# Counter for initial number of bait contigs
my $fished_counter = 0;		# Counter for number of contigs fished using the bait
my %contigs_hash;		# Hash to hold all contig names
my @headers_array;		# Array to hold all Fastg headers from the original input Fastg file
my %gc_hash;			# Hash to store GC values for each Fastg entry
my %length_hash;		# Hash to store length values for each Fastg entry


## USAGE OPTIONS

if (! @ARGV) { usage(); }
GetOptions (
	"input|i=s" => \$input_fastg,	# input Fastg file
	"bait|b=s" => \$input_bait_file,	# bait file for fishing
	"output|o=s" => \$output_prefix	# prefix for output files
)
or usage(); 


## MAIN

open (LOGFILE, "> $output_prefix\.log") or die ("Cannot open log file for writing:$!\n");	# Start the log file
	print LOGFILE "FASTG parsing tool by Brandon Seah 2014-10-15... \"Don't blame me!\"\n";
	print LOGFILE "Started FASTG parsing tool\n";
my $datetime = localtime();
	print LOGFILE $datetime, "\n";
	print LOGFILE "**********************************************************************\n";
	print LOGFILE "\n";
count_GC();
write_sif_output();
write_cytoscape_simple();
	print LOGFILE "Input FASTG file:       ", $input_fastg, "\n";
	print LOGFILE "  Contained ", $intra_nodes_counter, " contigs and ", $inter_nodes_counter, " connections (including self-connections)\n";
my %fishing_hash = %contigs_hash;	# create fishing hash from contigs hash
read_bait();				# read bait file, which converts values of all keys found in the initial bait file to 0 in the fishing hash
	print LOGFILE "Input bait file:        ", $input_bait_file, "\n";
	print LOGFILE "  Containing ", $baits_counter, " contigs of interest\n";
perform_fishing();
write_fastg_output();
	print LOGFILE "\n";
	print LOGFILE "Found ", $fished_counter, " contigs connected to the bait contigs (including bait contigs in total count)\n";
	print LOGFILE "Output files prefix:    ", $output_prefix, "\n";
$datetime = localtime();
	print LOGFILE "\n";
	print LOGFILE "Job finished\n";
	print LOGFILE $datetime, "\n";
	print LOGFILE "\n";
	print LOGFILE "**********************************************************************\n";
print STDERR "Job finished\n";
close (LOGFILE);


## SUBROUTINES

sub usage {				# Usage message - print and exit when script called without arguments
	print STDERR "\n";
	print STDERR " ******************************************************************************************************************************************\n";
	print STDERR "   FASTG format parsing tool ** Brandon Seah 2014-10-15 ** \"Don't blame me\"        \n";
	print STDERR "   Usage:                                                                           \n";
	print STDERR "   \$ perl fastg_parser_04.pl \\                                                    \n";
	print STDERR "      -i <input Fastg filename> \\                                                  \n";
	print STDERR "      -b <file with list of contigs as bait to retrieve connected contigs> \\       \n";
	print STDERR "      -o <output prefix for output of parser tool>                                  \n";
	print STDERR "   The parser will generate the following output files:                             \n";
	print STDERR "    <output_prefix>.simple.tsv        TSV format graph file of connectivity, each contig represented by one node\n";
        print STDERR "    <output_prefix>.attr.tsv          TSV format table with info for each contig, import as attribute table for simple graph in Cytoscape\n";
	print STDERR "    <output prefix>.sif               SIF format graph file of connectivity, each contig represented by two nodes (5\' and 3\' ends)\n";
	print STDERR "    <output prefix>.tsv               TSV format graph file of connectivity, each contig represented by two nodes (5\' and 3\' ends)\n";
	print STDERR "    <output prefix>.fasta             Fasta file with fished contigs                \n";
	print STDERR "    <output prefix>.fastg             Fastg file with fished contigs                \n";
	print STDERR "    <output prefix>.log               Log file with summary of the parser run       \n";
	print STDERR " *******************************************************************************************************************************************\n";
	print STDERR "\n";
	exit;
}

sub count_GC {
open(READFASTG, "< $input_fastg") or die ("Error! Cannot open Fastg input file: $input_fastg \: $!\n");
my $current_node;		# Current contig being counted
while(<READFASTG>) {
	chomp;
	if ($_ =~ /^\>/) {				# If this is a header line...
		my $stripped_line = $_;			# Parse the header
		$stripped_line =~ s/[\>;]//g;
		my @nodes_array = split (":", $stripped_line);
		$current_node=$nodes_array[0];		# Update the current contig
		$length_hash{$current_node}=0;		# Reset the counters
		$gc_hash{$current_node}=0;		# Reset the counters
	}
	else {						# Otherwise, keep counting
		$length_hash{$current_node} += length($_);	# Add base count to total length for this contig
		while ($_ =~ /[GCgc]/g) { $gc_hash{$current_node}++; }	# Count numbers of GC bases in this line
	}
}
close(READFASTG);
}

sub write_sif_output {						# Subroutine for writing network files in SIF and TSV format showing contig connectivity
open (SIFOUTPUT, "> $output_prefix\.sif") or die ("Error! Cannot open SIF file for output: $!\n");
open (TSVOUTPUT, "> $output_prefix\.tsv") or die ("Error! Cannot open TSV file for output: $!\n");
print TSVOUTPUT "node1\t", "interaction\t", "node2\t", "length\t", "coverage\t", "GC", "\n";		# header line for the TSV output file
open(FASTGINPUT, "< $input_fastg") or die("Error! Cannot open Fastg input file: $input_fastg \: $!\n");
while (<FASTGINPUT>) {
	chomp;
	if ($_ =~ /^\>/) { 
		my $stripped_line = $_;
		push @headers_array, $stripped_line;			# push header line verbatim to array of headers
		$stripped_line =~ s/[\>;]//g;				# remove extraneous characters from headers
		my @nodes_array = split(":", $stripped_line);		# split header into the component connected nodes
		$contigs_hash{$nodes_array[0]} = 1;			# save the node as a key in the contigs hash
		my @attributes_array = split("_", $nodes_array[0]);	# split the node name into its components to extract the length and coverage information
		my $gc_frac = $gc_hash{$nodes_array[0]}/$length_hash{$nodes_array[0]};	# calculate the GC%
		print SIFOUTPUT $nodes_array[0],"_5 ", "intra ", "$nodes_array[0]", "_3\n";		# print to graph SIF file, with contig represented as 5' and 3' nodes connected by an internal edge
		print TSVOUTPUT $nodes_array[0],"_5\t", "intra\t", $nodes_array[0], "_3\t", $length_hash{$nodes_array[0]}, "\t", $attributes_array[5],"\t", $gc_frac, "\t\n";	# print to graph TSV file
		$intra_nodes_counter++;
		if (scalar @nodes_array > 1) {
			for (my $i=1; $i < scalar @nodes_array; $i++) {
				if ($nodes_array[$i] =~ /'$/) {		# if the connected node is annotated with a ', meaning it is connected as reverse-complement to the first contig
					my $thisnode = $nodes_array[$i];
					$thisnode =~ s/\'//g;	# strip the ' character from node name
					print SIFOUTPUT $nodes_array[0],"_3 ", "inter ", $thisnode, "_3\n";	# print to graph SIF file, connection between 3' end of first contig and 3' end of this contig
					print TSVOUTPUT $nodes_array[0], "_3\t", "inter\t", $thisnode, "_3\t", "\t\n";
					$inter_nodes_counter++;
				}
				else {
					print SIFOUTPUT $nodes_array[0],"_3 ", "inter ", $nodes_array[$i], "_5\n";	# print to graph SIF file, connection between 3' end of first contig and 5' start of this contig
					print TSVOUTPUT $nodes_array[0], "_3\t", "inter\t", $nodes_array[$i], "_5\t", "\t\n";
					$inter_nodes_counter++;
				}
			}
		}
	}
	else { next; }
}
close(FASTGINPUT);
close (SIFOUTPUT);
close (TSVOUTPUT);
}

sub write_cytoscape_simple {						# Subroutine for writing "simple" network files in TSV format for Cytoscape import 
open (SIMOUTPUT, "> $output_prefix\.simple.tsv") or die ("Error! Cannot open TSV file for output: $!\n");
open (ATRTABLE, ">$output_prefix\.attr.tsv") or die ("Error! Cannot open TSV file for output: $!\n");	# TSV file containing node attributes for Cytoscape import
print ATRTABLE "node\t", "length\t", "coverage\t", "GC\n";						# header line for attributes table
print SIMOUTPUT "node1\t", "interaction\t", "node2\n";							# header line for the TSV output file
open(FASTGINPUT, "< $input_fastg") or die("Error! Cannot open Fastg input file: $input_fastg \: $!\n");
while (<FASTGINPUT>) {
	chomp;
	if ($_ =~ /^\>/) { 
		my $stripped_line = $_;
		push @headers_array, $stripped_line;			# push header line verbatim to array of headers
		$stripped_line =~ s/[\>;]//g;				# remove extraneous characters from headers
		my @nodes_array = split(":", $stripped_line);		# split header into the component connected nodes
		$contigs_hash{$nodes_array[0]} = 1;			# save the node as a key in the contigs hash
		my @attributes_array = split("_", $nodes_array[0]);	# split the node name into its components to extract the length and coverage information
		my $gc_frac = $gc_hash{$nodes_array[0]}/$length_hash{$nodes_array[0]};	# calculate the GC%
		print ATRTABLE $nodes_array[0], "\t", $length_hash{$nodes_array[0]}, "\t", $attributes_array[5], "\t", $gc_frac, "\n";	# print contig attributes to attributes table
		if (scalar @nodes_array == 1) {
			print SIMOUTPUT $nodes_array[0], "\t", "none\t", "\t\n";	# Do something when node is not connected to anything else
		}
		elsif (scalar @nodes_array > 1) {
			for (my $i=1; $i < scalar @nodes_array; $i++) {
				if ($nodes_array[$i] =~ /'$/) {				# if the connected node is annotated with a ', meaning it is connected as reverse-complement to the first contig
					my $thisnode = $nodes_array[$i];
					$thisnode =~ s/\'//g;				# strip the ' character from node name
					print SIMOUTPUT $nodes_array[0], "\t", "inter\t", $thisnode, "\n";
				}
				else {
					print SIMOUTPUT $nodes_array[0], "\t", "inter\t", $nodes_array[$i],"\n";	# print to graph SIF file, connection between 3' end of first contig and 5' start of this contig
				}
			}
		}
	}
	else { next; }
}
close(FASTGINPUT);
close (SIMOUTPUT);
close (ATRTABLE);
}


sub read_bait {				# Read the list of contig header names to use as bait for fishing
open (BAITINPUT, "< $input_bait_file") or die ("Error! Cannot open file containing list of headers for fishing $input_bait_file \: $!\n");
while(<BAITINPUT>) {
	chomp;
	my $stripped_line = $_;
	if ( $_ =~ /[\>;:]/ ) {
		$stripped_line =~ s/[\>;]//g;				# remove extraneous characters from headers
		if ( $stripped_line =~ /:/ ) {
			my @nodes_array = split (":", $stripped_line);
			$stripped_line = $nodes_array[0];
		}
	}
	$fishing_hash{$stripped_line} = 0;
	$baits_counter++;
}
close (BAITINPUT);
}

sub perform_fishing {						# extract contigs connected to the bait iteratively
	my $initial_count=0;					# define counter for the number of fished contigs, with dummy initial value
	my $current_count=1;					# define second counter for number of fished contigs, with dummy initial value
	while ($initial_count != $current_count) {				# While the iterative fishing is still yielding new returns...
		$initial_count=0;
		$current_count=0;
		foreach my $thiskey (keys %fishing_hash) {			# Count how many contigs are marked as fishing bait
			if ($fishing_hash{$thiskey} == 0) { $initial_count++; }
		}
		foreach my $elem (@headers_array) {				# For each header line in the Fastg file...
			$elem =~ s/[\>';]//g;					# strip extraneous characters
			my @split_elem = split(":", $elem);			# split the header into the component list of connected contigs
			my $flag = 0;						# define the flag for whether or not this contig is in the bait list
			foreach my $elem2 (@split_elem) {
				if ($fishing_hash{$elem2} == 0) { $flag = 1; }	# if this contig is in the bait list, set the flag to 1
			}
			if ($flag == 1) {					# if the flag is set to one
				foreach my $elem2 (@split_elem) {
					$fishing_hash{$elem2} = 0;		# add the connected contigs to the bait list by setting value in fishing_hash to 0
				}
			}
		}
		foreach my $thiskey2 (keys %fishing_hash) {			# Count how many contigs are marked as fishing bait after the fishing procedure
			if ($fishing_hash{$thiskey2} == 0) { $current_count++; }
		}
	}
}

sub write_fastg_output {							# Write fastg output containing all fished contigs
open (FASTGOUTPUT, "> $output_prefix\.fastg") or die ("Error! Cannot open Fastg file for output of fished contigs:$!\n");
open (FASTAOUTPUT, "> $output_prefix\.fasta") or die ("Error! Cannot open Fasta file for output of fished contigs:$!\n");
open (FASTGINPUT2, "< $input_fastg") or die ("Error! Cannot open Fastg input file: $input_fastg \:$!\n");
my $writeflag = 0;
while (<FASTGINPUT2>) {
	chomp;
	my $theheader;
	if ($_ =~ /^\>/) {
		my $dummy = $_;
		$dummy =~ s/[\>';]//g;					# remove special characters
		my @dummyarray = split (":",$dummy);			# split to an array
		$theheader = $dummyarray[0];				# get the name of the contig proper
		if ($fishing_hash{$theheader} == 0) { 			# if the contig is in the fished list, turn on write flag
			$writeflag = 1;
			$fished_counter++;			
		 }
		elsif ($fishing_hash{$theheader} == 1) {$writeflag = 0; }	# otherwise turn it off
	}
	if ($writeflag == 1) { 						# if writeflag is on, print the line to the output FASTG file
		print FASTGOUTPUT $_, "\n"; 
		if ($_ =~ /^\>/) { print FASTAOUTPUT ">",$theheader,"\n"; }
		else { print FASTAOUTPUT $_, "\n"; }
	}
}
close (FASTGINPUT2);
close (FASTGOUTPUT);
close (FASTAOUTPUT);
}
