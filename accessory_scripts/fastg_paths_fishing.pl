#!/usr/bin/env perl

# SPAdes-style Fastg and Scaffold.paths parsing tool
# Brandon Seah (kbseah@mpi-bremen.de)
# 2016-02-25

# Given output from SPAdes:
#   Fastg assembly graph,
#   scaffolds.paths (or contigs.paths) file (from SPAdes 3.6.2 onwards)
# And shortlist of scaffold IDs
# Return list of scaffolds connected to the initial list

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

## Global variables ##############################

my $version="2016-02-23";
my $fastg_file;
my $paths_file;
my $iter_depth = 0; # How many fishing iterations (0 = no limit)
my $iter_count = 0; # Counter for iterations of bait fishing
my $out="fastg_fishing"; # Output file prefix
my $bait_file;
my $rflag=0; # Flag if script called from within R
my %scaffolds_fullnames_hash;
my @bait_nodes_array;
my @bait_edges_array;
my %node_edge_hash; # Hash of nodes keyed by edges
my %edge_node_hash; # Hash of edges (arrayed) keyed by nodes
my %fastg_hash; # Hash of edges in Fastg file
my %edge_fishing_hash; # Hash of edges used for fishing
my %fished_nodes_hash; # Hash of nodes corresponding to fished edges

## Usage options #################################

if (! @ARGV) { usage(); } # Print usage statement if no arguments
GetOptions (
    "fastg|g=s" =>\$fastg_file, # Fastg file of assembly graph
    "paths|p=s" =>\$paths_file, # File scaffold.paths or contigs.paths
    "iter|i=i" =>\$iter_depth, # Number of iterations
    "output|o=s" =>\$out, # Output prefix
    "bait|b=s" =>\$bait_file, # List of scaffold names to fish
    "rflag|r" =>\$rflag # Report to stdout if called from R
)
or usage();

## MAIN ###########################################

my ($out_file, $out_path) = fileparse($out); # Parse output prefix supplied

open(my $outlog_fh, ">", $out_path.$out_file.".log") or die ("$!\n"); # Start output log
print $outlog_fh "Fastg fishing log\n";
print $outlog_fh "Script called: $0 \n";
print $outlog_fh "Version: $version\n";
print $outlog_fh scalar localtime() ."\n\n";

hash_nodes_edges();
read_bait_nodes();
#read_bait_edges(); # For checking

print $outlog_fh "Number of bait scaffolds: ". scalar @bait_nodes_array. "\n";
print $outlog_fh "Number of corresponding bait edges: ". scalar @bait_edges_array. "\n";

read_fastg();
perform_fishing_edges($iter_depth);
translate_fished_edges_to_nodes();

if ($rflag == 0) { # Save list of fished contigs to file by default
    open(my $outlist_fh, ">", $out_path.$out_file.".scafflist") or die ("$!\n");
    foreach my $thenode (sort {$a cmp $b} keys %fished_nodes_hash) {
        #print STDOUT $scaffolds_fullnames_hash{$thenode}."\t".$fished_nodes_hash{$thenode}."\n";
        print $outlist_fh $scaffolds_fullnames_hash{$thenode}."\n";
    }
    close ($outlist_fh);
} elsif ($rflag == 1) { # List of fished contigs to STDOUT if called from R function
    foreach my $thenode (sort {$a cmp $b} keys %fished_nodes_hash) {
        print STDOUT $scaffolds_fullnames_hash{$thenode}."\n";
    }
}

print $outlog_fh "Number of fishing iterations: $iter_count\n";
print $outlog_fh "Number of fished scaffolds: ". scalar (keys %fished_nodes_hash) ."\n";
print $outlog_fh "Total number of scaffolds: ". scalar (keys %scaffolds_fullnames_hash). "\n\n";
print $outlog_fh "Output files: \n";
print $outlog_fh "Log file: ".abs_path($out_path.$out_file).".log\n";
print $outlog_fh "List of fished scaffolds: ".abs_path($out_path.$out_file).".scafflist\n";

close($outlog_fh); # Close log file

## SUBROUTINES #####################################

sub usage {
    print STDERR "Fastg connectivity-based fishing\n\n";
    print STDERR "Version: $version\n\n";
    print STDERR "Finds contigs connected to a set of 'bait' contigs,\n";
    print STDERR "using Fastg-formatted assembly graph, and a corresponding\n";
    print STDERR "'paths' file linking scaffolds to graph edges.\n";
    print STDERR "These files are produced by SPAdes 3.6.2 onwards \n\n";
    print STDERR "Usage: perl $0 \\ \n";
    print STDERR "\t -g assembly_graph.fastg \\ \n";
    print STDERR "\t -p scaffolds.paths \\ \n";
    print STDERR "\t -o output_prefix \\ \n";
    print STDERR "\t -b list_of_bait_scaffolds \\ \n";
    print STDERR "\t -i [number of fishing iterations] \n";
    print STDERR "\t -r [Flag if called from within R] \n";
    print STDERR "\n";
    print STDERR "-i 0 means fishing iterates until no more contigs retrieved\n\n";
    print STDERR "Output: Logfile (prefix.log) and list of fished scaffolds (prefix.scafflist)\n\n";
    exit;
}

sub translate_fished_edges_to_nodes {
    # Convert list of fished edges to list of nodes
    foreach my $theedge (keys %edge_fishing_hash) {
        foreach my $thenode (@{$node_edge_hash{$theedge}}) {
            $fished_nodes_hash{$thenode}++;
        }
    }
}

sub perform_fishing_edges {
    my $depth = shift @_; # Iteration depth
    if ($depth == 0) { # I know this is inelegant
        $iter_depth = 2;
    }
    my $init_count=0; # Counter for no. of fished contigs
    my $curr_count=1; # Second counter
    while ($iter_count <= $iter_depth-1 # Loop until iteratively complete
           && $init_count != $curr_count) { # or when fishing is exhausted
        $iter_count++; # Iteration counter
        $init_count=0; # Reset counters for each iterative step
        $curr_count=0;
        if ($depth == 0) { # If exhaustive looping,
            $iter_depth++;  # keep shifting goalposts
        }
        # Count edges marked as bait
        foreach my $fishedge (keys %edge_fishing_hash) { 
            if (defined $edge_fishing_hash{$fishedge}
                && $edge_fishing_hash{$fishedge}  == 1) {
                $init_count++;
            }
        }
        # For all Fastg entries
        foreach my $fastg_entry (keys %fastg_hash) {
            # If edge is listed as a bait
            if (defined $edge_fishing_hash{$fastg_entry}
                && $edge_fishing_hash{$fastg_entry} == 1) { 
                foreach my $connected_edges (@{$fastg_hash{$fastg_entry}}) {
                    # Mark those connected edges as bait
                    $edge_fishing_hash{$connected_edges} = 1; 
                }
            }
        }
        # Count edges marked as bait after fishing
        foreach my $fishedge (keys %edge_fishing_hash) { 
            if (defined $edge_fishing_hash {$fishedge}
                && $edge_fishing_hash {$fishedge} == 1) {
                $curr_count++;
            }
        }
    }
    print $outlog_fh "Number of fished edges: ". $curr_count."\n";
}

sub read_fastg {
    open(FASTGIN, "<", $fastg_file) or die ("$!\n"); # Read in Fastg file
    while (<FASTGIN>) {
        chomp;
        # If a Fastg header line
        if ($_ =~ m/^>EDGE_(\d+)_.*:(.+);$/) { 
            my $current_node = $1;
            my $conn_line = $2;
            # Remove inverted commas, which mark revcomp
            $conn_line =~ s/'//g; 
            my @conn_line_array = split ",", $conn_line;
            foreach my $the_header (@conn_line_array) {
                if ($the_header =~ m/EDGE_(\d+)_/) {
                    push @{$fastg_hash {$current_node}}, $1;
                }
            }
        } else {
            next;
        }
    }
    close(FASTGIN);
}

sub read_bait_nodes {
    open(BAITIN, "< $bait_file") or die ("$!\n");
    # Put $bait_file within quotes otherwise - meaning STDIN not recognized (?!?)
    while (<BAITIN>) {
        chomp;
        if ($_ =~ m/NODE_(\d+)/) {
            push @bait_nodes_array, $1;
        }
        
        
    }
    close(BAITIN);
    foreach my $thenode (@bait_nodes_array) {
        foreach my $thebait (@{$edge_node_hash{$thenode}}) {
            push @bait_edges_array, $thebait;
            $edge_fishing_hash{$thebait} = 1; 
        }
    }
}

sub read_bait_edges { # For checking
    open(BAITIN, "<", $bait_file) or die ("$!\n");
    while (<BAITIN>) {
        chomp;
        push @bait_edges_array, $_;
        $edge_fishing_hash{$_} = 1; 
    }
    close(BAITIN);
}

sub hash_nodes_edges {
    open(PATHSIN, "<", $paths_file) or die ("Cannot open scaffold or contig paths file $!\n");
    my $skip_flag = 0;
    my $current_node;
    while (<PATHSIN>) {
        chomp;
        my $full_header = $_;
        # If header for reversed scaffold, ignore
        if ($full_header =~ m/NODE.*'$/) { 
            $skip_flag = 1;
            next;
            
        }
        # If header for forward scaffold...
        elsif ($full_header =~ m/NODE_(\d+)_.*\d+$/) { 
            $skip_flag = 0; # Turn off skip flag
            $current_node = $1;
            $scaffolds_fullnames_hash{$current_node} = $full_header; # Save full header name
            next;
        }
        if ($skip_flag == 1) {
            next;
        } else {
            my $edge_line = $_;
            # Remove unneeded chars, only interested in edge IDs
            $edge_line =~ s/[-\+;]//g; 
            my @edge_split = split ",", $edge_line;
            foreach my $the_edge (@edge_split) {
                push @{$node_edge_hash{$the_edge}}, $current_node; # Hash nodes with edges as key
                push @{$edge_node_hash{$current_node}}, $the_edge; # Hash edges with nodes as key
            }
        }
    }
    close (PATHSIN);
}

## Diagnostic functions #####################################

sub test_node_edge_hash {
     # Report hash contents for testing
    foreach my $the_edge (sort {$a cmp $b} keys %node_edge_hash) {
    print $the_edge ."\t";
    print join "\t", @{$node_edge_hash{$the_edge}};
    print "\n";
    }
}

sub test_fastg_hash {
    foreach my $theedge (sort {$a cmp $b} keys %fastg_hash) {
    print $theedge."\t";
    print join "\t", @{$fastg_hash{$theedge}};
    print "\n";
    }
}