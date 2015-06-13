#!/usr/bin/env perl

## Input format validator for gbtools

## Version 1 - 2015-06-13 - kbseah first version

use strict;
use warnings;
use Getopt::Long;

if (!defined @ARGV) {
    usage();
    exit;
}

my ($covstats_in, # Comma-separated list of covstats files
    $mark_in,     # Comma-separated list of marker taxonomy tables
    $ssu_in,      # SSU annotation table
    $trna_in,     # tRNA annotation table
    $user_in,     # Comma-separated list of user-supplied annotations
    $outdir);     # Output directory for putting cleaned up files

my $do_out = 0;
my @covstats_split;
my @mark_split;
my %scaffolds_hash;
my @mark_fields = ("scaffold", # Necessary fields for marker table files
                   "markerid",
                   "gene",
                   "Superkingdom",
                   "Phylum",
                   "Class",
                   "Order",
                   "Family",
                   "Genus",
                   "Species");
my @ssu_fields = ("scaffold", # Necessary fields for SSU table files
                  "SSUid",
                  "Superkingdom",
                  "Phylum",
                  "Class",
                  "Order",
                  "Family",
                  "Genus",
                  "Species");

GetOptions("covstats=s"=>\$covstats_in,
           "mark=s"=>\$mark_in,
           "ssu=s"=>\$ssu_in,
           "trna=s"=>\$trna_in,
           "user=s"=>\$user_in,
           "outdir|o=s"=>\$outdir) or die ("$!\n");

if (defined $outdir) {
    $do_out = 1;
}

## MAIN #######################################################################

print STDOUT "Validation of input files for gbtools\n\n";
if (!defined $covstats_in) {
    print STDERR "At least one covstats file is necessary to perform the input validation\n";
    print STDERR "Exiting...\n";
    exit;
}
## If only one covstats file supplied 
@covstats_split = split(",",$covstats_in);
if (scalar (@covstats_split) == 1) {
    check_covstats($covstats_in,1);
}
## Else if more than one covstats file supplied, the first one is used to
## hash the scaffold IDs, to check the other covstats files.
elsif (scalar (@covstats_split) > 1 ) {
    check_covstats($covstats_split[0],1);
    shift @covstats_split;
    foreach my $thecovstats (@covstats_split) {
        check_covstats($thecovstats,0);
    }
}
if (defined $mark_in) {
    @mark_split = split(",",$mark_in);
    foreach my $themark (@mark_split) {
        check_mark($themark, "Mark");
    }
}
if (defined $ssu_in) {
    check_mark($ssu_in, "SSU");
}
if (defined $trna_in) {
    check_trna($trna_in)
}
if (defined $user_in) {
    my @user_split = split (",",$user_in);
    foreach my $theuser (@user_split) {
        check_mark($theuser, "User");
    }
}

## SUBROUTINES ################################################################

sub usage {
    print STDERR "\n";
    print STDERR "Input format validator for gbtools package\n";
    print STDERR "\n";
    print STDERR "Usage: perl $0 --covstats FILE [--mark FILE] [--ssu FILE] [--trna FILE]\n\n";
    print STDERR "\t--covstats FILE Coverage statistics, output from BBtools pileup.sh, for example\n";
    print STDERR "\t--mark FILE Taxonomic classification of marker genes\n";
    print STDERR "\t--ssu FILE Taxonomic classification of SSU rRNA marker genes\n";
    print STDERR "\t--trna FILE tRNA markers and their type, output from tRNAscan-SE, for example\n";
    print STDERR "\t--user FILE User-supplied annotation file\n";
    print STDERR "\n";
    print STDERR "To be used with the gbtools package: https://github.com/kbseah/genome-bin-tools/\n\n";
}

sub check_covstats {
    my ($infile, $first) = @_;
    my $fatalerr = 0;
    my $linecount = 0;
    ## Open output file, if option to write given
    my $outfile = $infile."mod"; # Filename for output file - TODO : use rel paths
    if ($do_out==1) {
        open(OUT, ">", $outfile) or die ("$!\n");
    }
    my @necessary_fields = ("ID",
                            "Avg_fold",
                            "Length",
                            "Ref_GC");
    open(IN, "<", $infile) or die ("$1\n");
    ## Check the first line for comment character at beginning of line ########
    my $firstline = <IN>;
    $linecount++;
    chomp $firstline; # Strip EOL character
    if ($firstline =~ m/^#(.*)/) {
        print STDOUT "Covstats file $infile: first line has comment character \#\n";
        print STDOUT " Offending character will be removed in output\n";
        my $firstline = $1;
    }
    ## Check fields in first line (should be header line)
    my @firstline_split = split("\t",$firstline);
    #print join ("\t",@firstline_split),"\n";
    my %firstline_hash;
    for (my $i=0; $i < scalar(@firstline_split); $i++) {
        $firstline_hash{$firstline_split[$i]} = $i;
    }
    foreach my $necessary (@necessary_fields) {
        if (!exists $firstline_hash{$necessary}) {
            print STDOUT "Covstats file $infile: field $necessary missing from header line\n";
            $fatalerr++;
        }
    }
    if ($do_out==1 && $fatalerr ==0) {
        print OUT $firstline."\n";
    }
    ## Check that fields are numeric and appropriate ##########################
    while (<IN>) {
        if ($fatalerr > 0) { # Break and finish if fatal error detected.
            last;
        }
        $linecount++;
        chomp;
        my @linesplit = split ("\t", $_);
        ## Check that Avg_fold field is a number above zero ###################
        if (!is_number($linesplit[$firstline_hash{"Avg_fold"}])
            || $linesplit[$firstline_hash{"Avg_fold"}] < 0 ) {
            $fatalerr++;
            print STDOUT "Covstats file $infile line $linecount: Avg_fold field is not a positive number\n";
        }
        ## Check that Length field is a number above zero #####################
        if (!is_number($linesplit[$firstline_hash{"Length"}])
            || $linesplit[$firstline_hash{"Length"}] < 0) {
            $fatalerr++;
            print STDOUT "Covstats file $infile line $linecount: Length field is not a positive number\n";
        }
        ## Check that Ref_GC field is a number between 0, 1 ###################
        if (!is_number($linesplit[$firstline_hash{"Ref_GC"}])
            || $linesplit[$firstline_hash{"Ref_GC"}] > 1
            || $linesplit[$firstline_hash{"Ref_GC"}] < 0) {
            $fatalerr++;
            print STDOUT "Covstats file $infile line $linecount: Ref_GC field is not a number between 0, 1\n";
        }
        ## If this is the first covstats file, hash the scaffolds #############
        if ($first == 1) {
            $scaffolds_hash{$linesplit[$firstline_hash{"ID"}]} = 1;
        }
        ## Otherwise check that scaffold IDs  match first file ################
        elsif ($first != 1) {
            if (!defined $scaffolds_hash{$linesplit[$firstline_hash{"ID"}]} ) {
                $fatalerr++;
                print STDOUT "Covstats file $infile line $linecount: Scaffold ID doesn't match any IDs found in the first covstats file\n";
            }
        }
        ## Print validated line to output if no fatal errors found ############
        if ($do_out==1 && $fatalerr==0) {
            print OUT join("\t", @linesplit),"\n";
        }
        
    }
    close(IN);
    ## Close output file, if opened
    if ($do_out==1) {
        close(OUT);
    }
    if ($fatalerr == 0 ) {
        print STDOUT "Covstats file $infile: No fatal errors detected\n";
    }
}

sub check_mark { # Also works for SSU tables
    my ($infile, $type) = @_;
    my $fatalerr=0;
    my $linecount=0;
    my @necessary_fields;
    if ($type eq "Mark") {
        @necessary_fields = @mark_fields;
    } elsif ($type eq "SSU") {
        @necessary_fields = @ssu_fields;
    } elsif ($type eq "User") {
        @necessary_fields = ("scaffold");
    }
    my $outfile = $infile."mod"; # Filename for output file - TODO : use rel paths
    if ($do_out==1) {
        open(OUT, ">", $outfile) or die ("$!\n");
    }
    ## Open input file ########################################################
    open(IN, "<", $infile) or die ("$!\n");
    ## Check that firstline is a valid header #################################
    my $firstline = <IN>;
    $linecount++;
    chomp $firstline;
    my @firstline_split = split ("\t",$firstline);
    my %firstline_hash;
    for (my $i=0; $i < scalar(@firstline_split); $i++) {
        $firstline_hash{$firstline_split[$i]} = $i;
    }
    foreach my $necessary (@necessary_fields) {
        if (!exists $firstline_hash{$necessary}) {
            print STDOUT "Mark file $infile: field $necessary missing from header line\n";
            $fatalerr++;
        }
    }
    ## If no fatal errors, print header line to output ########################
    if ($do_out==1 && $fatalerr==0) {
        print OUT $firstline."\n";
    }
    
    ## Check subsequent lines #################################################
    while (<IN>) {
        if ($fatalerr > 0) { # Break loop if fatal error found
            last;
        }
        $linecount++;
        my $line = $_;
        chomp $line;
        ## Check for problematic characters ###################################
        if ($line =~ s/'//g) {
            print STDOUT "Mark file $infile line $linecount: Problematic character \' found. Stripping it from line...\n";
        }
        my @linesplit = split("\t",$line);
        ## Check that scaffold ID matches those in the first covstats file ####
        if (!defined $scaffolds_hash{$linesplit[$firstline_hash{"scaffold"}]}) {
            $fatalerr++;
            print STDOUT "Mark file $infile line $linecount: Scaffold ID doesn't match any ID in covstats file\n";
        }
        ## Print to output if validated #######################################
        if ($do_out == 1 && $fatalerr==0) {
            print OUT $line."\n";
        }
    }
    close(IN);
    ## Close output file, if opened ###########################################
    if ($do_out==1) {
        close(OUT);
    }
    ## Sound the all-clear ####################################################
    if ($fatalerr==0) {
        print STDOUT "$type file $infile: No fatal errors detected\n";
    }
}

sub check_trna {
    my ($infile) = @_;
    my $fatalerr=0;
    my $linecount=0;
    my $header;
    open(IN, "<", $infile) or die ("$!\n");
    ## Skip the first three lines (header)
    for (my $i=1; $i<=3; $i++) {
        $linecount++;
        $header=<IN>;
    }
    ## Check that scaffold IDs match
    while (<IN>) {
        if ($fatalerr>0) {
            last;
        }
        chomp;
        $linecount++;
        my @linesplit = split (/\s+/);
        if (!defined $scaffolds_hash{$linesplit[0]}) {
            $fatalerr++;
            print STDOUT "tRNA file $infile line $linecount: Scaffold ID doesn't match any ID in covstats file\n";
        }
        if (!defined $linesplit[4]) {
            print STDOUT "tRNA file $infile line $linecount: tRNA type appears to be missing \n";
        }
    }
    close(IN);
    if ($fatalerr==0) {
        print STDOUT "tRNA file $infile: No fatal errors detected\n";
    }
    
}

sub is_number {
    shift =~ /^\s*-?\d+\.?\d*\s*$/;
}