#!/usr/bin/env perl

=head1 NAME

input_validator.pl - Input format validator for gbtools package

=head1 SYNOPSIS

    perl input_validator.pl \
        --covstats I<file1>,I<file2> \
        --mark I<file1>,I<file2> \
        --ssu I<file> \
        --trna I<file> \
        --user I<file1>,I<file2> \
        --outdir I<string>

    perl input_validator.pl --help

=head1 DESCRIPTION

Validate input file formats for gbtools package in R. See Appendix to online
documentation at https://github.com/kbseah/genome-bin-tools for detailed
description of the file formats.

=head1 ARGUMENTS

=over 8

=item --covstats I<FILE1>,I<FILE2>,...

Comma-separated list of coverage statistics tables (no spaces between filenames).

=item --mark I<FILE1>,I<FILE2>,...

Comma-separated list of marker taxonomy tables (no spaces between filenames).

=item --ssu I<FILE>

SSU rRNA annotation table (output from I<get_ssu_for_genome_bin_tools.pl>)

=item --trna I<FILE>

tRNA annotation table (can directly use output of tRNAscan-SE)

=item --user I<FILE1>,I<FILE2>,...

Comma-separated list of user-defined annotation tables (no spaces between filenames).

=item --outdir|-o I<PATH>

Name for output folder. If supplied, corrected files for commonly-encountered
errors will be written to this folder. (NB: Will overwrite if folder and files
already exist!)

=item --log I<FILE>

Name for message log file.

Default: input_validator.log

=item --help|-h

Help message

=item --man

Full manual page

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
use Pod::Usage;
use File::Spec;
use File::Path qw(make_path);
use Text::Wrap;
use Time::Piece;

my @msg_log; # Log for messages

if (! @ARGV) {
    pod2usage(-message => "Insufficient options were supplied", -existatus => 2);
}

my ($covstats_in, # Comma-separated list of covstats files
    $mark_in,     # Comma-separated list of marker taxonomy tables
    $ssu_in,      # SSU annotation table
    $trna_in,     # tRNA annotation table
    $user_in,     # Comma-separated list of user-supplied annotations
    $outdir,      # Output directory for putting cleaned up files
    );

my $logfile = "input_validator.log"; # Default name for log file
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

GetOptions("covstats=s" => \$covstats_in,
           "mark=s" => \$mark_in,
           "ssu=s" => \$ssu_in,
           "trna=s" => \$trna_in,
           "user=s" => \$user_in,
           "outdir|o=s" => \$outdir,
           "log=s" => \$logfile,
           'help|h' => sub { pod2usage( -exitstatus => 2, -verbose => 2); },
           'man|m'=> sub { pod2usage ( -exitstatus => 0, -verbose => 2) }
           ) or pod2usage(-message => "Error in input arguments", -existatus => 2);

if (defined $outdir) {
    $do_out = 1;
}

## MAIN #######################################################################

msg ("Validation of input files for gbtools", \@msg_log);
if (!defined $covstats_in) {
    error ("At least one covstats file is necessary to perform the input validation", \@msg_log);
    exit;
}

# Make output directory if specified
if (defined $outdir) {
    $outdir = File::Spec->rel2abs ($outdir);
    msg ("Writing output to folder: $outdir", \@msg_log);
    make_path ($outdir);
}

## If only one covstats file supplied 
@covstats_split = split(",",$covstats_in);
if (scalar (@covstats_split) == 1) {
    check_covstats($covstats_in,1);
}
## Else if more than one covstats file supplied, the first one is used to hash
##  the scaffold IDs, to check the other covstats files.
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

msg ("Input validation complete. Thank you for using gbtools", \@msg_log);

# Dump message stack into log file
my $outpath = defined $outdir ? "$outdir/$logfile" : $logfile;
msg ("Writing log file to $outpath", \@msg_log);
open(LOGOUT, ">>", $outpath)
    or error ("Cannot open log file for writing. $!", \@msg_log);
print LOGOUT join "\n", @msg_log;
close(LOGOUT);

## SUBROUTINES ################################################################

sub check_covstats {
    my ($infile, $first) = @_;
    my ($invol,$indir,$infilename) = File::Spec->splitpath($infile);
    my $fatalerr = 0;
    my $linecount = 0;
    ## Open output file, if option to write given
    my $outfile;
    if ($do_out==1) {
        $outfile = $outdir."/".$infilename.".mod"; # Filename for output file
        msg ("Writing corrected covstats file to $outfile",\@msg_log);
        open(OUT, ">", $outfile) or error ("Cannot open $outfile for writing $!", \@msg_log);
    }
    my @necessary_fields = ("ID",
                            "Avg_fold",
                            "Length",
                            "Ref_GC");
    open(IN, "<", $infile) or error ("Cannot open file $infile: $!", \@msg_log);
    ## Check the first line for comment character at beginning of line ########
    my $firstline = <IN>;
    $linecount++;
    chomp $firstline; # Strip EOL character
    if ($firstline =~ m/^#(.*)/) {
        msg ("Covstats file $infile: first line has comment character \#", \@msg_log);
        if ($do_out == 1) {
            msg ("Offending character will be removed in output", \@msg_log);
        } else {
            msg ("Covstats file $infile: You must remove the comment character \# from the first line", \@msg_log);
        };
        $firstline = $1;
    }
    ## Check fields in first line (should be header line)
    my @firstline_split = split("\t",$firstline);
    my %firstline_hash;
    for (my $i=0; $i < scalar(@firstline_split); $i++) {
        $firstline_hash{$firstline_split[$i]} = $i;
    }
    foreach my $necessary (@necessary_fields) {
        if (!exists $firstline_hash{$necessary}) {
            error ("Covstats file $infile: compulsory field $necessary missing from header line", \@msg_log);
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
        ## Check that the line has same number of fields as header ###########
        if (scalar @linesplit != scalar (keys %firstline_hash)) {
            $fatalerr++;
            error ("Covstats file $infile line $linecount: Number of fields does not match header.", \@msg_log);
        }
        ## Check that Avg_fold field is a number above zero ###################
        if (!is_number($linesplit[$firstline_hash{"Avg_fold"}])
            || $linesplit[$firstline_hash{"Avg_fold"}] < 0 ) {
            $fatalerr++;
            error ("Covstats file $infile line $linecount: Avg_fold field is not a positive number", \@msg_log);
        }
        ## Check that Length field is a number above zero #####################
        if (!is_number($linesplit[$firstline_hash{"Length"}])
            || $linesplit[$firstline_hash{"Length"}] < 0) {
            $fatalerr++;
            error ("Covstats file $infile line $linecount: Length field is not a positive number", \@msg_log);
        }
        ## Check that Ref_GC field is a number between 0, 1 ###################
        if (!is_number($linesplit[$firstline_hash{"Ref_GC"}])
            || $linesplit[$firstline_hash{"Ref_GC"}] > 1
            || $linesplit[$firstline_hash{"Ref_GC"}] < 0) {
            $fatalerr++;
            error ("Covstats file $infile line $linecount: Ref_GC field is not a number between 0, 1", \@msg_log);
        }
        ## If this is the first covstats file, hash the scaffolds #############
        if ($first == 1) {
            $scaffolds_hash{$linesplit[$firstline_hash{"ID"}]} = 1;
        }
        ## Otherwise check that scaffold IDs  match first file ################
        elsif ($first != 1) {
            if (!defined $scaffolds_hash{$linesplit[$firstline_hash{"ID"}]} ) {
                $fatalerr++;
                error ("Covstats file $infile line $linecount: Scaffold ID doesn't match any IDs found in the first covstats file", \@msg_log);
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
        msg ("Covstats file $infile: No fatal errors detected", \@msg_log);
    }
}

sub check_mark { # Also works for SSU tables
    my ($infile, $type) = @_;
    my ($invol,$indir,$infilename) = File::Spec->splitpath($infile);
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
    # If output dir given, write corrected files #
    my $outfile;
    if ($do_out==1) {
        $outfile = $outdir."/".$infilename.".mod"; # Filename for output file
        msg ("Writing corrected marker table to $outfile", \@msg_log);
        open(OUT, ">", $outfile) or error ("Cannot open $outfile for writing: $!", \@msg_log);
    }
    ## Open input file ########################################################
    open(IN, "<", $infile) or error ("Cannot open file $infile: $!", \@msg_log);
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
            error ("Mark file $infile: field $necessary missing from header line", \@msg_log);
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
            msg ("Mark file $infile line $linecount: Problematic character \' found.", \@msg_log);
            msg ("Stripping it from line...", \@msg_log) unless $do_out != 1;
        }
        my @linesplit = split("\t",$line);
        ## Check that number of fields matches header #########################
        if (scalar @linesplit != scalar (keys %firstline_hash)) {
            $fatalerr++;
            error ("Mark file $infile line $linecount: Number of fields does not match header", \@msg_log);
        }
        ## Check that scaffold ID matches those in the first covstats file ####
        if (!defined $scaffolds_hash{$linesplit[$firstline_hash{"scaffold"}]}) {
            $fatalerr++;
            error ("Mark file $infile line $linecount: Scaffold ID doesn't match any ID in covstats file", \@msg_log);
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
        msg ("$type file $infile: No fatal errors detected", \@msg_log);
    }
}

sub check_trna {
    my ($infile) = @_;
    my $fatalerr=0;
    my $linecount=0;
    my $header;
    open(IN, "<", $infile) or error ("Cannot open file $infile: $!", \@msg_log);
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
            error ("tRNA file $infile line $linecount: Scaffold ID doesn't match any ID in covstats file" , \@msg_log);
        }
        if (!defined $linesplit[4]) {
            error ("tRNA file $infile line $linecount: tRNA type appears to be missing", \@msg_log);
        }
    }
    close(IN);
    if ($fatalerr==0) {
        msg ("tRNA file $infile: No fatal errors detected", \@msg_log);
    }
    
}

sub is_number {
    shift =~ /^\s*-?\d+\.?\d*\s*$/;
}

sub msg {
    my ($msg,   # Message to print
        $aref   # Message log
        ) = @_;
    my $time = localtime;
    my $line = "[".$time->hms."] $msg";
    push @$aref, $line;
    print STDERR "$line\n";
}

sub error {
    # Like msg except says "ERROR" in front
    my ($msg,
        $aref
       ) = @_;
    msg ("ERROR: $msg", $aref);
}