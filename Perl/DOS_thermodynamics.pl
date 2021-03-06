#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-07-21 15:26:01 mtw>
#
# Calculate thermodynamic properties from {sampled|reference} DOS
#
# usage: DOS_thermodynamics.pl -i <myDOS>
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2014 Michael Thomas Wolfinger <michael@wolfinger.eu>
# *  All rights reserved
# *
# *  This program is free software: you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation, either version 3 of the License, or
# *  (at your option) any later version.
# *
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program. If not, see <http://www.gnu.org/licenses/>.
# *
# *  This copyright notice MUST APPEAR in all copies of the script!
# ***********************************************************************

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($infile,$bn,$dir,$ext,$cmd);
my $K0 = 273.15;
my $T = 37+$K0;
my $gasconst =  1.98717;  # in [cal/K]
#my $kT = 0.00198717*4.16*$T;
my $kT = $T*$gasconst/1000.;
my $Z = 0.;
my $l = 0;
my ($debug,$fromsub,$fromdos) = (0)x3;
my %inDOS  = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("i=s"           => \$infile,
			   "s"             => sub{$fromsub = 1},
			   "d"             => sub{$fromdos = 1},
                           "-help"         => \&usage,
                           "v");
die "Cannot select both -s AND -d. Please choose one of them"
  if ($fromsub ==1 && $fromdos ==1);
die "ERROR in [DOS_thermodynamics]: no DOS input file provided"
  unless(defined $infile);
die "ERROR in [DOS_thermodynamics]: input DOS file not found"
  unless (-e $infile);

unless ($infile =~ /^\// || $infile =~ /^\.\//){$infile = "./".$infile;}

# parse input file (SUBOPT or DOS)
open(IN, "<", $infile) or die $!;
if ($fromdos == 1){ # parse DOS file from RNAwl
  while(<IN>){
    chomp;
    next if ($_ =~ /^#/);
    my $line = $_;
    die "error while parsing reference DOS file\n ---> $line <---\n"
      unless ($line =~ /(\-?\d+\.\d+)\s+(\d+)/);
    $l++;
    my $energy = $1;
    my $g_E = $2;
    $Z += $g_E * exp(-1*$energy/$kT);
  }
}
elsif ($fromsub == 1){ # parse RNAsubopt output
  while(<IN>){
    chomp;
    next if ($_ =~ /^[AUGC]/);
    my $line = $_;
    die "error while parsing subopt file\n ---> $line <---\n"
      unless ($line =~ /[\.\(\)]+\s+(-?\d+\.\d+)/);
    $l++;
    my $energy = $1;
    $Z += exp(-1*$energy/$kT);
  }
}

else{ die "Could not parse input file. Please choose -s OR -d exiting
...\n";}

print "Processed $l lines\n";
my $foo = -1*$kT*log($Z);
close(IN);
print "Z = $Z\n";
print "kTln(Z)= $foo\n";


sub usage {
  print <<EOF;
Calculate thermodynamic properties from subopt or {sampled|reference}
DOS

usage: $0 [options]
program specific options:                                default:
 -i      <file>   input file                             ($infile)
 -d               input file is {sampled|reference} DOS  ($fromdos)
 -s               input file is subopt                   ($fromsub)
 -d               debug output                           ($debug)
 -help            print this information

EOF
exit;
}
