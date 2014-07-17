#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-07-17 17:08:34 mtw>
#
# Evaluate (sampled) Density of States (DOS) vs. reference DOS
#
# usage: eval_sampledDOS.pl -s samples.dos -r reference.dos
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
my ($infile_s,$infile_r,$bn_r,$dir_r,$ext_r,$bn_s,$dir_s,$ext_s,$cmd);
my $debug         = 0;
my %refDOS        = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("s=s"           => \$infile_s,
			   "r=s"           => \$infile_r,
			   "d"             => sub{$debug=1},
                           "-help"         => \&usage,
                           "v");

die "ERROR in [eval_sampledDOS]: no sampled DOS input file provided"
  unless(defined $infile_s);
die "ERROR in [eval_sampledDOS]: no reference DOS input file provided"
  unless(defined $infile_r);
die "ERROR in [eval_sampledDOS]: input sampled DOS file not found"
  unless (-e $infile_s);
die "ERROR in [eval_sampledDOS]: input reference DOS file not found"
  unless (-e $infile_r);

unless ($infile_r =~ /^\// || $infile_r =~ /^\.\//){$infile_r = "./".$infile_r;}
unless ($infile_s =~ /^\// || $infile_s =~ /^\.\//){$infile_s = "./".$infile_s;}

print "#parsing reference DOS $infile_r\n";
print "#parsing sampled DOS $infile_s\n";

# 1) parse reference DOS
open(REF, "<", $infile_r) or die $!;
while(<REF>){
  chomp;
  next if ($_ =~ /^#/);
  my $line = $_;
  die "error while parsing reference DOS file\n ---> $line <---\n"
    unless ($line =~ /(\-?\d+\.\d+)\s+(\d+)/);
  $refDOS{$1}=log($2);
}
close(REF);

#print Dumper(\%refDOS);

# 2) parse sampled DOS
open(SAM, "<", $infile_s) or die $!;
while(<SAM>){
  my ($energy,$sDOS,$rDOS,$relErr);
  chomp;
  next if ($_ =~ /^#/);
  die "error while parsing sampled DOS file\n ---> $_ <---\n"
    unless ($_ =~ /(\-?\d+\.\d+)\s+(\d+\.\d)/);
  die ("Cannot take log of 0 at  $infile_s\n$!\n") if($2 == 0.);
  $energy = $1;
  $sDOS = log($2); # value in sampled DOS

  die "error: energy $energy not present in reference DOS\n"
    unless exists $refDOS{$energy};
  $rDOS = $refDOS{$energy}; # value in reference DOS

  # 3) evaluate statistics,ie compute relative error
  if($rDOS != 0){
    $relErr = abs(($sDOS-$rDOS)/$rDOS);
    print "$energy\t$relErr\n";# = abs(($sDOS-$rDOS)/$rDOS)\ \n";
  }
}

close(SAM);


sub usage {
  print <<EOF;
Evaluate (sampled) Density of States vs. reference DOS

usage: $0 [options]
program specific options:                                default:
 -r      <file>   reference DOS                          ($infile_r)
 -s      <file>   sampled DOS                            ($infile_s)
 -d               debug output                           ($debug)
 -help            print this information

EOF
exit;
}
