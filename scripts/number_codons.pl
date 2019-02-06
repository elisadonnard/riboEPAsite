#!/usr/bin/perl
#############################################################
# Program: number_codons.pl  -  Date: Sun Mar 26 16:33:13 EDT 2017
# Author: Elisa Donnard
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#############################################################

open (IN0,"<$ARGV[0]");

my %ref = ();
my %line = ();

while (<IN0>) {
    chomp $_;
    @tmp = split (/\t/, $_);
    $name = $tmp[4];
    $ref{$name} .= $tmp[1]."|";
    $line{$name}{$tmp[1]} = $_;
}

foreach $k (keys %ref) {
    $starts = $ref{$k};
    chomp $starts;
    @st = split (/\|/, $starts);
    @sorted = sort {$a <=> $b} @st;
    for ($i=0; $i<=$#sorted; $i++) {
	$start = $sorted[$i];
	print $i."\t".$line{$k}{$start}."\n";
    }
}
