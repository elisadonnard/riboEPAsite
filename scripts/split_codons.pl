#!/usr/bin/perl
#####################################
# Program: split_codons.pl  -  Date: Tue Jul 21 01:00:07 EDT 2015
# Autor: Elisa Donnard
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#####################################

open (IN0,"<$ARGV[0]"); # splitcodons.bed file intersect with cds.bed file intersect with bed12 file

while (<IN0>) {
    chomp $_;
    @tmp = split (/\t/, $_);
    $strand = $tmp[5];
    @blocks = split (/,/, $tmp[24]);
    @sizes = split (/,/, $tmp[23]);
    if ($strand eq "+") {
	$start1 = $tmp[1];
	$end1 = $tmp[8];
	$len = $end1 - $start1;
	for ($i=0; $i<=$#blocks; $i++) {
	    if ($tmp[7] == ($tmp[14] + $blocks[$i])) {
		$n = $i + 1;
	    }
	}
	$start2 = $tmp[14] + $blocks[$n];
	$end2 = $start2 + (3 - $len);
	print "$tmp[0]\t$start1\t$end1\t$start2\t$end2\t$strand\t$tmp[3]\n";
    }
    else {
	$start1 = $tmp[2];
	$end1 = $tmp[7];
	$len = $start1 - $end1;
	for ($i=0; $i<=$#blocks; $i++) {
	    if ($tmp[7] == ($tmp[14] + $blocks[$i])) {
		$n = $i - 1;
	    }
	}
	$start2 = $tmp[14] + $blocks[$n] + $sizes[$n];
	$end2 = $start2 - (3 - $len);
	print "$tmp[0]\t$end2\t$start2\t$end1\t$start1\t$strand\t$tmp[3]\n";
    }
}
