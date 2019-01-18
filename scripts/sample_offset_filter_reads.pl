#!/usr/bin/perl
#############################################################
# Program: sample_offset_filter_reads.pl  -  Date: Sat Feb  4 16:22:06 EST 2017
# Author: Elisa Donnard
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#############################################################

my %off=();
my %ref=();

open (IN0,"<$ARGV[0]"); # per sample offset for reads of each size

while (<IN0>) {
    chomp $_;
    @tmp = split (/\t/, $_);
    $sample = $tmp[0];
    $readsize = $tmp[1];
    $offset = $tmp[2];
    $off{$sample}{$readsize}=$offset;
}

open (IN1,"<$ARGV[1]"); # sample id to file name ref

while (<IN1>) {
    chomp $_;
    @tmp = split (/\t/, $_);
    $ref{$tmp[1]} = $tmp[0];
}

open (IN2,"<$ARGV[2]"); # sample bed file overlapped with codons
$filename = $ARGV[2];

while (<IN2>) {
    chomp $_;
    @tmp = split (/\t/, $_);
    $t = "\t";
    ###
    @sizes = split (/,/, $tmp[10]);
    @starts = split (/,/, $tmp[11]);
    $readsize = 0;
    for ($i=0; $i<=$#sizes; $i++) {
	$readsize += $sizes[$i];
    }
    $sample = $ref{$filename};
    $offset = $off{$sample}{$readsize};
    $strand = $tmp[5];
    $total = 0;
    $stai = 0;
    $stot = 0;
    $endi = 0;
    $etot = 0;
    $te = 0;
    $ns = 0;
    $ne = 0;
    if ($strand eq "+") {
	$readstart = $tmp[1];
	$codonstart = $tmp[13];
	$codonend = $tmp[14];
	###
	# For adding offset to split reads
	###
	for ($i=0; $i<=$#sizes; $i++) {
	    $total += $sizes[$i];
	    if ($total > $offset && $ns == 0) {
		$stai = $i;
		$stot = $total; # total up to here
		$ns = 1;
	    }
	    elsif ($total == $offset && $ns == 0) {
		$stai = $i + 1;
		$stot = $total + $sizes[$stai]; # total up to here
		$ns = 1;
	    }
	}
	$remain = $offset - ($stot - $sizes[$stai]); # how much we need to offset from the start of this piece
	##### 
	$psite = $readstart + $starts[$stai] + $remain;
	#print "$_ \n ******* $offset $remain $codonstart $psite \n ******** $starts[$stai] $sizes[$stai]\n";
	if ($codonstart == $psite) {
	    $pend = $codonend;
	    print $tmp[0].$t.$psite.$t.$pend.$t.$tmp[3].$t.$readsize.$t.$strand.$t.$tmp[15]."\n";
	}
    }
    else {
	$readstart = $tmp[1];
	$codonstart = $tmp[14];
	$codonend = $tmp[13];
	###
	# For adding offset to split reads
	###
	for ($i=$#sizes; $i>=0; $i--) {
	    $te += $sizes[$i];
	    if ($te > $offset && $ne == 0) {
		$endi = $i;
		$etot = $te;
		$ne = 1;
	    }
	    elsif ($te == $offset && $ne == 0) {
		$endi = $i - 1;
		$etot = $te + $sizes[$endi];
		$ne = 1;
	    }
	}
	$remain = $offset - ($etot - $sizes[$endi]);
        #####
	$psite = $readstart + $starts[$endi] + $sizes[$endi] - $remain;
	#print "$_ \n ******* $offset $remain $codonstart $psite \n ******** $etot $sizes[$endi] \n";
	if ($codonstart == $psite) {
	    $pend = $codonend;
	    print $tmp[0].$t.$pend.$t.$psite.$t.$tmp[3].$t.$readsize.$t.$strand.$t.$tmp[15]."\n";
	}
    }
}
