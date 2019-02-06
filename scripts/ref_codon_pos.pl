#!/usr/bin/perl
#####################################
# Program: codon_pos.pl  -  Date: Mon Dec 15 13:14:20 EST 2014
# Autor: Elisa Donnard
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#####################################
my $gene = "";
my $chr = "";
my $cds_start = "";
my $cds_end = "";
my $frame = "";
my $strand = "";
my $codon_start = "";
my $codon_end = "";
my %stop = ();
my %stopline = ();

open (IN0,"<$ARGV[0]");

while (<IN0>) {
    chomp $_;
    @tmp = split (/\t/, $_);
    $gene = $tmp[3];
    $chr = $tmp[0];
    $cds_start = $tmp[1];
    $cds_end = $tmp[2];
    $frame = $tmp[4];
    $strand = $tmp[5];
    if ($strand eq "+") {
	if ($frame == 0) {
	    $codon_start = $cds_start;
	    for ($i = $codon_start; $i < $cds_end; $i+=3) {
		$codon_end = $i + 1;
		$test_end = $i + 3;
		if ($test_end > $cds_end) {
		    print "$chr\t$i\t$codon_end\t$gene"."_split\t$frame\t$strand\n";
		}
		else {
		    print "$chr\t$i\t$codon_end\t$gene\t$frame\t$strand\n";
		}
	    }
	}
	else {
	    $codon_start = ($cds_start + $frame);
	    for ($i = $codon_start; $i < $cds_end; $i+=3) {
		$codon_end = $i + 1;
		$test_end = $i + 3;
		if ($test_end > $cds_end) {
		    print "$chr\t$i\t$codon_end\t$gene"."_split\t$frame\t$strand\n";
		}
		else {
		    print "$chr\t$i\t$codon_end\t$gene\t$frame\t$strand\n";
		}
	    }
	}
	#keep last pos CDS + genes
	if ($stop{$gene}) {
	    if ($stop{$gene} < $cds_end) {
		$stop{$gene} = $cds_end;
		$stop_end = $cds_end + 1;
		$stopline{$gene} = "$chr\t$cds_end\t$stop_end\t$gene\t$frame\t$strand\n";
	    }
	}
	else {
	    $stop{$gene} = $cds_end;
	    $stop_end = $cds_end + 1;
	    $stopline{$gene} = "$chr\t$cds_end\t$stop_end\t$gene\t$frame\t$strand\n";
	}
    }
    else {
	if ($frame == 0) {
	    $codon_start = $cds_end;
	    for ($j = $codon_start; $j > $cds_start; $j-=3) {
		$codon_end = $j - 1;
		$test_end = $j - 3;
		if ($test_end < $cds_start) {
		    print "$chr\t$codon_end\t$j\t$gene"."_split\t$frame\t$strand\n";
		}
		else {
		    print "$chr\t$codon_end\t$j\t$gene\t$frame\t$strand\n";
		}
	    }
	}
	else {
	    $codon_start = ($cds_end - $frame);
	    for ($j = $codon_start; $j > $cds_start; $j-=3) {
		$codon_end = $j - 1;
		$test_end = $j - 3;
		if ($test_end < $cds_start) {
		    print "$chr\t$codon_end\t$j\t$gene"."_split\t$frame\t$strand\n";
		}
		else {
		    print "$chr\t$codon_end\t$j\t$gene\t$frame\t$strand\n";
		}
	    }
	}
	#keep first pos CDS - genes
	if ($stop{$gene}) {
	    if ($stop{$gene} > $cds_start) {
		$stop{$gene} = $cds_start; 
		$stop_end = $cds_start - 1;
		$stopline{$gene} = "$chr\t$stop_end\t$cds_start\t$gene\t$frame\t$strand\n";
	    }
	}
	else {
	    $stop{$gene} = $cds_start;
	    $stop_end = $cds_start - 1;
	    $stopline{$gene} = "$chr\t$stop_end\t$cds_start\t$gene\t$frame\t$strand\n";
	}
    }
}

foreach $k (keys %stopline) {
    print "$stopline{$k}";
}
