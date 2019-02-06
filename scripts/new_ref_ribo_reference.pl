 #!/usr/bin/perl
#############################################################
# Program: ribo_reference.pl  -  Date: Tue Nov 24 15:54:05 EST 2015
# Author: Elisa Donnard
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
# creates a reference BED without the first 15 and last 5 codons of each transcript/gene
#
#############################################################

my %ref = ();

open (IN0,"<$ARGV[0]");

while (<IN0>) { #rsem transcript name and  gene name reference table
    chomp $_;
    @tmp = split (/\t/, $_);
    $gene = $tmp[0];
    $transc = $tmp[1];
    $ref{$transc} = $gene;
}
close (IN0);

open (IN1,"<$ARGV[1]");

my @sizes = ();
my @starts = ();
my $remain = "";
my $spos = "";
my $epos = "";
my $vnewsizes = "";
my $vnewstarts = "";
my @newsizes = ();
my @newstarts = ();
my $finalexons = "";
my $finalsizes = "";
my $finalstarts = "";
my $refpos = "";

while (<IN1>) { #CDS bed12 file
    chomp $_;
    my @tmp = split (/\t/, $_);
    my $id = $tmp[3];
    my $strand = $tmp[5];
    my $exons = $tmp[9];
    if ($strand eq "+") { # first the forward strand genes...
	@sizes = split (/,/, $tmp[10]);
	@starts = split (/,/, $tmp[11]);
	$total = 0;
	$stai = 0;
	$stot = 0;
	$endi = 0;
	$etot = 0;
	$te = 0;
	$ns = 0;
	$ne = 0;
	for ($i=0; $i<=$#sizes; $i++) {
	    $total += $sizes[$i];
	    if ($total > 45 && $ns == 0) {
		$stai = $i;
		$stot = $total; # total up to here
		$ns = 1;
	    }
	    elsif ($total == 45 && $ns == 0) {
		$stai = $i + 1;
		$stot = $total + $sizes[$stai]; # total up to here
		$ns = 1;
	    }
	}
	if ($total > 60) {
	    $remain = 45 - ($stot - $sizes[$stai]); # how much we need to remove from this exon  
	    ##### 
	    $spos = $tmp[1] + $starts[$stai] + $remain;
	    $vnewsizes = $sizes[$stai] - $remain;
	    $vnewsizes .= ",";
	    $vnewstarts = "0,";
	    $refpos = $starts[$stai] + $remain;
	    if ($stai < $#sizes)  {
		for ($i=($stai + 1); $i<=$#sizes; $i++) { # add the remaining exons
		    $vnewsizes .= $sizes[$i].",";
		    $n = $starts[$i] - $refpos;
		    $vnewstarts .= $n.",";
		}
	    }
	    # then comes the stop codon...
	    # use arrays made in start step
	    @newsizes = split (/,/, $vnewsizes);
	    @newstarts = split (/,/, $vnewstarts);
	    for ($i=$#newsizes; $i>=0; $i--) {
		$te += $newsizes[$i];
		if ($te > 15 && $ne == 0) {
		    $endi = $i;
		    $etot = $te;
		    $ne = 1;
		}
		elsif ($te == 15 && $ne == 0) {
		    $endi = $i - 1;
		    $etot = $te + $newsizes[$endi];
		    $ne = 1;
		}
	    }
	    $remain = 15 - ($etot - $newsizes[$endi]);
	    #####
	    $epos = $spos + $newstarts[$endi] + $newsizes[$endi] - $remain;
	    $finalsizes = "";
	    $finalstarts = "";
	    if ($endi > 0) {
		for ($i=0; $i<$endi; $i++) { # add the previous exons
		    $finalsizes .= $newsizes[$i].",";
		    $finalstarts .= $newstarts[$i].",";
		}
	    }
	    #new last exon
	    $finalsizes .= $newsizes[$endi] - $remain;
	    $finalstarts .= $newstarts[$endi];
	    @countexons = split (/,/, $finalstarts);
	    $finalexons = scalar(@countexons);
	    $finalsizes .= ",";
	    $finalstarts .= ",";
	}
    }
    else { # now the reverse strand genes...
	@sizes = split (/,/, $tmp[10]);
	@starts = split (/,/, $tmp[11]);
	$total = 0;
	$stai = 0;
	$stot = 0;
	$endi = 0;
	$etot = 0;
	$te = 0;
	$ns = 0;
	$ne = 0;
	for ($i=0; $i<=$#sizes; $i++) {
	    $total += $sizes[$i];
	    if ($total > 15 && $ns == 0) {
		$stai = $i;
		$stot = $total; # total up to here
		$ns = 1;
	    }
	    elsif ($total == 15 && $ns == 0) {
		$stai = $i + 1;
		$stot = $total + $sizes[$stai]; # total up to here
		$ns = 1;
	    }
	}
	if ($total > 60) {
	    $remain = 15 - ($stot - $sizes[$stai]); # what is missing after the $stai exon 
	    #####
	    $spos = $tmp[1] + $starts[$stai] + $remain;
	    $vnewsizes = $sizes[$stai] - $remain;
	    $vnewsizes .= ",";
	    $vnewstarts = "0,";
	    $refpos = $starts[$stai] + $remain;
	    if ($stai < $#sizes) {
		for ($i=($stai + 1); $i<=$#sizes; $i++) { # add the remaining exons
		    $vnewsizes .= $sizes[$i].",";
		    $n = $starts[$i] - $refpos;
		    $vnewstarts .= $n.",";
		}
	    }
	    # then comes the stop codon...
	    # use arrays made in start step
	    @newsizes = split (/,/, $vnewsizes);
	    @newstarts = split (/,/, $vnewstarts);
	    for ($i=$#newsizes; $i>=0; $i--) {
		$te += $newsizes[$i];
		if ($te > 45 && $ne == 0) {
		    $endi = $i;
		    $etot = $te;
		    $ne = 1;
		}
		elsif ($te == 45 && $ne == 0) {
		    $endi = $i - 1;
		    $etot = $te + $newsizes[$endi];
		    $ne = 1;
		}
	    }
	    $remain = 45 - ($etot - $newsizes[$endi]);
	    #####
	    $epos = $spos + $newstarts[$endi] + $newsizes[$endi] - $remain;
	    $finalsizes = "";
	    $finalstarts = "";
	    if ($endi > 0) {
		for ($i=0; $i<$endi; $i++) { # add the previous exons
		    $finalsizes .= $newsizes[$i].",";
		    $finalstarts .= $newstarts[$i].",";
		}
	    }
	    #new last exon
	    $finalsizes .= $newsizes[$endi] - $remain;
	    $finalstarts .= $newstarts[$endi];
	    @countexons = split (/,/, $finalstarts);
	    $finalexons = scalar(@countexons);
	    $finalsizes .= ",";
	    $finalstarts .= ",";
	}
    }
    #after all coordinate changes print line
    if ($total > 60) {
	print "$tmp[0]\t$spos\t$epos\t$id"."___"."$ref{$id}\t$tmp[4]\t$strand\t$spos\t$epos\t$tmp[8]\t$finalexons\t$finalsizes\t$finalstarts\n";
    }
}
