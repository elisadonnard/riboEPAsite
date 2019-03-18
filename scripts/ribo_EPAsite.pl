#!/usr/bin/perl
#####################################
# Program: ribo_EPAsite.pl  -  Date: Tue Jun 23 16:59:28 EDT 2015
# Autor: Elisa Donnard
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#####################################

my $sample = $ARGV[1];

open (OUTE,">codonEsite_$sample");
open (OUTP,">codonPsite_$sample");
open (OUTA,">codonAsite_$sample");
open (OUTN,">codonNsite_$sample");

my $lines = ();
my $sites = ();
my $codons = ();
my %num = ();
my %numid = ();
my %dna = ();
my %genetic_code = (
    'GCA'=>'A', #Alanine
    'GCC'=>'A', #Alanine
    'GCG'=>'A', #Alanine
    'GCT'=>'A', #Alanine
    'AGA'=>'R', #Arginine
    'AGG'=>'R', #Arginine
    'CGA'=>'R', #Arginine
    'CGC'=>'R', #Arginine
    'CGG'=>'R', #Arginine
    'CGT'=>'R', #Arginine
    'AAC'=>'N', #Asparagine
    'AAT'=>'N', #Asparagine
    'GAC'=>'D', #Aspartic acid
    'GAT'=>'D', #Aspartic acid
    'TGC'=>'C', #Cysteine
    'TGT'=>'C', #Cysteine
    'GAA'=>'E', #Glutamic acid
    'GAG'=>'E', #Glutamic acid
    'CAA'=>'Q', #Glutamine
    'CAG'=>'Q', #Glutamine
    'GGA'=>'G', #Glycine
    'GGC'=>'G', #Glycine
    'GGG'=>'G', #Glycine
    'GGT'=>'G', #Glycine
    'CAC'=>'H', #Histidine
    'CAT'=>'H', #Histidine
    'ATA'=>'I', #Isoleucine
    'ATC'=>'I', #Isoleucine
    'ATT'=>'I', #Isoleucine
    'TTA'=>'L', #Leucine
    'TTG'=>'L', #Leucine
    'CTA'=>'L', #Leucine
    'CTC'=>'L', #Leucine
    'CTG'=>'L', #Leucine
    'CTT'=>'L', #Leucine
    'AAA'=>'K', #Lysine
    'AAG'=>'K', #Lysine
    'ATG'=>'M', #Methionine
    'TTC'=>'F', #Phenylalanine
    'TTT'=>'F', #Phenylalanine
    'CCA'=>'P', #Proline
    'CCC'=>'P', #Proline
    'CCG'=>'P', #Proline
    'CCT'=>'P', #Proline
    'AGC'=>'S', #Serine
    'AGT'=>'S', #Serine
    'TCA'=>'S', #Serine
    'TCC'=>'S', #Serine
    'TCG'=>'S', #Serine
    'TCT'=>'S', #Serine
    'ACA'=>'T', #Threonine
    'ACC'=>'T', #Threonine
    'ACG'=>'T', #Threonine
    'ACT'=>'T', #Threonine
    'TGG'=>'W', #Tryptophan
    'TAC'=>'Y', #Tyrosine
    'TAT'=>'Y', #Tyrosine
    'GTA'=>'V', #Valine
    'GTC'=>'V', #Valine
    'GTG'=>'V', #Valine
    'GTT'=>'V', #Valine
    'TAA'=>'-', #STOP
    'TAG'=>'-', #STOP
    'TGA'=>'-', #STOP
    );

open (IN0,"<$ARGV[0]");

while (<IN0>) { #numbered.codons.aa.bed
    chomp $_;
    @tmp = split (/\t/, $_);
    $name = $tmp[5];
    $name =~ s/CDS_//;
    $name =~ s/_split//;
    $strand = $tmp[4];
    if ($strand eq "+") {
	$loc = $tmp[1]."_".$tmp[2]."_".$tmp[4]."_".$tmp[5];
    }
    else {
	$loc = $tmp[1]."_".$tmp[3]."_".$tmp[4]."_".$tmp[5];
    }
    $floc = $tmp[1]."_".$tmp[2]."_".$tmp[3]."_".$tmp[4]."_".$tmp[5];
    $dna{$floc} = $tmp[6]."\t".$tmp[7];
    $num{$name}{$tmp[0]} = $floc;
    $numid{$loc} = $tmp[0];
}
close (IN0);

open (IN1,"<$ARGV[1]");

while (<IN1>) {
    chomp $_;
    @tmp = split (/\t/, $_);
    $id = $tmp[3];
    $readsize = $tmp[4];
    $strand = $tmp[5];
    $start = $tmp[1];
    $end = $tmp[2];
    $chr = $tmp[0];
    $name = $tmp[6];
    $name =~ s/CDS_//;
    $name =~ s/_split//; 
    ##########
    # P site #
    ##########
    if ($strand eq "+"){
	$idloc = $chr."_".$start."_".$strand."_".$name;
    }
    else {
	$idloc = $chr."_".$end."_".$strand."_".$name;
    }
    if ($numid{$idloc}) {
	$numP = $numid{$idloc};
	$idloc = $num{$name}{$numP};
	@tmpP = split (/_/, $idloc);
	$chr = $tmpP[0];
	$start = $tmpP[1];
	$end = $tmpP[2];
	$strand = $tmpP[3];
	$codon = $dna{$idloc};
	print OUTP $chr."\t".$start."\t".$end."\t".$name."\t".$readsize."\t".$strand."\t".$name."\t".$codon."\n";
	################
	# A and E site #
	################
	$numE = $numP - 1;
	if ($num{$name}{$numE}) {
	    $idloc = $num{$name}{$numE};
	    @tmpE = split (/_/, $idloc);
	    $chr = $tmpE[0];
	    $start = $tmpE[1];
	    $end = $tmpE[2];
	    $strand = $tmpE[3];
	    $codon = $dna{$idloc};
	    if ($strand eq "+") {
		print OUTE $chr."\t".$start."\t".$end."\t".$name."\t".$readsize."\t".$strand."\t".$name."\t".$codon."\n";
	    }
	    else {
		print OUTA $chr."\t".$start."\t".$end."\t".$name."\t".$readsize."\t".$strand."\t".$name."\t".$codon."\n";
	    }
	}
	$numA = $numP + 1;
        if ($num{$name}{$numA}) {
            $idloc = $num{$name}{$numA};
            @tmpA = split (/_/, $idloc);
            $chr = $tmpA[0];
            $start = $tmpA[1];
            $end = $tmpA[2];
            $strand = $tmpA[3];
            $codon = $dna{$idloc};
            if ($strand eq "+") {
                print OUTA $chr."\t".$start."\t".$end."\t".$name."\t".$readsize."\t".$strand."\t".$name."\t".$codon."\n";
            }
            else {
                print OUTE $chr."\t".$start."\t".$end."\t".$name."\t".$readsize."\t".$strand."\t".$name."\t".$codon."\n";
            }
        }
	##########
	# Next 3 #
	##########
	if ($strand eq "+") {
	    for ($i=2; $i<=4; $i++) {
		$numN = $numP + $i;
		if ($num{$name}{$numN}) {
		    $idloc = $num{$name}{$numN};
		    @tmpN = split (/_/, $idloc);
		    $chr = $tmpN[0];
		    $start = $tmpN[1];
		    $end = $tmpN[2];
		    $strand = $tmpN[3];
		    $codon = $dna{$idloc};
		    print OUTN $chr."\t".$start."\t".$end."\t".$name."\t".$readsize."\t".$strand."\t".$name."\t".$codon."\n";
		}
	    }
	}
	else {
            for ($i=2; $i<=4; $i++) {
                $numN = $numP - $i;
                if ($num{$name}{$numN}) {
                    $idloc = $num{$name}{$numN};
                    @tmpN = split (/_/, $idloc);
                    $chr = $tmpN[0];
                    $start = $tmpN[1];
                    $end = $tmpN[2];
                    $strand = $tmpN[3];
                    $codon = $dna{$idloc};
                    print OUTN $chr."\t".$start."\t".$end."\t".$name."\t".$readsize."\t".$strand."\t".$name."\t".$codon."\n";
		}
	    }
	}
    }
}
