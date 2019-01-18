#!/usr/bin/perl
#####################################
# Program: get_splitcodon_seq.pl  -  Date: Wed Jul 22 20:41:03 EST 2015
# Autor: Elisa Donnard
#
# License: GPL - http://www.gnu.org/licenses/gpl.html
#
#####################################

my %dna = ();
my $chr = "";
my $start = "";
my $codon = "";
my $end = "";
my $strand = "";
my $toprint = "";
my %print = ();
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
$/ = ">";

while (<IN0>) { #genome fasta file
    chomp $_;
    my $seq = $_;
    if ($seq ne "") {
	if ($seq != /^>/) {
	    $seq =~ s/^/>/;
	}
	my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header print "$id\n";
	$seq =~ s/^>*.+\n//;  # remove FASTA header
	$seq =~ s/\n//g;  # remove endlines              
	if ($id !~ /Un|gl/) {
	    $dna{$id} = $seq;
	}
    }
}
close (IN0);

open (IN1,"<$ARGV[1]"); # txt file with splitcodon positions
$/ = "\n";

while (<IN1>) {
    chomp $_;
    @tmp = split (/\t/, $_);
    $gene = $tmp[6];
    $chr = $tmp[0];
    $start1 = $tmp[1];
    $end1 = $tmp[2];
    $start2 = $tmp[3];
    $a = $end1 - $start1;
    $end2 = $tmp[4];
    $strand = $tmp[5];
    $b = $end2 - $start2;
    $codona = substr($dna{$chr}, $start1, $a);
    $codonb = substr($dna{$chr}, $start2, $b);
    $codon = $codona.$codonb;
    $codon =~ tr/actg/ACTG/;
    if ($_ !~ /\+/) { # reverse complement for rev strand reads
	$codon =~ tr/ACGT/TGCA/;
	$codon = reverse($codon);
    }
    $aa = $genetic_code{$codon};
    $toprint = "$chr\t$start1\t$end2\t$strand\t$gene\t$codon\t$aa";
    $print{$toprint} = 1;
}

foreach $k (keys %print) {
    print "$k\n";
}
