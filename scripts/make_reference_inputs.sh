#!/bin/sh
egrep -w "CDS|start_codon|stop_codon" $1 > CDS_ref.ucsc.gtf
gtfToGenePred CDS_ref.ucsc.gtf CDS_ref.ucsc.genepred
genePredToBed CDS_ref.ucsc.genepred > CDS_ref.ucsc.bed
perl new_ref_ribo_reference.pl CDS_ref_gene_transc_ref.rsem CDS_ref.ucsc.bed > result.bed
sed "s/___//" result.bed > riboprofiling_ref_CDS.bed
bedToGenePred result.bed result.genepred
sed "s/\tchr/___chr/" result.genepred | awk -F "___" '{a++; print $1"\t"$3"\t"a"\t"$2}' > result.genepred.renamed
genePredToGtf file result.genepred.renamed riboprofiling_ref_CDS.gtf
grep CDS riboprofiling_ref_CDS.gtf > riboprofiling_ref_CDSonly.gtf
sed "s/;//;s/\"//g" riboprofiling_ref_CDSonly.gtf | awk '{OFS="\t"; print $1,$4-1,$5,$10,$8,$7}' > riboprofiling_ref_CDSonly.bed
perl ref_codon_pos.pl riboprofiling_ref_CDSonly.bed > nostart15cd_ref_ucsc.codon.starts.bed
grep split nostart15cd_ref_ucsc.codon.starts.bed > nostart15cd_ref_ucsc.splitcodons.bed
bedtools intersect -wo -a nostart15cd_ref_ucsc.splitcodons.bed -b nostart15_ref_ucsc.cds.bed > tmp
bedtools intersect -wo -a tmp -b riboprofiling_ref_CDS.bed | sed "s/CDS_//g;s/_split//" | awk '{if ($4==$10 && $4==$17) print $0}' > tmpp
perl split_codons.pl tmpp > nostart15cd_ref_split_codons_4pos.txt
grep -v split nostart15cd_ref_ucsc.codon.starts.bed | awk '{if ($6=="+") print $1,$2,$2+3,$4,$5,$6; else print $1,$2-2,$3,$4,$5,$6}' | sed "s/ /\t/g" > nostart15cd_ref_ucsc.wholecodons.bed
perl uniq_ucsc_codons.pl sacCer3/sacCer3.fa nostart15cd_ref_ucsc.wholecodons.bed > nostart15cd_ref_ucsc.codons.aa.bed
perl get_splitcodon_seq.pl sacCer3/sacCer3.fa nostart15cd_ref_split_codons_4pos.txt >> nostart15cd_ref_ucsc.codons.aa.bed
sed -i "s/CDS_//;s/_split//" nostart15cd_ref_ucsc.codons.aa.bed
