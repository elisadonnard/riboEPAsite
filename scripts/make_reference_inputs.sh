#!/bin/sh
mkdir tmp_files
egrep -w "CDS|start_codon|stop_codon" $1 > tmp_files/CDS_ref.ucsc.gtf
grep -w CDS tmp_files/CDS_ref.ucsc.gtf | sed "s/CDS/exon/" > tmp_files/oo
cat tmp_files/oo >> tmp_files/CDS_ref.ucsc.gtf
egrep "CDS|codon" tmp_files/CDS_ref.ucsc.gtf > tmp_files/tmp.gtf
gtfToGenePred tmp_files/tmp.gtf tmp_files/tmp.genepred
genePredToBed tmp_files/tmp.genepred > tmp_files/tmp.bed
perl new_ref_ribo_reference.pl tmp_files/CDS_ref_gene_transc_ref.rsem tmp_files/tmp.bed > tmp_files/result.bed
sed "s/___//" tmp_files/result.bed > tmp_files/riboprofiling_ref_CDS.bed
bedToGenePred tmp_files/result.bed tmp_files/result.genepred
sed "s/\tchr/___chr/" tmp_files/result.genepred | awk -F "___" '{a++; print $1"\t"$3"\t"a"\t"$2}' > tmp_files/result.genepred.renamed
genePredToGtf file tmp_files/result.genepred.renamed tmp_files/riboprofiling_ref_CDS.gtf
grep CDS tmp_files/riboprofiling_ref_CDS.gtf > tmp_files/riboprofiling_ref_CDSonly.gtf
sed "s/;//;s/\"//g" tmp_files/riboprofiling_ref_CDSonly.gtf | awk '{OFS="\t"; print $1,$4-1,$5,$10,$8,$7}' > tmp_files/riboprofiling_ref_CDSonly.bed
perl ref_codon_pos.pl tmp_files/riboprofiling_ref_CDSonly.bed > nostart15cd_ref_ucsc.codon.starts.bed
grep split nostart15cd_ref_ucsc.codon.starts.bed > tmp_files/nostart15cd_ref_ucsc.splitcodons.bed
grep CDS tmp_files/riboprofiling_ref_CDSonly.gtf | sed "s/;//;s/\"//g" | awk '{print $1,$4-1,$5,$3"_"$10,$8,$7}' | sed "s/ /\t/g" > tmp_files/nostart15_ref_ucsc.cds.bed
bedtools intersect -wo -a tmp_files/nostart15cd_ref_ucsc.splitcodons.bed -b tmp_files/nostart15_ref_ucsc.cds.bed > tmp_files/tmp
bedtools intersect -wo -a tmp_files/tmp -b tmp_files/riboprofiling_ref_CDS.bed | sed "s/CDS_//g;s/_split//" | awk '{if ($4==$10 && $4==$17) print $0}' > tmp_files/tmpp
perl split_codons.pl tmp_files/tmpp > tmp_files/nostart15cd_ref_split_codons_4pos.txt
grep -v split nostart15cd_ref_ucsc.codon.starts.bed | awk '{if ($6=="+") print $1,$2,$2+3,$4,$5,$6; else print $1,$2-2,$3,$4,$5,$6}' | sed "s/ /\t/g" > tmp_files/nostart15cd_ref_ucsc.wholecodons.bed
perl uniq_ucsc_codons.pl $2 tmp_files/nostart15cd_ref_ucsc.wholecodons.bed > nostart15cd_ref_ucsc.codons.aa.bed
perl get_splitcodon_seq.pl $2 tmp_files/nostart15cd_ref_split_codons_4pos.txt >> nostart15cd_ref_ucsc.codons.aa.bed
sed -i "s/CDS_//;s/_split//" nostart15cd_ref_ucsc.codons.aa.bed
perl number_codons.pl nostart15cd_ref_ucsc.codons.aa.bed > numbered_nostart15cd_ref_ucsc.codons.aa.bed
rm -rf tmp_files
