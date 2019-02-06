# riboEPAsite
Defining the codon in E/P/A sites from ribosome profiling data

## Getting Started

You will need alignment files from ribosome profiling of your interest

## Prerequisites

R;
Perl;
bedtools 2.25.0 or higher;
UCSC [Kent Utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)

## Step1: Generate the reference codon input files

Run the command below to generate the necessary reference input for your species of interest
Starting from the GTF annotation file for your species of interest and the genome fasta:

```
make_reference_inputs.sh ref_ucsc.gtf ref.fa
```

## Step2: Define the offset for each sample or batch of samples

Samples show only mild differences in the offset, we noticed it was more robust to determine offsets using all reads from a batch of ribosome profiling libraries prepared at the same time. To find the offsets, use the psite command from the plastid library for example described [here](https://plastid.readthedocs.io/en/latest/examples/p_site.html)

In the end create a tab separated file with the sample name, read size and offset per line

```
#sample		    #read size 	#offset
ribo2015_BY4741_rep1	29	13
ribo2015_BY4741_rep1	30	13
ribo2015_BY4741_rep1	31	13
ribo2015_BY4741_rep1	32	14
ribo2015_BY4741_rep1	33	14
ribo2015_trm4D_rep1	29	13
ribo2015_trm4D_rep1	30	13
ribo2015_trm4D_rep1	31	13
ribo2015_trm4D_rep1	32	14
ribo2015_trm4D_rep1	33	14
```

## Step3: Find the overlap between each read and the codon starts bed file

```
bedtools intersect -bed -s -split -wo -a ribo2015_BY4741_rep1.bam -b nostart15cd_sc3_ucsc.codon.starts.bed > ribo2015_BY4741_rep1.bed
```

## Step4: Use the calculated offsets to determine the correct overlapped reads

Create a tab separated file with your sample name and input bed file name (result from step3)
```
ribo2015_BY4741_rep1	ribo2015_BY4741_rep1.bed
ribo2015_BY4741_rep2	ribo2015_BY4741_rep2.bed
ribo2015_BY4741_rep3	ribo2015_BY4741_rep3.bed
ribo2015_BY4741_rep4	ribo2015_BY4741_rep4.bed
ribo2015_BY4741_rep5	ribo2015_BY4741_rep5.bed
ribo2015_BY4741_rep6	ribo2015_BY4741_rep6.bed
```

Then run:
```
perl sample_offset_filter_reads.pl sample_offset_long samples_bed ribo2015_BY4741_rep1.bed > offsetreads_ribo2015_BY4741_rep1.bed
```

## Step5: 

