# riboEPAsite
Defining the codon in E/P/A sites from ribosome profiling data

## Getting Started

You will need alignment files from ribosome profiling of your interest

## Prerequisites

R;
Perl;
bedtools 2.25.0 or higher;

## Step1: Generate the reference codon input files

The reference files for rat and yeast are available in the reference_files folder, for another species run the command below
Starting from the GTF annotation file for your species of interest:

```
make_reference_inputs.sh ref_ucsc.gtf
```

## Step2: Define the offset for each sample or batch of samples

Samples show only mild differences in the offset, we noticed it was more robust to determine offsets using all reads from a batch of ribosome profiling libraries prepared at the same time. To find the offsets, use the psite command from the plastid library for example described [here](https://plastid.readthedocs.io/en/latest/examples/p_site.html)

In the end create a file with the sample name, read size and offset per line

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

```
perl sample_offset_filter_reads.pl sample_offset_long samples_bed ribo2015_BY4741_rep1.bed > offsetreads_ribo2015_BY4741_rep1.bed
```

## Step5: 

