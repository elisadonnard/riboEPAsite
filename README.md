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

## Step3: 
