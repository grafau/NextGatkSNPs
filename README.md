# NextGatkSNPs
Nextflow pipelines for short read mapping and SNP calling - complies with GATK3.8 best practices. 

This is NextFlow pipeline that runs all processes to map short reads in fastq format to reference genome, 
filter and call SNPs. Branch princeV0.1 was optimized to work with slurm submission system. Please note 
that you will have to install bioinformatic packages (bioconda installations worked for me) and change 
configuration files (.conf) to reflect paths to your installations. In particular you will need BWA, 
PICARD, GATK, BGZIP and TABIX. This pipeline complies with GATK3.8 best practices and is unlikely to 
work with version 4.0+ without serious changes.

## PREPARE FILES FOR CALLING SNPS

The **process.nf** pipeline will run following processes:

  - map all fastq short reads to a reference genome independently,
  - sort all sam files and generate bam files,
  - merge mutliple bam files for the same indvidual,
  - remove PCR duplicates,
  - check if bam files are correct,
  - index merged bam file,
  - call haplotypes in g.vcf format,
  - check if vcf files are correct,
  - compare and index indvidual g.vcf files,
  - combine multiple g.vcf files.

This pipeline requires a configuration file, process.conf to find packages and modify cpu/time/memory usage. 
Additionally, it requires a comma separated "batch" file, which will include information about your fastq files. 
Finally, it will need a reference genome in fasta format indexed with bwa. 

Batch files format is as follows: 

  - column 1: individual names,
  - column 2: library name,
  - column 3: fastq files with P1 read,
  - column 4: fastq files with P2 reads.
  

If all files are present, you can run this pipeline as follows:

```
nextflow run -c process_3k_v0.3.conf process_3k_v0.3.nf --ref ref.fa --list batch.csv
```

## CALL SNPS

Now that you have (multiple) combined g.vcf files, you can use them to jointly call SNPs

The **call.nf** pipeline will run following processes:

  - split g.vcf files into chromosomes,
  - merge multiple combined g.vcf files for each chromosome,
  - genotype g.vcf haplotypes.
  

This pipeline requires a configuration file, call.conf to find packages and modify cpu/time/memory usage. Additionally, it requires a comma separated "chrom" file, which will include information about your chromosomes. Finally, it will need a reference genome in fasta format indexed with bwa.

Chrom files format is as follows: 

  - column 1: chromosome name,
  - column 2: bed file containing chromosome range,
  - column 3: combined g.vcf.
(you will need entries for all chromosomes for each combined g.vcf file)

If all files are present, you can run this pipeline as follows:

```
nextflow run -c call.conf call.nf --ref ref.fa --list chrom.txt
```

