
# Project Riz

### This project focuses on the analysis of genetic diversity in rice after irradiation.

### Folder structure

- **01.RawData/** : Raw sequencing data.
- **metadata/** : Sample information and morphological data 
- **qc/** : Quality control results.
- **fastqc/** : Quality control results.
- **multiqc/** : Quality control results.
- **cleaning/** : to eliminate poor quality bases and adapters if used
- **mapping/** : to align reads on the reference genome

##  FASTQC

### Create a “qc” directory in the working directory in which quality control results will be stored.

```bash
mkdir -p /home/rkouke/Project_Riz/qc
```
### Move to qc directory

```bash
cd /home/rkouke/Project_Riz/qc
```

### Create the “fastqc” directory for outputs

```bash
mkdir fastqc
```

### Move to "fastqc" directory

```bash
cd fastqc
```

### Load the “fastqc” module, taking the version into account

```bash
module load fastqc/0.12.1
```
### Display the help menu for the fastqc command

```bash
fastqc --help
```

### Run quality control command on fastq.gz sequencing data

```bash
fastqc /home/rkouke/01.RawData/*.fastq.gz -o /home/rkouke/Projet_Riz/qc/fastqc
```

### Check content
```bash
ls -lh
```

### Use generated html files to check read quality 

##  MULTIQC

## Create the “multiqc” directory for outputs 

```bash
mkdir -p /home/rkouke/Projet_Riz/qc/multiqc
```

## Move to the “multiqc” directory

```bash
cd /home/rkouke/Projet_Riz/qc/multiqc
```

## Load multiqc module to compress fastqc reports into a single report

```bash
module load multiqc/1.13
```

## Launch analysis

```bash
multiqc /home/rkouke/Projet_Riz/qc/fastqc -o /home/rkouke/Projet_Riz/qc/multiqc
```

### Check content
```bash
ls -lh
```

### Use generated html files to check read quality

## CLEANING

### cleaning depends on multiqc results

## MAPPING

### indexing of the reference genome

## Create a “ref” directory in the working directory in which reference genome will be stored

```bash
mkdir -p /home/rkouke/Projet_Riz/ref
```
## Move to the “ref” directory

```bash
cd /home/rkouke/Project_Riz/ref
```
## download reference genome
### ncbi-datasets-cli module loading

```bash
module load ncbi-datasets-cli/16.27.2
```
## download launch 
```bash
datasets download genome GCF_000147395.1 --include gff3,rna,cds,protein,genome,seq-report
```
### Check download 

```bash
ls -lh ncbi_dataset.zip
```
### decompress the ZIP file

```bash
unzip ncbi_dataset.zip
```

### REFseq indexation
### check the path (be in the directory which contains the reference genome)

```bash
pwd
```
## load mapping modules

```bash
module load bwa-mem2/2.2.1
```
 ### Display the help menu for the bwa-mem2 command

```bash
bwa-mem2 
```
## run indexing

```bash
bwa index reference_genome.fasta
```

### Check content

```bash
ls 
```

## Create a “mapping” directory in the working directory in which mapping results will be stored

```bash
mkdir -p /home/rkouke/Projet_Riz/mapping
```
## Move to the “mapping” directory

```bash
cd /home/rkouke/Project_Riz/mapping
```
## run a mapping test on one of the samples

```bash
bwa-mem2 mem reference_genome.fasta path_sample.R1.reads.fastq.gz path_sample.R2.reads.fastq.gz -t 8 -o sample.sam
```

## Convert sam file to bam with samtools view
###  load samtools modules

```bash
module load samtools/1.18 
```

```bash
samtools view -b -o sample.bam sample.sam 
```
## Mapping statistics with samtools flagstat

```bash
samtools flagstat sample.bam > sample.flagstat.txt
```
## view the file generated

```bash
cat sample.flagstat.txt 
```
## Extract reads mapped to the reference genome using samtools view

### Display the help menu for the samtools view command

```bash
samtooms view --help
```
## run extraction
```bash
samtools view -f 3 --threads 8 -h -b sample.bam -o sample.mapped_paired.bam
```
## view stats

```bash
samtools flagstat sample.mapped_paired.bam
```
## sort the Bam file

```bash
samtools sort sample.mapped_paired.bam -o sample.mapped_paired.sorted.bam --threads 8
```
### Check content

```bash
ls -lh sample.*
```

## Indexing bam file

```bash
samtools index sample.mapped_paired.sorted.bam
```
### Check content

```bash
ls -lrt sample.*
```




