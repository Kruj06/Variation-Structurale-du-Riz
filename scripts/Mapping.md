# Folder structure
# Mapping 
- **mapping /** : to align reads on the reference genome
- **mapping test/**: run a mapping test on one of the samples
- **mapping/** : run mapping on all samples using the for loop

  # Mapping test

#### indexing of the reference genome

#### Create a ref directory in the working directory in which reference genome will be stored

```bash
mkdir -p /home/rkouke/Projet_Riz/ref
```
#### Move to the ref directory

```bash
cd /home/rkouke/Project_Riz/ref
```
#### download reference genome
#### ncbi-datasets-cli module loading

```bash
module load ncbi-datasets-cli/16.27.2
```
#### download launch 
```bash
datasets download genome GCF_000147395.1 --include gff3,rna,cds,protein,genome,seq-report
```
#### Check download 

```bash
ls -lh ncbi_dataset.zip
```
#### decompress the ZIP file

```bash
unzip ncbi_dataset.zip
```

#### REFseq indexation
#### check the path (be in the directory which contains the reference genome)

```bash
pwd
```
#### load mapping modules

```bash
module load bwa-mem2/2.2.1
```
 #### Display the help menu for the bwa-mem2 command

```bash
bwa-mem2 
```
#### run indexing

```bash
bwa index reference_genome.fasta
```

#### Check content

```bash
ls 
```

#### Create a mapping/ directory in the working directory in which mapping results will be stored

```bash
mkdir -p /home/rkouke/Projet_Riz/mapping
```
#### Move to the mapping/ directory

```bash
cd /home/rkouke/Project_Riz/mapping
```
#### run a mapping test on one of the samples

```bash
bwa-mem2 mem reference_genome.fasta path_sample.R1.reads.fastq.gz path_sample.R2.reads.fastq.gz -t 8 -o sample.sam
```

#### Convert sam file to bam with samtools view
####  load samtools modules

```bash
module load samtools/1.18 
```

```bash
samtools view -b -o sample.bam sample.sam 
```
#### Mapping statistics with samtools flagstat

```bash
samtools flagstat sample.bam > sample.flagstat.txt
```
#### view the file generated

```bash
cat sample.flagstat.txt 
```
#### Extract reads mapped to the reference genome using samtools view

#### Display the help menu for the samtools view command

```bash
samtooms view --help
```
#### run extraction
```bash
samtools view -f 3 --threads 8 -h -b sample.bam -o sample.mapped_paired.bam
```
#### view stats

```bash
samtools flagstat sample.mapped_paired.bam
```
#### sort the Bam file

```bash
samtools sort sample.mapped_paired.bam -o sample.mapped_paired.sorted.bam --threads 8
```
#### Check content

```bash
ls -lh sample.*
```

#### Indexing bam file

```bash
samtools index sample.mapped_paired.sorted.bam
```
#### Check content

```bash
ls -lrt sample.*
```

