
# Project Riz

### This project focuses on the analysis of genetic diversity in rice after irradiation.

### Folder structure

- **01.RawData/** : Raw sequencing data.
- **metadata/** : Sample information and morphological data 
- **qc/** : Quality control results.
- **fastqc/** : Quality control results.
- **multiqc/** : Quality control results.
- **cleaning/** : to eliminate poor quality bases and adapters if used

##  FASTQC

#### Create a qc/ directory in the working directory in which quality control results will be stored.

```bash
mkdir -p /home/rkouke/Project_Riz/qc
```
#### Move to qc directory

```bash
cd /home/rkouke/Project_Riz/qc
```

#### Create the fastqc/ directory for outputs

```bash
mkdir fastqc
```

#### Move to fastqc/ directory

```bash
cd fastqc
```

#### Load the fastqc/ module, taking the version into account

```bash
module load fastqc/0.12.1
```
#### Display the help menu for the fastqc command

```bash
fastqc --help
```

#### Run quality control command on fastq.gz sequencing data

```bash
fastqc /home/rkouke/01.RawData/*.fastq.gz -o /home/rkouke/Projet_Riz/qc/fastqc
```

#### Check content
```bash
ls -lh
```

#### Use generated html files to check read quality 

##  MULTIQC

#### Create the multiqc/ directory for outputs 

```bash
mkdir -p /home/rkouke/Projet_Riz/qc/multiqc
```

#### Move to the multiqc/ directory

```bash
cd /home/rkouke/Projet_Riz/qc/multiqc
```

#### Load multiqc module to compress fastqc reports into a single report

```bash
module load multiqc/1.13
```

#### Launch analysis

```bash
multiqc /home/rkouke/Projet_Riz/qc/fastqc -o /home/rkouke/Projet_Riz/qc/multiqc
```

#### Check content
```bash
ls -lh
```

#### Use generated html files to check read quality

## CLEANING

#### cleaning depends on multiqc results





