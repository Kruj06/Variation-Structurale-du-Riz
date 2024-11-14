# Variation-Structurale-du-Riz
# __BIOINFORMATICS ANALYSIS STRATEGY - 

KOUKE R. Jaures, Broline AMURI

Jupyter inspired by the model created by C. Tranchant (DIADE-IRD), J. Orjuela (DIADE-IRD), F. Sabot (DIADE-IRD) and A. Dereeper (PHIM-IRD)
***

# <span style="color: #006E7F">Table of contents</span>
<a class="anchor" id="home"></a>


[PRACTICE I - Reads Quality Control](#quality) 
   * [Fastqc  `fastqc`](#quality)
   * [Multiqc `multiqc`](#quality)

[PRACTICE II - Mapping on all samples](#mapping) 
   * [Reference indexation  `bwa-mem2 index`](#refindex)
   * [Mapping with `bwa-mem2 mem`](#bwamem2-cmd)
   * [Convert sam to bam `samtools view`](#samtoolsview)
   * [Mapping statistics `samtools flagstat`](#flagstats)
   * [Filter correctly mapped reads `samtools view`](#corrmap)
   * [Sorting out filtered reads `samtools sort`](#sort)
   * [File indexing bam sorted  `samtools index`](#indexbam) 

[PRACTICE III- SNP calling](#SNP)
   * [Reference genome indexing `samtools faidx`](#refindex)
   * [Generate bcf file `bcftools mpileup`](#bcftools)
   * [Calling `bcftools call`](#calling)
   * [SNP statistiques `bcftools stats`](#stats)
   * [SNP annotation](#annotation)
   
   

***
