# dcomp
Sex chromosome dosage compensation across vertebrates

## Software versions

1. Ubuntu 18.04.4 LTS (Bionic Beaver)
2. R version 3.6.3 (2020-02-29)
3. Perl v5.26.1
2. [Subread package](http://subread.sourceforge.net) for RNAseq read alignments (version 2.0.1).
4. [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) for count data processing (version 3.28.1)
5. [tidyverse](https://www.tidyverse.org) for data cleanup and manipulation (version 1.3.0)


## Reference genomes and annotations
1. Chicken: GCF_000002315.6_GRCg6a_genomic.fna.gz, GCF_000002315.6_GRCg6a_genomic.gtf.gz
2. Platypus: GCF_004115215.1_mOrnAna1.p.v1_genomic.fna.gz, GCF_004115215.1_mOrnAna1.p.v1_genomic.gtf.gz
3. Opossum: GCF_000002295.2_MonDom5_genomic.fna.gz, GCF_000002295.2_MonDom5_genomic.gtf.gz
4. Mouse: GCF_000001635.27_GRCm39_genomic.fna.gz, GCF_000001635.27_GRCm39_genomic.gtf.gz

Referernce genomes were masked for the Y chromosome sequences in Platypus and Mouse and for the W chromosome in Chicken to reduce false positive mappings to these chromosomes for the sex lacking these chromosomes. Opossum assembly is from a female sample lacking the Y chromosome.

```
perl sexmask.pl GCF_000002315.6_GRCg6a_assembly_report.txt GCF_000002315.6_GRCg6a_genomic.fna.gz GCF_000002315.6_GRCg6a_genomic.zo.fna.gz W
perl sexmask.pl GCF_004115215.1_mOrnAna1.p.v1_assembly_report.txt GCF_004115215.1_mOrnAna1.p.v1_genomic.fna.gz GCF_004115215.1_mOrnAna1.p.v1_genomic.xo.fna.gz Y
perl sexmask.pl GCF_000001635.27_GRCm39_assembly_report.txt GCF_000001635.27_GRCm39_genomic.fna.gz GCF_000001635.27_GRCm39_genomic.xo.fna.gz Y
```
GTF files for platypus and opossum contained rRNA and tRNA annotations with empty value for the `gene_id` tag required by `subread` program. Therefore GTF files were modified to add `gene_id` tag value from the `product` tag.

```
zcat GCF_000002295.2_MonDom5_genomic.gtf.gz | perl -lne 'if ($_=~ /gene_id \"\";.*product \"(.*?)\"/){ $c++; $gid = "$1"; $_ =~ s/gene_id \"\";/gene_id \"$gid\";/; print $_} else{ print $_}' | gzip >GCF_000002295.2_MonDom5_genomic.mod.gtf.gz
zcat GCF_004115215.1_mOrnAna1.p.v1_genomic.gtf.gz | perl -lne 'if ($_=~ /gene_id \"\";.*product \"(.*?)\"/){ $c++; $gid = "$1"; $_ =~ s/gene_id \"\";/gene_id \"$gid\";/; print $_} else{ print $_ }' | gzip >GCF_004115215.1_mOrnAna1.p.v1_genomic.mod.gtf.gz
```

Genomes were indexed for read alignments. We chose parameters to create full index as a single block.
```
subread-buildindex -F -B -o reference_genome.fna reference_genome.fna.gz
```

## Read alignments

RNAseq reads were aligned against the species and **sex** specific reference genomes using the following command. Parameters were chosen such that (a) all possible subreads were chosen for mapping (`-n 300`) and (b) for multi-mapping reads, one random location was chosen for read assignments (`--multiMapping -B 1`). 

```
# Single-end format
subread-align -n 300 --multiMapping -B 1 -T 24 --rg-id [rgid] -a [XXX_genomic.gtf.gz] --sortReadsByCoordinates -i [XXX_genomic.fna] -r [inputreads.fq.gz] -t 0 -o [output.bam]
# Paired-end format
subread-align -n 300 --multiMapping -B 1 -T 24 --rg-id [rgid] -a [XXX_genomic.gtf.gz] --sortReadsByCoordinates -i [XXX_genomic.fna] -r [inputreads.R1.fq.gz] -R [inputreads.R2.fq.gz] -t 0 -o [output.bam]

```

## Feature counts

Features count was performed using the following command. Multiple files for a given species were processed simultaneously. Parameters were chosen such that (a) reads were assigned to all overlapping features if features were ovelapping (`-O`), (b) multi-mapping reads were counted (`-M`), and (c) strand-specificity of RNAseq was accounted for counting (`-s 0,1, or 2`). For paired-end data, fragments were counted instead of individual reads of a pair (`-p`).

```
featureCounts --byReadGroup -p -O -M -s 2 -a [XXX_genomic.gtf.gz] -o [outputcounts.txt] [Species-*.bam]
```

Strand parameter was `-s 2` for chicken and mouse samples, `-s 1` for opossum samples and `-s 0` for platypus samples.

## Counts data clean-up and processing

Counts data were processed to obtain normalised counts, merge specific metadata, counts, gene information using [dcanalysis.R](https://github.com/kango2/dcomp/blob/main/dcanalysis.R).

## Platypus PAR definitions

| chr   | par | start | end |
|:---:|:---:|---:|----:| 
| chrX1	| PAR1 | 1 | 45000000 |
| chrX2 | PAR1 | 1 | 5500000 |
| chrX2 | PAR2 | 20000000 | 29662106 |
| chrX3 | PAR1 | 1 | 18500000 |
| chrX3 | PAR2 | 31500000 | 33863336 |
| chrX4 | PAR1 | 1 | 5500000 |
| chrX5 | PAR1 | 1 | 8500000 |
| chrX5 | PAR2 | 69000000 | 70139320 |
