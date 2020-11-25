# dcomp
Sex chromosome dosage compensation across vertebrates

## External software used

1. [Subread package](http://subread.sourceforge.net) for RNAseq read alignments (version 2.0.1).
2. 

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



