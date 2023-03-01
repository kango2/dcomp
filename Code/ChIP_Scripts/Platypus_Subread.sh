#!/bin/bash

#PBS -N Subread
#PBS -l nodes=1:ppn=16
#PBS -l mem=100gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M myemailadress@generic.com
#PBS -m ae

module load subread/2.0.1

#F_H3K27me3
subread-align -P 6 -t 1 -i OrnAna1 \
-r F_H3K27me3.clean.trimmed.fq.gz \
-o F_H3K27me3.clean.subread.bam
#F_H3K9me2
subread-align -P 6 -t 1 -i OrnAna1 \
-r F_H3K9me2.clean.trimmed.fq.gz \
-o F_H3K9me2.clean.subread.bam
#F_H3K9me3
subread-align -P 6 -t 1 -i OrnAna1 \
-r F_H3K9me3.clean.trimmed.fq.gz \
-o F_H3K9me3.clean.subread.bam
#F_H4K20me1
subread-align -P 6 -t 1 -i OrnAna1 \
-r F_H4K20me1.clean.trimmed.fq.gz \
-o F_H4K20me1.clean.subread.bam
#F_RNA-polII
subread-align -P 6 -t 1 -i OrnAna1 \
-r F_RNA-polII.clean.trimmed.fq.gz \
-o F_RNA-polII.clean.subread.bam


#M_H3K27me3
subread-align -P 6 -t 1 -i OrnAna1 \
-r M_H3K27me3.clean.trimmed.fq.gz \
-o M_H3K27me3.clean.subread.bam
#M_H3K9me2
subread-align -P 6 -t 1 -i OrnAna1 \
-r M_H3K9me2.clean.trimmed.fq.gz \
-o M_H3K9me2.clean.subread.bam
#M_H3K9me3
subread-align -P 6 -t 1 -i OrnAna1 \
-r M_H3K9me3.clean.trimmed.fq.gz \
-o M_H3K9me3.clean.subread.bam
#M_H4K20me1
subread-align -P 6 -t 1 -i OrnAna1 \
-r M_H4K20me1.clean.trimmed.fq.gz \
-o M_H4K20me1.clean.subread.bam
#M_RNA-polII
subread-align -P 6 -t 1 -i OrnAna1 \
-r M_RNA-polII.clean.trimmed.fq.gz \
-o M_RNA-polII.clean.subread.bam
