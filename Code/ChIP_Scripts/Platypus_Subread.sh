#!/bin/bash

#PBS -N Subread
#PBS -l nodes=1:ppn=16
#PBS -l mem=100gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M a.milton@unsw.edu.au
#PBS -m ae

module load subread/2.0.1

#F_H3K27me3
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/F_H3K27me3.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/F_H3K27me3.clean.subread.bam
#F_H3K9me2
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/F_H3K9me2.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/F_H3K9me2.clean.subread.bam
#F_H3K9me3
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/F_H3K9me3.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/F_H3K9me3.clean.subread.bam
#F_H4K20me1
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/F_H4K20me1.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/F_H4K20me1.clean.subread.bam
#F_RNA-polII
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/F_RNA-polII.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/F_RNA-polII.clean.subread.bam


#M_H3K27me3
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/M_H3K27me3.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/M_H3K27me3.clean.subread.bam
#M_H3K9me2
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/M_H3K9me2.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/M_H3K9me2.clean.subread.bam
#M_H3K9me3
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/M_H3K9me3.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/M_H3K9me3.clean.subread.bam
#M_H4K20me1
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/M_H4K20me1.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/M_H4K20me1.clean.subread.bam
#M_RNA-polII
subread-align -P 6 -t 1 -i /srv/scratch/waters/Genomes/Platypus/SubreadIndex/OrnAna1 \
-r /srv/scratch/z5160329/OAN_ChIP_Ash/Trimmed_Reads/M_RNA-polII.clean.trimmed.fq.gz \
-o /srv/scratch/z5160329/OAN_ChIP_Ash/Subread/M_RNA-polII.clean.subread.bam
