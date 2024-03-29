#!/bin/bash

#PBS -N Subread_separate_HiC
#PBS -l nodes=1:ppn=16
#PBS -l mem=100gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M myemailadress@generic.com
#PBS -m ae

module load subread/2.0.1


#Platypus 5-1 (male)
subread-align -t 1 -i OrnAna1 \
-r mpimg_L20932-1_Platypus-5-1_S56_R1_001.fastq.gz \
--SAMoutput \
-o HiC_Platypus_5-1_male_forward.subread.sam

subread-align -t 1 -i OrnAna1 \
-r mpimg_L20932-1_Platypus-5-1_S56_R2_001.fastq.gz \
--SAMoutput \
-o HiC_Platypus_5-1_male_reverse.subread.sam

#Platypus 5-2 (male)
subread-align -t 1 -i OrnAna1 \
-r mpimg_L20933-1_Platypus-5-2_S57_R1_001.fastq.gz \
--SAMoutput \
-o HiC_Platypus_5-2_male_forward.subread.sam

subread-align -t 1 -i OrnAna1 \
-r mpimg_L20933-1_Platypus-5-2_S57_R2_001.fastq.gz \
--SAMoutput \
-o HiC_Platypus_5-2_male_reverse.subread.sam

#Platypus 5-3 (male)
subread-align -t 1 -i OrnAna1 \
-r mpimg_L20934-1_Platypus-5-3_S58_R1_001.fastq.gz \
--SAMoutput \
-o HiC_Platypus_5-3_male_forward.subread.sam

subread-align -t 1 -i OrnAna1 \
-r mpimg_L20934-1_Platypus-5-3_S58_R2_001.fastq.gz \
--SAMoutput \
-o HiC_Platypus_5-3_male_reverse.subread.sam

#Platypus 5-4 (male)
subread-align -t 1 -i OrnAna1 \
-r mpimg_L20935-1_Platypus-5-4_S59_R1_001.fastq.gz \
--SAMoutput \
-o HiC_Platypus_5-4_male_forward.subread.sam

subread-align -t 1 -i OrnAna1 \
-r mpimg_L20935-1_Platypus-5-4_S59_R2_001.fastq.gz \
--SAMoutput \
-o HiC_Platypus_5-4_male_reverse.subread.sam
