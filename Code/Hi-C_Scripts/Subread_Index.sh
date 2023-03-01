#!/bin/bash

#PBS -N Subjunc
#PBS -l nodes=1:ppn=16
#PBS -l mem=100gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M myemailadress@generic.com
#PBS -m ae

module load subread/2.0.1

subread-buildindex -o OrnAna1 GCF_004115215.2_mOrnAna1.pri.v4_genomic.fa
