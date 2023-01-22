#!/bin/bash -l 
#SBATCH -J star_1
#SBATCH -o star.output 
#SBATCH -e star.error 
#SBATCH --mail-user opalacios7@gmail.com 
#SBATCH --mail-type=ALL 
#SBATCH -t 00-04:00:00 
#SBATCH -A snic2019-3-671
#SBATCH -p core     
#SBATCH -n 10       
#SBATCH -M snowy   

# Load module
ml bioinfo-tools
ml star/2.7.2b

### This script will mapp the Emu RNA-seq reads to the droNov.fasta emu reference genome
reads=/proj/sllstore2017073/private/Octavio/Emu/RNA_seq_reads
brain=/proj/sllstore2017073/private/Octavio/Emu/RNA_seq_reads/female_brain
e_brain=/proj/sllstore2017073/private/Octavio/Emu/RNA_seq_reads/female_embryonic_brain

# Generating genome imdex
STAR --runThreadN 10 --runMode genomeGenerate \
 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping \
 --genomeFastaFiles /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/droNov.rename.fasta \
 --sjdbOverhang ReadLength-1

# Run STAR mappibng
STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $reads/droNov31_R1_val_1.fq.gz $reads/droNov31_R2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.ovary_1 --outSAMtype BAM SortedByCoordinate 

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $reads/droNov32_R1_val_1.fq.gz $reads/droNov32_R2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.ovary_2 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $reads/droNov33_R1_val_1.fq.gz $reads/droNov33_R2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.ovary_3 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $reads/droNov34_R1_val_1.fq.gz $reads/droNov34_R2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.ovary_4 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $reads/droNov35_R1_val_1.fq.gz $reads/droNov35_R2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.ovary_5 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $reads/droNov41_R1_val_1.fq.gz $reads/droNov41_R2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Male.testis_1 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $reads/droNov42_R1_val_1.fq.gz $reads/droNov42_R2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Male.testis_2 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $reads/droNov43_R1_val_1.fq.gz $reads/droNov43_R2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Male.testis_3 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $brain/SRR7026352_1_val_1.fq.gz $brain/SRR7026352_2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.brain_1 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $brain/SRR7026354_1_val_1.fq.gz $brain/SRR7026354_2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.brain_2 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $brain/SRR7026355_1_val_1.fq.gz $brain/SRR7026355_2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.brain_3 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $brain/SRR7026357_1_val_1.fq.gz $brain/SRR7026357_2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female.brain_4 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $e_brain/SRR787599_1_val_1.fq.gz $e_brain/SRR787599_2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female_embryonic_brain_1 --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping --readFilesIn $e_brain/SRR787601_1_val_1.fq.gz $e_brain/SRR787601_2_val_2.fq.gz --readFilesCommand zcat \
--outFileNamePrefix /proj/sllstore2017073/private/Octavio/Emu/STAR_mapping/Female_embryonic_brain_2 --outSAMtype BAM SortedByCoordinate
