#!/bin/bash -l
#SBATCH -J whatGene_taeGut
#SBATCH -o whatGene_taeGut.output
#SBATCH -e whatGene_taeGut.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 04:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 20

# whatGene pipeline

# copy files to SNIC_TMP
cp /proj/uppstore2018073/private/Valentina/2020Refugium/Intermediate/TaeGut/BWA/*.bam $SNIC_TMP

cd $SNIC_TMP

# set variables
DIR=/proj/uppstore2018073/private/Valentina/2020Refugium/Intermediate/TaeGut/SNPcalling/

# load modules
ml bioinfo-tools samtools/1.10

# FEMALES merge bam files and sort them
samtools merge -u - taeGut11_ERV.bam taeGut12_ERV.bam taeGut13_ERV.bam taeGut14_ERV.bam | samtools sort -T aln.sorted - -o gdna_plusb.bam -O BAM
samtools index gdna_plusb.bam

# MALES merge bam files and sort them
samtools merge -u - taeGut21_ERV.bam taeGut22_ERV.bam taeGut23_ERV.bam taeGut24_ERV.bam | samtools sort -T aln.sorted - -o gdna_zerob.bam -O BAM
samtools index gdna_zerob.bam

samtools merge -u - taeGut31_ERV.bam taeGut32_ERV.bam | samtools sort -T aln.sorted - -o rna_plusb.bam -O BAM
samtools index rna_plusb.bam
cp rna_plusb.b* $DIR

samtools sort -T aln.sorted -o rna_zerob.bam -O BAM taeGut41_ERV.bam
samtools index rna_zerob.bam
cp rna_zerob.b* $DIR

cp gdna_zerob* $DIR
cp gdna_plusb* $DIR

# load modules
conda deactivate
conda activate /proj/uppstore2018073/private/Valentina/Conda/SNPcalling

cd $DIR
REF=lycPyr2_rm2.1_merged.ERV.fasta

ls *.bam > ListOfBams

./bam_var_join.py $REF ListOfBams

./snp_calling_bchr.py baba.txt toico3.txt
