#!/bin/bash -l 
#SBATCH -J featureCounts_1
#SBATCH -o featureCounts.output 
#SBATCH -e featureCounts.error 
#SBATCH --mail-user opalacios7@gmail.com 
#SBATCH --mail-type=ALL 
#SBATCH -t 00-00:30:00 
#SBATCH -A snic2017-7-165
#SBATCH -p core     
#SBATCH -n 2       
#SBATCH -M snowy   

# Load module
ml bioinfo-tools
ml subread/2.0.0

# Summarize multiple paired-end datasets:
featureCounts -p -B -C -T 2 -A chrAliases.txt -t dispersed_repeat -g Target -a droNov.fasta.out.only.LTR_ERV.gff -o counts.only.LTR_ERV.txt Female.ovary_1Aligned.sortedByCoord.out.bam Female.ovary_2Aligned.sortedByCoord.out.bam Female.ovary_3Aligned.sortedByCoord.out.bam Female.ovary_4Aligned.sortedByCoord.out.bam Female.ovary_5Aligned.sortedByCoord.out.bam Female.brain_1Aligned.sortedByCoord.out.bam Female.brain_2Aligned.sortedByCoord.out.bam Female.brain_3Aligned.sortedByCoord.out.bam Female.brain_4Aligned.sortedByCoord.out.bam Female_embryonic_brain_1Aligned.sortedByCoord.out.bam Female_embryonic_brain_2Aligned.sortedByCoord.out.bam Male.testis_1Aligned.sortedByCoord.out.bam Male.testis_2Aligned.sortedByCoord.out.bam Male.testis_3Aligned.sortedByCoord.out.bam
