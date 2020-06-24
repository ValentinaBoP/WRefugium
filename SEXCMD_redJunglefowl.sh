#!/bin/bash -l 
#SBATCH -J SEXCMD_redJungle
#SBATCH -o SEXCMD_redJungle.output
#SBATCH -e SEXCMD_redJungle.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 12:00:00                           
#SBATCH -A snic2020-15-44   
#SBATCH -p core
#SBATCH -n 2

ml bioinfo-tools samtools bwa
conda deactivate
ml python/2.7.15
ml R_packages/3.6.0

SCRIPT=/proj/uppstore2018073/private/Valentina/2020Refugium/Code/SEXCMD
INDIR=/proj/uppstore2018073/private/Valentina/2020Refugium/Intermediate/GalGal/RedJungle
OUTDIR=/proj/uppstore2018073/private/Valentina/2020Refugium/Data/RedJungle/DNA
# sex_marker.galGal4.filtered.final.github.fasta downloaded from: https://github.com/lovemun/SEXCMD/tree/master/Examples/Chicken/sex_marker.galGal4.filtered.final.fasta
FINAL=$INDIR/sex_marker.galGal4.filtered.final.github.fasta

cd $OUTDIR

for INPUT in $( cat list_fastq )
do
 Rscript $SCRIPT/SEXCMD.R $FINAL 3 ZW $INPUT
done
