#!/bin/bash -l 
#SBATCH -J LTR_DIGEST_LycPyr
#SBATCH -o LTR_DIGEST_LycPyr.output
#SBATCH -e LTR_DIGEST_LycPyr.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 10:00:00                           
#SBATCH -A snic2020-15-57  
#SBATCH -p core
#SBATCH -n 4

ml bioinfo-tools GenomeTools/1.5.9 hmmer

SPECIES=LycPyr
GENOME=lycPyr7.4.fasta
NAME=${GENOME%.*}
HARVEST=$NAME.ltrharvest.out

mkdir /proj/uppstore2018073/private/Valentina/2020Refugium/Intermediate/FL/LTR_Digest/$SPECIES/
cd /proj/uppstore2018073/private/Valentina/2020Refugium/Intermediate/FL/LTR_Digest/$SPECIES/
ln -s /proj/sllstore2017073/private/BirdReferences/lycPyr7.4.fasta

#Index the genome file
gt suffixerator -db $GENOME -indexname $GENOME -tis -suf -lcp -des -ssp -sds -dna

#Run LTRharvest:
gt ltrharvest -index $GENOME -gff3 $GENOME.gff -out $GENOME.ltr.fa
gt ltrharvest -index $GENOME -seqids yes -tabout no > $HARVEST

#LTRdigest
gt gff3 -sortlines yes -retainids yes -tidy yes -fixregionboundaries yes -checkids yes $GENOME.gff > $GENOME.sorted.gff

mv $GENOME.sorted.gff $GENOME.gff

# HMM profiles downloaded from GyDB and pfam
gt ltrdigest -hmms /proj/uppstore2018073/private/Valentina/2020Refugium/Data/HMM/*hmm -aaout yes -outfileprefix ${GENOME}_ltrdigest $GENOME.gff $GENOME > ${GENOME}_ltrdigest_output_gff

# Remove false positive
gt select -rule_files /proj/uppstore2018073/private/Valentina/2020Refugium/Code/LTR_Digest/filter_protein_match.lua -- < ${GENOME}_ltrdigest_output_gff > ${GENOME}_ltrdigest_output_gff2


#GFF=lycPyr7.4.fasta_ltrdigest_output_gff2
#grep '#' $GFF | grep -v '##' | cut -c2- > temp
#awk '{a = NR-1; print $1, "seq"a}' OFS="\t" temp > conversion
#rm temp

#grep -v '#' $GFF > 2convert.gff2

#R

#conv = read.table("conversion", stringsAsFactors=F, header = F)
#data = read.table("2convert.gff2", stringsAsFactors=F, header =F)

#for(i in 1:nrow(conv)){
  
#  data$V1[data$V1 == conv$V2[i]] = conv$V1[i]
  
#}

#gff = "lycPyr7.4.fasta_ltrdigest_output_gff3"
## this gff3 table is given as Supplementary Table in the paper
#write.table(x = data, file = gff, sep = "\t", quote = F, col.names = F, row.names = F)

#boo = data$V3 == "LTR_retrotransposon"
#bed = data[boo, c(1,4,5,9,7)]
#bed$V9 = sub(pattern = ";.*", replacement = "", x = sub(pattern = "ID=", replacement = "", x = bed$V9))

#bedfile="lycPyr7.4.fasta_ltrdigest_output_bed"

#write.table(x = bed, file = bedfile, sep = "\t", quote = F, col.names = F, row.names = F)

#system("ml bioinfo-tools BEDTools && bedtools getfasta -fi lycPyr7.4.fasta -fo lycPyr7.4.fasta_ltrdigest_output_bed.fasta -bed lycPyr7.4.fasta_ltrdigest_output_bed -name -s")

