# Rpsblast to find ORFs in CR1 and DNA transposons

# working directory
cd /crex/proj/uppstore2018073/private/Valentina/2020Refugium/Intermediate/FL
mkdir DNA
mkdir CR1

# load modules
ml bioinfo-tools BEDTools R_packages/3.6.0 blast/2.10.1+

# CR1 identification
cd CR1

# link genomes
for SPECIES in galGal6a lycPyr7.4 droNov taeGut strHab calAnn
do
 ln -s /proj/uppstore2018073/private/Valentina/2020Refugium/Data/Genomes/$SPECIES.fasta ./$SPECIES.fasta
done

# link .out files
for SPECIES in galGal6a lycPyr7.4 droNov taeGut strHab calAnn
do
 ln -s /crex/proj/uppstore2018073/private/Valentina/2020Refugium/Intermediate/RMSK/bird_merged/$SPECIES.fasta.out ./$SPECIES.out
done

# convert .out in .bed, select for CR1 only, select for DNA transposons only and get sequences from the BED file
for SPECIES in galGal6a lycPyr7.4 droNov taeGut strHab calAnn
do
 /proj/sllstore2017073/private/scripts/RM2BED.sh -i $SPECIES.out -o $SPECIES.bed
 grep 'CR1' $SPECIES.bed > ${SPECIES}_CR1.bed
 bedtools getfasta -fi $SPECIES.fasta -fo ${SPECIES}_CR1.fasta -bed ${SPECIES}_CR1.bed -name
 grep 'DNA' $SPECIES.bed > ../DNA/${SPECIES}_DNA.bed
 bedtools getfasta -fi $SPECIES.fasta -fo ../DNA/${SPECIES}_DNA.fasta -bed ../DNA/${SPECIES}_DNA.bed -name
done

# make db for conserved domains
cd /proj/uppstore2018073/private/Valentina/2020Refugium/Data/
tar -zxvf transposon_domains.tar.gz
mv transposon_domains Transposon_domains
cd Transposon_domains
makeprofiledb -in transposase_domains.pn -out transposons -threshold 9.82 -scale 100.0 -dbtype rps -index true

cd -

Rscript --vanilla /proj/uppstore2018073/private/Valentina/2020Refugium/Code/orfCR1.R

cd ../DNA/

Rscript --vanilla /proj/uppstore2018073/private/Valentina/2020Refugium/Code/orfDNA.R
