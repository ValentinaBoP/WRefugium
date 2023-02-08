#! /usr/bin/env Rscript

# Code modified from James Galbraith's script orfChecker.R https://github.com/jamesdgalbraith/OrthologueRegions/blob/master/orfChecker.R

# it takes a fasta file with the DNA sequences of the transposons

# import libraries
library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(ORFik)

# get fasta names
filenames = list.files(pattern = "_DNA.fasta")

# get species name
species_name = sub(pattern = ".fasta", replacement = "", x = filenames)

for (i in 1:length(filenames)){

print(paste0("Analysing file: ", filenames[i]))

# read in sequences
raw_seq <- readDNAStringSet(filepath = filenames[i])

# create pseudoranges tibble of sequences
raw_tbl <- tibble(seqnames = names(raw_seq), start = 1, end = width(raw_seq)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of seqs with names
seq_group_tbl <- raw_tbl %>% dplyr::select(seqnames, group_name)

# find orfs over 1000 bp in sequences
#orfs_600 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 600) %>% as_tibble()

orfs_1000 = findORFsFasta(filenames[i], startCodon = startDefinition(1), minimumLength = 1000, is.circular = FALSE)

if(length(orfs_1000) == 0 ){
        print(paste0("No ORFs found in ", filenames[i]))
        next()
}

orfs_1000 = as_tibble(orfs_1000)
orfs_1000$idv_name = paste0(orfs_1000$seqnames, "#orf", 1:nrow(orfs_1000))

# make ranges object of orfs
orfs_1000_ranges <- GRanges(seqnames = orfs_1000$seqnames, ranges = IRanges(start = orfs_1000$start, end = orfs_1000$end))

# get seq of orfs and name orfs
orfs_seq <- Biostrings::getSeq(raw_seq, orfs_1000_ranges)
names(orfs_seq) <- orfs_1000$idv_name

# translate orf
orfs_aa_seq <- translate(orfs_seq, if.fuzzy.codon = "solve")
names(orfs_aa_seq) <- orfs_1000$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq, paste0(species_name[i], "_aa_orfs.fa"))
writeXStringSet(orfs_seq, paste0(species_name[i], "_nt_orfs.fa"))

print("Running rpsblast")

# search for transposase in orfs
orfs_aa_rps = character()
orfs_aa_rps = try(read_tsv(system(paste0("ml bioinfo-tools blast/2.10.1+ && rpsblast -query ", species_name[i], "_aa_orfs.fa -db /proj/uppstore2018073/private/Valentina/2020Refugium/Data/Transposon_domains/transposons -outfmt \"6 qseqid sseqid qstart qend sstart send length qlen slen pident stitle\" -evalue 0.01 -num_threads 12"), intern = T), col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle")), silent = TRUE)
if(class(orfs_aa_rps) == "try-error"){
  print(paste0("no domains found in ", filenames[i]))
  next()
}

# split rpsblast title toi make usable
orfs_aa_rps <- orfs_aa_rps %>% separate(stitle, into = c("code", "name", "description"), sep = ", ", extra = "drop")

# determine presence/absence of transposase domain
orfs_aa_rps_rt <- orfs_aa_rps %>% filter(grepl("transposase", description, ignore.case = T), (send - sstart + 1) >= 0.9 * slen)

if(nrow(orfs_aa_rps_rt) == 0 ){
        print(paste0("No intact DNA transposons in: ", filenames[i]))
        next()
}

# give seqnames their original names without the orf ids
orfs_aa_rps_rt = orfs_aa_rps_rt %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs <- orfs_aa_rps_rt %>% dplyr::select(seqnames) %>% base::unique()

if(nrow(intact_orfs) == 0){
        print(paste0("No intact DNA transposons in ", filenames[i]))
        next()
  } else {
        print(paste0("Saving: ", species_name[i], "_intact_orfs.bed"))
        intact_orfs_bed = intact_orfs %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
        intact_orfs_bed = intact_orfs_bed[,c(2:4,1)]
        write.table(x = as.data.frame(intact_orfs_bed), sep = "\t", quote = F, col.names = F, row.names = F, file = paste0(species_name[i], "_intact_orfs.bed"))
  }
}
