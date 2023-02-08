#! /usr/bin/env Rscript

# Code modified from James Galbraith's script orfChecker.R https://github.com/jamesdgalbraith/OrthologueRegions/blob/master/orfChecker.R

# import libraries
library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(ORFik)

# get fasta names
filenames = list.files(pattern = "_CR1.fasta")

# get species name
species_name = sub(pattern = ".fasta", replacement = "", x = filenames)

for (i in 1:length(filenames)){

print(paste0("Analysing file: ", filenames[i]))

# read in CR1 sequences
raw_seq <- readDNAStringSet(filepath = filenames[i])

# create pseudoranges tibble of CR1s
raw_tbl <- tibble(seqnames = names(raw_seq), start = 1, end = width(raw_seq)) %>% tidyr::separate(seqnames, into = c("seqnames", "group_name"), sep = "__")

# create table of CR1 seqs with names
seq_group_tbl <- raw_tbl %>% dplyr::select(seqnames, group_name)

# find orfs over 600bp in sequences
#orfs_600 <- ORFik::findORFs(raw_seq, startCodon = startDefinition(1), minimumLength = 600) %>% as_tibble()

orfs_600 = findORFsFasta(filenames[i], startCodon = startDefinition(1), minimumLength = 600, is.circular = FALSE)
orfs_600 = as_tibble(orfs_600)
orfs_600$idv_name = paste0(orfs_600$seqnames, "#orf", 1:nrow(orfs_600))

# make ranges object of orfs
orfs_600_ranges <- GRanges(seqnames = orfs_600$seqnames, ranges = IRanges(start = orfs_600$start, end = orfs_600$end))

# get seq of orfs and name orfs
orfs_seq <- Biostrings::getSeq(raw_seq, orfs_600_ranges)
names(orfs_seq) <- orfs_600$idv_name

# translate orf
orfs_aa_seq <- translate(orfs_seq, if.fuzzy.codon = "solve")
names(orfs_aa_seq) <- orfs_600$idv_name

# write nt and nt orfs to file
writeXStringSet(orfs_aa_seq, paste0(species_name[i], "_aa_orfs.fa"))
writeXStringSet(orfs_seq, paste0(species_name[i], "_nt_orfs.fa"))

print("Running rpsblast")

# search for EN and RT in orfs
orfs_aa_rps <- read_tsv(system(paste0("ml bioinfo-tools blast/2.10.1+ && rpsblast -query ", species_name[i], "_aa_orfs.fa -db /proj/uppstore2018073/private/Valentina/2020Refugium/Data/Transposon_domains/transposons -outfmt \"6 qseqid sseqid qstart qend sstart send length qlen slen pident stitle\" -evalue 0.01 -num_threads 12"), intern = T), col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle"))

# split rpsblast title toi make usable
orfs_aa_rps <- orfs_aa_rps %>% separate(stitle, into = c("code", "name", "description"), sep = ", ")

# determine presence/absence of RT and EN
orfs_aa_rps_rt <- orfs_aa_rps %>% filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.9 * slen)
orfs_aa_rps_en <- orfs_aa_rps %>% filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.9 * slen)

if(nrow(orfs_aa_rps_rt) == 0 | nrow(orfs_aa_rps_en) == 0){
        print(paste0("No intact CR1s in: ", filenames[i]))
        next()
}

# give seqnames their original names without the orf ids
orfs_aa_rps_rt = orfs_aa_rps_rt %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")
orfs_aa_rps_en = orfs_aa_rps_en %>% tidyr::separate(qseqid, into = c("seqnames", "orf_name"), sep = "#")

# select orfs with intact rt and en
intact_orfs <- orfs_aa_rps_rt %>% filter(seqnames %in% orfs_aa_rps_en$seqnames) %>% dplyr::select(seqnames) %>% base::unique()

if(nrow(intact_orfs) == 0){
        print(paste0("No intact CR1s with EN + RT in ", filenames[i]))
        next()
  } else {
        print(paste0("Saving: ", species_name[i], "_intact_orfs.bed"))
        intact_orfs_bed = intact_orfs %>% tidyr::separate(seqnames, into = c("element", "coordinates"), sep = "::") %>% tidyr::separate(coordinates, into = c("chromosome", "coordinates"), sep = ":") %>% tidyr::separate(coordinates, into = c("start", "end"), sep = "-")
        intact_orfs_bed = intact_orfs_bed[,c(2:4,1)]
        write.table(x = as.data.frame(intact_orfs_bed), sep = "\t", quote = F, col.names = F, row.names = F, file = paste0(species_name[i], "_intact_orfs.bed"))
  }
}
